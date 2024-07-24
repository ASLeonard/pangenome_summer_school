# Practical: Identification of small and structural variants.
Day 3: 4:00pm – 6:00pm

## Objectives
 - Call variation from genome-to-reference linear alignments
 - Call variation directly from pangenome alignments
 - Compare accuracy of variant calling

## Data and tools

Input
 - HiFi alignments
 - Assembled genomes
 - Pangenome

Tools
 - [vg](https://github.com/vgteam/vg)
 - [sniffles2](https://github.com/fritzsedlazeck/Sniffles)
 - [jasmine](https://github.com/mkirsche/Jasmine)
 - [bcftools](https://github.com/samtools/bcftools)

## Tasks

### Finding pangenome variation

We want to identify variation between the assemblies (and reference) directly from the pangenome.
We can do this through `vg deconstruct`, going from a pangenome _.gfa_ (with relevant path information) to a variant _.vcf_.
For this process, we need to specify **which** set of paths we consider to be the reference, as VCF requires a linear coordinate system.
We also want to specify the ploidy as 1, because even though we mostly work on diploid samples, each assembly represents only 1 copy of a genome.

```
# we may need to install some extra tools first
# however, vg is only available on linux
# can use conda instead of mamba if you prefer

mamba install -c conda-forge -c bioconda vg vcflib vcfbub


# we can also download the GFA from yesterday if we didn't reach that point
curl https://polybox.ethz.ch/index.php/s/rymYdwM8wxUfm1V/download > pangenome/bovines_with_P_lines.gfa

##LINUX ONLY
vg deconstruct -p "HER" --path-traversals --ploidy 1 -t 2 pangenome/bovines_with_P_lines.gfa > pangenome/minigraph.vcf
vcfbub -l 0 -a 10000 --input pangenome/minigraph.vcf | vcfwave -I 1000 -t 4 > pangenome/minigraph.decom.vcf


#alternatively download the file that it would have produced if vg won't run for you
curl https://polybox.ethz.ch/index.php/s/gMjAk0F3UZQP0UA/download > pangenome/minigraph.decom.vcf.gz
bgzip -d pangenome/minigraph.decom.vcf.gz
```

Currently, deconstructing tends to output too many variants, especially if the graph is not well formed.
One such example is
>HER     211583  >17>19  ACATATATATATGTATATATGTATATATATGTATATGTGTATATATATATATATATATATATATATAC    ATATAT  60      .       AC=2;AF=0.4;AN=5;AT=>17>18>19,>17>3240>19;NS=5;LV=0     GT      0       1       0       0       1

Where the reference starts "ACAT..." while the deletion alternate starts "ATAT...".
The initial "A" is common to both, while the next base is really a SNP (C→T), so this should be split into a SNP and a slightly smaller deletion.
However, the genotyping gets messed up, those "missing" samples **should** have the SNP.
All these tools are generally only semi-stable, so be careful!

>HER     211583  >17>19_1        ACATATATATATGTATATATGTATATATATGTATATGTGTATATATATATATATATATATATA A       60      .       AC=2;AF=0.400000;AN=5;AT=>17>18>19;NS=5;LV=0;ORIGIN=HER:211583;LEN=62;TYPE=del  GT      0       1       0       0       1  
HER     211650  >17>19_2        C       T       60      . AC=2;AF=0.400000;AN=5;AT=>17>18>19;NS=5;LV=0;ORIGIN=HER:211583;LEN=1;TYPE=snp   GT      0       .       0       0       .


For the moment, let's only consider structural variants, so we want to filter out small variants.

```
bcftools view -i 'abs(ILEN)>=50' -o pangenome/minigraph.SV.vcf pangenome/minigraph.decom.vcf
```

Remember as well that minigraph primarily will give us structural variants.
Rerunning these commands with a pggb graph would give us far more variants, but also much longer to run.

### Examining pangenome variation

In the VCF output of `vg deconstruct`, we also get some information on which nodes are involved in the variation.

We can pick some structural variants from the VCF and find the node IDs, and then visualise these regions in Bandage.
Can we confirm the genotypes make sense and follow the path information that we would expect?

We can also go in the opposite direction and find by visual inspection some interesting region in Bandage.
We can then search the pangenome VCF to see if that node appears and use that to follow up.

We also can look for regions with many missing genotypes (`.`) in the VCF.
How does the graph look in these regions?
Is it a simple mistake or does it look like an issue with our graph?

### Comparing to linear reference variation

We can compare against a more conventional approach, which we will take as the "truth" (although it has its own weaknesses).
We'll call SVs with sniffle2 using both the HiFi and ONT reads from the same individual.

```
sniffles --input data/OxO.HiFi.25.bam --vcf OxO.HiFi.sniffles.vcf --reference data/ARS-UCD1.2.fa.gz
sniffles --input data/OxO.ONT.25.bam --vcf OxO.ONT.sniffles.vcf --reference data/ARS-UCD1.2.fa.gz
```

Now we can check how many of the SVs we found through the graph and through the linear-reference approach are "the same".
We can also play around with parameters to determine how strict we want to be when discussing "the same" SV.
We'll use jasmine for this

```
jasmine --comma_filelist file_list=pangenome/minigraph.SV.vcf,OxO.HiFi.sniffles.vcf,OxO.ONT.sniffles.vcf threads=1 out_file=pangenome/SV_concordance.vcf genome_file=data/ARS-UCD1.2.fa.gz --pre_normalize --ignore_strand --allow_intrasample --normalize_type max_dist=50 max_dist_linear=.5 min_seq_id=0.5
```

We can then calculate the concordance with (you may want `grep -oE "SUPP_VEC=\d+"` instead on mac)

```
grep -oP "SUPP_VEC=\K\d+" pangenome/SV_concordance.vcf | sort | uniq -c > pangenome/SV_concordance.csv
```

which will tell us how many SVs are found across different samples.
Looking at regions where there is less agreement between the different variant sets may be useful to identify problematic regions in the graph, **or** regions where linear reference alignments of even long reads fail!

We can then run the following python code to make an UpSet plot.
We should be able to start an interactive python session by just typing `python` in the terminal.

```python
import matplotlib.pyplot as plt
import upsetplot

def parser(fname, names):
    data, memberships = [], []
    for line in open(fname):
        parts = line.split()
        data.append(int(parts[0]))
        memberships.append([names[i] for i,j in enumerate(parts[1]) if int(j)])
    return memberships, data

upsetplot.plot(upsetplot.from_memberships(*parser('pangenome/SV_concordance.csv',['graph','HiFi','ONT'])))
plt.show()
```

We can experiment with different SV merging "strictness", requiring more/less overlap of SVs, more/less sequence identity, etc.

```
#more strict
'max_dist=1 max_dist_linear=.05 min_seq_id=0.95'

#more lenient
'max_dist_linear=1 max_dist=5000'
```

Which will tell us how many SVs are common to both sets and how many are private to each.
The approaches are extremely different (as well as technical properties like alignment length/quality), but theoretically we want as high an agreement and as few privates as possible.
Note, this method only compares the REF/ALT status of a variant, without telling us anything about how correct the per-sample genotyping was. This is a much more complex problem addressed through tools like [truvari](https://github.com/ACEnglish/truvari).

We can do something similar for the small variants, again filtering the graph output

```
bcftools view -v snps --write-index=tbi -o graph.SNPs.vcf.gz graph.vcf
```

and comparing against a pre-made VCF of DeepVariant SNP calls using bcftools

```
bcftools isec -n +1 graph.SNPs.vcf.gz DV.SNPs.vcf.gz | awk '{++c[$5]} END {for (k in c) {print k,c[k]}}'
```
And again we expect as high overlap as possible.
