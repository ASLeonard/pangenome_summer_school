# Practical: Hands on pangenome assembly of a chromosome.
Day 2: 2:00pm – 6:00pm

## Objectives
We should be able to
 - Assemble genomes from long read sequencing data
 - Calculate quality metrics for assemblies
 - Build a graph pangenome from a collection of assemblies

## Data and tools

Input
 - "Whole genome" sequencing
   - Illumina (short)
   - PacBio CLR (long)
   - PacBio HiFi (long)
   - Nanopore (long)
 - Reference genome

Tools
 - [IGV](https://github.com/igvteam/igv)
 - [minimap2](https://github.com/lh3/minimap2)
 - [samtools](https://github.com/samtools/samtools)
 - [mdbg](https://github.com/ekimb/rust-mdbg)
 - [ragtag](https://github.com/malonge/RagTag)
 - [minigraph](https://github.com/lh3/minigraph)
 - [bandage](https://github.com/asl/BandageNG)

## Tasks

### Setting up the basics

We will install conda/mamba to manage most of the software installation and virtual environments.
After we install `mamba`, we will create an environment (called *ognigenoma* inside the file) and then activate it.
```
curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
bash Miniforge3-$(uname)-$(uname -m).sh

mamba env create --file environment.yml
mamba activate ognigenoma
```

We also want to download and install IGV to inspect alignments and bandage to visualise pangenomes.


https://igv.org/doc/desktop/#DownloadPage/
https://github.com/asl/BandageNG/releases

### Exploring the raw data

First we need to gather some of the data we will be working
```
#Let's make a directory to store the data in
mkdir data

#Download a neatly named cattle reference, and then index it
curl https://polybox.ethz.ch/index.php/s/d641EnjNf0TWw5f/download > data/ARS-UCD1.2.fa.gz
samtools faidx data/ARS-UCD1.2.fa.gz

#Download the HiFi reads from Original Braunvieh cow with sample accession SAMEA7759028.
#These are only from chromosome 25 for simplicity.

https://polybox.ethz.ch/index.php/s/C7fEDiQoDweJKFU/download > data/OxO.HiFi.25.fq.gz

#Download some pre-aligned short and ONT reads for the same animal and generate the csi index for it
curl https://polybox.ethz.ch/index.php/s/OrFlXd0G54Ntx6n/download > data/OxO.Illumina.25.bam
samtools index -@ 2 -c data/OxO.Illumina.25.bam

curl https://polybox.ethz.ch/index.php/s/saR530eBS8YimRp/download > data/OxO.ONT.25.bam
samtools index -@ 2 -c data/OxO.ONT.25.bam

```

Now we want to align the reads to the reference with minimap2
```
minimap2 -a -x map-hifi -t 4 data/ARS-UCD1.2.fa.gz data/OxO.HiFi.25.fastq.gz | samtools sort - -@ 4 --write-index -T $TMPDIR -o data/xO.HiFi.25.bam
```

Where we are piping the SAM output of minimap2 into samtools to directly create our BAM alignment file.
For minimap2, the parameters are `-a` to output SAM with base-level alignment, `-x map-hifi` to use mapping parameters appropriate for HiFi reads, and `-t 4` to use 4 threads.
For samtools, we are coordinate-sorting the output with 4 threads (`-@ 4`) and creating the ".csi" index file on the fly (`--write-index`).

---

We can now load the alignments into **IGV** to inspect them to get a feel for the data.
In particular, examine the alignments near the start/end of the chromosome as well as looking for large, complex variants.

Also load the Illumina and ONT r9 data and see how the alignments compare for different reads from the same sample!

---

At this stage, we can also generate some simple analyses for our sample.
We can generate a "consensus" sequence, effectively changing reference bases with variants from the alignments.
We can also directly call variants from the alignments and get a set of SNPs.
```
samtools consensus -@ 4 -X hifi -r 25 -a -l 0 -o OxO.HiFi.25.asm_consensus.fa data/OxO.HiFi.25.bam
bcftools mpileup --threads 4 -Ou -r 25 -f data/ARS-UCD1.2.fa.gz data/OxO.HiFi.25.bam | bcftools call --threads 4 -mv --ploidy 2 --skip-variants indels --write-index=tbi -o OxO.HiFi.25.vcf.gz
```

These are basic approaches, but might help give us a rough understanding of what to expect.


### _De novo_ assembly of a cattle chromosome

Now we have a basic overview of our data, we can generate a _de novo_ assembly.
This uses only the sequencing reads, so there currently is no _bias_ introduced by existing references which may have errors or true biological differences.

```
mkdir assemblies
#there is an unresolved bug in mdbg that prevents reading gzipped fastq, so just need to extract it first
gunzip -dc data/OxO.HiFi.25.fastq.gz > data/OxO.HiFi.25.fastq
rust-mdbg -k 31 -l 24 --density 0.003 --bf --minabund 2 --prefix assemblies/OxO.HiFi.25.mdbg --threads 4 data/OxO.HiFi.25.fastq

## magic simplify?

gfatools asm -u assemblies/OxO.HiFi.25.mdbg.gfa > assemblies/OxO.HiFi.25.mdbg.unitgs.gfa
to_basespace --gfa assemblies/OxO.HiFi.25.mdbg.unitgs.gfa --sequences assemblies/OxO.HiFi.25.mdbg
gfatools gfa2fa assemblies/OxO.HiFi.25.mdbg.unitgs.gfa.complete.gfa > assemblies/OxO.HiFi.25.asm_mdbg.fa
```

Now we have our genome assembly, we want to quantify some basic stats about it.
Let's go ahead and calculate the contiguity (N50).
If the assembly was perfect, we would expect 1 contig and a length of ~42 Mb.
We can recalculate the N**G**50 using the expected genome length.
Does this improve or worsen the N50?

```
mkdir tools
curl https://raw.githubusercontent.com/lh3/calN50/master/calN50.js > tools/calN50.js
k8 calN50.js assemblies/OxO.HiFi.25.asm_mdbg.fa
k8 calN50.js -L 42350435 assemblies/OxO.HiFi.25.asm_mdbg.fa

```

We can also calculate the conserved gene completeness based on the genes expected in this phylogeny.
We'll manually install this tool due to unnecessary conflicts in conda.

```
git clone https://github.com/huangnengCSU/compleasm.git
(cd compleasm; pip install .) #temporarily moves into compleasm directory and installs

compleasm download cetartiodactyla -L <dir>
compleasm run -a <asm.fa> -o completeness -l cetartiodactyla -L <dir> -t 4
```

Remember, we only assembled chromosome 25, so we should expect to find only a small fraction of the USCOs (around 460 for a good cattle genome).

We can also re-align our assembly to the reference and inspect in IGV
```
minimap2 -a -x asm5 -t 4 data/ARS-UCD1.2.fa.gz assemblies/OxO.HiFi.25.asm_mdbg.fa | samtools sort - -@ 4 --write-index -o assemblies/OxO.HiFi.25.asm_mdbg.ARS-UCD1.2.bam
```

Now using the "asm5" preset rather than "map-hifi", as we are mapping an assembly and not reads.
Where does the assembly look the best and where does it look the worst?
How does this relate to the choice of assembler?

We can also scaffold the contigs, so they have the same direction and layout as the reference.
Note, this **is** based on the reference, and so we are at risk of losing out on biological differences as we force our assembly to look more like the reference.
For high quality assemblies, this risk should be small.

```
ragtag.py scaffold data/ARS-UCD1.2.fa.gz assemblies/OxO.HiFi.25.asm_mdbg.fa -o assemblies/scaffolding -t 4
sed 's/_RagTag//g' assemblies/scaffolding/ragtag.scaffold.fasta > assemblies/OxO.HiFi.25.asm_mdbg.scaffolded.fa
```

We can now take a look at how well scaffolded our genome is:
 - how many gaps are there?
 - what regions of the reference go uncovered?

Since it is also so quick to generate these genome assemblies, what happens if we play around with the parameters (e.g., `-k`, `-l`, etc.)?
It is often difficult to understand how these parameters shape the outcome...


We can also run `hifiasm`, a more resource-intensive assembler but likely to produce better results.

Let's build a hifiasm assembly and then we can repeat the quality assessment of the assembly.
```
hifiasm -t 4 -o assemblies/OxO.HiFi.25.asm_hifiasm --primary data/OxO.HiFi.fastq.gz
gfatools gfa2fa assemblies/OxO.HiFi.25.asm_hifiasm.p_ctg.gfa > assemblies/OxO.HiFi.25.asm_hifiasm.fa
```

### Construct a chromosome pangenome

Now we have our assembly, we can gather some other publically available assemblies.
These have been named to follow [panSN](https://github.com/pangenome/PanSN-spec), so look like `sample#haplotype#contig`.
This allows multiple assemblies from the same chromosome to still be distinguished, as well as grouping samples with haplotype-resolved assemblies (e.g., `ZEB#1#15` and `ZEB#2#15` being the two haplotypes for chromosome 15 for a Zebu sample).

```
mkdir pangenome
curl https://polybox.ethz.ch/index.php/s/j0Fjt63F1mG2XkB/download > pangenome/OBV.fa.gz
curl https://polybox.ethz.ch/index.php/s/mNyJqihPEdWPNow/download > pangenome/BSW.fa.gz
curl https://polybox.ethz.ch/index.php/s/x4q6tWW3WRp1rBF/download > pangenome/SIM.fa.gz
curl https://polybox.ethz.ch/index.php/s/qSGYJVnOMIeCbW8/download > pangenome/NEL.fa.gz
curl https://polybox.ethz.ch/index.php/s/97V6opMeRhEt6Ek/download > pangenome/WIS.fa.gz
```

To get a quick overview, we can use minigraph to build a predominantly structural variation pangenome, using the reference genome as a "backbone" for variation.
First, we also rename the reference to follow panSN-spec (using the Hereford breed as the sample name).
We then run minigraph with the assemblies, starting with the reference and then in order of increasing in genomic distance.

```

samtools faidx data/ARS-UCD1.2.fa.gz 25 | awk '$1~/^>/ {gsub(/>/,">HER#0#",$1)}1' | bgzip -@ 2 -c > pangenome/HER.fa.gz
minigraph -cxggs -L 50 -j 0.01 pangenome/{HER,BSW,OBV,SIM,NEL,WIS}.fa.gz > pangenome/bovines.gfa
```

In particular, pay attention to the log output and consider the following questions:
 - How many events are added to the graph in each iteration?
 - Does that match expectations based on divergence?
 - Since minigraph is quick to run, what happens if change around the *order* of assemblies?

We can explore some other properties of this graph using
```
gfatools stat pangenome/bovines.gfa
gfatools bubble pangenome/bovines.gfa
```

We can also retrace the paths each assembly should take through the graph, again with minigraph
```
for i in HER BSW OBV SIM NEL WIS
do
  minigraph -cxasm --call -t 4 pangenome/bovines.gfa pangenome/${i}.fa.gz > pangenome/${i}.bed
done
```

We can also confirm that "non-reference" sequence can be lost with minigraph, by checking the starting coordinate of our assemblies compared to the initial coordinates of the reference.
```
curl https://raw.githubusercontent.com/lh3/minigraph/master/misc/mgutils.js > tools/mgutils.js
paste pangenome/{HER,BSW,OBV,SIM,NEL,WIS}.bed | k8 tools/mgutils.js merge - > pangenome/bovines.mg.variants
```

We can also approximately add P-lines to the minigraph gfa using a modified version, which will enable more general use of the file.

```
curl https://raw.githubusercontent.com/ASLeonard/minigraph/master/misc/mgutils.js > tools/mgutils_patched.js
{ cat pangenome/bovines.gfa ; paste pangenome/{HER,BSW,OBV,SIM,NEL,WIS}.bed | k8 tools/mgutils_patched.js path <(echo -e "HER\nBSW\nOBV\nSIM\nNEL\nWIS") - ; } > pangenome/bovines_with_P_lines.gfa
```
---

Since this graph is not too complicated, it should be relatively easy to load and explore in bandage.
Take a look around the graph and see which regions look "tangled" or very simple.

We can also play around with drawing the entire graph, subgregions by node ID, colouring specific paths by name, etc.

We can also query the graph for specific sequences using minimap2 through the "Graph search" function, and then colour the graph by hits.

Let's search for this sequence relating to the *PIGQ* gene.
> GCTACATCCACCTCATGTCCCCCTTCATCGAGCACATCCTGTGGCACTTGGGTCTGTCGGCCTGCCTGGGCCTAACGGTCGCCCTGTCCATCCTGTCGGACATCATCTCCCTTCTCACCTTCCACATCTACTGCTTCTACGTCTACGGCGCCAGGTGGGTGTGCCGCCCCCCCCCCTACCCCCCGGCGGGCTGGCGCGTGCAGCCCCGGGCCCTGGGCGCAGACAGGCTCCCGGGCCGGGCGGGTGTGGGGTGGCCCCCACCCTGACGGGAGTGGTCTGCAGGCTCTACTCCCTGAAGATCCACGGCCTGTCCTCGCTGTGGCGCCTGTTCCGAGGGAAGAAGTGGAACGTCCTGCGCCAGCGCGTGGACTCCTGCTCCTACGACCTGGACCAGGTACCAGCTCCGTCCAGCCGGCCTTCGGCCTGCGCCGGCCCATGCAGTGGCCCGGGCATCACAAGGGGCAGCCTCGCTCCTGCTCCGCAGGCCCCAGGGGTTCAGGGAGCTGCAAGGGGCAGGGTGTATACGGGGCCATAACGGGGCAGAACTCGCCGCAGGCGGGGACAGACCCTGCTCCGGGAGGCTGCAGGACAGCAGGCGCTGGAGGTGAGTGCGCATGGGGCCTGGCTCGTCCTGTTAGCTTGTCCACACCCGCCACAGGGACCCTGAGTGCTGCCCTTCAGGGTGTTCTCTTCCTGCCTAATTTTTACCCCCCTCCTGTCTAGTGTATGGGCTGCACACCTCTCAGCTGGTCCTGTCTCTACCCACTCTCTGCCCTTGTGTCACCAGGTCTGGAGGGCAGGCATCCCTCCCAGGACTGGGGAGGGCTGTCCCAACTGGCATCAGGCAGCGGTGGGCGCTCCACCTCTGGACTGGCGGAGGGGCTATGGTCAGTGCAGGTGGGGGCGCCCCAGCGTGGGGCTCTTGCTGCTGCTGATGCTGGCATCCTTTCTCCTCGTAGCTGTTCATCGGGACCTTGCTCTTCACCATCCTGATTTTCCTGCTTCCTACCACGGCCCTGTACTACCTGGTGTTCACCCTGGTGAGCCGAGTGCCACGTGCGCAGAGTGGTGACTGCCGGCTCTGCCTACCAGCGGGCAGCTGTGGGCAGACCACCCACCTGTGCCAGCCACAACAGAGGGTGGCTGTTCTGGACTCTGTGTACTGGGCCCTGGAGTGGACAGGGGTCCACCCTGGGGGTACCACCCCCTACTGCCCCGGCATACGCCGCTGGATGCATGGAGGGCAGGGCTCTCCGTAGTGCGGGGGGGCGGAGTGTGGCCAGCAGAGGGGCGTGGTACCCAGGTGGCCAGGCCGAGGGGCCAGTGGCGTGTGTCCCTGGGTTGGAGTGCCCCGTCTTCTGAGTGGGAGGGGGTTTCGTGACTCCTGTGCTCAGGCCTCAGGTTCAGGGCTGGGCCCCTGCACGCCCCGGGGACCTGTCAGGCCGCTCAGCTCAGGCATGTGCAGCATGGCTCTGGCCACGGTCCCAGGCCCAAGACCCTGAGCACTGTCTGTGCTGCAGCTCCGGCTCCTAGTGGTCACCGTGCAGGGCCTGGTCCATCTGCTCGTGGACCTCATCAACTCGCTGCCGCTGTACTCGCTTGGCCTCCGGCTGTGCCGGCCCTACAGGCTTGCGGGTAGGTCTGCACGCTGGCCAGGGAGCACGCTTGGGGCCCTGGAAGCCGCGGTTTGTGTGGCAGGCAGGCAGGGAGCACAGGGCCCTGGTTCTGGGGCGTGAGCTTCAGGGGGAGCCCGCCTTCAGGTTCAGAGGCTGCACGCGCAGCCCTCCTGCCTCCTGAGCCTTTCCTGGGATGAGTGTGCAGTTTGCTCGGAAGCGTGCCTGCAGTGTCGCATGGGCAGTGCCTGTGTGTGGCTCAGAAGAGAGACGTGTCCCCTGTGGCTGTGAGGCTGTGTGCCTGCTGGCATCCCTCCACTGGGGGCTGGCCCCCCAGAAGACTGCTGCCTGCTGCGGCCACACTGGCCAGCCCGGCCACCCCCACTGGCCCCCCTGCAGGCAGAGGCCCCTGTCTGCATGCTGACATCAGCTGCAGGTTGGCCAGATGCCCCTGACGGGAGGCCGTCCCCACTGCAGGCACAGGGAGGCCAGACCCGACTGGTCAGAGTTGTGGTCCAGCCTACCCAGTGCAGGGTCTTTGAACAGTCGGGGAAACTGGGGGTTCCTGGCGTGGCCCCCGGGTGACCAGTGTGCCCAGGAGGCCCTGGGCTGGTTTGCAGCCAGCAGCTCAGCCCACCCCAGAGCTTGTGGCGGCTGCTCGCAGCCAGTTTGGAGTCAGGAGACATGAATGTGCGTGTTTGTGGGCTTTAGGCCTTGAGTCTGCCCACCGGATGTCCTCAAACCCTCTTTCCGAGAGCTTCCTCCGCTGACCCCCAGGTGGCCCAGGCTCCCCCAGTCCCCACCACACCTCAGCCCCCGGACCTGGCACCTCCCTGGCCCCGCCTGCGCTCAGGCCTGCGCAGTCCTTTTAACTTCCCAGGAAAGGAAAGGCTGCAGATGCCGCCTTTTCTCTCCACTGCGTTGCCTTGGGGCCTCTTCATCAATTGATAATTAGCTCCTATTTCTCCTTTTAAGAACCTTCCACCTTGGGGATCCCCAGGCATGGGAGCCTACACTGCAGGCTAGAGTGTGGTGGGTCCCCGAGTCAGCTGTGCTGGAAGGCACTGTGATGAGCAGGGTGCCCACCCCCCAGTGGGTGCTGTGTCCACTCTAGGAAACCCCTGCCCTGCTCACAGAGCAGCACAGCCTGCTGGCTCTTCCCTGCCTGCAGACCCTGGGCTTTGGTCTCAGTGCGGTTCTGCTGGAGCTCTGGTGAGGGGCGGCTGCTGCCCACAACGTGCAGCCCTCCAAGGGGCGCACTGCGGGGACCGCCCCGCCTCACCGGGCCTCACTCTCCTCACAGCCGGCGTGAAGTTCCGGGTCCTGGAGCACGAGGCTGGCCGGCCCCTGCGCCTCCTGATGCAG

---

Now we can also build a more comprehensive pangenome, using pggb.
There are several tools required to run pggb (wfmash, seqwish, smoothxg, odgi, GFAffix), and so can be challenging to install without running through singularity.
These steps are also more compute-intensive, and so here we can just use a pre-built graph generated from the command

```
##LINUX ONLY
mkdir pangenome/pggb
cat pangenome/{HER,BSW,OBV,SIM,NEL,WIS}.fa.gz > pangenome/25.fa.gz
samtools faidx pangenome/25.fa.gz

pggb -i pangenome/25.fa.gz -o pangenome/pggb -t 4 -s 50k -p 97
cp pangenome/pggb/*smooth.final.gfa pangenome/25.pggb.gfa
cp pangenome/pggb/*smooth.final.og pangenome/25.pggb.og
```

We can calculate some similar stats using gfatools again, and compare to the minigraph graph generated from the **same** assemblies.
Similarly, we can load this graph in bandage and attempt to visualise the pangenome.

We can extract a subset of that graph using odgi (linux-only)

```
##LINUX ONLY
odgi extract -i pangenome/25.pggb.og -
```

and then much more manageably visualise the subgraph.


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
All these tools are work-in-progress, so be careful!

>HER     211583  >17>19_1        ACATATATATATGTATATATGTATATATATGTATATGTGTATATATATATATATATATATATA A       60      .       AC=2;AF=0.400000;AN=5;AT=>17>18>19;NS=5;LV=0;ORIGIN=HER:211583;LEN=62;TYPE=del  GT      0       1       0       0       1  
HER     211650  >17>19_2        C       T       60      . AC=2;AF=0.400000;AN=5;AT=>17>18>19;NS=5;LV=0;ORIGIN=HER:211583;LEN=1;TYPE=snp   GT      0       .       0       0       .


For the moment, let's only consider structural variants, so we want to filter out small variants.

```
bcftools view -i 'abs(ILEN)>=50' -o pangenome/minigraph.SV.vcf pangenome/minigraph.decom.vcf
```

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
```python
import matplotlib.pyplot as plt
import upsetplot

def parser(fname, names):
    data, memberships = [], []
    for line in open(fname):
        parts = line.split()
        data.append(int(parts[0])/norm)
        memberships.append([names[i] for i,j in enumerate(parts[1]) if int(j)])
    return memberships, data

upsetplot.plot(upsetplot.from_memberships(*parser('pangenome/SV_concordance.csv',['graph','HiFi','ONT'])))
plt.show()
```

We can experiment with different SV merging "strictness", requiring more/less overlap of SVs, more/less sequence identity, etc.
```
stricter: 'max_dist=1 max_dist_linear=.05 min_seq_id=0.95'
lenient: 'max_dist_linear=1 max_dist=5000'
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


# Practical: Using a pangenome to identify a known functional variant.
Day 5: 2:00pm – 6:00pm

## Objectives
 - Identify annotated genes overlapping/near pangenomic bubbles
 - Identify pangenomic regions associated with binary phenotypes

## Data and tools

Input
 - Pangenome
 - Reference genome annotation

Tools
 - [odgi](https://github.com/pangenome/odgi)
 - [impg](https://github.com/pangenome/impg)
 - [bedtools](https://github.com/arq5x/bedtools2)
 - [gafpack](https://github.com/ekg/gafpack)

## Tasks

To start, we want to get the ARS-UCD2.0 annotation in gff format, giving us the location of genes.

```
curl https://polybox.ethz.ch/index.php/s/gnReyfSopENjpxP/download > data/ARS-UCD1.2.exons.bed
```

### Intersecting pangenome bubbles and genes

We will then test for intersections between the annotation and pangenome "bubbles" (indicating some larger level of variation) from the minigraph gfa.

```
gfatools bubble pangenome/bovines_with_P_lines.gfa | cut -f -3,12 | sed 's/HER#0#//' > pangenome/minigraph.bubbles.bed
bedtools intersect -a data/ARS-UCD1.2.exons.bed -b pangenome/minigraph.bubbles.bed > pangenome/minigraph_SV_overlaps.bed
```

We can also require more overlap with e.g., `-f 0.5` requiring that at least 50% of the exon is overlapped by an SV to be reported.

We could also approach this from the SV VCF we created, and check for potential events that are private to a subset of assemblies.
We'll start by subsetting the VCF using `bcftools view` with the `-s <samples>` option within a [process substitution](https://en.wikipedia.org/wiki/Process_substitution) ("<(...)") to avoid creating intermediate files on disk.

```
bedtools intersect -a data/ARS-UCD1.2.exons.bed -b <(bcftools view -s <asm1,asm2> graph.SVs.vcf) > specific_genes.bed
```

We can also try and convert the annotation (via BED) into the pangenome coordindates, and then load that into bandage so we can visualise the events better.
```
curl https://github.com/AnimalGenomicsETH/KITKAT/blob/main/scripts/translate_bed_to_graph.py > tools/translate_bed_to_graph.py
python translate_bed_to_graph.py pangenome/bovines.gfa data/ARS-UCD1.2.exons.bed > pangenome/ARS-UCD1.2.exons.graph.bed
```

---

Let's see if we can find any interesting gene-SV regions in bandage!

---

### Pangenome association testing

We can also look for potential variation that could be associated with a binary phenotype.
This is a more "experimental" approach, but much more flexible in what we can identify.

Let's start by purely analysing what is already present in the graph.
We want to find regions where the graph clearly contains two (or more) "haplotypes", which could be associated with a phenotype.

We can split the graph up into small chunks (say 100 Kb) using `odgi extract`, and the find the [Jaccard similarity score](https://en.wikipedia.org/wiki/Jaccard_index) between each pair of assemblies within this chunk.

```
bedtools makewindows -g <(awk '$1==25' data/ARS-UCD1.2.fa.gz.fai) -w 100000 > pangenome/100Kb.windows.bed
odgi extract -i pangenome/25.pggb.og -b pangenome/100Kb.windows.bed -t 4 -s
for g in *.og
do
  odgi similarity -t 4 -i $g | sed 's/:[0-9]\+-[0-9]\+//g' | awk -v S=${{g}} 'NR>1 && $2>$1 {{print S,$1,$2,$6}}'
done | pigz -p 4 > pangenome/jaccard.100Kb.csv
```

We can then examine regions where the phenotype-grouped samples are all similar but are very different across groups.

Alternatively, we can make a presence/absence matrix, again chunking the genome into small windows.

```
 odgi pav -M -i pangenome/25.pggb.og -b <(echo -e "HER#0#25\t1000000\t2000000")
```

Now we can take this a step further, and integrate population level data from short read alignments to the pangenome.
We want to gather the node-level coverage per sample after aligning with `vg giraffe`.

```
for S in samples
do
  gafpack -g pangenome/25.pggb.gfa -a $S.gaf -l -c | awk -v S=$S '!/#/ {print S,$1,$2}' > $S.cov
done
```

We can then test to see if there is any statistically significant variation in node coverage between assigned phenotype groups.

---

Minigraph analysis version.

```python
import regex
from collections import defaultdict

paths = {line.split()[1]:{int(N) for N in regex.split(r'\D',line.split()[2])[:-1] if N} for line in open('pangenome/bovines_with_P_lines.gfa') if line[0] == 'P'}

nodes_as_keys = defaultdict(list)
for K,V in paths.items():
    for node in V:
        nodes_as_keys[node].append(K)

nodes_private_to_pairings = defaultdict(list)
for K,V in nodes_as_keys.items():
    nodes_private_to_pairings[tuple(sorted(V))].append(K)
nodes_private_to_pairings.keys()
```
