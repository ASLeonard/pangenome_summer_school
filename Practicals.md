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

### Exploring the raw data

First we need to gather some of the data we will be working
```
wget reference
wget reads
wget short reads
```

Now we want to align the reads to the reference with minimap2
```
minimap2 -a -x map-hifi -t 4 <reference> <reads> | samtools sort - -@ 4 --write-index -o F1.bam
```

Where we are piping the SAM output of minimap2 into samtools to directly create our BAM alignment file.
For minimap2, the parameters are `-a` to output SAM, `-x map-hifi` to use mapping parameters appropriate for HiFi reads, and `-t 4` to use 4 threads.
For samtools, we are coordinate-sorting the output with 4 threads (`-@ 4`) and creating the ".bai" index file on the fly (`--write-index`).

We can now load the alignments into IGV to inspect them to get a feel for the data.
In particular, examine the alignments near the start/end of the chromosome as well as looking for large, complex variants.

At this stage, we can also generate some basic output for our sample using 
```
samtools consensus
bcftools mpileup -Ou -f reference.fa alignments.bam | bcftools call -mv -Ob -o calls.bcf
```

These are basic approaches, but might help give us a rough understanding of what to expect.


### _De novo_ assembly of a cattle chromosome

Now we have a basic overview of our data, we can generate a _de novo_ assembly.
This uses only the sequencing reads, so there currently is no bias introduced by existing references which may have errors or true bioligical differences.

```
target/release/rust-mdbg example/reads-0.00.fa.gz -k 7 --density 0.0008 -l 10 --minabund 2 --prefix example
utils/magic_simplify example
```

Now we have our genome assembly, we want to quantify some basic stats about it.
Let's go ahead and calculate the contiguity (N50).
If the assembly was perfect, we would expect 1 contig and a length of ~50 Mb.
```
calN50.js <asm.fa>
```

and we can also calculate the conserved gene completeness
```
compleasm download cetartiodactyla -L <dir>
compleasm run -a <asm.fa> -o completeness -l cetartiodactyla -L <dir> -t 4
```

We can also re-align our assembly to the reference and inspect in IGV
```
minimap2 -a -x asm5 -t 4 <reference> <asm> | samtools sort - -@ 4 --write-index -o F1_asm.bam
```
Now using the "asm5" preset rather than "map-hifi".

We can also scaffold the contigs, so they have the same direction and layout as the reference.
Note, this **is** based on the reference, and so we are at risk of losing out on biological differences as we force our assembly to look more like the reference.
For high quality assemblies, this risk should be small.

```
ragtag
```

### Construct a chromosome pangenome

Now we have our assembly, we can gather some other publically available assemblies.
```
wget other
```

To get a quick overview, we can use minigraph to build a predominantly structural variation pangenome, using the reference genome as a "backbone" for variation.
```
minigraph -cxggs -c -L 50 -j 0.01 <reference> <asm> <other> > graph.gfa
```

We can explore properties of this graph using
```
gfatools stat graph.gfa
gfatools bubble graph.gfa
```

Since this graph is not too complicated, it should be relatively easy to load and explore in bandage.
Take a look around the graph and see which regions look "tangled" or very simple.

We can also retrace the paths each assembly should take through the graph, again with minigraph
```
minigraph -cxasm --call -t16 graph.gfa sample-asm.fa > sample.bed
```

We can also confirm that "non-reference" sequence can be lost with minigraph, by checking the starting coordinate of our assemblies compared to the initial coordinates of the reference.

Now let's build a more comprehensive pangenome, using pggb.
There are several tools required to run pggb (wfmash, seqwish, smoothxg, GFAffix), and so can be challenging to install without running through singularity.
These steps are also more compute-intensive, and so here we will just use a pre-built graph generated from the command
```
pggb ...
```

We can calculate some similar stats using gfatools again, and compare to the minigraph graph generated from the **same** assemblies.
Similarly, we can load this graph in bandage and attempt to visualise the pangenome.

We can extract a subset of that graph using odgi (linux-only)
```
odgi extract -i <...>
```

and then much more managably visualise the subgraph.


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
vg deconstruct -p "ARS_UCD2.0" --path-traversals --ploidy 1 -t 2 <graph.gfa>
```

Currently, deconstructing tends to output too many variants, especially if the graph is not well formed.
For the moment, let's only consider structural variants, so we want to filter out small variants.
```
bcftools view -i 'abs(ILEN)>=50' -o graph.SVs.vcf graph.vcf
```
### Comparing to liner reference variation

We can compare against a more conventional approach, which we will take as the "truth" (although it has its own weaknesses).

```
sniffles --input sample1.bam --snf sample1.snf
```

And then merge and sort across the samples
```
sniffles --input sample1.snf sample2.snf ... sampleN.snf --vcf multisample.vcf
```

Now we can check how many of the SVs we found through the graph and through the linear-reference approach are "the same".
We can also play around with parameters to determine how strict we want to be when discussing "the same" SV.

```
jasmine --comma_filelist file_list=graph.SVs.vcf,sniffles.SVs.vcf threads=1 out_file={output}  \
        genome_file={input.reference} --pre_normalize --ignore_strand --allow_intrasample --normalize_type \
        {params.settings}
```

We can then calculate the concordance with
```
grep -oP "SUPP_VEC=\K\d+" <jasmine_output>
```

Which will tell us how many SVs are common to both sets and how many are private to each.
The approaches are extremely different (as well as technical properties like alignment length/quality), but theoretically we want as high an agreement and as few privates as possible.
Note, this method only compares the REF/ALT status of a variant, without telling us anything about how correct the per-sample genotyping was. This is a much more complex problem addressed through tools like [truvari](https://github.com/ACEnglish/truvari).

We can do something similar for the small variants, again filtering the graph output
```
bcftools view -v snps --write-index -o graph.SNPs.vcf.gz graph.vcf
```

and comparing against a pre-made VCF of DeepVariant SNP calls using bcftools
```
bcftools isec -n +1 graph.SNPs.vcf.gz DV.SNPs.vcf.gz | awk '{++c[$5]} END {for (k in c) {print k,c[k]}}'
```
And again we expect as high overlap as possible.


# Practical: Using a pangenome to Identify a known functional variant.
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
wget annotation.gff
```

### Intersecting pangenome bubbles and genes

We will then test for intersections between the annotation and pangenome "bubbles" (indicating some larger level of variation) from the minigraph gfa.
```
bedtools intersect -a annotation.gff -b minigraph.bubbles.bed > 
```

We could also approach this from the SV VCF we created, and check for potential events that are private to a subset of assemblies.
We'll start by subsetting the VCF using `bcftools view` with the `-s <samples>` option within a [process substitution](https://en.wikipedia.org/wiki/Process_substitution) ("<(...)") to avoid creating intermediate files on disk.

```
bedtools intersect -a annotation.gff -b <(bcftools view -s <asm1,asm2> graph.SVs.vcf) > specific_genes.bed
```

We can also try and convert the gff (via BED) into the pangenome coordindates, and then load that into bandage so we can visualise the events better.
```
awk gff > bed
https://github.com/AnimalGenomicsETH/KITKAT/blob/main/scripts/translate_bed_to_graph.py gfa bed
```

<_image on how to load BED in bandage>...

### Pangenome association testing

We can also look for potential variation that could be associated with a binary phenotype.
This is a more "experimental" approach, but much more flexible in what we can identify.

Let's start by purely analysing what is already present in the graph.
We want to find regions where the graph clearly contains two (or more) "haplotypes", which could be associated with a phenotype.

We can split the graph up into small chunks (say 1 Kb) using `odgi extract`, and the find the [Jaccard similarity score](https://en.wikipedia.org/wiki/Jaccard_index) between each pair of assemblies within this chunk.

```
odgi extract -i {params.og} -b {params.bed} -t {threads} -s
for g in *.og
do
  odgi similarity -t {threads} -i $g | sed 's/:[0-9]\+-[0-9]\+//g' | awk -v S=${{g}} 'NR>1 && $2>$1 {{print S,$1,$2,$6}}'
done | pigz -p {threads} > {params._output}
```

We can then examine regions where the phenotype-grouped samples are all similar but are very different across groups.

Alternatively, we can make a presence/absence matrix, again chunking the genome into small windows.

```
 odgi pav -M -i 26.pggb.og -b <(echo -e "ARS_UCD1.2#0#26\t1000000\t2000000")
```

Now we can take this a step further, and integrate population level data from short read alignments to the pangenome.
We want to gather the node-level coverage per sample after aligning with `vg giraffe`.

```
for S in samples
do
  gafpack -g pggb.gfa -a $S.gaf -l -c | awk -v S=$S '!/#/ {print S,$1,$2}' > $S.cov
done
```

We can then test to see if there is any statistically significant variation in node coverage between assigned phenotype groups.


