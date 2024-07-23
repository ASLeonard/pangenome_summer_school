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

mamba env create --file environment.yaml
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

curl https://polybox.ethz.ch/index.php/s/C7fEDiQoDweJKFU/download > data/OxO.HiFi.25.fq.gz

#Download some pre-aligned short and ONT reads for the same animal and generate the csi index for it
curl https://polybox.ethz.ch/index.php/s/OrFlXd0G54Ntx6n/download > data/OxO.Illumina.25.bam
samtools index -@ 2 -c data/OxO.Illumina.25.bam

curl https://polybox.ethz.ch/index.php/s/saR530eBS8YimRp/download > data/OxO.ONT.25.bam
samtools index -@ 2 -c data/OxO.ONT.25.bam
```

Now we want to align the reads to the reference with minimap2
```
minimap2 -a -x map-hifi -t 4 data/ARS-UCD1.2.fa.gz data/OxO.HiFi.25.fq.gz | samtools sort - -@ 4 --write-index -T $TMPDIR -o data/OxO.HiFi.25.bam
```

Where we are piping the SAM output of minimap2 into samtools to directly create our BAM alignment file.
For minimap2, the parameters are `-a` to output SAM with base-level alignment, `-x map-hifi` to use mapping parameters appropriate for HiFi reads, and `-t 4` to use 4 threads.
For samtools, we are coordinate-sorting the output with 4 threads (`-@ 4`) and creating the ".csi" index file on the fly (`--write-index`).

We can download the file if you have any issues from
```
https://polybox.ethz.ch/index.php/s/tva6Jt7hbooIjg8`/download data/OxO.HiFi.25.bam
```

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
gunzip -dc data/OxO.HiFi.25.fq.gz > data/OxO.HiFi.25.fq
rust-mdbg -k 31 -l 24 --density 0.003 --bf --minabund 2 --prefix assemblies/OxO.HiFi.25.mdbg --threads 4 data/OxO.HiFi.25.fq

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
hifiasm -t 4 -o assemblies/OxO.HiFi.25.asm_hifiasm --primary data/OxO.HiFi.fq.gz
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
>GCTACATCCACCTCATGTCCCCCTTCATCGAGCACATCCTGTGGCACTTGGGTCTGTCGGCCTGCCTGGGCCTAACGGTCGCCCTGTCCATCCTGTCGGACATCATCTCCCTTCTCACCTTCCACATCTACTGCTTCTACGTCTACGGCGCCAGGTGGGTGTGCCGCCCCCCCCCCTACCCCCCGGCGGGCTGGCGCGTGCAGCCCCGGGCCCTGGGCGCAGACAGGCTCCCGGGCCGGGCGGGTGTGGGGTGGCCCCCACCCTGACGGGAGTGGTCTGCAGGCTCTACTCCCTGAAGATCCACGGCCTGTCCTCGCTGTGGCGCCTGTTCCGAGGGAAGAAGTGGAACGTCCTGCGCCAGCGCGTGGACTCCTGCTCCTACGACCTGGACCAGGTACCAGCTCCGTCCAGCCGGCCTTCGGCCTGCGCCGGCCCATGCAGTGGCCCGGGCATCACAAGGGGCAGCCTCGCTCCTGCTCCGCAGGCCCCAGGGGTTCAGGGAGCTGCAAGGGGCAGGGTGTATACGGGGCCATAACGGGGCAGAACTCGCCGCAGGCGGGGACAGACCCTGCTCCGGGAGGCTGCAGGACAGCAGGCGCTGGAGGTGAGTGCGCATGGGGCCTGGCTCGTCCTGTTAGCTTGTCCACACCCGCCACAGGGACCCTGAGTGCTGCCCTTCAGGGTGTTCTCTTCCTGCCTAATTTTTACCCCCCTCCTGTCTAGTGTATGGGCTGCACACCTCTCAGCTGGTCCTGTCTCTACCCACTCTCTGCCCTTGTGTCACCAGGTCTGGAGGGCAGGCATCCCTCCCAGGACTGGGGAGGGCTGTCCCAACTGGCATCAGGCAGCGGTGGGCGCTCCACCTCTGGACTGGCGGAGGGGCTATGGTCAGTGCAGGTGGGGGCGCCCCAGCGTGGGGCTCTTGCTGCTGCTGATGCTGGCATCCTTTCTCCTCGTAGCTGTTCATCGGGACCTTGCTCTTCACCATCCTGATTTTCCTGCTTCCTACCACGGCCCTGTACTACCTGGTGTTCACCCTGGTGAGCCGAGTGCCACGTGCGCAGAGTGGTGACTGCCGGCTCTGCCTACCAGCGGGCAGCTGTGGGCAGACCACCCACCTGTGCCAGCCACAACAGAGGGTGGCTGTTCTGGACTCTGTGTACTGGGCCCTGGAGTGGACAGGGGTCCACCCTGGGGGTACCACCCCCTACTGCCCCGGCATACGCCGCTGGATGCATGGAGGGCAGGGCTCTCCGTAGTGCGGGGGGGCGGAGTGTGGCCAGCAGAGGGGCGTGGTACCCAGGTGGCCAGGCCGAGGGGCCAGTGGCGTGTGTCCCTGGGTTGGAGTGCCCCGTCTTCTGAGTGGGAGGGGGTTTCGTGACTCCTGTGCTCAGGCCTCAGGTTCAGGGCTGGGCCCCTGCACGCCCCGGGGACCTGTCAGGCCGCTCAGCTCAGGCATGTGCAGCATGGCTCTGGCCACGGTCCCAGGCCCAAGACCCTGAGCACTGTCTGTGCTGCAGCTCCGGCTCCTAGTGGTCACCGTGCAGGGCCTGGTCCATCTGCTCGTGGACCTCATCAACTCGCTGCCGCTGTACTCGCTTGGCCTCCGGCTGTGCCGGCCCTACAGGCTTGCGGGTAGGTCTGCACGCTGGCCAGGGAGCACGCTTGGGGCCCTGGAAGCCGCGGTTTGTGTGGCAGGCAGGCAGGGAGCACAGGGCCCTGGTTCTGGGGCGTGAGCTTCAGGGGGAGCCCGCCTTCAGGTTCAGAGGCTGCACGCGCAGCCCTCCTGCCTCCTGAGCCTTTCCTGGGATGAGTGTGCAGTTTGCTCGGAAGCGTGCCTGCAGTGTCGCATGGGCAGTGCCTGTGTGTGGCTCAGAAGAGAGACGTGTCCCCTGTGGCTGTGAGGCTGTGTGCCTGCTGGCATCCCTCCACTGGGGGCTGGCCCCCCAGAAGACTGCTGCCTGCTGCGGCCACACTGGCCAGCCCGGCCACCCCCACTGGCCCCCCTGCAGGCAGAGGCCCCTGTCTGCATGCTGACATCAGCTGCAGGTTGGCCAGATGCCCCTGACGGGAGGCCGTCCCCACTGCAGGCACAGGGAGGCCAGACCCGACTGGTCAGAGTTGTGGTCCAGCCTACCCAGTGCAGGGTCTTTGAACAGTCGGGGAAACTGGGGGTTCCTGGCGTGGCCCCCGGGTGACCAGTGTGCCCAGGAGGCCCTGGGCTGGTTTGCAGCCAGCAGCTCAGCCCACCCCAGAGCTTGTGGCGGCTGCTCGCAGCCAGTTTGGAGTCAGGAGACATGAATGTGCGTGTTTGTGGGCTTTAGGCCTTGAGTCTGCCCACCGGATGTCCTCAAACCCTCTTTCCGAGAGCTTCCTCCGCTGACCCCCAGGTGGCCCAGGCTCCCCCAGTCCCCACCACACCTCAGCCCCCGGACCTGGCACCTCCCTGGCCCCGCCTGCGCTCAGGCCTGCGCAGTCCTTTTAACTTCCCAGGAAAGGAAAGGCTGCAGATGCCGCCTTTTCTCTCCACTGCGTTGCCTTGGGGCCTCTTCATCAATTGATAATTAGCTCCTATTTCTCCTTTTAAGAACCTTCCACCTTGGGGATCCCCAGGCATGGGAGCCTACACTGCAGGCTAGAGTGTGGTGGGTCCCCGAGTCAGCTGTGCTGGAAGGCACTGTGATGAGCAGGGTGCCCACCCCCCAGTGGGTGCTGTGTCCACTCTAGGAAACCCCTGCCCTGCTCACAGAGCAGCACAGCCTGCTGGCTCTTCCCTGCCTGCAGACCCTGGGCTTTGGTCTCAGTGCGGTTCTGCTGGAGCTCTGGTGAGGGGCGGCTGCTGCCCACAACGTGCAGCCCTCCAAGGGGCGCACTGCGGGGACCGCCCCGCCTCACCGGGCCTCACTCTCCTCACAGCCGGCGTGAAGTTCCGGGTCCTGGAGCACGAGGCTGGCCGGCCCCTGCGCCTCCTGATGCAG

### Constructing a baselevel pangenome

Now we can also build a more comprehensive pangenome, using `pggb`.
There are several tools required to run `pggb` (`wfmash`, `seqwish`, `smoothxg`, `odgi`, `GFAffix`), and so can be challenging to install without running through singularity.
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

We can extract a subset of that graph using odgi (linux-only).

```
##LINUX ONLY
odgi extract -i pangenome/25.pggb.og -t 2 -r "HER:30000000-31000000" -o /dev/stdout |\
odgi view -t 2 -i - -g > pangenome/25.subgraph.gfa
```

and then much more manageably visualise the subgraph for the region we are interested in.
