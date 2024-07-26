# Practical: Using a pangenome to identify a known functional variant.
Day 5: 2:00pm – 6:00pm

## Objectives
 - Identify annotated genes overlapping/near pangenomic bubbles
 - Identify pangenomic regions associated with binary phenotypes

## Data and tools

Input
 - Pangenome
 - Reference genome annotation
 - Short read GAF alignments

Tools
 - [odgi](https://github.com/pangenome/odgi)
 - [bedtools](https://github.com/arq5x/bedtools2)
 - [gafpack](https://github.com/ekg/gafpack)

## Tasks

To start, we want to get the ARS-UCD1.2 annotation in gff format, giving us the location of genes.

```
curl https://polybox.ethz.ch/index.php/s/gnReyfSopENjpxP/download > data/ARS-UCD1.2.exons.bed
```

### Intersecting pangenome bubbles and genes

We will then test for intersections between the annotation and pangenome "bubbles" (indicating some larger level of variation) from the minigraph gfa.

```
gfatools bubble pangenome/bovines_with_P_lines.gfa | cut -f -3,12 | sed 's/HER#0#//' > pangenome/minigraph.bubbles.bed
bedtools intersect -a data/ARS-UCD1.2.exons.bed -b pangenome/minigraph.bubbles.bed > pangenome/minigraph_SV_overlaps.bed
```

We can also require more overlap with `bedtools intersect` using e.g., `-f 0.5` requiring that at least 50% of the exon is overlapped by an SV to be reported (or `-F 0.5` that 50% of the SV is overlapping an exon).

We could also approach this from the SV VCF we created (from Day 3 afternoon), and check for potential events that are private to a subset of assemblies.
We'll start by subsetting the VCF using `bcftools view` with the `-s <sample names>` option within a [process substitution](https://en.wikipedia.org/wiki/Process_substitution) ("<(...)") to avoid creating intermediate files on disk.

For SVs private to Wisent

```
bedtools intersect -a data/ARS-UCD1.2.exons.bed -b <(bcftools view -s WIS -x pangenome/minigraph.SV.vcf) > Wisent_specific_genes.bed
```

or for ones segregating in Swiss breeds

```
bedtools intersect -a data/ARS-UCD1.2.exons.bed -b <(bcftools view -s OBV,SIM,BSW -x pangenome/minigraph.SV.vcf) > Swiss_specific_genes.bed
```

We can also try and convert the annotation (via BED) into the pangenome coordinates, and then load that into bandage so we can visualise the events better (requires BandageNG).

```
curl https://raw.githubusercontent.com/AnimalGenomicsETH/pangenome_KIT/main/scripts/translate_bed_to_graph.py > tools/translate_bed_to_graph.py
python tools/translate_bed_to_graph.py pangenome/bovines.gfa data/ARS-UCD1.2.exons.bed > pangenome/ARS-UCD1.2.exons.graph.bed
```

---

Let's see if we can find any interesting gene-SV regions in Bandage!
We could search again for complex regions as we did before and inspect those regions for genes.

```
gfatools bubble pangenome/bovines_with_P_lines.gfa | sort -k4,4nr | head
```

---

### Pangenome association testing

We can also look for potential variation that could be associated with a binary phenotype.
This is a more "experimental" approach, but much more flexible in what we can identify.

Let's start by purely analysing what is already present in the graph.
We want to find regions where the graph clearly contains two (or more) "haplotypes", which could be associated with a phenotype.

We want to install a python module called regex with pip into our conda environment.

```
pip install regex
```

We should then be able to run something like the following.
For whatever reason, we need to copy and paste this in 3 steps.

```python
import regex
from collections import defaultdict

paths = {line.split()[1]:{int(N) for N in regex.split(r'\D',line.split()[2])[:-1] if N} for line in open('pangenome/bovines_with_P_lines.gfa') if line[0] == 'P'}

nodes_as_keys = defaultdict(list)
for K,V in paths.items():
    for node in V:
        nodes_as_keys[node].append(K)

#COPY AND PASTE ABOVE FIRST

nodes_private_to_pairings = defaultdict(list)
for K,V in nodes_as_keys.items():
    nodes_private_to_pairings[tuple(sorted(V))].append(K)

#COPY AND PASTE ABOVE FIRST

#we can see all combinations of samples that span nodes
print(nodes_private_to_pairings.keys())

#this will find us which nodes are only spanned by BSW, OBV, and SIM
print(nodes_private_to_pairings[('BSW','OBV','SIM')])
```

We can see there are quite a few nodes that are private to those three assemblies.
We can investigate those nodes in Bandage and see if the regions look interesting or if they are "by chance".

### Better pangenome testing

You can download the pggb graph here if you didn't make your own.

```
curl https://polybox.ethz.ch/index.php/s/dOPVst0FFGlrixa/download > pangenome/25.pggb.gfa.gz
unzip pangenome/25.pggb.gfa.gz
```

We can split the graph up into small chunks (say 1000 Kb) using `odgi extract`, and the find the [Jaccard similarity score](https://en.wikipedia.org/wiki/Jaccard_index) between each pair of assemblies within this chunk.

```
bedtools makewindows -g <(awk '$1~/HER/' pangenome/25.fa.gz.fai) -w 1000000 > pangenome/1000Kb.windows.bed

mkdir pangenome_subgraphs
cd pangenome_subgraphs

odgi extract -i ../pangenome/25.pggb.gfa -b ../pangenome/1000Kb.windows.bed -t 4 -s
for g in *.og
do
  odgi similarity -t 4 -i $g | sed 's/:[0-9]\+-[0-9]\+//g' | awk -v S=$g 'NR>1 && $2>$1 {print S,$1,$2,$6}'
done > ../pangenome/jaccard.100Kb.csv
```

We can then examine regions where the phenotype-grouped samples are all similar but are very different across groups.
Making this into a simple ratio to find peaks requires a bit more detailed work shown [here](https://github.com/AnimalGenomicsETH/pangenome_KIT/blob/main/snakepit/pangenome_analysis.smk#L174).

Alternatively, we can make an averaged presence/absence matrix, again chunking the genome into small windows.

```
odgi pav -M -i pangenome/25.pggb.gfa -b pangenome/1000Kb.windows.bed
```

The output also effectively shows how similar each sample is to the others by measuring the "Presence/absence variants ratio".
The bins are very large (and so not so sensitive), but we can already see the most diverged assembly (WIS) tends to have a lower ratio, supporting that the pangenome reflects that sample is the most different.

### Short read coverage analysis

Now we could take this a step further, and integrate population level data from short read alignments to the pangenome.
We want to gather the node-level coverage per sample after aligning with `vg giraffe`.

We have to install rust (a programming language) and compile `gafpack`.
If you are using the docker image, this should already be available

```
mamba create -n rust -c conda-forge rust

git clone https://github.com/ekg/gafpack.git
cd gafpack
cargo build --release
cd ..

# if you installed gafpack this way, be sure to use the absolute or relative path rather than just calling gafpack
# from this directory it should be "gafpack/target/release/gafpack", so we can check it is there with

ls $PWD/gafpack/target/release/gafpack
```


We can then calculate the coverage of short read sequencing over the pangenome nodes.

```
# if you aligned to the pggb graph or renamed the file, be sure to replace these values

gafpack -g pangenome/bovines_with_P_lines.gfa -a pangenome/sample.gaf -l -c | awk -v S=$S '!/#/ {print S,$1,$2}' > pangenome/summed_coverage.txt

```

With many additional samples, we could then test to see if there is any statistically significant variation in node coverage between the assigned phenotype groups.
