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
