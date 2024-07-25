# Practical: Aligning sequences to a pangenome
Day 4: 11:00am – 12:30pm

## Objectives
We should be able to
 - align short reads to a pangenome graph
 - align long reads to a pangenome graph
 - understand the GAF alignment output

## Data and tools

Input
 - short reads
 - long reads
 - Pangenome

Tools
 - [vg](https://github.com/vgteam/vg)
 - [GraphAligner](https://github.com/maickrau/GraphAligner)

## Tasks

### Setting up

We should all have been able to work with `vg` yesterday, so we should already be able to run `vg giraffe`.

We can try to install GraphAligner, which seemed to work with

```
mamba create -n GraphAligner -c conda-forge -c bioconda graphaligner
```

GraphAligner has some different dependencies to earlier tools we've used, so we are creating a new environment (the `-n GraphAligner`) to allow a stress-free install.

### Aligning short reads to a pangenome

Currently, `vg giraffe` is one of the best pangenome aligners, but the long read mode is not quite ready yet.
If we are able to run `vg`, we can prepare the giraffe indices with

```
vg autoindex -t 4 --workflow giraffe -g pangenome/25.pggb.gfa -r "HER"  -p pangenome/vg_index
# this step is not necessary in most cases, but a strange filesystem feature that bottlenecks jobs on a typical high performance computer
chmod 444 pangenome/vg_index.dist
```

We can extract the short reads from the original bam with

```
samtools collate -Ou data/OxO.Illumina.SR.bam | samtools fastq -N -o data/OxO.Illumina.fq.gz
```

And then run the giraffe mapping for short read samples, with the output in GAF.

```
vg giraffe -t 4 -Z pangenome/vg_index.giraffe.gbz -m pangenome/vg_index.min -d pangenome/vg_index.dist --named-coordinates -o gaf --interleaved -f data/OxO.Illumina.fq.gz > pangenome/sample.gaf
```

Note as well, running without `--named-coordinates` leads to very different output, as `vg` internally tends to work with "chopped up" graphs with 32bp nodes.

Always check the nodes of your alignments match what you expect from your graph!
For example, getting alignments to node ">51251" when you expect your graph to only have 4000 nodes means something has gone wrong.

There are also some strange bugs in `vg` that haven't been fully resolved in this mode, so be careful for anything that looks strange.
Many of these bugs have been fixed through feedback, but these tools are still relatively new.
Check out some more discussion [here](https://github.com/vgteam/vg/issues/3996) and [here](https://github.com/vgteam/vg/issues/4249) if you are interested.

### Aligning long reads to a pangenome

We can align the HiFi data used to create the Original Braunvieh haplotype back to the pangenome graph.

```
# we can activate the env for GraphAligner
mamba activate GraphAligner

# run the alignment
GraphAligner -g pangenome/bovines_with_P_lines.gfa  -f data/OxO.HiFi.25.fastq.gz -a pangenome/OxO.HiFi.25.gaf -x vg -t 4
```

GraphAligner can also align to de Bruijn graphs (more k-mer based than our sequence-based graphs), so the `-x vg` indicates we want to use preset parameters relevant for a sequence graph.

### Quantifying the alignments

This is just a simple example, so the effect may not be so large, but let's try and quantify how much better the alignments are to the graph compared to the linear reference.

We can count the number of mismatches from the alignment records, typically encoded as the "NM:i:" tag.

From the bam files, we can print all the records, only keep the "NM" tag (`grep -oE "NM:i:\d+"`), only keep the actual number of mismatches (`cut -d':' -f 3`), and then sum these up (`awk '{++n;c+=$1} END {print c,n}'`).

```
# we need samtools back, so reactivate your original environment (use the correct line for your name)

mamba activate ognigenoma
#OR
mamba activate summerschool2

samtools view data/OxO.HiFi.25.bam | grep -oE "NM:i:\d+" | cut -d':' -f 3 | awk '{++n;c+=$1} END {print "number of alignments: "n"\nnumber of mismatches: ",c}'
```

And then we will do the same thing for the graph alignment file
from the gaf files

```
grep -oE "NM:i:\d+" pangenome/OxO.HiFi.25.gaf | cut -d':' -f 3 | awk '{++n;c+=$1} END {print "number of alignments: "n"\nnumber of mismatches: ",c}'
```

Some reads might align in many parts, so let's look for the top 10 reads with the most number of multiple alignments.

```
awk '{++R[$1]} END {for (r in R) {print r,R[r]}}' pangenome/OxO.HiFi.25.gaf | sort -k1,1nr | head

# or maybe a bit easier to read but should be identical output

cut -f 1 pangenome/OxO.HiFi.25.gaf | sort | uniq -c | sort -k1,1nr | head
```

There is not a great way to visualise these alignments similar to how we can view BAM alignments in IGV unfortunately.
But we can explore where they map to in the graph.

We can also try and convert the alignment format into a BED file which we can load with bandage.
This code is quite fragile and would need more robust parsing to handle reads which align across multiple nodes rather than just ignoring them, but it should work.
We first sort by the 10th column for aligned length, so we can see longer matches.

The BED format required by Bandage is also a bit ugly and not well documented after the first 3 standard columns, but these appears to work.
We should be able to load this in Bandage and see where the reads aligned!

```
sort -k10,10nr pangenome/OxO.HiFi.25.gaf | awk -v OFS='\t' '$6~/>/ {n=1gsub(">","",$6); if (n==1) {print $6,$8,$9,$1,$10,"+",$8,$9,"150,150,0";next}} {n=gsub("<","",$6); if (n==1) {print $6,$8,$9,$1,$10,"-",$8,$9,"0,150,150"}}' > pangenome/OxO.HiFi.25.GA.bed
```
