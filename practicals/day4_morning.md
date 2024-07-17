# Practical: Aligning sequences to a pangenome
Day 4: 11:00am – 12:30pm

## Objectives
We should be able to
 - A
 - B
 - C

## Data and tools

Input


Tools


## Tasks

### Setting up

Need to install GraphAligner, which seemed to work with
```
mamba create -n GraphAligner -c conda-forge -c bioconda graphaligner
```

### Aligning reads to a pangenome

We can align the HiFi data used to create the Original Braunvieh haplotype back to the pangenome graph.
```
GraphAligner -g pangenome/bovines_with_P_lines.gfa  -f data/OxO.HiFi.25.fastq.gz -a test.gaf -x vg -t 4
```

Currently, `vg giraffe` is one of the best pangenome aligners, but the long read mode is not quite ready yet.
If we are able to run `vg`, we can prepare the giraffe indices with
```
vg autoindex -t 4--workflow giraffe -g pangenome/25.pggb.gfa -r "HER"  -p pangenome/vg_index
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

There are also some strange bugs in `vg` that haven't been fully resolved in this mode, so be careful for anything that looks strange.
Check out some more discussion [here](https://github.com/vgteam/vg/issues/3996) and [here](https://github.com/vgteam/vg/issues/4249) if you are interested.

### Quantifying the alignments

There is not a great way to visualise these alignments unfortunately.

We can count the number of mismatches from the alignments

From the bam files
```
samtools view OxO.HiFi.bam | grep -oE "NM:i:\d+" | cut -d':' -f 3 | awk '{++n;c+=$1} END {print c,n}'
```

from the gaf files
```
grep -oE "NM:i:\d+" test.gaf | cut -d':' -f 3 | awk '{c+=$1;++n} END {print c,n}'
```

Some reads align in many parts. \
We can find the top 10 reads with
```
cut -f 1 test.gaf | sort | uniq -c | sort -k1,1nr | head
```

And then explore where in the graph they map, what does it look like?
