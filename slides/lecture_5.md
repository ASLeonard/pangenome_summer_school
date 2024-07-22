---
# (Day 2: 8.30am – 10.30am)
title: Genome assembly and pangenome graphs
author: Alexander Leonard
institute: ETH Zürich
date: Day 2
output:
  beamer_presentation:
    theme: Boadilla
    keep_tex: true  
---

# Caveat emptor

There are not any pangenome "curriculums". \
These are ideas *we think* are useful to help apply pangenomics to your own research.

Get involved and discuss any questions or ideas of your own!

---

# Objectives

By the end of the lecture, we should be able to:

 - assemble genomes from long reads and evaluate them
 - understand common pangenome terminology and formats
 - build pangenomes graphs from input assemblies

---

# Overview

- Long read sequencing and assembly
- Some pangenome basics
- Building a pangenome

---

# Sequencing data

Genome sequencing falls in several "generations":

 - Sanger (manual)
 - "Next generation" (high throughput)
 - Third generation (long reads)

Perhaps we are due for a fourth generation designation soon...

---

# Sequencing data

Improvements in sequencing happen at a breakneck pace.

![sequencing improvements](img/sequencing_technologies.svg){ width=50% }

---

# Genome assembly

We most likely work with chromosomes that are too large to sequence directly. \
Genome sequencing provides us with many smaller fragments.

. . .

We need to *reassemble* all the reads to reconstruct the original genome sequence.

---

# Genome assembly — Theory

Consider the simple case

```ruby
TTAGGCAA  
    GCAAGTCCCA  
         TCCCATTAA
```

. . .

The assembled sequence would be "TTAGGCAAGTCCCATTAA".\
We call this a **contig** (refering to a **contig**uous region of the genome).

---

# Genome assembly — Theory

However, one possible type of issue

```ruby
TTAGGCAA  
    GCAAGTCATCAT  
             CATCATCATCCC  
          CATCATCATCCC
```

. . .

It is already ambiguous which read (the 3rd or 4th) is better, and so we have to "guess" the genome sequence.

---

# Genome assembly — Theory

Many genomes are unfortunately full of complex repeats that even with *perfect* reads cannot be resolved.

```ruby
TTAGGCAA  
    GCAAGTCCCA  
         TCACATTAA
           ^
```

Now, we add in sequencing errors and assembly gets even harder.

---

# Genome assembly — Theory

Various theoretical arguments for how much sequencing coverage is needed to "guarantee" we can reassemble the genome correctly.

. . .

In short, long+accurate+high coverage is the best.

---

# Genome assembly — Short reads

Confidently placing 150 (or smaller!) basepair sequencing reads was incredibly challenging. \
Absolutely no hope of spanning complex repeats.

. . .

Highly fragmented assemblies, but each piece is generally correct.

---

# Genome assembly — Long reads

Long reads solved some of the issues, but initially were extremely error-prone (>80% accuracy).

. . .

This required different assembly algorithms (de Bruijn graphs versus Overlap-Layout-Consensus).

---

# Genome assembly — Accurate long reads

With long and accurate reads, many problems disappeared.\
Assemblers like `hifiasm` (https://github.com/chhylp123/hifiasm) enabled drastic improvements in both:

 - genome quality
 - compute resources
 - data required

. . .

Jarvis et al. **Semi-automated assembly of high-quality diploid human reference genomes**. *Nature* (2022)

---

# Genome assembly — Accurate long reads

Even more improvements with additional sequencing:

 - ONT ultralong reads
 - Hi-C
 - SBB/mod short reads

. . .

Rautiainen et al. **Telomere-to-telomere assembly of diploid chromosomes with Verkko**. *Nature Biotechnology* (2023)

---

# Genome assembly — Triobinning

Accurate long reads also enabled distinguishing genomic variation from errors **on each read**.

Many samples of interest are not haploid, we want to assemble each haplotype.

![triobinning](img/triobinning.svg){ width=50% }

---

# Genome assembly — Triobinning

For diploids, this is easier. \
We can use parental reads to assign phases *over heretozygous regions*.

. . .

For higher ploidy, this is challenging but feasible in theory.

---

# Genome assessment

How can we investigate how *good* our genome assemblies are?

. . .

Is the assembly:

 - **contiguous**?
 - **complete**?
 - **correct**?

---

# Genome assessment — Contiguity

The easiest metric has the hardest definition:

>Given a set of contigs, the N50 is defined as the sequence length of the shortest contig such that 50% of the entire assembly is contained in contigs or scaffolds equal to or larger than this value.

---

# Genome assessment — Contiguity

Basically, are the chromosomes mostly in one piece each.

![NGx](img/N50.svg){ width=40% }

---

# Genome assessment — Completeness

We can exploit evolutionary conservation!

Search the assembly for the set of **U**niversal **S**ingle-**C**opy **O**rthologs (USCOs).

. . .

Since these genes are conserved, anything that is missing contributes to "incompleteness".

---

# Genome assessment — Completeness

We can also generate *k*-mers from short read sequencing.

. . .

Any sequencing *k*-mer missing from the assembly contributes to "incompleteness".

---

# Genome assessment — Correctness

We can align the genome to the reference and call SNPs. \
The more SNPs, the more presumed errors.

. . .

**†** REFERENCE BIAS **†**

More diverged samples should have more SNPs.

---

# Genome assessment — Correctness

We again can make use of *k*-mers from short read sequencing.

. . .

Assembly *k*-mers not found in the sequencing are more confidently errors.

. . .

Assembly *k*-mers found too often or not often enough are potential false collapses/duplications.

---

# Genome assessment — Correctness

![k-mer spectrum](img/kmer_spectrum.svg){ width=60% }

Can calculate a "QV" score from these.

---

# Genome assessment

There are many limitations to any "summary" metrics:

. . .

:::incremental
 - Many of these are global measures, but many misassemblies are local.
 - How do we pick the "best" assembly if the metrics are not strictly better in one option?
 - They were developed based on an outdated "state of the art".
:::

---

# Genome assessment

Some of these metrics are quickly becoming pointless with routine, high-quality assemblies.

. . .

How will we assess genomes in the near future? \
Will we even *need* to assess them?

---

# A quick break

And then we move into pangenomes!

---

# Pangenome complications

Comparing assemblies

TODO: coordinates not the same?


---


# Pangenomes


TODO: purpose of reference vs analysis

We now have a lot of genome assemblies, what are we going to do?

We broadly want to:

 - collapse similar regions into one sequence
 - represent variation as related regions

---

# Pangenomes

How "compact"/"strict" we want can vary.

Consider a region with a SNP every other base, should that be one big bubble or lots of small ones?

---

# Pangenomes

Ideally some happy intermediate between nucleotide-level and redundant sequence.

![extreme pangenome types](img/extreme_graph_types.svg){ width=90% }

---

# Pangenomes

Something like

![extreme pangenome types](img/ideal_graph_types.svg){ width=90% }

---

# Pangenome terminology

**Pangenome**: a collection of assemblies

**Graph**: a type of pangenome representation with nodes and edges

**Nodes**: some sequence

**Edges**: connection between contiguous sequences

**Bubbles**: regions of variation

---

# Genome file formats

Most sequencing data (or anything representing genomes) are in *fasta*/*q*.

Sequence alignments are generally in *SAM*/*BAM*.

Other miscellaneous files like *BED*, *GFF*, etc.

---

# Pangenome file formats

New file formats!

**GFA**: Graphical Fragment assembly \
Three main components:

 - S-lines: the sequence of the nodes
 - L-lines: how the graph is connected with edges
 - P-lines: how a "sample" traverses the graph (*optional*)

```ruby
S s1  AATTTACC
S s2  GGTAT
S s3  T
L s1  + s2  + 0M
L s1  + s3  + 0M
L s2  + s3  - 0M
P ME  s1+,s2+,s3+x
P YOU s1+,s3+
```

There are other, less used lines (**W**alk/**J**ump).

---

# Pangenome file formats

Most downstream tools have their own "efficient" representation:

 - `.og`
 - `.vg`
 - `.xg`
 - `.gbz`

These graphs contain a lot of information. \
GFA is human-readable and can be stored better for computer operations.

---

# Pangenome file formats

**GAF**: Graph Alignment Format \
A graph "superset" of **PAF** (Pairwise mApping Format).

Similar to SAM/BAM, broadly capturing:

 - which read
 - aligns to where
 - and how good it was

 Likewise, this is human-readable, and so some tools prefer the binary version `.gam`.

---

# Pangenomes

TODO: we focus on graph/sequence pangenoems
Overview on
 - pangenome types

 TODO: Create pangenome from VCF (but all the biases remain)

---

# Building pangenomes — tools

Several types:

 - `minigraph` (https://github.com/lh3/minigraph)
 - `cactus` (https://github.com/ComparativeGenomicsToolkit/cactus)
 - `pggb` (https://github.com/pangenome/pggb)
 - `pgr-tk` (https://github.com/cschin/pgr-tk)

Some specialised types:

 - `pangene` (https://github.com/lh3/pangene)

---

# Building pangenomes — tools

|             | Structural variation | Small variation | Reference-based | N+1      | Compute needed |
|-------------|:--------------------:|:---------------:|:---------------:|:--------:|:--------------:|
| `minigraph` | Yes                  | No              | Yes             | Easy     | Laptop         |
| `cactus`    | Yes                  | Yes             | No-ish          | Easy-ish | Cluster        |
| `pggb`      | Yes                  | Yes             | No              | Rebuild  | Big cluster    |

---

 # Building pangenomes — tools

`pggb` is lossless, so is `cactus`?

Can perfectly reconstruct any assembly from the graph, this is not true for minigraph.

---

# Building pangenomes — steps

Some form of alignment

Some form of collapsing homology

TODO: pggb params?

---

# Building bigger pangenomes

Many large scale efforts in progress:

 - Human Pangenome Research Council
 - Vertebrate Genome Project
 - Darwin Tree of Life

Several hundreds of genomes. \
Millions of core hours!

---

# Building bigger pangenomes

The "N+1" problem is a **big** problem for consortia:

 - rebuild with annual data freezes?
 - "stable" graph and "development" graph?

---

# Building bigger pangenomes

How do these problems scale for compute resources? \
What will be bottlenecks in the near future?

. . .

 - `wfmash` is $\mathcal{O}(n^2)$
 - `seqwish` is memory/disk hungry
 - `GFAffix` is almost IO bound

. . .

*Will rate of development beat rate of sequencing?*

---

# Pangenomes

Overview on
 - pangenome types
 - pros/cons
 - different pipelines for "variation graphs"

---

# Assembly pangenomes

Some genome assemblers start by building graphs representing variation in the sequencing reads.

. . .

*Isn't that a type of pangenome?*

. . .

Could we represent `n=2+` genomes as graphs?

Could we represent somatic mutations as graphs?

---

# Assembly pangenomes

The answer is: *maybe?*

"Graph-to-graph" alignment is a harder problem, but an intriguing idea to consider.

. . .

We will almost certainly never work with **less** variation, so what will that future look like?

---

# Summary

. . .

:::incremental
 - Accurate long reads have *effectively* "solved" genome assembly
   - Some complex organisms (plants especially) still have challenges
 - Graph pangenomes are great for representing all this newly accessible variation
 - Different pangenomes have different (dis)advantages
:::

---

# Questions?

And then coffee
