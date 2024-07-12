# Lecture 1

Introduction to pangenomics

Dr. Alexander Leonard  \
ETH Zürich

[alleonard@ethz.ch](alleonard@ethz.ch)

---

[//]: # (Day 1: 8.30am - 9.15am)

# What is a "genome"?

Something obvious

What is "my genome"

Pangenome simple introduction
What is a genome?
Reference genome
My genome(s)
What is a pangenome
Multiple assemblies
Graph structure?
Tradeoff between complexity and potential

---

# What is a "genome"?

Something obvious

What ~is "my genome"~ are my genomes


---

# What is a "pangenome"?

First used in connection to bacteria -- many E. coli genomes share ~50% of genes, so the "set" of genes was the "pangenome".

Now we refer to a collection of *many* genomes as a pangenome.

In reality there are many types of pangenomes.
We will only discuss one flavour.

---

# So why just now?

Pangenomes may seem very obvious, but

 - Assemblies *were* extremely expensive/hard to produce
 - The scale of diversity (in humans at least) was underestimated
 - Alignment/graph algorithms weren't possible

---

# Some guiding steps

The goals of this course are to
 - perform analyses with long read data
 - construct and assess pangenomes
 - use pangenomes to find meaningful variation

---

[//]: # (Day 2: 8.30am – 10.30am)

# Lecture 5

Genome assembly and pangenome graphs

Dr. Alexander Leonard  \
ETH Zürich

[alleonard@ethz.ch](alleonard@ethz.ch)

---

# Objectives

We want to cover
 - X
 - Y

---

# Genome assembly

Genome assembly x several slides
 - assembly theory
 - assembly QC

---

# Genome assembly

Sequencing has come a long way in recent years, but we still only get fragments still.
(Maybe not for long, ONT yeast project...)

We need to *reassemble* all the reads to reconstruct the original genome sequence.

---

# Genome assembly -- theory

Consider the simple case

```
TTAGGCAA  
    GCAAGTCCCA  
         TCCCATTAA
```

The assembled sequence would be "TTAGGCAAGTCCCATTAA"

---

# Genome assembly

However, one possible type of issue

```
TTAGGCAA  
    GCAAGTCATCAT  
             CATCATCATCCC  
          CATCATCATCCC
```

It is already ambiguous which read (the 3rd or 4th) is better, and so we have to "guess" the genome sequence.

---

# Genome assembly

Many genomes are unfortunately full of complex repeats that even with *perfect* reads cannot be resolved.
Now, we add in sequencing errors and assembly gets even harder.

```
TTAGGCAA  
    GCAAGTCCCA  
         TCACATTAA
           ^
```

---

# Genome assembly -- short reads

Confidently placing 150 (or smaller!) basepair sequencing reads was incredibly challenging.

Long reads solved some of the issues, but initially were extremely error-prone (>80% accuracy).

This required different assembly algorithms.

---

# Genome assembly -- long reads

HiFi assembly is easier (not solved)

Hint at T2T recipe + verkko + graphs (motivate single-sample pangenome)

---

# Genome assembly

Genome assembly x several slides
 - assembly theory
 - assembly QC

---

# Pangenomes

Overview on
 - pangenome types
 - pros/cons
 - different pipelines for "variation graphs"

---

# Pangenomes -- variation graphs

Overview on
 - pangenome types
 - pros/cons
 - different pipelines for "variation graphs"

---

# Pangenomes

Overview on
 - pangenome types
 - pros/cons
 - different pipelines for "variation graphs"

---

# Pangenomes

Overview on
 - pangenome types
 - pros/cons
 - different pipelines for "variation graphs"

---

# Pangenomes

Overview on
 - pangenome types
 - pros/cons
 - different pipelines for "variation graphs"

---

# Some more stuff

Pangenome overview
Pangenome assembly methods and graphical representation.
Description of software available for pangenome assembly, running
assembly software visualization tools.

# HPRC/VPG/dToL updates

Consider big projects, with millions of core hours

Overview on
 - pangenome types
 - pros/cons
 - different pipelines for "variation graphs"

---

# Additional thoughts

Some genome assemblers start by building graphs representing variation in the sequencing reads.

*Isn't that a type of pangenome?*

Could we represent `n=2+` genomes as graphs?

Could we represent somatic mutations as graphs?

---

# Additional thoughts

The answer is **maybe?**

"Graph-to-graph" alignment is a harder problem, but an intriguing idea to consider.

We will almost certainly never work with **less** variation, so what will that future look like?

---

# Summary

 - Accurate long reads have *effectively* "solved" genome assembly
   - Some complex organisms (plants especially) still have challenges
 - Graph pangenomes are great for representing all this newly accessible variation
 - Different pangenomes have different (dis)advantages

---

# Coffee

---

[//]: # (Day 2: 11.00am – 12.30pm)

# Lecture 6

Pangenome Visualisation

Dr. Alexander Leonard  \
ETH Zürich

---

# Pangenomes

IGV is useful for understanding alignments/bed
Is there a pangenome equivalent? (not really)

---

# Pangenomes

Multiple sequence comparisons can be confusing

---

# Pangenomes

Some nice example of a graph bubble.

But graphs can be confusing as well.

Some graph hairball.

---

# Pangenomes

Visualisation of a pangenome
Want to see broad synteny (odgi viz)
Want to see specific loci (bandage)

Hint at many other tools
Gfateaus
waragraph
pantools

---

# Practical: Hands on pangenome assembly of a chromosome.

Goals of this afternoon.
Part 1:
 - Assemble an entire chromosome from long reads
 - make a chromosome

Part 2:
 - A
 - B

---

[//]: # (Day 3: 8.30am – 10.30am)


# Lecture 7

Assessment of the pangenome

Dr. Alexander Leonard  \
ETH Zürich

[alleonard@ethz.ch](alleonard@ethz.ch)

---

# Okay now for really

Assessment of the pangenome:
Identification of assembly errors, finding rogue sequences, quality
metrics, assessment of completeness of the pangenome

Layering information into visualisation
Gene annotations
Per-node values
paths

Here we include a figure
![Figure](./bandage_issue.svg)

---

[//]: # (Day 3: 11.00am – 12.30am)

# Lecture 8

Assessment of the pangenome part 2

Dr. Alexander Leonard  \
ETH Zürich

[alleonard@ethz.ch](alleonard@ethz.ch)

---

# Pangenome quality checking

Quality check, error identification and correction of pangenomes

---

# What makes a good pangenome?

Not sure yet

---

[//]: # (Day 5: 8.30am – 10.15am)

# Lecture 13

Finding associations within the pangenome

Dr. Alexander Leonard  \
ETH Zürich

[alleonard@ethz.ch](alleonard@ethz.ch)

---

# Functional variant hunting examples in livestock

Primarily lift over IGGSy examples

---

# Gaur deletion of *TAS2R46*

Total deletion of two coding genes
 - *TAS2R46*
 - *OR6AA1*

Easy to find as total overlap

---

# Wisent deletion of *THRSP*

Partial deletion of gene (complete deletion of exon 1).
Other exon is noncoding, so effectively a full knockout.
Short reads support deletion (and ancestral reads).

---

# Some observations

So far we've examined only deletions, why?
 - Nonreference sequence is rarely annotated
 - Coordinates for insertions are not helpful
 - Duplications are also not obvious in graphs

Can we find QTL outside of annotated elements

---

# More general approaches

How to query graph

---

# White headed phenotype

Jaccard metric slides x2-3

Short read alignments (and issues)

---

# Practical: Using a pangenome to Identify a known functional variant

What we gonna do
