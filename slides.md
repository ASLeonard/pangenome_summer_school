# Lecture 0

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

# Genome assembly — theory

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

# Recap

We learned about XX

---

# Objectives

By the end of the lecture, we should be able to:
 - Create figures of pangenomic bubbles
 - Interactively explore pangenomes with different layers of data
 - Identify large-scale homology between pangenome assemblies

---

# Visualising genomic data

What types of genomic data do we normally try and visualise?

[//]: # (Interactive question)

---


# Visualising genomic data

IGV (**Integrative Genomics Viewer**) is a useful tool for visualising different formats of genomic data:
 - read alignments
 - bed files
 - gene annotations

*Seeing* the data can often influence later analyses.

---

# Visualising genomic data

TODO: discuss more about how we can do asm/eval stuff
TODO: add IGV screenshot

---

# Visualising genomic data

There are many other ways to visualise genomic data, such as:
 - JBrowse
 - Ribbon
 - USCS Genome Browser

---

# Visualising pangenomic data

Is there a pangenomic equivalent?

**Sadly, not really...**

But it depends what are we interested in?
 - Viewing relationship between many assemblies ✓
 - Viewing alignments/annotations on pangenome graphs ✖

---

# Visualising pangenomic data

How do we visualise the *GFA* output of pangenome construction?

One of the most common tools is `Bandage` (https://github.com/asl/BandageNG). \
It has several advantages:
 - easy to install
 - quick to load

---

# Bandage visualisations

A relatively easy example of a minigraph bubble.

TODO: add image

---

# Bandage visualisations

A moderate example of a pggb region.

TODO: add image

---

# Bandage visualisations

A complex example of a pggb tangle.

TODO: add image

---

# Using Bandage interactively

Beyond making images, we can use Bandage for:
 - searching sequence
 - annotating paths
 - other stuff

 TODO: expand this into multiple slides

---

# Beyond bubbles

Bandage is a powerful tool for working on a *local* scale. \
How can we look at pangenomes (and the relationships between assemblies) on a *global* scale?


Something like synteny Circos plots?
TODO: add image

---

# Pangenomic synteny

We can easily construct "multi-assembly" synteny plots, but are they helpful?

TODO: add wisent synteny plot

# Pangenomic synteny

Another *critical* pangenome tool is `odgi` (https://github.com/pangenome/odgi).

We will use `odgi viz`.

# Odgi visualisations

TODO: wisent example for viz?

Nodes are "ordered" left to right, but what does that mean?
How do we interpret the links (graph topology)?

---

# Odgi visualisations

TODO: color by depth
TODO: expand a bit more with inject etc?

---

# Other pangenome visualisers

Visualisation of a pangenome
Want to see broad synteny (odgi viz)
Want to see specific loci (bandage)

Hint at many other tools
 - Waragraph (and earlier gfaestus; depends on odgi layout)
 - sequenceTubeMap
 - Panache
 - Panagram

Many other "Pan" puns at https://github.com/colindaven/awesome-pangenomes.

---

# Combinations

Adapted from https://www.annualreviews.org/content/journals/10.1146/annurev-genom-120219-080406

---

# Summary

TODO: add summary

---

# Practical: Hands on pangenome assembly of a chromosome.

Goals of this afternoon. \
Part 1:
 - Assemble an entire chromosome from long reads
 - Build a chromosome pangenome from six assemblies

Part 2:
 - Visually explore minigraph and pggb pangenomes
 - Examine sequences and annotated features in pangenomes

 ---

# Questions?

And then lunch



---

[//]: # (Day 3: 8.30am – 10.30am)


# Lecture 7

Assessment of the pangenome

Dr. Alexander Leonard  \
ETH Zürich

[alleonard@ethz.ch](alleonard@ethz.ch)

---

# Recap

Let's recap that lecture
TODO: add recap

---

# Objectives

We want to achieve X
TODO: add objectives

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

# Summary

TODO: add summary

---

# Questions?

And then coffee

---

[//]: # (Day 3: 11.00am – 12.30am)

# Lecture 8

Assessment of the pangenome part 2

Dr. Alexander Leonard  \
ETH Zürich

[alleonard@ethz.ch](alleonard@ethz.ch)

---
# Recap

Let's recap that lecture
TODO: add recap

---

# Objectives

We want to achieve X
TODO: add objectives

---

# Pangenome quality checking

Quality check, error identification and correction of pangenomes

---

# What makes a good pangenome?

Not sure yet

---

# Summary

TODO: add summary

---

# Questions?

And then lunch

---

[//]: # (Day 5: 8.30am – 10.15am)

# Lecture 13

Finding functional associations within livestock pangenomes

Dr. Alexander Leonard  \
ETH Zürich

[alleonard@ethz.ch](alleonard@ethz.ch)

---

# Recap

Let's recap that lecture
TODO: add recap about short read alignment?

---

# Objectives

You should be able to
 - Find SV disrupted genes in a pangenome
 - TODO: check these
 - Identify possible QTL responsible for binary traits

TODO: add objectives

---

# Functional variant hunting examples in livestock

Pangenomes are big and contain a lot of variation. \
Can we prioritise variants that are *more likely* to be meaningful?

What can we do with:
 - pangenomes
 - reference annotations
 - breed/group-specific traits
 - short read sequencing

---

# Functional variant hunting examples in livestock

Important to remember, this is just **prioritising**, not proving.

In the future, these may be the "first step" of a research project.

Some details might be simplified in the following examples.

---

# Gaur deletion of *TAS2R46*

Based on work from https://www.nature.com/articles/s41467-022-30680-2.

Pangenome containing five cattle breeds (*B. t. taurus* and *B. t. indicus*) and one wild bovine (*Bos gaurus*).

Gaur diverged from cattle ~5M years ago, so what evolved in that time? \
Lots of mutations obviously, but which ones matter...

TODO: image of gaur?

---

# Gaur deletion of *TAS2R46*


We start with the reference annotation, converting GFF to BED.

```
1	ensembl	gene	339070	350389	.	-	.	gene_id "ENSBTAG00000006648"
1	ensembl	CDS		350267	350389	.	-	0	gene_id "ENSBTAG00000006648"

```
We then check structural variants unique to gaur/cattle that totally overlap genes.
Manually curate the handful of "affected" genes.

We'll do this later in the practical session!

---

# Gaur deletion of *TAS2R46*

BTA5 has a 17 Kb deletion unique to gaur that spans
 - *TAS2R46*
 - *OR6AA1*

TODO: add gaur pangenome images

Easy to find as total overlap

---

# eVNTR for *LOC112449094* in cattle

Based on https://link.springer.com/article/10.1186/s13059-023-02969-y.

TODO: How did we find the VNTR in the graph?


---

# eVNTR for *LOC112449094* in cattle

TODO: How did we genotype it, and do the functional part


---

# Wisent deletion of *THRSP*

Based on the work in https://www.biorxiv.org/content/10.1101/2024.04.08.588592v1.

European (and American) Bison diverged ~1.5M years ago, and also have distinct habitats/adaptions.

We can be a bit more precise, and look for clade-specific SVs that overlap genes.

TODO: add upsetplot

---

# Wisent deletion of *THRSP*

We find a partial deletion of *THRSP* (complete deletion of exon 1). \
The other exon is noncoding, so this is effectively a full knockout.

Short reads support deletion (and ancestral reads).

TODO: add bandage

---

# Wisent deletion of *THRSP*

We have a limited sample size (*n*=2 for bison and *n*=5 for non-bison). \
SRA/ENA have many more short read sequencing datasets for bison/cattle.

> Short read sequencing with a pangenome is a poweful tool!!

TODO: add short read plot

---

# The good, the bad, and the unknown

So far we've mostly examined only deletions, why?
 - nonreference sequence is rarely annotated
 - coordinates for insertions are not helpful
 - duplications/inversions are also not always obvious in graphs


Can we find QTL outside of annotated elements?

---

# Functional nonreference sequence

Based on the work in https://www.pnas.org/doi/abs/10.1073/pnas.2101056118.

We can identify tens of megabases of nonreference sequence, but then what? \
Since that sequence is *nonreference*, we likely have limited knowledge about it.

Unmapped reads now can align to the reference, so we can locally assembly a nonreference transcriptome. \
These putative genes can show differential expression, indicating potential functional consequences.

TODO: add image


---

# Functional nonreference sequence

"Reference" genomes being annotated is practical, not ideal. \
We now can have many assemblies per population/breed/species. \
Annotation is still computationally expensive, but Ensembl and others are working on pangenome annotation.

As **annotated** genomes become more readily available, "nonreference" will fade away.

---

# More general approaches

Approach the same problem but from the opposite end.


Rather than looking for genes that are different amongst the assemblies, we want to look for where in the pangenome do assemblies separate into groups.

How can we query the graph for such information?

---

# More general approaches

Linear equivalents to this question are:
 - GWAS
 - signatures of selection
 - Fst
 - XP-EHH

 Many of which are possible, but have not yet been implemented efficiently for pangenomes.

---

# White headed phenotype

Based on the work in https://www.biorxiv.org/content/10.1101/2024.02.02.578587v1.

Let's use a simple and easy-to-record phenotype: head color.

This trait is generally breed-defined as well.

TODO: image of hereford, brown swiss, simmental, chianina

---

# White headed phenotype

We built a large pangenome with multiple assemblies per trait group (white- or colour-headed).

We now need to "scan" the pangenome to find potential QTL. \
But how do we scan?

Broadly, we want to find where in the graph
 - W=W
 - C=C
 - W≠C

Looking for segregating nodes (sequence) is a start.

# White headed phenotype

Define the Jaccard similarity metric.

How do we slide across the genome?

Identify peaks.

Limitations
> related to trait consistency/divergence \
Size of QTL (a SNP would get lost)

---

# White headed phenotype

The identified region is extremely interesting, but our sample size is still limited. \
Testing hundreds of assemblies is (currently) impossible (for every step), so we again turn to large short read sequencing datasets.

We can also infer which samples cover (or don't) which paths through the bubble of interest.

---

# Summary

Pangenomes are a powerful *tool* for exploring functionally relevant sequence.

More annotated genomes (and more genomes containing functionally different sequence) will further empower these approaches.

Short read sequencing is still incredibly useful **when** combined with a pangenome graph.

Many of these methods are still "experimental" and rarely straightforward.

---

# Practical: Using a pangenome to identify functional variants

Goals of this afternoon:
 - Identify annotated genes overlapping/near pangenomic bubbles
 - Identify pangenomic regions associated with binary phenotypes
 - TODO: this is the last practical, so some summary?

---

# Questions?

And then coffee
