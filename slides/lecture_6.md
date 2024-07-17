[//]: # (Day 2: 11.00am – 12.30pm)

# Lecture 6

Pangenome Visualisation

Dr. Alexander Leonard  \
ETH Zürich

[alleonard@ethz.ch](alleonard@ethz.ch)

---

# Recap

We have learned about:
 - how accurate long reads have reshaped genome assembly
 - pangenomes are a natural step to encapsulating all the new genomes
 - There is no "one true pangenome", but each has their own strengths

---

# Objectives

By the end of the lecture, we should be able to:
 - create figures of pangenomic bubbles
 - interactively explore pangenomes with different layers of data
 - identify large-scale homology between pangenome assemblies

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

*Seeing* the data can often influence later analyses. \
Too many/few reads where we expect them, overlap of variants and complex annotations, etc.

---

# Visualising genomic data

TODO: add IGV screenshot with data

---

# Visualising genomic data

There are many other ways to visualise genomic data, such as:
 - JBrowse
 - Ribbon
 - USCS Genome Browser

---

# Visualising pangenomic data

Is there a pangenomic equivalent?

**Sadly, not really...** \
Everything is more complicated in the pangenomic world.

But it depends what are we interested in.
 - Viewing relationship between many assemblies ✓
 - Viewing alignments/annotations on pangenome graphs ✖

---

# Visualising pangenomic data

How do we visualise the *GFA* output of pangenome construction?

One of the most common tools is `Bandage` (https://github.com/asl/BandageNG). \
It has several advantages:
 - easy to install
 - quick to load small-moderate graphs
 - lots of extra functionality

---

# `Bandage` visualisations

A relatively easy example of a minigraph bubble.

TODO: add image

---

# `Bandage` visualisations

A moderate example of a pggb region.

TODO: add image

---

# `Bandage` visualisations

A complex example of a pggb tangle.

TODO: add image

---

# Using `Bandage` interactively

Beyond view graphs, we can use Bandage for:
 - searching for sequence hits (`blastn`, `minimap2`, etc.)
 - annotating paths
 - loading BED files

 TODO: expand this into multiple slides

---

# Using `Bandage` interactively

[//]: # (Some interactive examples)

---

# Beyond bubbles

Bandage is a powerful tool for working on a *local* scale. \
How can we look at pangenomes (and the relationships between assemblies) on a *global* scale?

Something like synteny Circos plots?

---

# Beyond bubbles

TODO: add circos image

---

# Pangenomic synteny

We can easily construct "multi-assembly" synteny plots, but are they helpful?

TODO: add wisent synteny plot

---

# Pangenomic synteny

Many one-to-one alignments is **not** the same as many-to-many alignments.

Very easy to misinterpret or even miss key relationships. \
But, this can be a helpful *stepping stone* to transition to pangenomic concepts.

---

# Pangenomic synteny

Viewing too much information can be just as unhelpful as viewing too little.

TODO: All vs all paf plot

---

# Pangenomic synteny

Another *critical* pangenome tool is `odgi` (https://github.com/pangenome/odgi).

> odgi is a play on the Italian word "oggi" (/ˈɔd.dʒi/), which means "today". \
>As of 2019, a standard refrain in genomics is that genome graphs will be useful in x years. \
>But, if we make them efficient and scalable, they will be useful today.

We can use `odgi viz` to get something in between 1D linear (easy) and pangenome (informative) views. \
This bins the pangenome and produces a linear, static visualisation of the graph.

---

# `odgi` visualisations

TODO: wisent example for viz?

Nodes are "ordered" left to right, but what does that mean?
How do we interpret the links (graph topology)?

---

# `odgi` visualisations

There are many additional layers of information we can use:
 - inversions
 - traversal depth
 - any `odgi inject` BED files

We can also plot a "compressed" mode, and see which regions are variable.

TODO: color by depth
TODO: expand a bit more with inject etc?

---

# `odgi` visualisations

We can also plot a "compressed" mode, and see which regions are variable.

TODO: add compressed plot

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

# Pangenome visualisations

These tools are also likely best used in combination, so we can understand the graph at different scales.

Adapted from https://www.annualreviews.org/content/journals/10.1146/annurev-genom-120219-080406
TODO: adapt this

---

# Summary

- Visualising your data is **critical**, even more so for pangenomes! TODO: probably add an example for this.
- `Bandage` is powerful for interactively exploring "local" regions.
  - Increasing graph complexity will be *impossible* to responsively display.
- `odgi` is powerful for statically visualising entire graphs.

---

# Practical: Hands on pangenome assembly of a chromosome.

Goals of this afternoon. \
**Part 1:**
 - Assemble an entire chromosome from long reads
 - Build a chromosome pangenome from six assemblies

**Part 2:**
 - Visually explore minigraph and pggb pangenomes
 - Examine sequences and annotated features in pangenomes

 ---

# Questions?

And then 🍽️

---
