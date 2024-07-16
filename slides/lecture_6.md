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
