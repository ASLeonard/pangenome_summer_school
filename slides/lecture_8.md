---
# (Day 3: 11.00am – 12.30am)
title: Working with pangenomes
author: Alexander Leonard
institute: ETH Zürich
date: Day 3
output:
  beamer_presentation:
    theme: Boadilla
    keep_tex: true  
---

# Recap

We have learned about:
 - assessing the general graph metrics of a pangenome
 - determining the core and accessory components of a pangenome
 - TODO

 Assessment of the pangenome:
Identification of assembly errors, finding rogue sequences, quality
metrics, assessment of completeness of the pangenome
Alex Leonard
10.30am – 11.00am Coffee Break
11.00am – 12.30am Quality check, error identification and correction of pangenomes


---

# Objectives

By the end of the lecture, we should be able to:
 - identify poorly constructed regions
 - edit GFA files to resolve pangenome errors
 - how to rerun?

  (https://timkahlke.github.io/LongRead_tutorials/BAN.html)

---

# Overview

something

 - A
 - B

 Trimming out low quality nodes etc (vg clip, odgi depth?)

---

# Pangenome quality checking

Quality check, error identification and correction of pangenomes



---

# Pangenome correction

"Bad idea", but we can manually edit graphs
(some semblance to Hi-C scaffolding)

---

# A quick break

And then some aspirational methods!

---

# Bonus topics

# Black box pangenomes

We know pangenomes improve analyses, but some fields are hesistant to change. \
So what can we do?

. . .

How much "the algorithm" do we understand for accepted programs? \
Some degree of input → black box → output...

---

# Black box pangenomes

`vg surject` pipeline can:

 - align reads to a pangenome
 - "surject" reads back to **a** linear reference
 - variant call with `DeepVariant` from a typical BAM file

. . .

Lose some of the benefit by discarding non-reference alignments, but net positive.

---

# Personalised pangenomes

Imagine a population graph containing **all** variation.

. . .

Currently, too much variation:

 - is prohibitively slow to work with
 - can be detrimental to accuracy (too many combinations)

---

# Personalised pangenomes

---

# TODO

Where should these go?
 - `vg surject`
 - personalised pangenomes
 -  "pangenome communities"

---

# Summary

. . .

:::incremental:::
 - A
 - B
 - C
:::

---

# Questions?

And then lunch
