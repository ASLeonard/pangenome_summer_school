
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
