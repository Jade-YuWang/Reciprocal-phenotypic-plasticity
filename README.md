📘 Reciprocal-phenotypic-plasticity

Code Repository for “A computational ecological genetic model of phenotypic plasticity in species interactions”

This repository contains all analysis scripts, datasets, simulation code, and mapping functions used in the article:

Wang et al., 2025 – A computational ecological genetic model of phenotypic plasticity in species interactions

The purpose of this repository is to:

Reproduce the results of the manuscript

Provide open-source implementation of the ecological–genetic mapping framework

Facilitate reuse and extension of the model by other researchers

Enable transparency through publicly available data and code


🧬 Overview

This repo includes all four modules necessary for full reproduction of the study:

Estimation of phenotypic plasticity from growth curves

Mapping of reciprocal genetic effects

Bi-dimensional inter-genomic SNP interaction analysis

Variance dissection and simulation validation


📁 Repository Structure

Reciprocal-phenotypic-plasticity/
│
├── Reciprocal_phenotypic_plasticity/
│   ├── Reciprocal_phenotypic_plasticity.R
│   ├── function.R
│   ├── E_mo.csv
│   ├── S_mo.csv
│   ├── ES_E_co.csv
│   ├── ES_S_co.csv
│   ├── E-SNP.txt
│   ├── S-SNP.txt
│
├── phenotypin_plasticity_par/
│   ├── E.coli_phenotypic_plasticity.R
│   ├── S.aureus_phenotypic_plasticity.R
│   ├── E_mo.csv
│   ├── S_mo.csv
│   ├── ES_E_co.csv
│   ├── ES_S_co.csv
│   ├── E-SNP.txt
│   ├── S-SNP.txt
│
├── Bi-dimensional mapping of plasticity genes/
│   ├── Bi-dimensional&Ven.mapping.R
│   ├── function.R
│   ├── all.e.sig.csv
│   ├── all.s.sig.csv
│   ├── loci.ee.0.01.csv
│   ├── loci.s.0.01.csv
│   ├── E-SNP.txt
│   ├── S-SNP.txt
│
├── Dissecting genetic architecture/
│   ├── Dissecting_genetic_architecture.R
│   ├── e.var.fun.RData
│   ├── all.e.sig.csv
│   ├── all.s.sig.csv
│   ├── loci.ee.0.01.csv
│   ├── loci.s.0.01.csv
│   ├── E-SNP.txt
│   ├── S-SNP.txt
│
├── simulation/
│   ├── simulation0.05.R
│   ├── simulation0.1.R
│   ├── simulationp=q.R
│   ├── simulationp!q.R
│
└── README.md


📄 Citation

If you use this repository, please cite:

Wang Y., He X., Yang D., Zhao J., Jin Y., Wu R.
A computational ecological genetic model of phenotypic plasticity in species interactions. 2025.

🤝 Contact

Corresponding author: Prof. Rongling Wu
GitHub maintainer: Yu Wang
