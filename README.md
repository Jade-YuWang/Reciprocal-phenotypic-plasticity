📘 Reciprocal-phenotypic-plasticity
Code Repository for
“A computational ecological genetic model of phenotypic plasticity in species interactions”

This repository contains all analysis scripts, datasets, simulation modules, and mapping functions accompanying:

Wang et al., 2025 — A computational ecological genetic model of phenotypic plasticity in species interactions

The repository is designed to:

Reproduce all analyses and figures reported in the manuscript

Provide an open-source implementation of the ecological–genetic mapping framework

Enable reuse and extension by researchers working on community ecology, microbial interactions, or quantitative genetics

Ensure full transparency through publicly available data and executable scripts

🧬 Overview

The complete analysis pipeline consists of four major modules that collectively implement reciprocal phenotypic plasticity modeling, inter-genomic genetic mapping, and genetic architecture dissection:

Estimation of phenotypic plasticity from monoculture and co-culture growth curves

Functional mapping of reciprocal genetic effects using coFunMap

Bi-dimensional SNP interaction mapping between E. coli and S. aureus

Variance decomposition of direct, indirect, and epistatic components, with simulation-based validation

Each component is organized into a dedicated folder for clarity and reproducibility.
```
Reciprocal-phenotypic-plasticity/
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

```

🔬 Pipeline Overview (Data → Results)
```
                ┌───────────────────────────┐
                │   Raw Data (CSV / SNP)    │
                │  - Monoculture growth     │
                │  - Co-culture growth      │
                │  - E- & S-SNP genotypes   │
                └──────────────┬────────────┘
                               │
                               ▼
   ┌──────────────────────────────────────────────────────┐
   │ 1. Phenotypic Processing (phenotypin_plasticity_par) │
   │  - Fit logistic growth for monoculture & co-culture  │
   │  - Compute plasticity curves y(t) = x(t) – z(t)      │
   └───────────────────────────┬──────────────────────────┘
                               │
                               ▼
     ┌──────────────────────────────────────────────────┐
     │ 2. Functional Mapping (Reciprocal_phenotypic...) │
     │  - Fit coFunMap model                            │
     │  - Detect defensive / offensive loci             │
     │  - Detect horizontal epistasis                   │
     │  - Output significant SNPs                       │
     └─────────────────────────┬────────────────────────┘
                               │  
                               ▼  
     ┌──────────────────────────────────────────────────┐
     │ 3. Bi-dimensional SNP Scan                       │
     │   (Bi-dimensional mapping of plasticity genes)   │
     │  - Evaluate E × S SNP combinations               │
     │  - Produce Venn partitions & E-S interaction map │
     └─────────────────────────┬────────────────────────┘
                               │
                               ▼
       ┌──────────────────────────────────────────────┐
       │ 4. Genetic Architecture Dissection           │
       │   (Dissecting genetic architecture)          │  
       │  - Decompose variance: direct / indirect /   │
       │    epistatic components                      │
       │  - Produce plots and summaries               │
       └───────────────────────┬──────────────────────┘
                               │
                               ▼
         ┌────────────────────────────────────────────┐
         │ 5. Simulation Validation (simulation/)     │      
         │  - Check identifiability and mapping power │
         │  - Validate model stability                │
         └────────────────────────────────────────────┘

```

📄 Citation

If you use this repository, please cite:

Wang Y., He X., Yang D., Zhao J., Jin Y., Wu R.
A computational ecological genetic model of phenotypic plasticity in species interactions. 2025.

🤝 Contact

Corresponding author: Prof. Rongling Wu
Repository maintainer: Yu Wang
