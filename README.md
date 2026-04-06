# FlyDeconv

FlyDeconv is an R package for large-scale cell-type deconvolution analysis across multiple *Drosophila melanogaster* tissues.  
It provides a unified, standardized, and extensible framework integrating **52 deconvolution methods** for systematic analysis and benchmarking.

---

## 🔬 Overview

Cell-type deconvolution aims to infer cellular composition from bulk transcriptomic data.  
Despite rapid methodological development, applications in *Drosophila* systems remain fragmented, with limited standardization across methods and datasets.

FlyDeconv addresses this challenge by providing:

- A **unified execution framework** for 52 deconvolution methods  
- **Curated multi-tissue reference datasets** derived from Fly Cell Atlas  
- Standardized input/output formats for reproducibility  
- Integrated support for benchmarking, comparison, and downstream analysis  

---

## ✨ Key Features

- **Comprehensive method integration**
  - 52 deconvolution methods implemented in a unified pipeline
  - Support for both:
    - Reference-based (RB) methods (e.g., DWLS, SECRET)
    - Semi-reference-free / marker-based (SMF) methods (e.g., GSVA, MCPcounter)

- **Multi-tissue support**
  - Built-in resources covering 15 *Drosophila* tissues

- **Standardized workflow**
  - Unified function interface: `flydeconv()`
  - Harmonized data structures across methods

- **Benchmark-ready**
  - Automatic recording of:
    - Estimated cell proportions
    - Runtime (time)
    - Memory usage

- **Reproducibility**
  - Fully reproducible pipelines
  - Optional Docker-based execution environment

---

## 📦 Installation

### Install from GitHub

```r
# install.packages("remotes")
remotes::install_github("srr0214/FlyDeconv")