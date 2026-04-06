# FlyDeconv

FlyDeconv is an R package for cell-type deconvolution analysis across multiple *Drosophila melanogaster* tissues.  
It provides a unified, standardized, and extensible framework for running diverse deconvolution methods with built-in internal reference resources.

---

## 🔬 Overview

Cell-type deconvolution aims to estimate cellular composition from bulk transcriptomic data.  
While numerous deconvolution methods have been developed, their application to *Drosophila* systems remains fragmented and lacks standardized workflows.

FlyDeconv addresses this gap by:

- Integrating multiple deconvolution methods under a unified interface
- Providing curated internal reference datasets for *Drosophila*
- Standardizing input/output formats across methods
- Supporting benchmarking and comparative analysis
- Enabling seamless downstream integration (e.g., Shiny apps)

---

## ✨ Key Features

- **Unified interface**: Run multiple deconvolution methods via a single function (`flydeconv()`)
- **Multi-tissue support**: Built-in resources for 15 *Drosophila* tissues
- **Method diversity**:
  - Reference-based methods (e.g., DWLS, SECRET)
  - Marker-based methods (e.g., GSVA, MCPcounter)
- **Automatic output management**:
  - Cell proportion estimates
  - Runtime benchmarking (time and memory)
- **Reproducible framework**: Standardized data structures and execution pipelines

---

## 📦 Installation

### Install from GitHub

```r
# install.packages("remotes")
remotes::install_github("YOUR_GITHUB_USERNAME/FlyDeconv")