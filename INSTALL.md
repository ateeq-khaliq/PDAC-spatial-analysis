# Installation Guide

## System Requirements

### Hardware
- **CPU:** Multi-core processor (8+ cores recommended)
- **RAM:** 32GB minimum, 64GB recommended
- **Storage:** 100GB free disk space
- **OS:** Linux, macOS, or Windows (with WSL2)

### Software
- **R:** Version 4.3.0 or higher
- **Python:** Version 3.8 or higher
- **RStudio:** Recommended for interactive analysis

---

## Installation Steps

### 1. Install R and RStudio

**Linux (Ubuntu/Debian):**
```bash
sudo apt update
sudo apt install r-base r-base-dev
```

**macOS:**
```bash
brew install r
```

**Windows:**
Download from [CRAN](https://cran.r-project.org/)

**RStudio:**
Download from [RStudio website](https://posit.co/download/rstudio-desktop/)

### 2. Install R Packages

Open R or RStudio and run:
```R
# Set CRAN mirror
options(repos = c(CRAN = "https://cloud.r-project.org/"))

# Install basic packages
install.packages(c(
  "Seurat",           # v4.4.0
  "ggplot2",          # v3.4.2
  "dplyr",            # v1.1.2
  "tidyr",            # v1.3.0
  "patchwork",        # v1.1.2
  "viridis",          # v0.6.3
  "RColorBrewer",     # v1.1-3
  "pheatmap",         # v1.0.12
  "scales",           # v1.2.1
  "cowplot",          # v1.1.1
  "devtools",         # v2.4.5
  "BiocManager"       # v1.30.21
))

# Install Bioconductor packages
BiocManager::install(c(
  "GSVA",                    # v1.50.0
  "ComplexHeatmap",          # v2.16.0
  "SingleCellExperiment"     # v1.22.0
))

# Install spatial analysis packages
install.packages("remotes")
remotes::install_github("dmcable/spacexr", build_vignettes = FALSE)  # RCTD v2.2.0

# Install custom tools
devtools::install_github("ateeq-khaliq/SpatialCoherence")

# Install other GitHub packages
devtools::install_github("saeyslab/nichenetr")
devtools::install_github("navinlabcode/copykat")
BiocManager::install("mistyR")
```

### 3. Install Python and Packages

**Install Python (if not already installed):**

**Linux:**
```bash
sudo apt install python3 python3-pip python3-venv
```

**macOS:**
```bash
brew install python3
```

**Windows:**
Download from [python.org](https://www.python.org/)

**Create virtual environment:**
```bash
cd PDAC-spatial-analysis
python3 -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate
```

**Install Python packages:**
```bash
pip install --upgrade pip
pip install -r requirements.txt
```

---

## Verification

Test your installation:
```R
# Test R packages
library(Seurat)
library(SpatialCoherence)
library(GSVA)
library(spacexr)
print("R packages loaded successfully!")
```
```python
# Test Python packages
import scanpy as sc
import pandas as pd
import numpy as np
print("Python packages loaded successfully!")
```

---

## Common Installation Issues

### Issue 1: Cannot install Seurat
**Solution:**
```bash
# Install system dependencies (Ubuntu/Debian)
sudo apt install libhdf5-dev libcurl4-openssl-dev libssl-dev libxml2-dev
```

### Issue 2: RCTD/spacexr installation fails
**Solution:**
```R
# Install with specific options
remotes::install_github("dmcable/spacexr", 
                        build_vignettes = FALSE, 
                        force = TRUE)
```

### Issue 3: Python package conflicts
**Solution:**
```bash
# Create fresh environment
python3 -m venv fresh_env
source fresh_env/bin/activate
pip install -r requirements.txt
```

### Issue 4: Out of memory errors
**Solution:**
- Increase system RAM or use swap space
- Process samples in batches
- Use data.table for large files

---

## Getting Help

If you encounter issues:

1. Check [Troubleshooting Guide](docs/troubleshooting.md)
2. Search [GitHub Issues](https://github.com/ateeq-khaliq/PDAC-spatial-analysis/issues)
3. Contact: [your.email@institution.edu]

---

## Next Steps

After installation, proceed to [TUTORIAL.md](TUTORIAL.md) for analysis walkthrough.
