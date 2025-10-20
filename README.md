# Spatial Transcriptomics Reveals Chemotherapy-Reinforced Hierarchical Resistance Architecture in Pancreatic Cancer

[![DOI](https://img.shields.io/badge/DOI-10.1016%2Fj.cell.XXXX-blue)](https://doi.org/XXXXX)
[![GEO](https://img.shields.io/badge/GEO-GSE######-orange)](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE######)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

**Authors:** Ateeq Khaliq   
**Institution:** Indiana University school of Medicine   
**Correspondence:** akhaliq@iu.edu.   

---

## 📄 Citation

If you use this code or data, please cite:

```
Khaliq et al. (2025). Spatial Transcriptomics Reveals Chemotherapy-Reinforced 
Hierarchical Resistance Architecture in Pancreatic Cancer. Cell. 
DOI: [Add after publication]
```

---

## 🔬 Overview

This repository contains the complete computational analysis pipeline for our spatial transcriptomics study of pancreatic ductal adenocarcinoma (PDAC). We analyzed **66 tissue sections** from **36 patients** (treatment-naïve and post-chemotherapy) using 10X Visium spatial transcriptomics.

### Key Findings

1. **Conserved Hierarchical Architecture**: PDAC exhibits a four-layer spatial organization:
   - Tumor Core (immune-depleted, hypoxic)
   - Tumor-Stromal Interface (myCAF-enriched barriers)
   - Stromal-Immune Regulatory Belt (T cell accumulation boundary)
   - Peripheral Zone (normal tissue, vasculature)

2. **Treatment Paradox**: Chemotherapy reinforces resistance architecture despite achieving cytoreduction:
   - Consolidates stromal barriers through myCAF and SPP1+ macrophage recruitment
   - Creates metabolic stress sanctuaries with coordinated hypoxic adaptation
   - Expands immunosuppressive fields via dispersed iCAF proliferation

3. **Therapeutic Implications**: Territory-centric strategies needed to disrupt coordinated multi-layered resistance networks rather than targeting individual pathways

---

## 📊 Dataset

### Spatial Transcriptomics (10X Visium)
- **66 tissue sections** (216,180 spots after QC)
  - 40 sections from 20 treatment-naïve patients
  - 26 sections from 16 post-neoadjuvant chemotherapy patients
  - Multiregional sampling (tumor cores, invasive margins, stromal interfaces)
- **Median per spot**: 8,542 UMIs, 3,287 genes

### Multiplexed Proteomics (COMET)
- **4 representative samples** (124,889 cells)
  - 37-antibody panel (epithelial, stromal, immune markers)
  - Single-cell spatial validation
  - Cross-platform correlation (R = 0.71, p < 0.001)

### Single-Cell RNA-seq Validation
- **Independent dataset**: GSE205013
  - 27 samples, 141,186 cells
  - Treatment comparison validation
  - Cell-intrinsic transcriptional reprogramming confirmed

### Data Availability

- **Spatial transcriptomics:** [GEO GSE######](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE######)
- **COMET proteomics:** [Zenodo DOI: 10.5281/zenodo.######](https://doi.org/10.5281/zenodo.######)
- **Integrated Seurat objects:** [Zenodo DOI: 10.5281/zenodo.######](https://doi.org/10.5281/zenodo.######)

---

## 🛠️ Installation

See [INSTALL.md](INSTALL.md) for detailed installation instructions.

### Quick Start

**System Requirements:**
- R ≥ 4.3.0
- Python ≥ 3.8
- 32GB RAM minimum (64GB recommended)
- ~100GB free disk space

**Install R packages:**
```R
# Core packages
install.packages(c("Seurat", "ggplot2", "dplyr", "tidyr", "patchwork"))

# Bioconductor packages
BiocManager::install(c("spacexr", "GSVA", "ComplexHeatmap"))

# Custom tools
devtools::install_github("ateeq-khaliq/SpatialCoherence")
```

**Install Python packages:**
```bash
pip install -r requirements.txt
```

---

## 📖 Tutorial

See [TUTORIAL.md](TUTORIAL.md) for step-by-step analysis guide with example data.

### Analysis Pipeline Overview

```
Raw Visium Data (Space Ranger)
    ↓
01_spatial_preprocessing.R ──────→ QC, normalization, filtering
    ↓
02_cell_type_deconvolution.R ────→ RCTD cell type inference (3 tiers)
    ↓
03_spatial_ecotype_identification.R → Identify 10 spatial ecotypes (SE1-SE10)
    ↓
04_STARFysh_hub_identification.R → Deep learning spatial hub analysis
    ↓
05_spatial_coherence_analysis.R ─→ Quantify territorial organization
    ↓
06_hierarchical_network.R ───────→ Define 4-layer architecture
    ↓
07_pathway_enrichment.R ─────────→ MSigDB Hallmark pathway GSEA
    ↓
08_MISTy_spatial_modeling.R ─────→ Multi-scale interaction analysis
    ↓
09_treatment_comparison.R ───────→ Naive vs treated differential analysis
    ↓
10_geospatial_hotspot.R ─────────→ Hypoxic/proliferative domain mapping
    ↓
11_CNA_analysis.R ───────────────→ CopyKAT malignancy validation
    ↓
12_figure_generation.R ──────────→ Publication-quality figures
```

---

## 📁 Repository Structure

```
PDAC-spatial-analysis/
├── README.md                          # This file
├── INSTALL.md                         # Installation guide
├── TUTORIAL.md                        # Step-by-step tutorial
├── LICENSE                            # MIT License
├── requirements.txt                   # Python dependencies
├── environment.yml                    # Conda environment
├── renv.lock                         # R package versions
│
├── scripts/                          # Main analysis scripts (12 scripts)
│   ├── 01_spatial_preprocessing.R
│   ├── 02_cell_type_deconvolution.R
│   ├── 03_spatial_ecotype_identification.R
│   ├── 04_STARFysh_hub_identification.R
│   ├── 05_spatial_coherence_analysis.R
│   ├── 06_hierarchical_network_construction.R
│   ├── 07_pathway_enrichment_analysis.R
│   ├── 08_MISTy_spatial_modeling.R
│   ├── 09_treatment_comparison.R
│   ├── 10_geospatial_hotspot_analysis.R
│   ├── 11_CNA_analysis.R
│   └── 12_figure_generation.R
│
├── comet_analysis/                   # COMET validation pipeline
│   ├── 01_pseudospot_generation.R
│   ├── 02_STalign_alignment.R
│   ├── 03_validation_analysis.R
│   └── 04_cross_platform_correlation.R
│
├── utils/                            # Utility functions
│   ├── pathology_annotation.R
│   ├── visualization_functions.R
│   ├── statistical_tests.R
│   └── helper_functions.R
│
├── config/                           # Configuration files
│   ├── file_paths.R
│   └── analysis_parameters.R
│
├── data/                             # Example data (small files only)
│   └── example_data/
│       ├── example_visium_subset.rds
│       └── example_metadata.csv
│
├── outputs/                          # Example outputs
│   └── example_outputs/
│       ├── example_ecotype_plot.pdf
│       └── example_coherence_results.csv
│
└── docs/                             # Additional documentation
    ├── analysis_workflow.md
    ├── data_format_specifications.md
    └── troubleshooting.md
```

---

## 🚀 Quick Start Example

```R
# Load required libraries
library(Seurat)
library(SpatialCoherence)
library(ggplot2)

# Load example data
seurat_obj <- readRDS("data/example_data/example_visium_subset.rds")

# Run spatial coherence analysis
coherence_results <- calculate_advanced_spatial_coherence(
  seurat_object = seurat_obj,
  ecotype_column = "spatial_ecotype",
  n_permutations = 100,
  seed = 123
)

# Visualize spatial ecotypes
SpatialDimPlot(seurat_obj, group.by = "spatial_ecotype")

# Visualize coherence scores
plot_coherence_results(coherence_results)
```

---

## 🎯 Key Computational Methods

### Spatial Ecotype Identification
- **Method**: ISCHIA clustering on RCTD-deconvolved cell type proportions
- **Optimal K determination**: Elbow method, gap statistics, Calinski-Harabasz index
- **Validation**: Intra-ecotype vs inter-ecotype correlation (p < 0.001)

### Spatial Coherence Analysis
- **Tool**: [SpatialCoherence](https://github.com/ateeq-khaliq/SpatialCoherence) v1.0.0
- **Metric**: Normalized coherence score (0-1 scale)
- **Classification threshold**: 0.47 (organized vs disorganized)
- **Statistical validation**: 100 permutations

### Hierarchical Architecture Definition
- **Network analysis**: Ecotype connectivity across organizational scales
- **Spatial context networks**: Binary through quaternary combinations
- **Layer assignment**: Based on connectivity patterns and functional roles

### Pathway Analysis
- **Database**: MSigDB Hallmark gene sets (50 pathways)
- **Method**: GSVA single-sample GSEA
- **Validation**: Distance-dependent gradient analysis (R = -0.71 to 0.68)

### Treatment Response Analysis
- **Approach**: Multi-dimensional (coherence + composition + geospatial positioning)
- **Hotspot detection**: SpottedPy scale-independent analysis (10-300 spots)
- **Validation**: Independent scRNA-seq dataset (GSE205013)

---

## 📊 Spatial Ecotypes (SE1-SE10)

| Ecotype | Description | Layer | Key Cell Types |
|---------|-------------|-------|----------------|
| SE1 | Myeloid-enriched | Regulatory Belt | SPP1+ macrophages, C1Q+ macrophages, T cells |
| SE2 | Lymph node structures | Peripheral | B cells, T cells |
| SE3 | Tumor-stromal interface | Interface | Tumor cells, myCAFs |
| SE4 | Normal pancreatic tissue | Peripheral | Normal epithelial cells |
| SE5 | Vascular/endothelial | Peripheral | Endothelial cells, pericytes |
| SE6 | myCAF-enriched stromal | Interface/Belt | myCAFs (α-SMA+, TAGLN+) |
| SE7 | Malignant tumor core | Core | Tumor epithelial cells, proliferative |
| SE8 | Lymphoid-enriched | Regulatory Belt | T cells, NK cells, anti-tumor signatures |
| SE9 | iCAF-enriched stromal | Regulatory Belt | iCAFs (IL-6+, CXCL12+) |
| SE10 | Hypoxic tumor core | Core | Tumor cells, hypoxic/glycolytic |

---

## 🔬 Key Findings Summary

### 1. Conserved Four-Layer Architecture

**Evidence:**
- Consistent spatial organization across 66 specimens
- Hierarchical network analysis reveals multi-scale organization
- Spatial coherence varies by layer (Core: 0.35±0.15, Peripheral: 0.83±0.11)

**Functional Specialization:**
- **Core**: Proliferation, hypoxia, glycolysis (immune-depleted)
- **Interface**: Matrix remodeling, angiogenesis, EMT (physical barriers)
- **Regulatory Belt**: Interferon response, immune checkpoints (immune confinement)
- **Peripheral**: Oxidative phosphorylation, normal tissue functions

### 2. Treatment-Reinforced Resistance

**Compositional Changes:**
- Tumor depletion: 52.8% → 31.5% (p = 0.003)
- Stromal expansion: → 68.5% (myCAF p = 0.046, iCAF p = 0.028)

**Spatial Reorganization:**
- Enhanced coherence in stromal ecotypes (SE1 p = 0.002, SE6 p = 0.028, SE5 p = 0.041)
- No coherence change in tumor cores (SE7 p = 0.131, SE10 p = 0.387)

**Metabolic Reprogramming:**
- Treated cores: ↑ Hypoxia, EMT, PI3K-AKT, Angiogenesis
- Naive cores: ↑ OXPHOS, G2M checkpoint, Inflammatory signaling
- Validated in independent scRNA-seq (all p < 0.001)

### 3. Multi-Layered Defense Architecture

**Proximal Defense (Interface Layer):**
- myCAF accumulation around hypoxic cores (mean distance: 89.3 μm, p < 0.01)
- SPP1+ macrophage recruitment (p < 0.01)
- Angiogenic support for metabolically stressed zones (p < 0.01)

**Distal Defense (Regulatory Belt):**
- iCAF expansion without clustering (composition p = 0.028, coherence p = 0.322)
- Long-range immunosuppression via diffusible factors (IL-6, CXCL12, LIF)
- Systematic immune exclusion from hypoxic sanctuaries (p < 0.05)

---

## 🔗 Related Repositories

- **SpatialCoherence Tool:** [github.com/ateeq-khaliq/SpatialCoherence](https://github.com/ateeq-khaliq/SpatialCoherence)
  - Quantifies territorial organization in spatial transcriptomics data
  - Calculates normalized coherence scores with permutation testing
  
- **scRNA-seq Validation Analysis:** [github.com/ateeq-khaliq/NYU_Singlecell_analysis](https://github.com/ateeq-khaliq/NYU_Singlecell_analysis)
  - Independent validation using GSE205013 dataset
  - Treatment-induced transcriptional reprogramming analysis

---

## 💻 Software Dependencies

**R Packages:**
- Seurat (4.4.0) - Spatial transcriptomics analysis
- spacexr (2.2.0) - RCTD cell type deconvolution
- GSVA (1.50.0) - Pathway enrichment analysis
- SpatialCoherence (1.0.0) - Spatial organization quantification
- CopyKAT (1.1.0) - Copy number aberration analysis
- mistyR (1.0.0) - Multi-view spatial modeling

**Python Packages:**
- scanpy (1.9.1) - Single-cell analysis
- SpottedPy (1.0.0) - Geospatial hotspot analysis
- STARFysh (1.0.0) - Deep learning spatial deconvolution
- STalign (1.0.0) - Spatial alignment

See [INSTALL.md](INSTALL.md) for complete list and installation instructions.

---

## 📈 Expected Runtime

**Per sample (on 64GB RAM, 16-core machine):**
- Preprocessing & QC: ~5-10 minutes
- Cell type deconvolution: ~15-30 minutes
- Ecotype identification: ~5-10 minutes
- Spatial coherence: ~10-15 minutes
- Pathway enrichment: ~5-10 minutes
- Complete pipeline (single sample): ~1-2 hours

**Full cohort (66 samples):**
- Sequential processing: ~60-120 hours
- Parallel processing (10 samples): ~6-12 hours

---

## 🐛 Troubleshooting

**Out of memory errors:**
- Increase RAM or use swap space
- Process samples in smaller batches
- Use `gc()` to free memory between steps

**RCTD fails:**
- Verify reference dataset format
- Check for NA values in expression matrix
- Ensure cell type annotations are present

**Spatial coherence calculation slow:**
- Reduce permutations for testing (e.g., 10 instead of 100)
- Process ecotypes in parallel

See [docs/troubleshooting.md](docs/troubleshooting.md) for detailed solutions.

---

## 📧 Contact & Support

- **Primary Contact:** Ateeq Khaliq - akhaliq@iu.edu
- **Lab Website:** [URL]
- **Report Issues:** [GitHub Issues](https://github.com/ateeq-khaliq/PDAC-spatial-analysis/issues)
- **Questions:** Use GitHub Discussions or contact directly

---

## 🙏 Acknowledgments

We thank:
- The patients who generously contributed tissue samples
- [Sequencing Core Facility] for technical support
- [Computational Core] for computational resources
- Funding agencies: [List specific grants with numbers]

---

## 📜 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

You are free to:
- Use the code for academic and commercial purposes
- Modify and distribute the code
- Include it in proprietary software

Under the condition that:
- You include the original license and copyright notice

---

## 🔄 Version History

- **v1.0.0** (2025-XX-XX): Initial release with Cell publication
  - Complete analysis pipeline for 66 PDAC samples
  - Spatial ecotype identification
  - Spatial coherence analysis
  - Treatment comparison analysis
  - COMET validation pipeline

---

## 📚 Additional Resources

- **Paper:** [Cell DOI - Add after publication]
- **Preprint:** [bioRxiv DOI - if applicable]
- **Interactive Atlas:** [URL - if created]
- **Video Tutorial:** [URL - if created]
- **Presentation Slides:** [URL - if available]

---

## 🌟 Citation

If this work has been useful for your research, please cite:

```bibtex
@article{khaliq2025spatial,
  title={Spatial Transcriptomics Reveals Chemotherapy-Reinforced Hierarchical Resistance Architecture in Pancreatic Cancer},
  author={Khaliq, Ateeq and [Co-authors]},
  journal={Cell},
  year={2025},
  doi={[Add after publication]}
}
```

---

**Last Updated:** October 17, 2025

**Maintained by:** Ateeq Khaliq ([@ateeq-khaliq](https://github.com/ateeq-khaliq))

**Repository:** https://github.com/ateeq-khaliq/PDAC-spatial-analysis
