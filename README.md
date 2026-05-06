# GRaFT – Grain-Resolved Fabric and Texture Analysis

Grain-Resolved Fabric and Texture analysis workflow for EBSD datasets processed using the MAPClean ecosystem.

**GRaFT** is the third stage of the MAPClean ecosystem and is designed for process-agnostic quantitative texture and microstructural analysis of reconstructed EBSD grain datasets. The workflow operates on the final grain outputs produced by **GRaMC** and separates computation into two linked environments:

- **MATLAB / MTEX** for all EBSD-derived calculations, map generation, checkpointing, and CSV export.
- **Python** for downstream statistical visualisation, cross-sample comparison, and interpretative analysis.

The workflow automatically processes all phase-specific `*_finalGrains.mat` files exported from GRaMC and produces:
- publication-ready EBSD maps
- grain-resolved quantitative datasets
- contact-network datasets
- pixel-scale deformation datasets
- sample-level texture statistics

---

# Key Features

## Core Texture Analysis
- Orientation Distribution Function (ODF) calculation
- Pole Figure (PF) generation
- Inverse Pole Figure (IPF) generation
- Texture-strength quantification
- Fabric tensor analysis

## Grain-Scale Quantification
- Equivalent circular diameter (ECD)
- Aspect ratio
- Long-axis orientation
- Grain Orientation Spread (GOS)
- Ellipse fitting
- Crystal-shape visualisation

## Pixel-Scale Deformation Analysis
- Grain Reference Orientation Deviation (GROD)
- Kernel Average Misorientation (KAM)
- GROD/GOS ratio mapping

## Contact and Pair Analysis
- Grain-pair misorientation analysis
- Long-axis angular comparison
- Crystallographic axis comparison

## Multi-Phase Support
Currently configured for:
- Anorthite
- Forsterite
- Diopside

---

# Workflow Structure

## MATLAB / MTEX Stage
The MATLAB script performs:
- EBSD loading
- grain-object analysis
- texture calculations
- map generation
- CSV export
- checkpoint creation

## Python Notebook Stage
The companion Python notebook performs:
- statistical plotting
- population comparison
- distribution analysis
- cross-sample synthesis
- interpretative visualisation

The MATLAB stage exports all required CSV datasets for direct use in the Python notebook.

---

# Maps Produced

Each phase exports the following maps:

| Map ID | Description |
| :--- | :--- |
| 01 | Pole Figure (PF) |
| 02 | Inverse Pole Figure (IPF) |
| 03 | Aspect Ratio Map |
| 04 | Grain Size (ECD) Map |
| 05 | Long-Axis Orientation Map |
| 06 | GOS Map |
| 07 | GROD Map |
| 08 | GROD/GOS Ratio Map |
| 09 | Crystal Shape Map |
| 10 | KAM Map |
| 11 | Ellipse Fit Map |

---

# CSV Exports

## All Phases
- `*_grains.csv`

## Anorthite Only
- `*_contacts.csv`
- `*_pixels.csv`

## Summary
- `AllSamples_TextureStats.csv`

---

# Requirements

- MATLAB (tested on R2024b)
- MTEX Toolbox (tested on v6.0.0)

Python notebook dependencies depend on the downstream statistical workflow.

---

# Installation

1. Place `GRaFT.m` in the same workspace as:
   - `MAPClean.m`
   - `GRaMC.m`

2. Ensure the `checkpoints/` folder contains:
   ```text
   *_finalGrains.mat
   ```

3. Open MATLAB and run:
   ```matlab
   GRaFT
   ```

---

# Directory Structure

```text
My_EBSD_Project/
├── DataFiles/
├── checkpoints/
├── MAPClean/
├── GRaMC/
├── GRaFT/
├── MAPClean.m
├── GRaMC.m
├── GRaFT.m
└── GRaFT.ipynb
```

---

# Stage Control Flags

The workflow can be selectively controlled using:

```matlab
run_statsCalc        = true;
run_samplePlots      = true;
run_exportCSVs       = true;
run_contactAnalysis  = true;
run_pixelExport      = true;
```

---

# Parameters

| Parameter | Description |
| :--- | :--- |
| `halfwidth` | ODF halfwidth |
| `resolution` | Pole figure grid resolution |
| `exportRes` | PNG export resolution |
| `mapLineWidth` | Boundary linewidth |
| `minGrainsForODF` | Minimum grains required for ODF calculation |
| `saveIntermediateMAT` | Save intermediate checkpoint MAT files |

---

# The MAPClean Ecosystem

## Stage 1 — MAPClean
Pixel-level EBSD restoration and adaptive cleaning.

## Stage 2 — GRaMC
Grain reconstruction and microstructural restoration.

## Stage 3 — GRaFT
Grain-resolved fabric and texture analysis.

---

# Citation

If you use this workflow, please cite the associated publication once available.

---

# License

This project is licensed under the GPL v3 License.
