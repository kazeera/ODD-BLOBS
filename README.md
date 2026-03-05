# R-ODD-BLOBS
### (One Dimensional Data – Boolean Logic Binning System)

ODD-BLOBS is a pipeline for modeling DNA replication structures using quantitative chromatin fiber data.
During DNA replication, the **replication fork** forms at the boundary between replicated and unreplicated DNA. While fork activity has traditionally been inferred using genetic, molecular, and sequencing approaches, these methods do not directly visualize fork structure along individual chromatin fibers.
ODD-BLOBS analyzes chromatin fiber intensity data to identify:

- replicated DNA regions  
- replication forks  
- unreplicated DNA  
- protein localization and co-localization along the fiber  

This enables modeling of replication structures and protein behavior along individual DNA fibers.

<img src="visual description/3_Qualitative To Quantitative.JPG?raw=true" width="600"></img>
<img src="visual description/4_ODD-BLOBS_Logic.JPG?raw=true" width="600"></img>
<img src="visual description/8_Application.JPG?raw=true" width="600"></img>

---

# Credits

Conceptualized by:
**Dr. Sarah Sabatinos and Marc Green**

Published in: Sabatinos, S. A., & Green, M. D. (2018). *A Chromatin Fiber Analysis Pipeline to Model DNA Synthesis and Structures in Fission Yeast.* In **Genome Instability** (pp. 509-526)

Implemented in **R** by:
**Kazeera Aliar and Kerenza Cheng**

---
# Methods of Running ODD-BLOBS

## Method 1 — RShiny Visualization App (Recommended)

An interactive **RShiny application** is included in this repository for running ODD-BLOBS and visualizing fiber data.

The Shiny interface allows users to:

- upload fiber intensity tables  
- map imaging channels (DNA, BrdU, proteins)  
- adjust analysis thresholds  
- visualize fibers as stacked heatmaps  
- summarize protein distribution across replication regions  
- export figures as PDF  

### Running the Shiny app


Launch the app: 
The interface will open in your browser.

---

### App Interface

#### Tab 1 — Fiber Visualization

Displays stacked intensity tracks for each channel:

* DNA control
* BrdU (replication signal)
* Protein 1
* Protein 2 (optional)

Features:

* zoomable fiber visualization
* customizable lane labels
* customizable lane colors
* channel mapping flexibility
* export figure as PDF

---

#### Tab 2 — Region Summary

Displays a bar plot summarizing protein localization across fiber regions:

* replicated DNA
* forks
* unreplicated DNA

Uses a color-blind friendly **viridis palette**.

---

### User Guide

Detailed instructions for using the Shiny app are available here:

**[ODD-BLOBS Shiny User Guide](USER_GUIDE_ODDBLOBS_SHINY.pdf)**

---

## Method 2 — Running the R Scripts Directly

The original ODD-BLOBS analysis can also be run directly using the R scripts.

This approach produces JSON output files and may be useful for automated analysis pipelines.

---

# Files

### oddblobs_.R

Located in:

```
scripts/r/
```

Main script used to:

* read fiber data and user-defined parameters
* threshold intensity arrays
* identify replication tracts
* define fork boundaries
* detect protein localization

---

### functions_.R

Contains helper functions used to process and reformat data during the analysis pipeline.

---

# Input

## Fiber data table (.txt)

Input files contain fluorescence intensity arrays measured along chromatin fibers.

Each channel represents intensity along a single fiber.

Typical columns include:

```
Channel 1
Channel 2
Channel 3
Channel 4
X (pixel)
Y (pixel)
X (microns)
Y (microns)
```

Only the **Channel columns** are used by ODD-BLOBS.

---

## Command-line arguments

The original script requires eight parameters:

1. experiment name
2. file name
3. tract threshold
4. protein 1 threshold
5. protein 2 threshold
6. fork pixels into replicated region (pR)
7. fork pixels into unreplicated region (pU)
8. smoothing parameter (close gaps ≤ X pixels)

---

# Output

## table1.json

Contains sequential regions along the fiber trace.

Each object corresponds to a region:

* forkOpen
* replicated
* forkClose
* unreplicated

Example:

```
{
  "Region": "Replicated",
  "Start": 100,
  "End": 111,
  "Protein1": "[ ]",
  "Protein2": "[101,102,103,104,110,111]"
}
```

---

## table2.json

Summarizes protein distribution across region types.

Example:

```
{
  "Size": 466,
  "Prot1Percents": 44.44,
  "Prot2Percents": 46.82,
  "_row": "Tracts"
}
```

---

# Dependencies

R ≥ 3.4.4
Required R packages:

```
jsonlite
shiny
ggplot2
dplyr
tidyr
bslib
viridis
```

---

# Recommendation

For most users, the **RShiny interface is the preferred way to run ODD-BLOBS**, as it provides:

* interactive visualization
* easier parameter tuning
* direct figure export for publication.
