# ODD-BLOBS
### (One Dimensional Data - Boolean Logic Binning System)

In biology, the replication fork is a structure that is formed by DNA helicase during DNA replication between the areas of "unreplicated" and replicated DNA.
Fork initiation, structure, and progression has been inferred by genetic and molecular methods such as DNA combing, chip-ChIP, and sequencing, but there is a critical gap in knowledge as to what the fork actually looks like.
  
The localization of proteins, such as Replication Protein A (RPA) and Histone H2A, during replication can be used to model fork structures and dissociation. 
  <br /> 

### Credits
* conceptualized by Dr. Sarah Sabatinos and Marc Green [published in Sabatinos, S. A., & Green, M. D. (2018). A Chromatin Fiber Analysis Pipeline to Model DNA Synthesis and Structures in Fission Yeast. In Genome Instability (pp. 509-526)]
* translated from VBA and written in R by Kazeera Aliar

### Purpose: To model DNA replication using quantitative chromatin fiber data.
* ODD-BLOBS first defines areas of 1) DNA replication and consequently, 2) replication fork and 3) unreplicated regions.
* It then checks for protein localization and co-localization along the fiber as well as in these regions.

<img src="visual description/3_Qualitative To Quantitative.JPG?raw=true" width="600"></img>
<img src="visual description/4_ODD-BLOBS_Logic.JPG?raw=true" width="600"></img>
<img src="visual description/8_Application.JPG?raw=true" width="600"></img>

### Dependencies (R statistical environment)
* R >= 3.4.4 (Mar 2018)
* jsonlite > 1.5.9
   
## Files
#### A) oddblobs_.R 
located in scripts/r is the main script used to:
* read in fiber data and user-defined arguments (thresholds, etc.)
* threshold intensity arrays
* find location of tracts (areas of replication)
* define forks at start and end of tracts
* find location of proteins
 
#### B) functions_.R 
located in scripts/r contain functions that process and reformat data  
  
## Input
#### A) fiber data table (.txt)
- includes 3 main color channels of pixel intensities (BrdU for replication and 2 proteins)
- each channel is an array, 1 x n pixels, where n is the length of the fiber (and number of rows)
 
#### B) 8 command-line arguments:
1) experiment name
2) file name
3) tract threshold - channel 2 values higher than this value will be considered "replicated"
4) protein 1 threshold
5) protein 2 threshold
6) fork - number of pixels into replicated region (pR)
7) fork - number of pixels into unreplicated region (pU)
8) smooth it - close gaps of less than or equal to this number (pixels) in array x
 
## Output
#### 1) table1.json = positions of all regions in sequential order along the fiber trace 
               - one object {} is one region (either forkOpen, replication, forkClose, unreplicated)
               - one index each for start and end
               - array indices for locations of protein 1 and 2
               - e.g. of one object: 
                {
                    "Region": "Replicated",
                    "Start": 100,
                    "End": 111,
                    "Protein1": "[  ]", 
                    "Protein2": "[ 101,102,103,104,110,111 ]"
                },
 
#### 2) table2.json = summarizes percent of protein amount across regions (percents add up to a 100)
               - one object is one region (either fork, replicated unreplicated)
               - total size in pixels, percents of protein 1 and protein 2
               - e.g of one object:
                {
                    "Size": 466,
                    "Prot1Percents": 44.4444,
                    "Prot2Percents": 46.8208,
                    "_row": "Tracts"
                },
 
 
### Update log:
* oddblobs6, functions4 (latest, Nov 6, 2020) - added annotations, removed unnecessary code<br /> 
* oddblobs5 - differentiates between data files with protein data and those without using a flag "prot_exists"<br />

![](README_visual.pdf?raw=true)
