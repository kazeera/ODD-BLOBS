#' Last updated: Nov 06, 2020 by Kazeera Aliar
#'
#' oddblobs_.R is the main script used to:
#' * read in fiber data and user-defined arguments (thresholds, etc.)
#' * threshold intensity arrays 
#' * find location of tracts (areas of replication)
#' * define forks at start and end of tracts
#' * find location of proteins

#' Notes: 
#' tract = replicated area
#' fork = replication forks/tip windows at the end of each tract, where replicated and unreplicated areas of DNA meet

# Install necessary package if not in library
if(! "jsonlite" %in% installed.packages())
  install.packages("jsonlite")
# Load library
library(jsonlite) 

# Read in array of 8 strings passed in as command-line arguments (see ReadMe for more info)
args <- commandArgs(trailingOnly = TRUE)

# Import functions from a separate file
source(file.path(getwd(),'..','scripts','r','functions4.R'))

# # Test (Note: Select block of code and control + shift + C to uncomment or comment multiple lines in RStudio)
# source("C:/Users/Uzair/Documents/GitHub/ODD-BLOBS/scripts/r/functions3.R") # change path
# # Socket server
# sock <- make.socket("10.0.1.29", 3000)
# on.exit(close.socket(sock))

# Define codes for binary off/on, fork open and close 
OFF <- 0L # also represents unreplicated region
ON <- 1L  # also represents replicated region
FORK_OPEN <- 2L
FORK_CLOSE <- 3L

setwd(file.path(getwd(),'..','appdata','fibers'))
# setwd(file.path('C:','Users','Uzair','Documents','GitHub','Tangible-Chromatin','appdata','fibers'))
# setwd('C:/Users/Uzair/Documents/GitHub/ODD-BLOBS/appdata/') # change path

# Specify experiment name and file name
EXPTname <- args[1]
FILENAME <- args[2] 

# Read in file and convert to data frame
fp <- file(file.path(getwd(), EXPTname, FILENAME), 'r')
# df <- read.csv(fp, sep = ',', skipNul= TRUE)
df <- read.table(fp, sep = '\t',header = TRUE, fileEncoding = "UTF-16LE", skipNul = TRUE)
close(fp)

# Remove unnecessary rows (always first 5)
df <- df[-c(1:5), ] 

# Check if there is info about proteins 1 and 2
# Number of columns in file without protein channels:
num_col_without_prot <- 6L

# Initialize and change flag if the number of columns exceed columns without proteins
prot_exists <-  FALSE
if (ncol(df) > num_col_without_prot) 
  prot_exists <- TRUE

# Define columns depending on whether protein info exists
if (prot_exists) {
  tract_col <- 3L
  prot1_col <- 4L
  prot2_col <- 5L
} else {
  tract_col <- 2L
}

# Read in respective pixel intensity channels/arrays 
tract_data <- df[,tract_col]
if (prot_exists) {
  prot1_data <- df[,prot1_col]
  prot2_data <- df[,prot2_col]
}

# Define smooth it value
smooth_val <- as.integer(args[8])   

# Note: 
# First threshold an array = go along the array and assign "ON" if above threshold and "OFF" otherwise
# Next "smooth" the array (close gaps of x pixels) based on the smooth it value, x

# PARAMETER TRACT THRESHOLD CHANGED
tract_thres <- as.integer(args[3])  
# array1 is thresholded array
tract_array1 <- threshold_array(tract_data, tract_thres)
# array2 is after smooth it 
tract_array2 <- smooth_it(tract_array1, smooth_val)

# Find tracts data frame with tract locations and lengths
tracts_df <- find_tracts(tract_array2)
tract_starts <- tracts_df[, "Starts At"]
tract_ends <- tracts_df[, "Ends At"]
tract_lengths <- tracts_df$'Length'

# FORK WINDOW CHANGED
pixels_in_rep <- as.integer(args[6])     
pixels_in_unrep <- as.integer(args[7])   
# array3 is after finding forks
tract_array3 <- add_forks(tract_array2, tract_starts, tract_ends, tract_lengths, pixels_in_unrep, pixels_in_rep, length(tract_array2))

# Make a data frame with region start and end indices
table1 <- find_regions(tract_array3)

# PARAMETER PROT 1 THRESHOLD CHANGED
if (prot_exists) {
  # Command line arg for protein 1 threshold
  prot1_thres <- as.integer(args[4])      
  # Threshold and smooth intensity array
  prot1_array1 <- threshold_array(prot1_data, prot1_thres)
  prot1_array2 <- smooth_it(prot1_array1, smooth_val)
} else {
  # If there is no protein info, make a placeholder array with all OFF values
  prot1_array2 <- rep(OFF, length(tract_array1))
}
# Append protein 1 column to table1
table1 <- cbind(table1, 
                Protein1 = prot_indices_in_region(prot1_array2, nrow(table1), table1$"Start", table1$"End"))

# PARAMETER PROT 2 THRESHOLD CHANGED
# Repeat with protein 2
if (prot_exists) { 
  prot2_thres <- as.integer(args[5])
  prot2_array1 <- threshold_array(prot2_data, prot2_thres)
  prot2_array2 <- smooth_it(prot2_array1, smooth_val)
} else {
  prot2_array2 <- rep(OFF, length(tract_array1))
}
table1 <- cbind(table1, 
                Protein2 = prot_indices_in_region(prot2_array2, nrow(table1), table1$"Start", table1$"End"))

# Add co-localization indices to table1 (where protein 1 and 2 are at the same pixel/position)
colocalized_array <- prot1_array2 == ON & prot2_array2 == ON
table1 <- cbind(table1, 
                Colocalization = prot_indices_in_region(colocalized_array, nrow(table1), table1$"Start", table1$"End"))

# Find co-localization of protein and replace all TRUE with ON (or 1) and FALSE with OFF (or 0)
colocalized_prot_array <- as.integer(as.logical(prot1_array2 == ON & prot2_array2 == ON))

# This table has information about region sizes and counts and pixels of protein in those regions
table2 <- get_region_sizes(tract_array3)

# # Get counts of proteins (number of pixels)
# table2 <- cbind(table2, 
#                 Prot1Count = get_prot_counts(prot1_array2, tract_array3),
#                 Prot2Count = get_prot_counts(prot2_array2, tract_array3),
#                 ColocalizationCount = get_prot_counts(colocalized_prot_array, tract_array3))

# Get relative percent of proteins across regions
table2 <- cbind(table2, 
                Prot1Percent = get_prot_percents(prot1_array2, tract_array3),
                Prot2Percent = get_prot_percents(prot2_array2, tract_array3),
                ColocalizationPercent = get_prot_percents(colocalized_prot_array, tract_array3),
                Filename = "")

# Clean up
rm(tracts_df)
table1$Code = NULL

#  Print output to console 
out <- c(toJSON(EXPTname), toJSON(FILENAME), toJSON(table1), toJSON(table2))
print(toJSON(out))

# # Test 
# stream_out(tracts_df[1:5,], con=socketConnection("10.0.1.29",3000), pagesize=5, verbose=TRUE,prefix="")
# In output folder: 
# # Make tables into json format
# write(file='output1.json', toJSON(table1))
# write.csv(table1, file = "table1.csv")
# out1 <- toJSON(table1)
# write.socket(sock,out1)
# write(file='output2.json', toJSON(table2))
# write.csv(table2, file = "table2.csv")
# out2 <- toJSON(table2)
# write.socket(sock,out2)

# Clear R environment - remove all variables
rm(list=ls())