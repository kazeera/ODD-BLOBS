#Last updated: Mar 21, 2018 by KA
#This version is compatible with fiber analysis data from Excel spreadsheets
#USER INPUT required

#### 1. Set up R Workspace ####
#--------------------------
#load libraries and functions
library(readxl)
library(openxlsx)
source("~/functions.R")

#USER INPUT in your data spreadsheets, which row does the table (with headers) start?
xlstart_row <- 6

#USER INPUT parameters
#smooth_it
smooth_value <- 3# readinteger("smooth it value: ")

#fork_window
pixels_in_rep <- 5#readinteger("pixels into replicated: ")
pixels_in_unrep <- 5# readinteger("pixels into unreplicated: ")

#thresholds
tract_thres <- 200#readinteger("tract threshold: ")
prot1_thres <- 4000#readinteger("protein 1 threshold: ")
prot2_thres <- 400#readinteger("threshold: ", 400)

#USER INPUT output filename
OUT_FILE <- "outfile.xlsx"

#make user parameter table
Threshold <- c(tract_thres, prot1_thres, prot2_thres)
user_parameters <- data.frame(Threshold, Smooth=c(smooth_value), "Fork Window: Replicated"=c(pixels_in_rep," ", " "), "Unreplicated"=c(pixels_in_unrep," ", " "))
row.names(user_parameters) = c("Tract", "Protein 1", "Protein 2")

#define constants 
OFF <- 0L
ON <- 1L
FORK_OPEN <- 2L
FORK_CLOSE <- 3L

#### 2. Read in file #### 
path <- file.choose()
sheet_names <- excel_sheets(path)

#### 3. Perform Analysis #### 

#Create an Excel Workbook and get sheet names from file
wb <- createWorkbook()
#addWorksheet(wb, paste("overall analysis"))

#for each sheet, run the analysis
for (i in seq_along(sheet_names)) {
  #make a new worksheet in output file for each analysis
  addWorksheet(wb, paste(sheet_names[i]," analysis"))
  df <- read_excel(path, sheet=i, skip=xlstart_row-1)
  
  #define pixel intensity arrays from data file *must be in this order, or you can change it here
  tract_data <- as.integer(unlist(df[,2]))
  prot1_data <- as.integer(unlist(df[,3]))
  prot2_data <- as.integer(unlist(df[,4]))
  
#~~~~ FINDING TRACTS (without forks)~~~~ 
  #array 1 is thresholded array (ON = above threshold, OFF = below) 
  tract_array1 <- threshold_array(tract_data, tract_thres)
  
  #array 2 is after smooth it
  tract_array2 <- smooth_it(tract_array1, smooth_value)
  
  #make dataframe with tract locations and lengths
  tracts_table <- find_tracts(tract_array2)
  #get mean and standard deviation of tracts before adding forks
  mean_before <- ifelse(is.null(tracts_table), 0, round(mean(tracts_table$Length)))
  sd_before <- ifelse(is.null(tracts_table), 0, round(sd(tracts_table$Length)))
 
 #if no tracts are found, record that no tracts are found and skip the analysis for the sheet
  if (is.null(tracts_table)) {
    writeData(wb, i, paste("ODD-BLOBS run on", date(),"No tracts found.."), startCol = 1, startRow = 1)
    writeData(wb, i, df, startCol = 10, startRow = 1, colNames=TRUE, rowNames = TRUE)
    saveWorkbook(wb,file=OUT_FILE,overwrite = TRUE)
    next #move on to next sheet in datafile
  }
  
  #otherwise, make a table of tract information
  tracts_table1 <- tracts_table
  tract_starts <- tracts_table[, "Starts At"]
  tract_ends <- tracts_table[, "Ends At"]
  tract_lengths <- tracts_table$'Length'

#~~~~ FINDING TRACTS (with forks)~~~~ 
  #array 3 is after adding forks in array
  tract_array3 <- add_forks(tract_array2, tract_starts, tract_ends, tract_lengths, pixels_in_unrep, pixels_in_rep, length(tract_array2))
  
  #remake tracts table with fork windows added in
  tracts_table2 <- find_tracts(tract_array3)
  #get mean and standard deviation of tracts
  mean_after <- ifelse(is.null(tracts_table2), 0, round(mean(tracts_table2$Length)))
  sd_after <- ifelse(is.null(tracts_table2), 0, round(sd(tracts_table2$Length)))

  #WRITE TO OVERALL DATA (in first sheet)
  avg_len_after[i] <- mean_after
  avg_sd_after[i] <- sd_after

#~~~~ FINDING PROTEINS IN REGIONS~~~~ 
#--- table1 has region start and end indices of forkopen, replicated, forkclose and unreplicated
  table1 <- find_regions(tract_array3)
#PROTEIN1
  #threshold and smooth protein 1 channel data
  prot1_array1 <- threshold_array(prot1_data, prot1_thres)
  prot1_array2 <- smooth_it(prot1_array1, smooth_value)
  
  #find indices where prot1 is ON and add to table1
  prot1_on_indices <-which(prot1_array2 == ON)
  prot1_counts <- count_occurrences(prot1_on_indices, table1$"Start", table1$"End")
  table1 <- cbind(table1, Protein1Count=ifelse(prot1_counts>0, prot1_counts, " "))
 
#PROTEIN2

  #threshold and smooth protein 2 channel data  
  prot2_array1 <- threshold_array(prot2_data, prot2_thres)
  prot2_array2 <- smooth_it(prot2_array1, smooth_value)
  
  #find indices where prot2 is ON and add to table1
  prot2_on_indices <-which(prot2_array2 == ON)
  prot2_counts <- count_occurrences(prot2_on_indices, table1$"Start", table1$"End")
  table1 <- cbind(table1, Protein2Count=ifelse(prot2_counts>0, prot2_counts, " "))
  #uncomment line below to show the indices of protein occurences
  #table1 <- cbind(table1, Protein2Indices=prot_indices_in_region(prot2_array2, nrow(table1), table1$"Code", table1$"Start", table1$"End"))

#COLOCALIZATION
  #add colocalization data to the table
  colocalized_indices <- which (prot1_array2==ON & prot2_array2==ON)
  coloc_counts <- count_occurrences(colocalized_indices, table1$"Start", table1$"End")
  table1 <- cbind(table1, Colocalization = as.logical(coloc_counts>0))
  #uncomment the next line to see number of colocalization events
  #table1 <- cbind(table1, ColocalizationCount=ifelse(coloc_counts>0, coloc_counts, " "))

#--- table2 has information about region sizes and percents of protein in those regions
  table2 <- get_region_sizes(tract_array3, table1$Code)
  table2 <- cbind(table2, Prot1Percents=get_prot_percents(prot1_array2, tract_array3))
  table2 <- cbind(table2, Prot2Percents=get_prot_percents(prot2_array2, tract_array3))
  
  ##find colocalization of protein and replace all TRUE with ON (or 1) and FALSE with OFF (or 0)
  #colocalized_prot_array <- as.integer(as.logical(prot1_array2==ON & prot2_array2==ON))
  #table2 <- cbind(table2, Colocalization=get_prot_percents(colocalized_prot_array, tract_array3))
  
  #clean up 
  table1$Code = NULL
  
#### 4. Print tables to file #### 
  j <-1 #row index
  #1. date and time of run
  writeData(wb, i+1, c(paste("ODD-BLOBS run on", date()),"Notes: ",""), startCol = 1, startRow = j) 
  j <- j+3
  
  #2. parameters (i.e. input variables by user)
  writeData(wb, i+1, user_parameters, startCol = 1, startRow = j, colNames=TRUE, rowNames = TRUE) 
  j <- j+nrow(user_parameters)+2

  #3. tract length and individual tract info without/with forks
  writeData(wb, i+1, paste("Before adding fork windows:  avg tract length: ",mean_before, " pixels +/- ",sd_before), startCol = 1, startRow = j)
  writeData(wb, i+1, paste("After adding fork windows, avg tract length: ", mean_after, " pixels +/- ", sd_after), startCol = 7, startRow = j)
  j <- j+2
  
  writeData(wb, i+1, tracts_table1, startCol = 1, startRow = j, colNames=TRUE, rowNames = TRUE)
  writeData(wb, i+1, table2, startCol = 7, startRow = j, colNames=TRUE, rowNames = TRUE) 
  j <- j+nrow(table2) + 2
  writeData(wb, i+1, tracts_table2, startCol = 7, startRow = j, colNames=TRUE, rowNames = TRUE) 
  
  ##4. tract length and aggregate prot1,prot2,colocalization in regions (forks, replicated, unrepicated)
  # writeData(wb, i+1, paste("After adding fork windows, avg tract length: ", mean_after, " pixels +/- ", sd_after), startCol = 1, startRow = j)
  # j <- j+2

  ##5. individual tract info (with forks)
  # writeData(wb, i+1, tracts_table2, startCol = 1, startRow = j, colNames=TRUE, rowNames = TRUE) 
  # j <- j+ifelse(is.null(tracts_table2), 0, nrow(tracts_table2)) + 2

  #6. table1 - protein info in each region
  writeData(wb, i+1, table1, startCol = 13, startRow = 2, colNames=TRUE, rowNames = TRUE) 

  #save the workbook with new sheet 
  saveWorkbook(wb,file=OUT_FILE,overwrite = TRUE)
  
} #closing bracket for for loop