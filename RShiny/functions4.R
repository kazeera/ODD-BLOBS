#' Functions for oddblobs_.R
#' Last updated: Nov 06, 2020 by Kazeera Aliar

#' Assigns 1 if intensity is greater or equal to threshold and 0 if it is less
#' @param x A numeric vector.
#' @param threshold A number.
#' @return A vector of 0s (below threshold) and 1s (above).
#' 
threshold_array <- function(x, threshold) {
   #  When x is less than the threshold, assign 0, else 0
   pmax(x>= threshold, 0)
}


#' Closes gaps of less than or equal to smooth_val (pixels) in array x
#' @param x A numeric vector.
#' @param smooth_val A number.
#' @return A numeric vector where gaps of 1 to smooth_val have been "filled in".
#' 
smooth_it <- function(x, smooth_val) {
  # Converts NA to zero
  x[is.na(x)] <- 0
  
  # Return original array if the smooth value is zero to reduce computation
  if (smooth_val == 0)  { return(x) }
  
  # Finds runs of OFF that are less than equal to the value of smooth it
  rle(x == OFF) -> y
  runs <- which(y$value==TRUE & y$lengths <= smooth_val)
  
  # Find the starts and ends of these indices
  cumsum(y$lengths) -> indices
  ends <- indices[runs]
  new <- ifelse(runs>1, runs-1,0)
  starts <- indices[ new] +1
  if (0 %in% new ) starts <- c(1, starts)
  
  # Find the sequences of these indices (ie. all indices in "gap" that need to be smoothed)
  gap_indices <- NULL
  for (i in seq_along(starts)) { 
    gap_indices <- c(gap_indices, starts[i]:ends[i])
  }
  # Replace all binary OFFs to ONs
  x[gap_indices] <- ON
  
  return (x)
}


#' Finds the areas of replication in a tract array (ie. runs of ONs or 1s)
#' @param x A numeric vector, tract array 0s/1s.
#' @return A data frame with basic tract information (number, starts, ends, length in pixels).
#' 
find_tracts <- function(x) {
    min_length <- 2
   # Find continuous sequence of ones
   tract_logic <- rle(x == ON)
   tract_runs <- which (tract_logic$values == TRUE & tract_logic$lengths >= min_length)
   
   if (any(tract_runs)){
     # Convert to indices
     tracts_indices <- cumsum(tract_logic$lengths)
     # Find tract end positions
     ends <- tracts_indices[tract_runs]
     # Find tract starts
     new <- ifelse(tract_runs>1, tract_runs-1,0)
     starts <- tracts_indices[ new ] +1
     if (0 %in% new ) starts <- c(1, starts)
   }
   
   # Make and return a table with tract's start and end on same row
   data.frame("Tract No." = 1:length(starts), 
              "Starts At" = starts,
              "Ends At" = ends, 
              "Length" = ends - starts+1)
}


#' Replaces "fork window" (area at tip or ends of tract) in tract array with fork open and close numeric codes.
#' @param x A numeric vector, specify tract array 0s/1s. 
#' @param starts, ends, lengths  Numeric vectors specifying where tracts occur on x. 
#' @param pU, pR Integers indicating the number of pixels into Unreplicated and Replicated regions of each tract. 
#' @return A numeric vector like x but where there are fork region codes at the end of each tract (where 0 and 1s meet)
#' 
add_forks <- function (x, starts, ends, lengths, pU, pR, MAX_L){
  # Find indices of fork based on tract locations
  fork1_starts <- starts - pU
  fork1_ends <- starts + pR -1
  fork2_starts <- ends - pR +1
  fork2_ends <- ends + pU
  
  # Correct for outside zones
  MIN_L <- 1 # minimum length of tract
  fork1_starts[fork1_starts < MIN_L] <- MIN_L
  fork2_starts[fork2_starts < MIN_L] <- MIN_L
  fork1_ends[fork1_ends > MAX_L] <- MAX_L
  fork2_ends[fork2_ends > MAX_L] <- MAX_L
  
  # Correct for fork zone overlap for each tract
  for (i in seq_along(lengths)) 
  {
    if (lengths[i] <= pU+pR){
      mid <- starts[i]+ floor(lengths[i]/2)
      fork1_ends[i] <- mid
      fork2_starts[i] <- mid +1 
    }
    # Replace fork indices with "fork" values to differentiate from unrep/rep areas
    x[ c(fork1_starts[i]:fork1_ends[i]) ] <- FORK_OPEN
    x[ c(fork2_starts[i]:fork2_ends[i]) ] <- FORK_CLOSE
  }
  return (x)
}


# Finds unreplicated, replicated, fork regions and put start and end indices into a table.
#' @param x A numeric array indicating a tract array with forks 0s/1s/2s/3s.
#' @return A data frame with start/end indices of unreplicated, replicated, fork regions.
#' 
find_regions <- function(x) {
  # Find start and ends of regions, as well as lengths
  equal_runs <- rle(x)
  ends <- cumsum(equal_runs$lengths)
  starts <- c(1,ends[-length(ends)]+1)
  
  # Helper function:
  # Depending on codes defined in oddblobs_.R, match each code to its region name (e.g. OFF = 0 = "Unreplicated")
  translate_region_code <- function(element) {
    switch(element,
           OFF = "Unreplicated", # represented by 0 in x
           ON = "Replicated", #1
           FORK_OPEN = "ForkOpen", #2
           FORK_CLOSE = "ForkClose" #3
    )
  }
  # Translate numeric codes of cumulative sums to region name
  region_codes <- sapply(equal_runs$values+1, translate_region_code)
  
  # Return data frame with region info
  data.frame(Region = region_codes,
             Code = equal_runs$values, 
             Start = starts, 
             End = ends)
}


#' Finds all the locations of protein in each region.
#' @param x A numeric protein array 0s/1s.
#' @param starts, ends  Numeric vectors specifying where regions (fork, unrep, rep) occur. 
#' @param nrows The number of regions.
#' @return Arrays of string with 0 to multiple protein indices in each region, e.g. "", "5 6", "9", ""...
#' 
prot_indices_in_region <- function (x, nrows, starts, ends){
  # Initialize a vector with protein indices
  prot_indices <- vector(length = nrows) 
  
  # Assign protein on indices to an array
  prot_ons <- which(x == ON) 
  
  # Go through each row in data frame and see where a protein ON occurs between the start/end indices of a tract
  for(i in 1:nrows) {
  # l <- lapply(1:nrows, function(i){
    x <- prot_ons[which(prot_ons >= starts[i] & prot_ons <= ends[i])]
    # add this vector as a string to the correct pos'n in the vector
    # need to optimize because this is super slow
    # x_as_array <- paste("[",  paste(x, collapse=","), "]", collapse="")
    x_as_array <- paste(x, collapse=",")
    prot_indices[i] <- ifelse(length(x_as_array)==0,"",x_as_array)
  }
  return (prot_indices)
}
  

#' Calculates sizes of tracts, forks, and unreplicated regions, and total respectively.
#' @param x A numeric tract array after finding forks.
#' @return A data frame with one column indicating the sizes of each region.
#' 
get_region_sizes <- function(x) {
  # Calculate sizes of each regions depending on numeric region codes defined in oddblobs_.R 
  rep_size <- sum(x == ON)
  unrep_size <- sum(x == OFF)
  fork_size <- sum(x == FORK_CLOSE) + sum(x == FORK_OPEN)
  total_size <- rep_size + fork_size + unrep_size # OR length(x)
  
  # Make and return data frame
  data.frame(Size=c("Replicated" = rep_size, 
                    "Forks" = fork_size, 
                    "Unreplicated" = unrep_size, 
                    "Total" = total_size))
}


#' Summarizes amount of protein in tracts, forks, and unreplicated regions, and total respectively.
#' @param prot_array A numeric vector of indicating the locations of protein.
#' @param tract_array A numeric vector of indicating the locations of protein.
#' @return A data frame indicating percents of protein in each regions.
#' 
get_prot_percents <- function(prot_array, tract_array) {
  # Get the number of proteins in each region
  prot_in_tracts <- count_matches(prot_array, tract_array, ON)
  prot_in_forks  <- count_matches(prot_array, tract_array, FORK_OPEN) + count_matches(prot_array, tract_array, FORK_OPEN)
  prot_in_unrep <- count_matches(prot_array, tract_array, OFF)
  total_prot_ONS <- sum(prot_array == ON)
  
  # Combine
  prot_counts <- c(prot_in_tracts, prot_in_forks, prot_in_unrep, total_prot_ONS)
  
  # Return protein percents
  prot_counts*100/prot_counts[4] # value at index 4 is total prot ons
}


#' Counts the number of overlapping indices of protein ON and region.
#' @param p, tr Numeric vectors indicating protein and tract arrays.
#' @param region_code Number of region code  (0 in unreplicated, 1 in tracts, 2 in forks - defined in oddblobs_.R).
#' @return number of matches between protein ONs and respective region.
#' 
count_matches <- function(p, tr, region_code) {
  sum(p == ON & tr == region_code)
}

