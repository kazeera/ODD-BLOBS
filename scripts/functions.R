readinteger <- function(msg)
{ 
  num <- readline(prompt=msg)
  return(as.integer(num))
}

makeIntensityPlot<- function(y,title,colour) {
  plot(1:length(y), y, main=title,  xlab="Pixel", ylab="Intensity",type="h",col=colour)
}

#binning - assigning 1 if intesity is greater or equal to threshold, or 0 if it is less
threshold_array <- function(x, threshold) {
  pmax(x>= threshold, 0) 
}

#helps smooth_it
is_smooth_it_needed <- function(x, i, smooth_value) {
  left <- FALSE 
  right <- FALSE
  for (j in i:i+smooth_value) {
    if (j == ON) {
      right == TRUE
      break
    }
  } 
  for (j in i-smooth_value:i) {
    if (j == ON) {
      left == TRUE
      break
    }
  } 
  return (right & left)
}


#smooth_it function closes gaps of less than or equal to smooth_value (pixels) in array x
#need to optimize in R to work on whole arrays (without using packages)

smooth_it <- function(x, smooth_value) {
  #return original array if the smooth value is zero to reduce computation
  if (smooth_value == 0)  { return(x) }
  
  #accounts for gaps - turns NA to zero
  x[is.na(x)] <- 0
  
  size_of_x <- length(x)
  start_of_middle <- smooth_value
  end_of_middle <- size_of_x - smooth_value
  
  for (i in seq_along(x)) {
    if (x[i] == OFF) {
      #if we need to close a gap in the middle of the array
      if (i > start_of_middle &  i < end_of_middle) {
        need_to_smooth_left <- is_smooth_it_needed(x, i-smooth_value, i) 
        need_to_smooth_right <- is_smooth_it_needed(x, i, i+smooth_value)
        
        #if (x[i - smooth_value] == ON  & x[i + smooth_value] == ON &
        if (x[i-1] == ON & need_to_smooth_left & need_to_smooth_right) {
          j <- i
          while(j <= i+smooth_value & j < end_of_middle) {
            if (x[j] == OFF) {
              x[j] <- ON
            }
            j <- j + 1
          }
          i <- j
        }
      }
      #else if we need to close the gaps at the beginning of the array
      else if (i <= start_of_middle) {
        j <- i
        need_to_smooth <- is_smooth_it_needed(x, 1, smooth_value+1)  #1 = first index
        
        while (j <= start_of_middle & need_to_smooth) {
          if (x[j] == OFF) { #assuming data length is larger than smoothvalue
            x[j] <- ON
          }
          j <- j + 1
        }
      }
      #else if we need to close the gaps at the end of the array
      else if (i >= end_of_middle) {
        j <- i
        if (x[size_of_x] == 1){ #TODO Something wrong when last value of array is 1
          #need_to_smooth_right <- help_smooth_it(x, j, size_of_x)
          #need_to_smooth_left <- help_smooth_it(x, end_of_middle, j)
          #need_to_smooth <- need_to_smooth_right & need_to_smooth_left
          need_to_smooth <- is_smooth_it_needed(x, j, size_of_x)
        }
        else {
          need_to_smooth <- is_smooth_it_needed(x, end_of_middle, size_of_x)
        }
        
        if (x[j] == OFF) {
          while(j <= size_of_x  ) {
            if (x[j] == OFF & need_to_smooth) {
              x[j] <- ON
            }
            j <- j + 1
          }
        }
      }
      else {
        print("error doing smooth_it")
        break;
      }
    }
    #if x[i] is 1, go to the next element
  }
  return (x)
}


#takes in the brdu array info and returns a dataframe with basic tract information (number, starts, ends, length in pixels) 
find_tracts <- function(x) {
  min_length <- 2
  #find continuous sequence of ones
  tract_logic <- rle(x == ON)
  tract_runs <- which (tract_logic$values == TRUE & tract_logic$lengths >= min_length)
  
 if (any(tract_runs)){
    #convert to indices
    tracts_indices <-cumsum(tract_logic$lengths)
    #find tract ends
    ends <- tracts_indices[tract_runs]
    #find tract starts
    new <- ifelse(tract_runs>1, tract_runs-1,0)
    starts <- tracts_indices[ new ] +1
    if (0 %in% new ) starts <- c(1, starts)
 }
  else {
    return(NULL)
  }
  
  tract_num <- 1:length(starts)
  tract_lengths <- ends - starts+1
  #populate a matrix with start and end data on same row
  tract_df <- data.frame(tract_num, starts,ends, tract_lengths)
  colnames(tract_df) <- c("Tract No.", "Starts At", "Ends At", "Length")
  return (tract_df)
  
}

#add_forks replaces "fork window" ON/OFF in tract array with fork open and close values 
add_forks <-function (x, starts, ends, lengths, pU, pR, MAX_LENGTH){
  MIN_LENGTH <- 1
  fork1_starts <- starts - pU
  fork1_ends <- starts + pR -1
  fork2_starts <- ends - pR +1
  fork2_ends <- ends + pU
  
  #correct for outside zones
  fork1_starts[fork1_starts < MIN_LENGTH] <- MIN_LENGTH
  fork2_starts[fork2_starts < MIN_LENGTH] <- MIN_LENGTH
  fork1_ends[fork1_ends > MAX_LENGTH] <- MAX_LENGTH
  fork2_ends[fork2_ends > MAX_LENGTH] <- MAX_LENGTH
  
  #correct for fork zone overlap
  for (i in seq_along(lengths)) 
  {
    if (lengths[i] <= pU+pR){
      mid <- starts[i]+ floor(lengths[i]/2)
      fork1_ends[i] <- mid
      fork2_starts[i] <- mid +1 
    }
    #replace fork indices with "fork" values to differentiate from unrep/rep areas
    x[ c(fork1_starts[i]:fork1_ends[i]) ] <- FORK_OPEN
    x[ c(fork2_starts[i]:fork2_ends[i]) ] <- FORK_CLOSE
  }
  return(x)
}

#find unrep, replicated, fork regions and put start and end indices into a table
find_regions <- function(x) {
  equal_runs<- rle(x)
  ends <- cumsum(equal_runs$lengths)
  starts <- c(1,ends[-length(ends)]+1)
  
  translate_region_code <- function(element) {
    switch(element,
           OFF = "Unreplicated",
           ON = "Replicated",
           FORK_OPEN= "ForkOpen",
           FORK_CLOSE= "ForkClose"
    )
  }
  region_codes <- sapply(equal_runs$values+1, translate_region_code)
  region_df <- data.frame(Region=region_codes,Code=equal_runs$values, Start=starts, End=ends)
  return (region_df)
}

#prot_indices_in_region returns vector with strings of protein "ON" indices that occur between start and end of a region(tracts or fork)
prot_indices_in_region <- function (prot_array, nrows, codes, starts, ends){
  prot_indices <- vector(length = nrows) 
  #assign protein on indices to an array
  prot_ons <- which(prot_array == ON) 
  #go through each row in data frame and see where a protein ON occurs between the start/end indices of a tract
  for(i in 1:nrows) {
    x <- prot_ons[which(prot_ons >= starts[i] & prot_ons <= ends[i])]
    #add this vector as a string to the correct pos'n in the vector
    #need to optimize because this is super slow
    prot_indices[i] <- paste(x, collapse=" ")
  }
  return (prot_indices)
}

check_colocalization <- function(prot1_array, prot2_array){
  return ()
}

#returns vector with sizes of tracts, forks, and unreplicated regions, and total respectively
get_region_sizes <- function(tract_array, all_regions) {
  
  #compute total ons in each region (fork, tract, unreplicated)
  rep_size <- sum(tract_array==ON)
  unrep_size <- sum(tract_array==OFF)
  fork_size <- sum(tract_array==FORK_CLOSE | tract_array==FORK_OPEN)
  total_size <- rep_size + fork_size + unrep_size #OR length(x)
  
  #make data frame with this info
  size_df <-data.frame("Total Ons"=c(rep_size, fork_size, unrep_size, total_size))
  
  #compute average sizes of regions
  rep_size <- rep_size/num_of_regions(ON, all_regions)
  unrep_size <- unrep_size/num_of_regions(OFF, all_regions)
  total_forks <- num_of_regions(FORK_CLOSE, all_regions) + num_of_regions(FORK_OPEN, all_regions)
  fork_size <- unrep_size/total_forks
  total_size <- NA
  
  #add to dataframe
  size_df <- cbind(size_df, "Average Size"=c(rep_size, fork_size, unrep_size, total_size))
  #sapply to make this integer
  #sapply(size_df$"Average Size", round)
  
  #add rownames and return
  rownames(size_df) <- c("Tracts", "Forks", "Unreplicated", "Total")
  return(size_df) 
}

num_of_regions <- function(code_key, all_regions) {
  return (sum(all_regions==code_key))
}

#returns vector with info about amount of protein in tracts, forks, and unreplicated regions, and total respectively
get_prot_percents <- function(prot_array, tract_array) {
  prot_in_tracts <- count_matches(prot_array, tract_array, ON)
  prot_in_forks  <- count_matches(prot_array, tract_array, FORK_OPEN) + count_matches(prot_array, tract_array, FORK_OPEN)
  prot_in_unrep <- count_matches(prot_array, tract_array, OFF)
  total_prot_ONS <- sum(prot_array == ON)
  
  prot_counts <- c(prot_in_tracts, prot_in_forks, prot_in_unrep, total_prot_ONS) 
  #conditional so that 0 is not in denominator
  prot_counts[4]<- ifelse(prot_counts[4]!= 0,total_prot_ONS, 1)
  prot_percents <- prot_counts*100/prot_counts[4] #value at index 4 is total prot ons
  
  return (prot_percents)
}

#returns the number of matches between prot ONs and respective region (region_val = 0 in unreplicated, 1 in tracts, 2 in forks)
count_matches <- function(prot_array, tract_array, region_code) {
  return (sum(prot_array == ON & tract_array == region_code))
}

#prot_indices_in_region returns vector with strings of protein "ON" indices that occur between start and end of a region(tracts or fork)
count_occurrences <- function (on_indices, starts, ends){
  #count number of occurences in each range
  counts <- sapply(1:length(starts), function(x) sum(on_indices>=starts[x] & on_indices<=ends[x]))
  return (counts)
}



