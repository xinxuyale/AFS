# Boundary Correction

# Input variables:

# order1: the left order. An n1 by 2 matrix.
# order2: the right order. An n2 by 2 matrix.

# Define the BC function 
BC <- function(order1, order2) {
  # Change the column names and format of the dataset. 
  colnames(order1) <- c("wv", "intens") 
  order <- data.frame(order1)
  colnames(order2) <- c("wv", "intens") 
  order <- data.frame(order2)
  # n1, n2 records the number2 of pixels of order1 and order2.
  n1 <- dim(order1)[1]
  n2 <- dim(order2)[1]  
  
  # start: the start index of the overlap region in the left order.
  # end: the end index of the overlap region in the right order.
  start <- which(order1$wv>=order2$wv[1])[1]
  end <- which(order2$wv>order1$wv[n1])[1]-1
  # temp: records the intensities of the overlap region.
  temp <- order1$intens[start:n1]

  # This chunk of code does the weighted average. 
  n_temp <- length(temp)
  for (i in 1:n_temp){
    w <- i/n_temp
    temp[i] <- (1-w)*order1$intens[start+i-1]+w*order2$intens[i]
  }  

  # Revise the original orders with the corrected region.
  order1$intens <- c(order1$intens[1:(start-1)], temp)
  order2$intens <- c(temp, order2$intens[(end+1):n2])

  # Return the two corrected orders in a list variable.
  return(list(order1, order2))
}
