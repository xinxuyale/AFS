# Alpha-shape and Lab Source Fitting to Spectrum algorithm (ALSFS)

# Load in the 'alphahull' package.
library("alphahull")

# Input variables:

# order: the order of spectrum to remove blaze function. It is an n by 2 matrix,
#   where n is the number of pixels. Each row is the wavelength and intensity at 
#   each pixel.
# led: the corresponding order of lab source spectrum. It is also an n by 2 matrix.
# q: the parameter q, uppder q quantile within each window will be used to 
#   do linear transformation on the lab source spectrum. 
# d: the smoothing parameter for local polynomial regression, which is the 
#   proportion of neighboring points to be used when fitting at one point. 

# Define the ALSFS function 
ALSFS <- function(order, led, q = 0.95, d=0.25) {
  # Default value of q and d are 0.95 and 0.25.
  # Change the column names and format of the dataset. 
  colnames(order) <- c("wv", "intens") 
  order <- data.frame(order)
  # n records the number of pixels.
  n <- dim(order)[1]
  # Variable u is the parameter u in the step 1 of ALSFS algorithm. It scales the intensity vector.
  u <- (range(order$wv)[2]-range(order$wv)[1])/10/max(order$intens)
  order$intens <- order$intens*u 
  
  # Let alpha be 1/6 of the wavelength range of the whole order. 
  alpha <- (range(order$wv)[2]-range(order$wv)[1])/6
  # Use ahull() function from package "alphahull" to run alpha shape on the order.
  obj<- ahull(order$wv, order$intens, alpha)
  # Variable as is a k by 2 matrix. Each row records the index of the two vertices of 
  # a segment in the alpha shape. k is the number of segments in this alpha shape.
  as <- obj$arcs[obj$arcs[,3]>0,7:8]

  
  # This chunk of code detects loops in the boundary of the alpha shape. 
  # Ususally there is only one loop and so the indices just one column of variable as.
  # Variable loop is a list.
  # The indices of the k-th loop are recorded in the k-th element of variable loop. 
  loop <- list()
  k <- 0
  while(length(dim(as))>0){
    k <- k+1
    loop[[k]] <- as[1,]
    as <- as[-1,]
    while(length(dim(as))>0 & loop[[k]][length(loop[[k]])]!=loop[[k]][1]){
      if (loop[[k]][length(loop[[k]])] == as[1,1]) {
        loop[[k]] <- c(loop[[k]], as[1,2])
        as <- as[-1,]
      } else {
        ind <- which(as == loop[[k]][length(loop[[k]])])
        if (ind <= dim(as)[1]) {
          loop[[k]] <- c(loop[[k]], as[ind,2])
          as <- as[-ind,]
        } else {
          ind <- ind-dim(as)[1]
          loop[[k]] <- c(loop[[k]], as[ind,1])
          as <- as[-ind,]
        }
      }
    }
  }
  loop[[k]] <- c(loop[[k]], loop[[k]][1])
  
  # Use the loops to get the set W_alpha. 
  # Variable Wa is a vector recording the indices of points in W_alpha.
  Wa <- 0
  for (k in 1:length(loop)){
    temp <- loop[[k]]
    temp <- temp[-length(temp)]
    temp <- temp[temp<n]
    max_k <- max(temp)
    min_k <- min(temp)
    len_k <- length(temp)
    as_k <- temp
    if (!(as_k[1] == min_k & as_k[len_k] == max_k)) {
      if (which(as_k==min_k) < which(as_k==max_k)) as_k <- as_k[which(as_k==min_k):which(as_k==max_k)]
      if (which(as_k==min_k) > which(as_k==max_k)) as_k <- c(as_k[which(as_k==min_k):len_k], as_k[1:which(as_k==max_k)])
    }
    Wa <- c(Wa, as_k)
  }
  Wa <- sort(Wa)
  Wa <- Wa[-1]
  
  # AS is an n by 2 matrix recording tilde(AS_alpha). Each row is the wavelength and intensity of one pixel.
  AS <- order
  for (i in 1:(n-1)){
    a <- Wa[which(i<Wa)[1]-1]
    b <- Wa[which(i<Wa)[1]]
    AS$intens[i] <- AS$intens[a]+(AS$intens[b]-AS$intens[a])*((AS$wv[i]-AS$wv[a])/(AS$wv[b]-AS$wv[a]))
  }
  # Run a local polynomial on tilde(AS_alpha), as described in step 3 of the AFS algorithm.
  # Use the function loess() to run a second order local polynomial.
  # Variable AS_fit is the local polynomial regression object.
  AS_fit <- loess(intens ~ wv, data = AS, degree = 2, span = d, control = loess.control(surface = "direct"))
  # B1 records hat(B_1), which is the primary estimate of the blaze function.
  B1 <- predict(AS_fit, data.frame(wv = order$wv))
  
  # Add a new column called select to the matrix order. 
  # order$select records hat(y^(1)).
  order$select <- order$intens/B1
  # Calculate Q_2q-1 in step 3 of the ALSFS algorithm.
  Q <- quantile(order$select, 1-(1-q)*2)
  
  # Make indices in Wa to the format of small windows. 
  # Each row of the variable window is a pair of neighboring indices in Wa.
  window <- cbind(Wa[1:(length(Wa)-1)], Wa[2:length(Wa)])

  # This chunk of code select the points whose intensities are in the top q quantile 
  # within each window and also in the overall top Q_2q-1 quantile.
  # The point indices are recorded in variable index, which is S_alpha, q in step 4
  # of the ALSFS algorithm.
  index <- 0
  for (i in 1:(dim(window)[1])) {
    loc_window <- window[i,]
    temp <- order[loc_window[1]:loc_window[2],]
    index.i <- which(temp$select>=quantile(temp$select, q) & temp$select >= Q)
    index.add <- loc_window[1]+index.i-1 
    index <- c(index, index.add)
  }
  index <- unique(index[-which(index==0)])
  index <- sort(index)
  
  # The following chunk of code does step 5 of the ALSFS algorithm.
  # The function optim() is used to calculate the optimization of the three 
  # linear transformation parameters.
  # The final estimate is in variable B2.
  m <- length(index)
  led[,2] <- led[,2]/max(led[,2])*max(order$intens)
  Xnew <- cbind(rep(1, m), led[index,2], led[index,1])
  beta <- c(0, 1, 0)
  fn <- function(x) sum((order$intens[index]/(Xnew%*%x)-1)^2)
  beta <- optim(beta, fn)$par
  B2 <- beta[1]+beta[2]*led[,2]+beta[3]*led[,1]
  
  # Return the blaze-removed spectrum.
  return(order$intens/B2)
}

