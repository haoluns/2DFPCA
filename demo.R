# Demo on MNIST data

# Installing the required packages
required_packages <- c("dplyr", "fda", "lsei", "raster")
(missed_packages <- setdiff(required_packages, rownames(installed.packages())))
if (length(missed_packages)) {
  sapply(missed_packages, install.packages)
}
library(dplyr)
library(fda)
library(lsei)

# Load the data and source code
load("Mnist.RData")
source("funs_2DFPCA.R")




# Only use the figure 0 in the MNIST data
fulldata <- observed
selectindex <- which(label == 0)
observed <- fulldata[selectindex]

# Center the data 
meanobserve <- Reduce("+", observed)/length(selectindex)
observed <- lapply(observed, function(x) {
  x-meanobserve
})
timepoints1 <- timepoints1[selectindex]
timepoints2 <- timepoints2[selectindex]
label <- label[selectindex]


# Initialize the spline parameters 
maxpixelsize <- 28
nbasis <- 12
library(fda)
beta1 <- matrix(1, nrow = nbasis, ncol = nbasis)
for (i in 1:nbasis) {
  beta1[, i] <- seq(0.01, 0.12, by = 0.01)
}

basis1 <- create.bspline.basis(rangeval = c(1, maxpixelsize), nbasis = nbasis, norder = 4)
basis2 <- create.bspline.basis(rangeval = c(1, maxpixelsize), nbasis = nbasis, norder = 4)

initializeGlobalXmat(timepoints1, timepoints2, basis1, basis2)
	
previous_beta <- list()
pc_list <- list()
result_list <- list()




# Fit the first 3 FPCs
numFPC <- 3
for (i in 1:numFPC) {
  if (i == 1) {
	
	 # First FPC
    res_first <- first_FPC_2d_image(beta1, observed, timepoints1, timepoints2, basis1, basis2, threshold = 1e-4)
    result_list[[1]] <- res_first
    previous_beta[[1]] <- res_first$beta
    pc_list[[1]] <- res_first$pc_fit
  } else {
    
	 # Higher order FPCs
    res_higherorder <- second_FPC_conditional_2d_image(beta1, pc_index = i, observed, timepoints1, timepoints2, basis1, basis2, betalist = previous_beta, threshold = 1e-4)
    result_list[[i]] <- res_higherorder
    previous_beta[[i]] <- res_higherorder$beta
    pc_list[[i]] <- res_higherorder$pc_fit
  }
}



# Plot the first 3 FPCs	
par(mfrow = c(1,3))
for (i in 1:numFPC){

  plotindex <- i
  res <- result_list[[plotindex]]
  sfit <- result_list[[length(result_list)]]$sfit 	

  plotmat = matrix(0, ncol = 28, nrow = 28)	
  for (x in 1:28){
    for (y in 1:28){
	    plotmat[x, y] <- (eval.fd2d(x, y, res$pc_fit) * 1000)
    }
  }
  
  library(raster)
  plot(raster(plotmat), col = grey.colors(10, start = 0, end = 1), main = paste0("FPC ", plotindex))
}

