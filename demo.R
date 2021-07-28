# Demo on MNIST data

required_packages <- c("dplyr", "fda", "lsei")
(missed_packages <- setdiff(required_packages, rownames(installed.packages())))
if (length(missed_packages)) {
  sapply(missed_packages, install.packages)
}



load("Mnist.RData")
source("funs_2DFPCA.R")


library(dplyr)
library(fda)
library(lsei)


fulldata <- observed
selectindex <- which(label == 0)
observed <- fulldata[selectindex]

meanobserve <- Reduce("+", observed)/length(selectindex)
observed <- lapply(observed, function(x) {
  x-meanobserve
})
timepoints1 <- timepoints1[selectindex]
timepoints2 <- timepoints2[selectindex]
label <- label[selectindex]



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


# Find the first 2 FPCs
numFPC <- 2
for (i in 1:numFPC) {
  if (i == 1) {
    res_first <- first_FPC_2d_image(beta1, observed, timepoints1, timepoints2, basis1, basis2, threshold = 1e-4)
    result_list[[1]] <- res_first
    previous_beta[[1]] <- res_first$beta
    pc_list[[1]] <- res_first$pc_fit
  } else {
    res_second <- second_FPC_conditional_2d_image(beta1, pc_index = i, observed, timepoints1, timepoints2, basis1, basis2, betalist = previous_beta, threshold = 1e-4)
    result_list[[i]] <- res_second
    previous_beta[[i]] <- res_second$beta
    pc_list[[i]] <- res_second$pc_fit
  }
}
