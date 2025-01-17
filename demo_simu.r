# Demo on Simulation studies

# Installing the required packages
required_packages <- c("dplyr", "fda", "lsei")
(missed_packages <- setdiff(required_packages, rownames(installed.packages())))
if (length(missed_packages)) {
  sapply(missed_packages, install.packages)
}
library(dplyr)
library(fda)
library(lsei)

# Load the data and source code
load("fpc_0and1x.RData")
source("funs_2DFPCA.R")


# Initialize the parameters and read in the true surfaces
result_list_global <<- result_list

values_global <- list()
for (pcindex in 1:3) {
  valuev <- sapply(1:length(timepoints1[[1]]), function(index) {
    tempvalue2 <- eval.fd2d.image(index, result_list_global[[pcindex]]$pc_fit)
    tempvalue2
  })
  values_global[[pcindex]] <- valuev
}



# Function for simulating data surfaces where the error follows normal distribution
simu_gauss <- function(numofsimu = 100) {
  simulatedobserved <- list()

  for (iter in 1:numofsimu) {
    outputv <- c()

    var1 <- rnorm(1, 0, 100)
    var2 <- rnorm(1, 0, 50)
    var3 <- rnorm(1, 0, 25)
    varlist <- c(var1, var2, var3)
    error <- rnorm(length(timepoints1[[1]]), 0, 0.25)
    for (index in 1:length(timepoints1[[1]])) {
      e <- error[index]
      tempvalue <- 0
      
      for (pcindex in 1:numofsimupc) {
        tempvalue <- tempvalue + values_global[[pcindex]][index] * varlist[pcindex] + e
      }

      outputv <- c(outputv, tempvalue)
    }
    simulatedobserved[[iter]] <- outputv
  }
  return(simulatedobserved)
}


# Function for simulating data surfaces where the error follows non-normal distribution
simu_nongauss <- function(numofsimu = 100) {
  simulatedobserved <- list()

  for (iter in 1:numofsimu) {
    outputv <- c()

    var1 <- rgamma(1, 1, 0.0100) - 100
    var2 <- rgamma(1, 1, 0.0250) - 40
    var3 <- rgamma(1, 1, 0.05) - 20
    varlist <- c(var1, var2, var3)
    error <- rnorm(length(timepoints1[[1]]), 0, 0.25)
    for (index in 1:length(timepoints1[[1]])) {
      e <- error[index]
      tempvalue <- 0
      
      for (pcindex in 1:numofsimupc) {
        tempvalue <- tempvalue + values_global[[pcindex]][index] * varlist[pcindex] + e
      }

      outputv <- c(outputv, tempvalue)
    }
    simulatedobserved[[iter]] <- outputv
  }
  return(simulatedobserved)
}




# Number of FPCs
numofsimupc <- 3

# Number of replications
numreplica <- 50

# Number of trajectories
ncurvelist <- c(100, 400)

# Number of basis splines 
nbasislist <- c(8, 12, 16)
ind <- 1

# Simulating from normal or non-normal error distribution 
gauss <- TRUE

outputlist <- list()
for (iii in 1:length(ncurvelist)) {
  for (jjj in 1:length(nbasislist)) {
    nbasis <- nbasislist[jjj]
    numofcurve <- ncurvelist[iii]
    nb <- nbasis
    nc <- numofcurve

	# Initialize the spline parameters 
    maxpixelsize <- 28
    library(fda)
    beta1 <- matrix(1, nrow = nbasis, ncol = nbasis)
    for (i in 1:nbasis) {
      beta1[, i] <- seq(0.01, 0.01 * nbasis, by = 0.01)
    }
    basis1 <- create.bspline.basis(rangeval = c(1, maxpixelsize), nbasis = nbasis, norder = 4)
    basis2 <- create.bspline.basis(rangeval = c(1, maxpixelsize), nbasis = nbasis, norder = 4)
    initializeGlobalXmat(timepoints1, timepoints2, basis1, basis2)


    monteresult <- list()
    for (monte in 1:numreplica) {
      
      previous_beta <- list()
      pc_list <- list()
      result_list <- list()

		# Simulating the data surfaces
      if (gauss == TRUE) {
        sobserved <- simu_gauss(numofcurve)
      }
      if (gauss == FALSE) {
        sobserved <- simu_nongauss(numofcurve)
      }

      for (i in 1:numofsimupc) {
        if (i == 1) {
			# First FPC
          res_first <- first_FPC_2d_image(beta1, sobserved, timepoints1, timepoints2, basis1, basis2, threshold = 1e-5, minit = 3)
          result_list[[1]] <- res_first
          previous_beta[[1]] <- res_first$beta
          pc_list[[1]] <- res_first$pc_fit
        } else {
			# Higher order FPCs
          res_higherorder <- second_FPC_conditional_2d_image(beta1, pc_index = i, sobserved, timepoints1, timepoints2, basis1, basis2, betalist = previous_beta, threshold = 1e-5, minit = 1)
          result_list[[i]] <- res_higherorder
          previous_beta[[i]] <- res_higherorder$beta
          pc_list[[i]] <- res_higherorder$pc_fit
        }
      }

		# Evaluating the mean squared error of the FPCs
      result <- c()
      for (pcindex in 1:numofsimupc) {
        outputv <- c()
        o2 <- c()
        o1 <- c()
        for (index in 1:length(timepoints1[[1]])) {
          tempvalue <- eval.fd2d.image(index, result_list[[pcindex]]$pc_fit)
          tempvalue2 <- values_global[[pcindex]][index]

          outputv <- c(outputv, min((tempvalue - tempvalue2)^2, (tempvalue + tempvalue2)^2))
          o2 <- c(o2, tempvalue2)
          o1 <- c(o1, tempvalue)
        }
        result <- c(result, mean(outputv))
      }

		# Evaluating the integrated mean squared error of the FPCs
      outputv2 <- c()
      for (pcindex in 1:numofsimupc) {
        xy <- 2 * inprod.fd2dx(result_list[[pcindex]]$pc_fit, result_list_global[[pcindex]]$pc_fit)
        x2 <- inprod.fd2d(result_list[[pcindex]]$pc_fit, result_list[[pcindex]]$pc_fit)
        y2 <- inprod.fd2d(result_list_global[[pcindex]]$pc_fit, result_list_global[[pcindex]]$pc_fit)

        outputv2 <- c(outputv2, min(x2 + y2 - xy, x2 + y2 + xy))
      }

      monteresult[[monte]] <- c(result, outputv2)
    }

	# Store the results
    tt <- monteresult
    result <- apply(tt %>% do.call(rbind, .), 2, mean)

    outputlist[[ind]] <- c(nb, nc, result)
    ind <- ind + 1
  }
}


outputmat <- Reduce(rbind, outputlist)
outputmat[,3:5] <- outputmat[,3:5] * 1000
row.names(outputmat) <- rep("", nrow(outputmat))
colnames(outputmat) <- c('Number of splines', 'Number of curves', 'MSE1', 'MSE2', 'MSE3',"IMSE1", 'IMSE2', 'IMSE3' )

outputmat