library(dplyr)
library(fda)
library(lsei)

create.2d.basis = function(b1,b2){
	#spline_basis1 =create.bspline.basis(rangeval=dat$Month_from_baseline%>%range,nbasis=10,norder=4)
	#spline_basis2 =create.bspline.basis(rangeval=dat$Month_from_baseline%>%range,nbasis=10,norder=4)
	return (list(basis1=b1,basis2=b2))
}

fd2d = function(betas,basis2d){
	#betas is a matrix of coefficients, beta_ij = coefficient for the ith basis1 and jth basis2
	prod1 = inprod(basis2d$basis1, basis2d$basis1,0,0)
	prod2 = inprod(basis2d$basis2, basis2d$basis2,0,0)
	
	return(list(betas=betas,basis1=basis2d$basis1,basis2=basis2d$basis2,prod1=prod1,prod2=prod2))
}


fd2d = function(betas,basis1,basis2){
	#betas is a matrix of coefficients, beta_ij = coefficient for the ith basis1 and jth basis2
	prod1 = inprod(basis1, basis1,0,0)
	prod2 = inprod(basis2, basis2,0,0)
	
	return(list(betas=betas,basis1=basis1,basis2=basis2,prod1=prod1,prod2=prod2))
}



norm.fd2d=function(fd2dobject){

betas = fd2dobject$betas
basis1 = fd2dobject$basis1
basis2 = fd2dobject$basis2
prod1 = fd2dobject$prod1
prod2 = fd2dobject$prod2

numbasis1=basis1$nbasis
numbasis2=basis2$nbasis



temp=0
 for(i in 1:numbasis1){
 for(j in 1:numbasis2){
 for(ip in 1:numbasis1){
 for(jp in 1:numbasis2){
 temp = temp +betas[i,j]*betas[ip,jp]*prod1[i,ip]*prod2[j,jp]
}
 }
 }
 }
 
return (temp)
}





inprod.fd2d=function(fd2dobject1,fd2dobject2){
#Both objects should share the same basis
betas1 = fd2dobject1$betas
basis1 = fd2dobject1$basis1
basis2 = fd2dobject1$basis2
prod1 = fd2dobject1$prod1
prod2 = fd2dobject1$prod2

betas2 = fd2dobject2$betas

numbasis1=basis1$nbasis
numbasis2=basis2$nbasis



temp=0
 for(i in 1:numbasis1){
 for(j in 1:numbasis2){
 for(ip in 1:numbasis1){
 for(jp in 1:numbasis2){
 temp = temp +betas1[i,j]*betas2[ip,jp]*prod1[i,ip]*prod2[j,jp]
}
 }
 }
 }
 
return (temp)
}



inprod.fd2dx=function(fd2dobject1,fd2dobject2){
#Both objects no need share the same basis
betas1 = fd2dobject1$betas
basis11 = fd2dobject1$basis1
basis21 = fd2dobject1$basis2

basis12 = fd2dobject2$basis1
basis22 = fd2dobject2$basis2


prod1 = inprod(basis11, basis12,0,0)
prod2 = prod1

betas2 = fd2dobject2$betas

numbasis11=basis11$nbasis
numbasis21=basis21$nbasis
numbasis12=basis12$nbasis
numbasis22=basis22$nbasis



temp=0
 for(i in 1:numbasis11){
 for(j in 1:numbasis21){
 for(ip in 1:numbasis12){
 for(jp in 1:numbasis22){
 temp = temp +betas1[i,j]*betas2[ip,jp]*prod1[i,ip]*prod2[j,jp]
}
 }
 }
 }
 
return (temp)
}



constraint.fd2d=function(fd2dobject1){
#Both objects should share the same basis

betas1 = fd2dobject1$betas

prod1 = fd2dobject1$prod1
prod2 = fd2dobject1$prod2

const = betas1

numbasis1=basis1$nbasis
numbasis2=basis2$nbasis

for (i in 1:numbasis1){
for (j in 1:numbasis2){
	const[i,j] = sum( betas1*prod1[,i] %*% t(prod2[,j]))
}}


return (c(const))
}



eval.fd2d=function(timei1,timei2, fd2dobject){
# sum_i  sum_ j   beta_ij, phi_i(t1)   phi_j(t2)

betas = fd2dobject$betas
basis1 = fd2dobject$basis1
basis2 = fd2dobject$basis2

T = length(timei1)
value = rep(0,T)
for (t in 1:T){

	betanew = sapply(1:dim(betas)[1],function(r){
	temppc = fd(betas[r,],basis1)
	eval.fd(timei1[t],temppc)
	}
	)
	value[t] = eval.fd(timei2[t],fd(betanew,basis2))

}
return (value)

}



eval.fd2d.image<-function(indextimei, fd2dobject){

sapply(indextimei,function(x){sum(fd2dobject$betas * globalxmatlist[[x]])})
}



eval.fd2d.image.nonstandard <-function(indextimei, fd2dobject,timei1,timei2){

xmatlist  <- lapply(1:length(timei1),function(t){
	
	evalb1 = eval.basis(timei1[t], basis1)
	evalb2 = eval.basis(timei2[t], basis2)
	(t(evalb1) %*% evalb2)				
	}
)

sapply(indextimei,function(x){sum(fd2dobject$betas *xmatlist[[x]])})	

}


eval.fd2d.image.nonstandard2 <-function(indextimei, fd2dobject,xmatlist){

sapply(indextimei,function(x){sum(fd2dobject$betas *xmatlist[[x]])})	

}

	
	

initializeGlobalXmat<-function(timepoints1,timepoints2,basis1,basis2){
timei1 =timepoints1[[1]]
timei2 =timepoints2[[1]]	

globalxmat  <<- lapply(1:length(timei1),function(t){
	
	evalb1 = eval.basis(timei1[t], basis1)
	evalb2 = eval.basis(timei2[t], basis2)
	c(t(evalb1) %*% evalb2)				 #Need check
	}
	)%>%do.call(rbind,.)


globalxmatlist  <<- lapply(1:length(timei1),function(t){
	
	evalb1 = eval.basis(timei1[t], basis1)
	evalb2 = eval.basis(timei2[t], basis2)
	(t(evalb1) %*% evalb2)				
	}
	)
	
}

	

first_FPC_2d_image = function(beta1,observed,timepoints1,timepoints2,basis1,basis2,threshold=1e-3,minit=3){
	thresh = 1
	it = 1

	value = -1e14
	pc_fit = fd2d(beta1, basis1,basis2)
	pc_fit$betas = 1/sqrt(norm.fd2d(pc_fit))*pc_fit$betas
	
	beta1 = pc_fit$betas
	
	tempsum = t(globalxmat%>%as.matrix)%*%(as.matrix(globalxmat))
		
		
	while (thresh >  threshold|it<minit){
	
		beta1_before = beta1
		value_before = value
		
		timei1 =timepoints1[[1]]
		timei2 =timepoints2[[1]]
		#iterationxmat = eval.fd2d(timei1,timei2,pc_fit)
		iterationxmat = eval.fd2d.image(1:length(timei1),pc_fit)
		
		alpha_fit = function(subj_index){
			timei1 =timepoints1[[subj_index]]
			timei2 =timepoints2[[subj_index]]
			
			xmati = iterationxmat
			model = lm.fit(x=as.matrix(xmati),y=as.matrix(observed[[subj_index]]))			
			
			#model = lm(observed[[subj_index]]~0+xmati)
			sf = model$coefficient%>%as.numeric
			rf = ((model%>%residuals)^2)%>%mean
			c(sf,rf)			
		}
	
		result = lapply((1:length(observed)),alpha_fit)%>%do.call(rbind,.)
		sfit = result[,1]


		rfit=  result[,2]
		value = mean(rfit^2)
		#print(value)
		if(abs(value_before - value/value_before) < threshold) break

				
		yem = do.call(c,observed)
		
		
		
		if (FALSE){
			temp = matrix(0,nrow=144,ncol=144)
			
			runlist=c()
			for(ii in 1:144){
			for(jj in 1:144){
			runlist = rbind(runlist,c(ii,jj))
			}}
			
			sfitlong = lapply(1:length(observed),function(subj_index){
					rep(sfit[subj_index],dim(globalxmat)[1])
				})%>%do.call(c,.)
				
			temp = lapply(1:dim(runlist)[1],function(index){
				ii=runlist[index,1]
				jj=runlist[index,2]
											
				temp1=rep(globalxmat[,ii],length(observed))*sfitlong
				temp2=rep(globalxmat[,jj],length(observed))*sfitlong
				
				return(sum(temp1*temp2))
			})
			
			for(index in 1:dim(runlist)[1]){
				ii=runlist[index,1]
				jj=runlist[index,2]
											
				temp1=rep(globalxmat[,ii],length(observed))*sfitlong
				temp2=rep(globalxmat[,jj],length(observed))*sfitlong
				
				temp[ii,jj]=(sum(temp1*temp2))
				#print(index)
			
			}
			
			
			for(ii in 1:144){
			for(jj in 1:144){
				xalpha = lapply(1:length(observed),function(subj_index){
					xmati = cbind(globalxmat[,ii],globalxmat[,jj])
					
					xmati*sfit[subj_index]

				})%>%do.call(rbind,.)
				
				temp[ii,jj]=sum(xalpha[,1]*xalpha[,2])
				
			}
			}
		}
		
		
		
		
		#memory overflow issue
		if(FALSE){
		xalpha = lapply(1:length(observed),function(subj_index){
			xmati = globalxmat
			
			xmati*sfit[subj_index]

		})%>%do.call(rbind,.)
		}
		
		xalpha2 = lapply(1:length(observed),function(subj_index){
										
					tempsum*sfit[subj_index]^2

				})%>%Reduce('+',.)	
		
				
		
		
		globalxmatt= t(globalxmat %>% as.matrix)
		xyalpha = lapply(1:length(observed),function(index){
			(globalxmatt %*% observed[[index]])*sfit[index]
		})%>%Reduce('+',.)	
			
		#beta1 = solve(t(xalpha%>%as.matrix)%*%(as.matrix(xalpha)) , t(xalpha%>%as.matrix)%*%as.matrix(yem))
		beta1 = solve(xalpha2 ,xyalpha)
		beta1 = matrix(beta1,nrow = basis1$nbasis,ncol=basis2$nbasis)
		pc_fit = fd2d(beta1, basis1,basis2)
		pc_fit$betas = 1/sqrt(norm.fd2d(pc_fit))*pc_fit$betas
		beta1 = pc_fit$betas
		
		thresh =  max(abs(c(beta1_before)-c(beta1)))
		it = it+1
		 {print(paste0("Iteration ",it));print(paste0("Tol: ",as.numeric(thresh)));}


	}##while end
	
	
	return(list(beta=beta1,pc_fit = pc_fit,sfit=sfit,thresh = thresh,it=it,value=value))
	
}



third_FPC_conditional_2d_image  = function(beta3, pc_index, observed,timepoints1,timepoints2,basis1,basis2, betalist  , threshold=1e-4,minit=1){


start_time <- Sys.time()


thresh = 1
it = 1
value = -1e14

pc_fit = fd2d(beta3, basis1,basis2)
pc_fit$betas = 1/sqrt(norm.fd2d(pc_fit))*pc_fit$betas
beta3 = pc_fit$betas
	

pc_fits_previous = lapply(1:length(betalist), function(x){
pc_fit1 = fd2d(betalist[[x]], basis1,basis2)
pc_fit1$betas = 1/sqrt(norm.fd2d(pc_fit1))*pc_fit1$betas
pc_fit1
})


observed2 =observed
observed2[which(sapply(observed,length)<=pc_index)]=NULL
timepoints1.2 = timepoints1
timepoints1.2[which(sapply(observed,length)<=pc_index)]=NULL
timepoints2.2 = timepoints2
timepoints2.2[which(sapply(observed,length)<=pc_index)]=NULL


timei1 =timepoints1.2[[1]]
timei2 =timepoints2.2[[1]]
xmat_previous_list = lapply(pc_fits_previous, function(x){
		#eval.fd2d(timei1,timei2,x)
		eval.fd2d.image(1:length(timei1),x)
		})
		
xmat_previous = xmat_previous_list%>%do.call(cbind,.)


tempsum = t(globalxmat%>%as.matrix)%*%(as.matrix(globalxmat))
globalxmatt= t(globalxmat %>% as.matrix)
	
cmat = lapply(pc_fits_previous,constraint.fd2d)%>%do.call(rbind,.)

end_time <- Sys.time()
#print(end_time - start_time	)

convcount = 0		
while (thresh >  threshold|it<minit){
	
	
start_time <- Sys.time()


	beta3_before = beta3
	value_before = value
	#iterationxmat = eval.fd2d(timei1,timei2,pc_fit)
	iterationxmat = eval.fd2d.image(1:length(timei1),pc_fit)
	
	alpha_fit = function(subj_index){
		
		timei1 =timepoints1.2[[subj_index]]
		timei2 =timepoints2.2[[subj_index]]
						
		xmati = iterationxmat
		model = lm.fit(y=as.matrix(observed2[[subj_index]]), x=cbind(xmat_previous,xmati))
		
		#model = lm(observed2[[subj_index]]~xmat_previous+xmati+0)
		sf = model$coefficient%>%as.numeric
		rf = ((model%>%residuals)^2)%>%mean
		c(sf,rf)
	}

	result= lapply(1:length(observed2),alpha_fit)%>%do.call(rbind,.)
	sfit = result[,1:pc_index]
	rfit=  result[,ncol(result)]
	value = as.numeric(mean(rfit^2))

	#if(abs(value_before - value/value_before) < threshold) {print('ResidualBreak');break;}
	
	
	N = sapply(observed2,length)%>%sum

	
	if (FALSE){
	yem = lapply(1:length(observed2), function(index){
		
	yfits_previous = lapply(1:length(pc_fits_previous), function(x){
		sfit[index,x]*xmat_previous_list[[x]]%>%as.numeric
		})%>%do.call(rbind,.)%>%colSums

	(observed2[[index]] - yfits_previous)/sqrt(N)
	})%>%do.call(c,.)

			
	xalpha = lapply(1:length(observed2),function(subj_index){

			
			xmati = globalxmat

			(xmati*sfit[subj_index,ncol(sfit)])/sqrt(N)

	})%>%do.call(rbind,.)
	A = xalpha
	qmat = 2*(t(A)%*%A)
	pmat = as.numeric(-2*t(yem)%*%A)
	}
		
	###new 
	
	xalpha2 = lapply(1:length(observed2),function(subj_index){
									
				tempsum*sfit[subj_index,ncol(sfit)]^2

			})%>%Reduce('+',.)/N
				
	yemobs = lapply(1:length(observed2), function(index){
		
	yfits_previous = lapply(1:length(pc_fits_previous), function(x){
		sfit[index,x]*xmat_previous_list[[x]]%>%as.numeric
		})%>%do.call(rbind,.)%>%colSums

	(observed2[[index]] - yfits_previous)
	})
	
	
	xyalpha = lapply(1:length(observed),function(index){
		(globalxmatt %*% yemobs[[index]])*sfit[index,ncol(sfit)]/N
	})%>%Reduce('+',.)	
		
	pmat = as.numeric(-2*xyalpha)		
	qmat = 2*xalpha2
	
	
	
	
	beta3  = lsei::qp(qmat, pmat, cmat, d=rep(0,length(betalist)))
	#sum(as.numeric(cmat)*c(beta3))

	beta3 = matrix(beta3,nrow = basis1$nbasis,ncol=basis2$nbasis)
	
	if (abs(beta3_before[1,1]) -abs(beta3[1,1]) < 0.001){
		if(beta3_before[1,1]*beta3[1,1] < 0 ){
			beta3= (-1)* beta3
		}
	}
	
	
	pc_fit = fd2d(beta3, basis1,basis2)
	pc_fit$betas = 1/sqrt(norm.fd2d(pc_fit))*pc_fit$betas
	beta3 = pc_fit$betas
    #inprod.fd2d(pc_fits_previous[[1]],pc_fit)
	
	threshbefore = thresh	 
	thresh =  max(abs(c(beta3_before))-abs(c(beta3)))
	if (thresh > threshbefore){convcount = convcount + 1}
	
	if(convcount > 5){
		break;
	}
	
	loss_value= 0  ##mean((yem - A%*%c(beta3))^2)*10000  ##cannot compute the loss memory too much
	
	it = it+1
	#templist[[it]]=beta3

	
	{print(paste0("FPC ", pc_index, ",Iteration ",it));print(paste0("Tol: ",as.numeric(thresh)));}
	
end_time <- Sys.time()
#print(end_time - start_time	)
}#while end

return(list(beta=beta3,pc_fit = pc_fit, sfit = sfit, previous_beta =betalist, thresh = thresh,it=it,value=value,gamma=gamma))

}




pred_soap_2d_fast= function(ylist,pc_list,previous_beta,timepoints1,timepoints2,basis1,basis2){


timei1 =timepoints1[[1]]
timei2 =timepoints2[[1]]


xmat = lapply(1:length(previous_beta), function(i){
	pc_fit  = pc_list[[i]]
	
	eval.fd2d.image(1:length(timei1),pc_fit)
	
	})%>%do.call(cbind,.)
	
	
cfits = lapply(1:length(ylist), function(x){
	
	index_pc = length(previous_beta)
	
	model = lm.fit(x=as.matrix(xmat[,1:index_pc]),y=as.matrix(ylist[[x]]))
	
	cfit= model%>%coef%>%as.numeric
	cfit
})%>%do.call(rbind,.)

list(sfit = cfits)

}














pred_soap_2d = function(ylist,pc_list,previous_beta,timepoints1,timepoints2,basis1,basis2){


timei1 =timepoints1[[1]]
timei2 =timepoints2[[1]]

cfits = lapply(1:length(ylist), function(x){
	
	xmat = lapply(1:length(previous_beta), function(i){
		pc_fit  = pc_list[[i]]
		
		eval.fd2d.image(1:length(timei1),pc_fit)
		
	})%>%do.call(cbind,.)
	
	index_pc = length(previous_beta)
	
	cfit = lm(ylist[[x]]~0+xmat[,1:index_pc])%>%coef%>%as.numeric
		
	cfit
})%>%do.call(rbind,.)

yfits = lapply(1:nrow(cfits), function(i){
	cfit = cfits[i,]
	
	yfit = lapply(1:length(pc_list),function(index){
		eval.fd2d.image(1:length(timei1),pc_list[[index]])*cfit[index]
	})%>%Reduce('+',.)
	yfit
}
)%>%do.call(cbind,.)


residuals = sapply(1:length(ylist), function(x){
	
	resid = ylist[[x]] - yfits[,x]
	if (all(resid==0)) {
		return(NA)
	} else {
	return(mean(resid^2))
	}
	
})
sigmahat = mean(residuals[!is.na(residuals)])

list(sigmahat = sigmahat, yfit=  yfits,sfit = cfits)

}

pred_2d<-function(sfits,pc_list,previous_beta,timepoints1,timepoints2,basis1,basis2){

	cfits = sfits 

	timei1 = timepoints1[[1]]

	cfit = cfits[1,]
	
	yfit = lapply(1:length(pc_list),function(index){
		eval.fd2d.image(1:length(timei1),pc_list[[index]])*cfit[index]
	})%>%Reduce('+',.)		
	
	return (yfit)

}

vectomatrix<-function(vec,picturesize=28){

	firstmat = matrix(0,ncol=picturesize,nrow=picturesize)	
	for (x in 1:picturesize){
	for (y in 1:picturesize){
		
		index= intersect(which(timepoints1[[1]]==x),which(timepoints2[[1]]==y))
		firstmat[x,y] = as.integer(vec[index])
		if (firstmat[x,y]<0){
			firstmat[x,y]=0
		}
		if (firstmat[x,y]>255){
			firstmat[x,y]=255
		}
	}
	}
return (firstmat)
}


pred_soap_2d.nonstandard = function(ylist,pc_list,previous_beta,timepoints1,timepoints2,basis1,basis2){

timei1 =timepoints1
timei2 =timepoints2

xmatlist  <- lapply(1:length(timei1),function(t){
	
	evalb1 = eval.basis(timei1[t], basis1)
	evalb2 = eval.basis(timei2[t], basis2)
	(t(evalb1) %*% evalb2)				
	}
)



cfits = lapply(1:length(ylist), function(x){
	
	xmat = lapply(1:length(previous_beta), function(i){
		pc_fit  = pc_list[[i]]
		
		eval.fd2d.image.nonstandard2(1:length(timei1),pc_fit,xmatlist)
		
	})%>%do.call(cbind,.)
	
	index_pc = length(previous_beta)
	
	cfit = lm(ylist[[x]]~0+xmat[,1:index_pc])%>%coef%>%as.numeric
		
	cfit
})%>%do.call(rbind,.)

yfits = lapply(1:nrow(cfits), function(i){
	cfit = cfits[i,]
	
	yfit = lapply(1:length(pc_list),function(index){
		eval.fd2d.image.nonstandard2(1:length(timei1),pc_list[[index]],xmatlist)*cfit[index]
	})%>%Reduce('+',.)		
	yfit
}
)%>%do.call(cbind,.)


residuals = sapply(1:length(ylist), function(x){
	
	resid = ylist[[x]] - yfits[,x]
	if (all(resid==0)) {
		return(NA)
	} else {
	return(mean(resid^2))
	}
	
})
sigmahat = mean(residuals[!is.na(residuals)])

list(sigmahat = sigmahat, yfit=  yfits,sfit = cfits)

}


