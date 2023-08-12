library(lubridate)
library(dplyr)
library(sp)
library(geosphere) # for lat-long distance in m
library(VineCopula) # for bivariate copula
library(survival)
library(pbv)
library(mvtnorm)
library(nloptr)
library(numDeriv)
library(zeallot)
library(statmod) 
library(distrEx)
library(distr)
# optional R packages
# library(scoringRules)
# library(riskRegression) # baseline hazard of Cox PH
# library(pbivnorm) # for bivariate normal CDF
# library(ggplot2)

# All 2^p combinations of 1/2
binary.combination <- function (p) {
    retval <- matrix(0, nrow = 2^p, ncol = p)
    for (n in 1:p) {
        retval[, n] <- rep(c(rep(1, (2^p/2^n)), rep(2, (2^p/2^n))), 
            length = 2^p)
    }
    retval
}

# Given a group of record, return matrix of pairs within a hard distance threshold (in km)
dist.binary <- function(longlat, row_ind, cutoff=50, fun=geosphere::distGeo) {
# convert to m
cutoff.m <- cutoff* 1e3

n <- NROW(longlat)

ans <- NULL

if (n >=2) {

for (i in 2:n) {
j <- 1:(i-1)
tmp <- data.frame(i, j, distance=fun(longlat[i,], longlat[j,])) %>% 
		filter(distance <= cutoff.m)

if (NROW(tmp) > 0) {ans <- c(ans, list(tmp))}
}

if (is.list(ans)) {
	tmp <- do.call(rbind,ans) 
}else {
	tmp <- ans
}

# return row index and distance (in km)
ans <- cbind(row_ind[tmp[,1]], row_ind[tmp[,2]], tmp[,3]/1e3)
}

ans
}

# Faster than v1
dist.binary.v2 <- function(longlat, row_ind, cutoff=50, fun=geosphere::distGeo) {

n <- NROW(longlat)

ans <- tmp <- NULL

if (n >=2) {
tmp <- data.frame(expand.grid(ind1 = 1:n, ind2 = 1:n)) %>% filter(ind1 > ind2)

longlat1 <- longlat[tmp$ind1,] ; longlat2 <- longlat[tmp$ind2,]

tmp <- tmp %>% mutate(distance=fun(longlat1, longlat2)/1e3) %>% filter(distance <= cutoff)
}

if (NROW(tmp) > 0) ans <- cbind(row_ind[tmp[,1]], row_ind[tmp[,2]], tmp[,3])

ans
}



# For GROUPED record in dplyr, return list of matrices for pairs within a hard distance threshold (chosen in km)
dist.group <- function(dat.tmp, cutoff=50, fun=geosphere::distGeo) {

dist.binary(dat.tmp %>% select(LON, LAT), dat.tmp$row_ind, 
	cutoff, fun=fun)
}

dist.group.v2 <- function(dat.tmp, cutoff=50, fun=geosphere::distGeo) {

dist.binary.v2(dat.tmp %>% select(LON, LAT), dat.tmp$row_ind, 
	cutoff, fun=fun)
}

# Given a group of record, return matrix of pairs of r nearest neighbor within a hard distance threshold (default: Inf) in km
dist.neighbor <- function(longlat, row_ind, r=10, fun=geosphere::distGeo) {


n <- NROW(longlat)
r <- min(r, n)

order.r <- function(x, r) {order(x)[1:r]}


ans.nn <- tmp <- tmp.matrix <- pairwise.distance <- NULL

if (n >=2) {
tmp <- data.frame(expand.grid(ind1 = 1:n, ind2 = 1:n)) %>% filter(ind1 > ind2) # only lower triangular element without diagonal

longlat1 <- longlat[tmp$ind1,] ; longlat2 <- longlat[tmp$ind2,]

tmp <- tmp %>% mutate(distance=fun(longlat1, longlat2)/1e3)

tmp.matrix <- matrix(0, nrow=n, ncol=n)
tmp.matrix[(tmp[,2]-1)*n + tmp[,1]] <- tmp[,3]
tmp.matrix <- tmp.matrix + t(tmp.matrix) + diag(Inf,n) # get full distance matrix and set self-distance as Inf


tmp.matrix <- cbind(1:n, t(apply(tmp.matrix, 1, order.r, r=r))) # get indices for r nearest neighbor

}

ans.nn <- tmp.matrix

if (NROW(tmp) > 0) {
for (j in 1:(r+1)) {
ans.nn[,j] <- row_ind[tmp.matrix[,j]] # convert back to row_ind in full dataset
}

pairwise.distance <- tmp
pairwise.distance <- pairwise.distance %>% mutate(ind1 = row_ind[ind1], ind2 = row_ind[ind2])

ans.nn <- data.frame(ans.nn) ; names(ans.nn) <- c("row_ind", paste0("nn",1:r))

}



list(nearest_neighbor=ans.nn, pairwise_distance=pairwise.distance)
}


dist.group.neighbor <- function(dat.tmp, r=10, fun=geosphere::distGeo) {
dist.neighbor(dat.tmp %>% select(LON, LAT), dat.tmp$row_ind, 
	r=r, fun=fun)
}


###### Function for calculating bivariate Gaussian copula #####
# given parameter (rho, psi) and distance

# input:
# u1, u2: vector of probability
# distance: vector of distance of pairs
# rho: reparameterization of factor [0,1]
# psi: practical range, > 0
# kap: nugget effect [0,1]

# C(u1, u2) : Copula cdf
Cop.cdf <- function(u1, u2, distance, rho, psi, kap) {

delta <- (1 - rho)*kap*exp(-3*distance/psi) + rho

pbv::pbv_rcpp_pbvnorm(x=qnorm(u1), y=qnorm(u2), rho=delta)
}

# c(u1, u2): Copula pdf
Cop.pdf <- function(u1, u2, distance, rho, psi, kap) {

x <- qnorm(u1) ; y <- qnorm(u2)
delta <- (1 - rho)*kap*exp(-3*distance/psi) + rho

pbv::pbv_rcpp_dbvnorm(x, y, rho=delta, use_log=FALSE)/
pbv::pbv_rcpp_dbvnorm(x, y, rho=0*delta, use_log=FALSE)
}


C.fun <- function(u1,u2,distance, para) {
tmp <- data.frame(u1=u1, u2=u2) %>% 
		mutate(	u1 = ifelse( 1 - u1 < 10*.Machine$double.eps, u1 - 10*.Machine$double.eps, u1),
				u2 = ifelse( 1 - u2 < 10*.Machine$double.eps, u2 - 10*.Machine$double.eps, u2)) %>%
		mutate(	u1 = ifelse( u1 < 10*.Machine$double.eps, 10*.Machine$double.eps, u1),
				u2 = ifelse( u2 < 10*.Machine$double.eps, 10*.Machine$double.eps, u2))
				
Cop.cdf(tmp$u1, tmp$u2, distance, para[1], para[2], para[3])
}


c.fun <- function(u1,u2, distance,para) {

tmp <- data.frame(u1=u1, u2=u2) %>% 
		mutate(	u1 = ifelse( 1 - u1 < 10*.Machine$double.eps, u1 - 10*.Machine$double.eps, u1),
				u2 = ifelse( 1 - u2 < 10*.Machine$double.eps, u2 - 10*.Machine$double.eps, u2)) %>%
		mutate(	u1 = ifelse( u1 < 10*.Machine$double.eps, 10*.Machine$double.eps, u1),
				u2 = ifelse( u2 < 10*.Machine$double.eps, 10*.Machine$double.eps, u2))
				
Cop.pdf(tmp$u1, tmp$u2, distance, para[1], para[2], para[3])
}


###### Gaussian copula (p-dimension version mainly for prediction)
# C(u1, u2,...up|R_matrix)

# By Quasi-Monte-Carlo procedure
Cop.cdf.general <- function(F_matrix, R_matrix, type, abseps=1e-20, maxpts=1e5, steps=4096) {
if (NCOL(F_matrix) == 1) F_matrix <- matrix(c(F_matrix), nrow=1)

U_matrix <- qnorm(F_matrix)

p.tmp.Bretz <- function(u, R_matrix=R_matrix, maxpts=maxpts, abseps=abseps) {
mvtnorm::pmvnorm(lower=rep(-Inf,length(u)), upper=u,
	sigma=R_matrix, 
	algorithm=mvtnorm::GenzBretz(maxpts=maxpts, abseps=abseps))
}

p.tmp.Miwa <- function(u, R_matrix=R_matrix, steps=steps) {
mvtnorm::pmvnorm(lower=rep(-Inf,length(u)), upper=u,
	sigma=R_matrix, 
	algorithm=mvtnorm::Miwa(steps=steps))
}

ans <- NA

if (type == "Miwa") {
	U_matrix[U_matrix == Inf] <- 1e20 # need to replace Inf to large number for numerical stability
	ans <- apply(U_matrix,1, p.tmp.Miwa, R_matrix=R_matrix,  steps=steps)
} else {
	ans <- apply(U_matrix,1, p.tmp.Bretz, R_matrix=R_matrix, abseps=abseps, maxpts=maxpts)
}

ans
}


# c(u1,...,up| R_matrix)
Cop.pdf.general <- function(F_matrix, R_matrix) {

if (NCOL(F_matrix) == 1) F_matrix <- matrix(c(F_matrix), nrow=1)

U_matrix <- qnorm(F_matrix)

eig_obj <- eigen(R_matrix)
eigen_value <- eig_obj$values
eigen_matrix <- eig_obj$vectors

det_R_minus_half <- 1/sqrt(prod(eigen_value))
R_inv_minus_I <- eigen_matrix%*%diag(1/eigen_value)%*%t(eigen_matrix) -
			diag(1, NCOL(F_matrix))

pdf.tmp <- function(q_u, det_R_minus_half, R_inv_minus_I) {
exp(-0.5*t(q_u)%*%R_inv_minus_I%*%q_u)
}
det_R_minus_half*apply(U_matrix,1, pdf.tmp , det_R_minus_half=det_R_minus_half, R_inv_minus_I=R_inv_minus_I)
}

# function to calculate the correlation matrix
R_calculation <- function(index, para_copula, pair_dist_tmp) {
c(rho, psi, kappa) %<-% para_copula
tmp <- length(index)
tau_matrix <- matrix(Inf, nrow=tmp, ncol=tmp)

for (i in 1:(tmp-1)) {
	for (j in (i+1):tmp) {
		tau_matrix[i,j] <- as.numeric(pair_dist_tmp %>% 
				filter(ind1 == index[i] & ind2 == index[j]) %>% 
				select(distance))[1]
	}
}

tau_matrix <- kappa*exp(-3*tau_matrix/psi)
B_matrix <- tau_matrix + t(tau_matrix) + diag(1, tmp)

(1 - rho)*B_matrix + rho
}


# function to drop one particular event
drop.one.event <- function(dat, pair.dat, event_index) {
	
	dat <- dat %>% ungroup()
	unique_event_ind <- unique(dat$event_ind) # get unique event index
	
	event_contribute <- data.frame(drop_event_ind = event_index, 
			drop_evt_date=as.Date("1900-12-03", "%Y-%m-%d"), drop_COUNTY=NA,
			drop_stage1=NA, drop_stage2=NA)
			
	dat_drop_one_event 	<- 	dat %>% filter(event_ind != event_index) %>% 
		mutate(row_ind_new = row_number())
		
	dat_tmp <- dat_drop_one_event %>% select(row_ind, row_ind_new)
		
	pair.dat_drop_one_event <- pair.dat %>% filter(ind1 %in% dat_drop_one_event$row_ind)
		
	tmp1 <- pair.dat_drop_one_event %>% select(ind1) ; names(tmp1) <- "row_ind"
	suppressMessages(tmp1 <- tmp1 %>% left_join(., dat_tmp))
		
	tmp2 <- pair.dat_drop_one_event %>% select(ind2) ; names(tmp2) <- "row_ind"
	suppressMessages(tmp2 <- tmp2 %>% left_join(., dat_tmp))
		
		
	pair.dat_drop_one_event <- pair.dat_drop_one_event %>% 
				mutate(ind1 = tmp1$row_ind_new, ind2 = tmp2$row_ind_new)
		
	event_contribute[1,c(1,4,5)] <- c(event_index,
			c(NROW(dat) - NROW(dat_drop_one_event), 
			NROW(pair.dat) - NROW(pair.dat_drop_one_event)))
		
	event_contribute[1,2:3] <- as.data.frame(dat %>% 
			filter(event_ind == event_index) %>% select(evt_date, N_COUNTY) %>% head(1))
			
		
list(dat_drop_one_event=dat_drop_one_event, 
	pair.dat_drop_one_event=pair.dat_drop_one_event,
	event_contribute=event_contribute)
}

# function to filter single event
single.event <- function(dat, pair.dat, event_index) {
	
	dat <- dat %>% ungroup()
	unique_event_ind <- unique(dat$event_ind) # get unique event index
	
	dat_single_event 	<- 	dat %>% filter(event_ind == event_index) %>% 
		mutate(row_ind_new = row_number())
		
	dat_tmp <- dat_single_event %>% select(row_ind, row_ind_new)
		
	pair.dat_single_event <- pair.dat %>% filter(ind1 %in% dat_single_event$row_ind)
		
	tmp1 <- pair.dat_single_event %>% select(ind1) ; names(tmp1) <- "row_ind"
	suppressMessages(tmp1 <- tmp1 %>% left_join(., dat_tmp))
		
	tmp2 <- pair.dat_single_event %>% select(ind2) ; names(tmp2) <- "row_ind"
	suppressMessages(tmp2 <- tmp2 %>% left_join(., dat_tmp))
		
		
	pair.dat_single_event <- pair.dat_single_event %>% 
				mutate(ind1 = tmp1$row_ind_new, ind2 = tmp2$row_ind_new)
		
	event_contribute[1,c(1,4,5)] <- c(event_index,
			c(NROW(dat) - NROW(dat_single_event), 
			NROW(pair.dat) - NROW(pair.dat_single_event)))
		
	event_contribute[1,2:3] <- as.data.frame(dat %>% 
			filter(event_ind == event_index) %>% select(evt_date, N_COUNTY) %>% head(1))
			
		
list(dat_single_event=dat_single_event, 
	pair.dat_single_event=pair.dat_single_event
	)
}




# function to drop by a logic vector for rows
drop.by.vector <- function(dat, pair.dat, logic.vector, drop.label="drop") {
	
	if (is.data.frame(logic.vector)) {
		logic.vector <- unname(logic.vector %>% unlist())
	}
	
	dat <- dat %>% ungroup() %>% mutate(tmp_group = logic.vector)
	
		
	dat_drop 	<- 	dat %>% filter(!(tmp_group %in% drop.label)) %>% 
		mutate(row_ind_new = row_number())
		
	dat_tmp <- dat_drop %>% select(row_ind, row_ind_new)
		
	pair.dat_drop <- pair.dat %>% filter(ind1 %in% dat_drop$row_ind)
		
	tmp1 <- pair.dat_drop %>% select(ind1) ; names(tmp1) <- "row_ind"
	suppressMessages(tmp1 <- tmp1 %>% left_join(., dat_tmp))
		
	tmp2 <- pair.dat_drop %>% select(ind2) ; names(tmp2) <- "row_ind"
	suppressMessages(tmp2 <- tmp2 %>% left_join(., dat_tmp))
		
		
	pair.dat_drop <- pair.dat_drop %>% 
				mutate(ind1 = tmp1$row_ind_new, ind2 = tmp2$row_ind_new)
	
	event_contribute <- c(drop_stage1=c(NROW(dat) - NROW(dat_drop), 
		drop_stage2=NROW(pair.dat) - NROW(pair.dat_drop)))


list(dat_drop=dat_drop, 
	pair.dat_drop=pair.dat_drop,
	event_contribute=event_contribute)
}


# Functions for updating the density k=3

# f(z_{ij}|t_{ij}): just use the bivariate copula
new_f3 <- function(z_tmp, F_T, nu, mu, cop) {
	
	n_tmp <- length(z_tmp)
	f_Z <- dgamma(z_tmp, shape = nu, scale= mu/nu)
	F_Z <- pgamma(z_tmp, shape = nu, scale= mu/nu)
	F_T <- rep(F_T, n_tmp)
	h_Z <- VineCopula::BiCopHfunc1(F_T, F_Z, cop)
	h 	<- VineCopula::BiCopPDF(F_T, F_Z, cop)
		
	f3 <- exp(log(f_Z) + log(h))
	ifelse(f3 < 0 | is.nan(f3), .Machine$double.xmin, f3)
}

# z_{ij}*f(z_{ij}|t_{ij})
Z_times_new_f3 <- function(z_tmp, F_T, nu, mu, cop)  {
	z_tmp*new_f3(z_tmp, F_T, nu, mu, cop)
}

# f(z_{ij} | t_{ij}, {z_ij'}s): use both bivariate and Gaussian copula
new_f_Z <- function(z_tmp, F_T, nu, mu, cop, F_matrix, R_matrix, r1) {
	n_tmp <- length(z_tmp)
	f_Z <- dgamma(z_tmp, shape = nu, scale= mu/nu)
	F_Z <- pgamma(z_tmp, shape = nu, scale= mu/nu)
	F_T <- rep(F_T, n_tmp)
	h_Z <- VineCopula::BiCopHfunc1(F_T, F_Z, cop)
	h 	<- VineCopula::BiCopPDF(F_T, F_Z, cop)
	
	F3 <- h_Z; f3 <- exp(log(f_Z) + log(h))
		
	F3 <- ifelse(F3 <= 0 | is.nan(F3), .Machine$double.xmin, F3)
	f3 <- ifelse(f3 < 0 | is.nan(f3), .Machine$double.xmin, f3)
		
	F_matrix <-  matrix(F_matrix, nrow=n_tmp, ncol=length(F_matrix), byrow=TRUE)
	F_matrix[,1] <- F3
		
	numerator <- Cop.pdf.general(F_matrix, R_matrix)
		
	if (r1 > 0) { # at least one neighbor in denominator
		denominator <- Cop.pdf.general(F_matrix[1,-1], R_matrix[-1,-1])
	} else { # only one neighbor in denominator c = dC(u)/du = du/du = 1
		denominator <- 1
	}
	
	f3*numerator/denominator
}

# z_{ij}*f(z_{ij} | t_{ij}, {z_ij'}s)	
Z_times_new_f_Z <- function(z_tmp, F_T, nu, mu, cop, F_matrix, R_matrix, r1)  {
	z_tmp*new_f_Z(z_tmp, F_T, nu, mu, cop, F_matrix, R_matrix, r1)
}



# Two stage estimation for k = 1
# dat		: raw dataset contains all variables (WITH index of event group and row index)
# pair.dat	: matrix of incies of pair and distance
# form		: formula type for the first stage GLM fitting
# grad		: logical indicator if  require APPROXIMATE gradient for two stages
# hessian	: locgical indicator if require hessian matrix of two stages
# ini		: initial value of stage 2
# tol.bound	: Small bound for the parameter space
# factr		: relative improment for declaring convergence for L-BFGS-B (default: 1e8)
# progress	: logical indicator if require printing progress

est1.twostage <- function(dat, pair.dat, form, grad=TRUE, hessian=FALSE,
	ini=c(0.5, 1, 0.3), tol.bound=1e-3, factr=1e7, progress=TRUE) {
	
F.ind <- function(ind) {F0[ind]}
y.ind <- function(ind,m=m.glm) {unname(m$y[ind])}


# Stage 1
m.glm <- glm(form, family=binomial(link="logit"), data=dat, 
			x = TRUE, y = TRUE, model=FALSE)
F0 <- 1 - fitted(m.glm)

beta.est <- coef(m.glm)

# stage 1 log-likelihood as a data.frame
loglik.stage1 <- function(par.glm, X=X, y=y) {
upper.b <-  log(.Machine$double.xmax)
lower.b <- -upper.b

sum(data.frame(X.glm = X%*%par.glm) %>% 
		mutate(X.glm = case_when(
		X.glm > upper.b ~ upper.b,
		X.glm < lower.b ~ lower.b,
		TRUE ~ X.glm
		)) %>% mutate(eXb = exp(X.glm)) %>% mutate(F0 =1/(1+eXb)
		, loglik = y*X.glm  - log(1 + eXb)) %>% select(loglik))
}


if (progress) cat("Stage 1 done, coef of glm:", beta.est,"\n")

pair.dat <- pair.dat %>% mutate(
		F01 = F.ind(ind1), F02 = F.ind(ind2), 
		y1=y.ind(ind1), y2=y.ind(ind2)) %>%
		mutate(y12 = case_when(
		y1 == 0 & y2 == 1 ~ 1,			# for 01
		y1 == 1 & y2 == 0 ~ 2,			# for 10
		y1 == 1 & y2 == 1 ~ 3, 			# for 11
		TRUE ~ 0 						# for 00
		)) %>% select(-c(y1,y2))


# function for second state
logCL <- function(para, dat=pair.dat, progress) {


ans <- sum(dat %>% mutate(C1 = C.fun(F01, F02, distance, para=para)) %>%
		mutate(f1 = case_when(
		y12 == 1 ~ F01 - C1,					# for 01
		y12 == 2 ~ F02 - C1,					# for 10
		y12 == 3 ~ 1 - F01 - F02 + C1, 			# for 11
		TRUE ~ C1 								# for 00
		)) %>% 
		mutate(f1 = ifelse(
		f1 <= 0, .Machine$double.xmin,f1)) %>%
		mutate(logf = log(f1)) %>% select(logf)) 


if (progress) cat("logCL:", ans, "; para:", para,"\n")

ans
}

# Stage 2
obj <- optim(par=ini, fn=logCL, method="L-BFGS-B",
			dat=pair.dat, progress=progress,
			lower=rep(tol.bound, 3), 
			upper=c(1 - tol.bound, Inf, 1 - tol.bound),
			control=list(fnscale=-1, factr=factr))


est <- c(beta.est,obj$par)
names(est) <- c(names(beta.est), "rho", "psi", "kappa")

result <- list(optim.obj=obj, est=est)

if (grad) {

if (progress) cat("Estiamtion done, calculating numerical gradient\n")

result$score <- c(
			numDeriv::grad(loglik.stage1, beta.est, 
				X=model.matrix(m.glm), y=m.glm$y),
			numDeriv::grad(logCL, obj$par, dat=pair.dat, progress=FALSE)	
			)
names(result$score) <- c(names(beta.est), "rho", "psi", "kappa")

}


if (hessian) {

if (progress) cat("Estiamtion done, calculating numerical hessian\n")

result$hessian <- list(stage1 = numDeriv::hessian(loglik.stage1, beta.est, 
				X=model.matrix(m.glm), y=m.glm$y), # -solve(vcov(m.glm))
				stage2 = numDeriv::hessian(logCL, obj$par, dat=pair.dat, progress=FALSE)	
			)
}

result
}


# Two stage estimation for k = 1 without A0
# dat		: raw dataset contains all variables (WITH index of event group and row index)
# pair.dat	: matrix of incies of pair and distance
# form		: formula type for the first stage GLM fitting
# grad		: logical indicator if  require APPROXIMATE gradient for two stages
# hessian	: locgical indicator if require hessian matrix of two stages
# ini		: initial value of stage 2 (Note: rho is fixed at 0) only psi, kappi here
# tol.bound	: Small bound for the parameter space
# factr		: relative improment for declaring convergence for L-BFGS-B (default: 1e8)
# progress	: logical indicator if require printing progress
est1.dropA0.twostage <- function(dat, pair.dat, form, grad=TRUE, hessian=FALSE,
	ini=c(1, 0.3), tol.bound=1e-3, factr=1e7, progress=TRUE) {
	
F.ind <- function(ind) {F0[ind]}
y.ind <- function(ind,m=m.glm) {unname(m$y[ind])}


# Stage 1
m.glm <- glm(form, family=binomial(link="logit"), data=dat, 
			x = TRUE, y = TRUE, model=FALSE)
F0 <- 1 - fitted(m.glm)

beta.est <- coef(m.glm)

# stage 1 log-likelihood as a data.frame
loglik.stage1 <- function(par.glm, X=X, y=y) {
upper.b <-  log(.Machine$double.xmax)
lower.b <- -upper.b

sum(data.frame(X.glm = X%*%par.glm) %>% 
		mutate(X.glm = case_when(
		X.glm > upper.b ~ upper.b,
		X.glm < lower.b ~ lower.b,
		TRUE ~ X.glm
		)) %>% mutate(eXb = exp(X.glm)) %>% mutate(F0 =1/(1+eXb)
		, loglik = y*X.glm  - log(1 + eXb)) %>% select(loglik))
}


if (progress) cat("Stage 1 done, coef of glm:", beta.est,"\n")

pair.dat <- pair.dat %>% mutate(
		F01 = F.ind(ind1), F02 = F.ind(ind2), 
		y1=y.ind(ind1), y2=y.ind(ind2)) %>%
		mutate(y12 = case_when(
		y1 == 0 & y2 == 1 ~ 1,			# for 01
		y1 == 1 & y2 == 0 ~ 2,			# for 10
		y1 == 1 & y2 == 1 ~ 3, 			# for 11
		TRUE ~ 0 						# for 00
		)) %>% select(-c(y1,y2))


# function for second state
logCL <- function(para, dat=pair.dat, progress) {


ans <- sum(dat %>% mutate(C1 = C.fun(F01, F02, distance, para=para)) %>%
		mutate(f1 = case_when(
		y12 == 1 ~ F01 - C1,					# for 01
		y12 == 2 ~ F02 - C1,					# for 10
		y12 == 3 ~ 1 - F01 - F02 + C1, 			# for 11
		TRUE ~ C1 								# for 00
		)) %>% 
		mutate(f1 = ifelse(
		f1 <= 0, .Machine$double.xmin,f1)) %>%
		mutate(logf = log(f1)) %>% select(logf)) 


if (progress) cat("logCL:", ans, "; para:", para,"\n")

ans
}

# Stage 2
obj <- optim(par=ini, fn=logCL, method="L-BFGS-B",
			dat=pair.dat, progress=progress,
			lower=rep(tol.bound, 2), 
			upper=c(Inf, 1 - tol.bound),
			control=list(fnscale=-1, factr=factr))


est <- c(beta.est,0,obj$par)
names(est) <- c(names(beta.est), "rho", "psi", "kappa")

result <- list(optim.obj=obj, est=est)

if (grad) {

if (progress) cat("Estiamtion done, calculating numerical gradient\n")

result$score <- c(
			numDeriv::grad(loglik.stage1, beta.est, 
				X=model.matrix(m.glm), y=m.glm$y),
			0,
			numDeriv::grad(logCL, obj$par, dat=pair.dat, progress=FALSE)	
			)
names(result$score) <- c(names(beta.est), "rho", "psi", "kappa")

}


if (hessian) {

if (progress) cat("Estiamtion done, calculating numerical hessian\n")

result$hessian <- list(stage1 = numDeriv::hessian(loglik.stage1, beta.est, 
				X=model.matrix(m.glm), y=m.glm$y), # -solve(vcov(m.glm))
				stage2 = numDeriv::hessian(logCL, obj$par, dat=pair.dat, progress=FALSE)	
			)
}

result
}


# Two stage estimation for k = 1 at event level at once (example only, probably should be parallelized)
# dat		: raw dataset contains all variables (WITH index of event group and row index)
# pair.dat	: matrix of incies of pair and distance
# form		: formula type for the first stage GLM fitting
# grad		: logical indicator if  require APPROXIMATE gradient for two stages
# hessian	: locgical indicator if require hessian matrix of two stages
# ini		: initial value of stage 2 (probably from est1.twostage())
# tol.bound	: Small bound for the parameter space
# factr		: relative improment for declaring convergence for L-BFGS-B (default: 1e8)
# progress	: logical indicator if require printing progress in term of event
# save.name	: save file as the assigned name at working directory if not NULL
est1.twostage.event <- function(dat, pair.dat, form, grad=TRUE, hessian=FALSE,
	ini=c(0.4, 2, 0.22), tol.bound=1e-3, factr=1e7, progress=TRUE, save.name=NULL) {
	
	dat <- dat %>% ungroup()
	unique_event_ind <- unique(dat$event_ind) # get unique event index
	
	event_contribute <- data.frame(drop_event_ind = unique_event_ind, 
			drop_evt_date=as.Date("1900-12-03", "%Y-%m-%d"), drop_COUNTY=NA,
			drop_stage1=NA, drop_stage2=NA)
	
	cnt <- 1
	
	result <- NULL
	
	for (i in 1:length(unique_event_ind)) {
		# prepare drop one event dataset and pair.dat and reindexing
	
		c(dat_drop_one_event,
		pair.dat_drop_one_event,
		event_contribute) %<-% drop.one.event(dat.binary, pair.dat, event_index=i)
		
		
		event_contribute[i,c(1,4,5)] <- c(unique.event[i],
				c(NROW(dat) - NROW(dat_drop_one_event), 
				NROW(pair.dat) - NROW(pair.dat_drop_one_event)))
		
		event_contribute[i,2:3] <- as.data.frame(dat %>% filter(event_ind == unique.event[i]) %>% 
				select(evt_date, N_COUNTY) %>% head(1))
		
		# checking
		#tmp <- cbind(dat_drop_one_event[pair.dat_drop_one_event$ind1,c("evt_date", "N_COUNTY")],
		#		dat_drop_one_event[pair.dat_drop_one_event$ind2,c("evt_date", "N_COUNTY")])
		# any(tmp[,1] != tmp[,3]) ; any(tmp[,2] != tmp[,4]) # should be FALSE FALSE
				
		# 1e-3*geosphere::distGeo(
		# as.numeric(dat_drop_one_event[pair.dat_drop_one_event$ind1[1],c("LON", "LAT")]),
		# as.numeric(dat_drop_one_event[pair.dat_drop_one_event$ind2[1],c("LON", "LAT")]))
		# should match with pair.dat_drop_one_event[1,"distance"]
		
		# Core run after dropping one event
		result[[i]] <- est1.twostage(dat_drop_one_event, pair.dat_drop_one_event, form, grad=grad, hessian=hessian,
				ini=ini, tol.bound=tol.bound, factr=factr, progress=progress)
		
		if (progress) {	cat(i, "out of", length(unique_event_ind),"!\n")	}
		if (!is.null(save.name)) {
			ans <- list(result=result, event_contribute=event_contribute) 
			save(ans, file=save.name)
		}
	}
	
list(result=result, event_contribute=event_contribute)
}



# Full estimation log-likelihood when k = 1
# pair.dat	: matrix of incies of pair and distance
# form		: formula type for the first stage GLM fitting
# grad		: logical indicator if  require APPROXIMATE gradient for two stages
# hessian	: locgical indicator if require hessian matrix of two stages
# ini		: initial value for all parameter (probably from two-stage approach)
# tol.bound	: Small bound for the parameter space
# factr		: relative improment for declaring convergence for L-BFGS-B (default: 1e8)
# progress	: logical indicator if require printing progress
	
est1.full <- function(dat, pair.dat, form, grad=TRUE, hessian=FALSE,
	ini, tol.bound=1e-3, factr=1e7, progress=TRUE) {
	
	
m.glm <- glm(form, family=binomial(link="logit"), data=dat, 
			x = TRUE, y = TRUE, model=FALSE)

X <- model.matrix(m.glm)
y <- m.glm$y

upper.b <-  log(.Machine$double.xmax)
lower.b <- -upper.b

y.ind <- function(ind) {unname(y[ind])}

pair.dat <- pair.dat %>% mutate( 
		y1=y.ind(ind1), y2=y.ind(ind2)) %>%
		mutate(y12 = case_when(
		y1 == 0 & y2 == 1 ~ 1,			# for 01
		y1 == 1 & y2 == 0 ~ 2,			# for 10
		y1 == 1 & y2 == 1 ~ 3, 			# for 11
		TRUE ~ 0 						# for 00
		)) %>% select(-c(y1,y2))
		
# stage 1 log-likelihood as a data.frame
loglik.stage1 <- function(par.glm, X=X, y=y) {
upper.b <-  log(.Machine$double.xmax)
lower.b <- -upper.b

df.tmp <- data.frame(X.glm = X%*%par.glm) %>% 
		mutate(X.glm = case_when(
		X.glm > upper.b ~ upper.b,
		X.glm < lower.b ~ lower.b,
		TRUE ~ X.glm
		)) %>% mutate(eXb = exp(X.glm)) %>% mutate(F0 =1/(1+eXb)
		, loglik = y*X.glm  - log(1 + eXb)) %>% select(loglik, F0)
		
		
list(loglik= sum(df.tmp$loglik), F0 = df.tmp$F0)
}

# stage 2 pairwise likelihood given par.glm
logCL <- function(par.copula, pair.dat=pair.dat, F0=F0) {

F.ind <- function(ind) {F0[ind]}

sum(pair.dat %>% mutate(F01 = F.ind(ind1), F02 = F.ind(ind2)) %>%
		mutate(C1 = C.fun(F01, F02, distance, para=par.copula)) %>%
		mutate(f1 = case_when(
		y12 == 1 ~ F01 - C1,					# for 01
		y12 == 2 ~ F02 - C1,					# for 10
		y12 == 3 ~ 1 - F01 - F02 + C1, 			# for 11
		TRUE ~ C1 								# for 00
		)) %>% 
		mutate(f1 = ifelse(
		f1 <= 0, .Machine$double.xmin,f1)) %>%
		mutate(logf = log(f1)) %>% select(logf))
}

para <- ini

loglik.full <- function(para, X, y, pair.dat, progress) {

len.para <- length(para)

par.glm <- para[1:(len.para - 3)]
par.copula <- para[(len.para - 2):len.para]

# calculate marginal model
stage1 <- loglik.stage1(par.glm, X=X, y=y) 

# combine Copula model given par.glm (and hence F_1)
ans <- stage1$loglik + logCL(par.copula, pair.dat=pair.dat, F0=stage1$F0)

if (progress) cat("logLik:", ans,"; para:", c(round(para,3)),"\n")

ans
}

# trick to make sure that regression coefficient won't cause overflow in stage 1
beta.lb <- -log(.Machine$double.xmax)/apply(abs(X),2,max)
beta.ub <- -beta.lb

obj <- optim(par=ini, fn=loglik.full, X=X, y=y, pair.dat=pair.dat, progress=progress,
		method="L-BFGS-B",
		lower=c(beta.lb, rep(tol.bound, 3)), 
		upper=c(beta.ub, 1 - tol.bound, Inf, 1 - tol.bound),
		control=list(fnscale=-1, factr=factr))


est <- c(obj$par)

names(est) <- c(names(coef(m.glm)), "rho", "psi", "kappa")

result <- list(obj=obj, est=est)

if (grad) {

if (progress) cat("Estimation done, calculating numerical gradient\n")

result$score <- numDeriv::grad(loglik.full, est, X=X, y=y, pair.dat=pair.dat, progress=FALSE)
}

if (hessian) {

if (progress) cat("Estimation done, calculating numerical hessian\n")

result$hessian <- numDeriv::hessian(loglik.full, est, X=X, y=y, pair.dat=pair.dat, progress=FALSE)
}

result
}


# Full estimation weighted log-likelihood when k = 1
# pair.dat	: matrix of incies of pair and distance
# form		: formula type for the first stage GLM fitting
# grad		: logical indicator if  require APPROXIMATE gradient for two stages
# hessian	: locgical indicator if require hessian matrix of two stages
# ini		: initial value for all parameter (probably from two-stage approach)
# tol.bound	: Small bound for the parameter space
# factr		: relative improment for declaring convergence for L-BFGS-B (default: 1e8)
# progress	: logical indicator if require printing progress
	
est1.full.weight <- function(dat, pair.dat, form, grad=TRUE, hessian=FALSE,
	ini, tol.bound=1e-3, factr=1e7, progress=TRUE) {

W1 <- NROW(dat) ; W2 <- NROW(pair.dat)	

m.glm <- glm(form, family=binomial(link="logit"), data=dat, 
			x = TRUE, y = TRUE, model=FALSE)

X <- model.matrix(m.glm)
y <- m.glm$y

upper.b <-  log(.Machine$double.xmax)
lower.b <- -upper.b

y.ind <- function(ind) {unname(y[ind])}

pair.dat <- pair.dat %>% mutate( 
		y1=y.ind(ind1), y2=y.ind(ind2)) %>%
		mutate(y12 = case_when(
		y1 == 0 & y2 == 1 ~ 1,			# for 01
		y1 == 1 & y2 == 0 ~ 2,			# for 10
		y1 == 1 & y2 == 1 ~ 3, 			# for 11
		TRUE ~ 0 						# for 00
		)) %>% select(-c(y1,y2))
		
# stage 1 log-likelihood as a data.frame
loglik.stage1 <- function(par.glm, X=X, y=y) {
upper.b <-  log(.Machine$double.xmax)
lower.b <- -upper.b

df.tmp <- data.frame(X.glm = X%*%par.glm) %>% 
		mutate(X.glm = case_when(
		X.glm > upper.b ~ upper.b,
		X.glm < lower.b ~ lower.b,
		TRUE ~ X.glm
		)) %>% mutate(eXb = exp(X.glm)) %>% mutate(F0 =1/(1+eXb)
		, loglik = y*X.glm  - log(1 + eXb)) %>% select(loglik, F0)
		
		
list(loglik= sum(df.tmp$loglik), F0 = df.tmp$F0)
}

# stage 2 pairwise likelihood given par.glm
logCL <- function(par.copula, pair.dat=pair.dat, F0=F0) {

F.ind <- function(ind) {F0[ind]}

sum(pair.dat %>% mutate(F01 = F.ind(ind1), F02 = F.ind(ind2)) %>%
		mutate(C1 = C.fun(F01, F02, distance, para=par.copula)) %>%
		mutate(f1 = case_when(
		y12 == 1 ~ F01 - C1,					# for 01
		y12 == 2 ~ F02 - C1,					# for 10
		y12 == 3 ~ 1 - F01 - F02 + C1, 			# for 11
		TRUE ~ C1 								# for 00
		)) %>% 
		mutate(f1 = ifelse(
		f1 <= 0, .Machine$double.xmin,f1)) %>%
		mutate(logf = log(f1)) %>% select(logf))
}

para <- ini

loglik.full.weight <- function(para, X, y, pair.dat, W1, W2, progress) {

len.para <- length(para)

par.glm <- para[1:(len.para - 3)]
par.copula <- para[(len.para - 2):len.para]

# calculate marginal model
stage1 <- loglik.stage1(par.glm, X=X, y=y) 

# combine Copula model given par.glm (and hence F_1)
ans <- stage1$loglik/W1 + logCL(par.copula, pair.dat=pair.dat, F0=stage1$F0)/W2

if (progress) cat("weighted logLik:", ans,"; para:", c(round(para,3)),"\n")

ans
}

# trick to make sure that regression coefficient won't cause overflow in stage 1
beta.lb <- -log(.Machine$double.xmax)/apply(abs(X),2,max)
beta.ub <- -beta.lb

obj <- optim(par=ini, fn=loglik.full.weight, X=X, y=y, pair.dat=pair.dat, W1=W1, W2=W2, progress=progress,
		method="L-BFGS-B",
		lower=c(beta.lb, rep(tol.bound, 3)), 
		upper=c(beta.ub, 1 - tol.bound, Inf, 1 - tol.bound),
		control=list(fnscale=-1, factr=factr))


est <- c(obj$par)

names(est) <- c(names(coef(m.glm)), "rho", "psi", "kappa")

result <- list(obj=obj, est=est)

if (grad) {

if (progress) cat("Estimation done, calculating numerical gradient\n")

result$score <- numDeriv::grad(loglik.full.weight, est, X=X, y=y, pair.dat=pair.dat, W1=W1, W2=W2, progress=FALSE)
}

if (hessian) {

if (progress) cat("Estimation done, calculating numerical hessian\n")

result$hessian <- numDeriv::hessian(loglik.full.weight, est, X=X, y=y, pair.dat=pair.dat, W1=W1, W2=W2, progress=FALSE)
}

result
}


# Predictive probability for k = 1 based on (at most) r nearest neighbor only with claim = 1
# for SINGLE hail event
pred1.nn.single.event <- function(X_event, evt_ind, y_event, 
		para_est, longlat, r_max=5, type=c("Miwa", "Bretz"), abseps=1e-20, maxpts=1e5, steps=4096,
		fun=geosphere::distGeo,  progress=TRUE) {

if (r_max >= 10 && type == "Miwa") {
cat("Warning: Gaussian Copula CDF is numerical unstable and slow with Miwa when order (r_max) is large, switch to Bretz approximation!\n")

type <- "Bretz"
}


# order.r <- function(x, r) {order(x)[1:r]}

pair_to_matrix <- function(pairwise_distance, n) {
pairwise_matrix <- matrix(0, nrow=n, ncol=n)
pairwise_matrix[(pairwise_distance$ind2 - 1)*n + pairwise_distance$ind1] <- pairwise_distance$distance
pairwise_matrix + t(pairwise_matrix) + diag(Inf,n)
}

F_combination <- function(F_tmp, all_comb) {
result <- all_comb
for (j in 1:NCOL(result)) {
F_TRUE_FALSE <- c(F_tmp[j], 1)

result[,j] <- F_TRUE_FALSE[all_comb[,j]]
}

result
}



n <- NROW(X_event)
n_claim <- sum(y_event)

claim_ind <- (1:n)[y_event == 1]

if (length(unique(evt_ind)) != 1) { 
	cat("Not single event data is supported!")
	return(rep(NA, n))
}

para_beta <- para_est[1:(length(para_est)-3)]
para_copula <- tail(para_est,3)
c(rho, psi, kappa) %<-% para_copula

# Given the model, calculate marginal F1 = 1 - P(y_kj|x_j) = 1/(1+exp(X_j beta)
result <- F1 <- 1/( 1 + exp(c(X_event%*%para_beta)))


if (r_max == 0) {
cat("r_max equals 0, just use marginal probability!\n")
return (1- F1)
}

if (n == 1 || n_claim == 0) {
result <- 1 - F1 # 3 cases here: only one data point or no claim or (one claim + y = 1)
} else {

# another case with multiple records in event


# Case 0: When y is 0, can use all 1 neighbor
r0 <- min(r_max, n_claim)

all_comb_y0 <- binary.combination(r0 + 1)
sign_C_y0 <- (-1)^apply(all_comb_y0,1,sum)


# Case 1: When y = 1
r1 <- min(r_max, n_claim - 1) 

# Case 1a: at least two 1
if (r1 > 0) {
all_comb_y1 <- binary.combination(r1 + 1)
sign_C_y1 <- (-1)^apply(all_comb_y1,1,sum)
}


# loop here, maybe able to make it faster by lapply()?
for (i in 1:n) {


if (n_claim  > 1) {

	if (y_event[i] == 0) {
	# all pair about claim = 1 and the index
	pairwise_distance <- data.frame(expand.grid(ind1 = c(i,claim_ind), ind2 = claim_ind)) %>% filter (ind1 != ind2) %>%
		mutate(distance=fun(longlat[ind1,], longlat[ind2,])/1e3) 
		
	ind_tmp <- c(i, unname(unlist(
				pairwise_distance %>% filter(ind1 == i) %>% slice_min(distance, n=r0, with_ties=FALSE) %>% select(ind2)))
				)
	
	F_tmp <- F1[ind_tmp]
	F_matrix <- F_combination(F_tmp, all_comb_y0)
	
	R_matrix <- R_calculation(ind_tmp, para_copula, pairwise_distance)
	numerator <- sum(Cop.cdf.general(F_matrix, R_matrix=R_matrix, 
					type=type, abseps=abseps, maxpts=maxpts, steps=steps)*sign_C_y0)
	
	if (r0 > 1) {
	denominator <- sum(Cop.cdf.general(F_matrix[-c(1:2^r0),-1], R_matrix=R_matrix[-1,-1], 
						type=type, abseps=abseps, maxpts=maxpts, steps=steps)*sign_C_y0[-c(1:2^r0)])
	} else { denominator <- 1 - F_tmp[2] }
	result[i] <- numerator/denominator
	
	
	
	} else {
	pairwise_distance <- data.frame(expand.grid(ind1 = claim_ind, ind2 = claim_ind)) %>% filter (ind1 != ind2) %>%
		mutate(distance=fun(longlat[ind1,], longlat[ind2,])/1e3)
	
	ind_tmp <- c(i, unname(unlist(
				pairwise_distance %>% filter(ind1 == i) %>% slice_min(distance, n=r1, with_ties=FALSE) %>% select(ind2)))
				)
				
	F_tmp <- F1[ind_tmp]
	F_matrix <- F_combination(F_tmp, all_comb_y1)


	R_matrix <- R_calculation(ind_tmp, para_copula, pairwise_distance)
	numerator <- sum(Cop.cdf.general(F_matrix, R_matrix=R_matrix, 
				type=type, abseps=abseps, maxpts=maxpts, steps=steps)*sign_C_y1)
	
	if (r1 > 1) {
	denominator <- sum(Cop.cdf.general(F_matrix[-c(1:2^r1),-1], R_matrix=R_matrix[-1,-1], 
				type=type, abseps=abseps, maxpts=maxpts, steps=steps)*sign_C_y1[-c(1:2^r1)])
	} else {denominator <- 1 - F_tmp[2] }
				
	result[i] <- numerator/denominator
	
	}

} else { # i.e. only one claim
	if (y_event[i] == 0) {
	# all pair about claim = 1 and the index
	pairwise_distance <- data.frame(expand.grid(ind1 = c(i,claim_ind), ind2 = claim_ind)) %>% filter (ind1 != ind2) %>%
		mutate(distance=fun(longlat[ind1,], longlat[ind2,])/1e3) 
		
	ind_tmp <- c(i, unname(unlist(
				pairwise_distance %>% filter(ind1 == i) %>% slice_min(distance, n=r0, with_ties=FALSE) %>% select(ind2)))
				)
	
	F_tmp <- F1[ind_tmp]
	F_matrix <- F_combination(F_tmp, all_comb_y0)
	
	R_matrix <- R_calculation(ind_tmp, para_copula, pairwise_distance)
	numerator <- sum(Cop.cdf.general(F_matrix, R_matrix=R_matrix, 
				type=type, abseps=abseps, maxpts=maxpts, steps=steps)*sign_C_y0)
				
	denominator <- 1 - F_tmp[2] # onyly one claim left, (-1)^0*C(F(1)) + (-1)*C(F(0))
	result[i] <- numerator/denominator
	
	} else {
	result[i] <- 1 - F1[i]
	}



}
if (progress) cat(r0, r1, round(result[i],4), i,"out of", n, "\n")
}
}

result	
}


# Predictive probability for k = 1 based on (at most) r nearest claim neighbor within a small distance
# for SINGLE hail event
pred1.nnwithin.single.event <- function(X_event, evt_ind, y_event, within_dist=0.5,  
		para_est, longlat, r_max=5, type=c("Miwa", "Bretz"), abseps=1e-20, maxpts=1e5, steps=4096,
		fun=geosphere::distGeo,  progress=TRUE) {

if (r_max >= 10 && type == "Miwa") {
cat("Warning: Gaussian Copula CDF is numerical unstable and slow with Miwa when order (r_max) is large, switch to Bretz approximation!\n")

type <- "Bretz"
}


order.r <- function(x, r) {order(x)[1:r]}

pair_to_matrix <- function(pairwise_distance, n) {
pairwise_matrix <- matrix(0, nrow=n, ncol=n)
pairwise_matrix[(pairwise_distance$ind2 - 1)*n + pairwise_distance$ind1] <- pairwise_distance$distance
pairwise_matrix + t(pairwise_matrix) + diag(Inf,n)
}

F_combination <- function(F_tmp, all_comb) {
result <- all_comb
for (j in 1:NCOL(result)) {
F_TRUE_FALSE <- c(F_tmp[j], 1)

result[,j] <- F_TRUE_FALSE[all_comb[,j]]
}

result
}



n <- NROW(X_event)
n_claim <- sum(y_event)

claim_ind <- (1:n)[y_event == 1]

if (length(unique(evt_ind)) != 1) { 
	cat("Not single event data is supported!")
	return(rep(NA, n))
}

para_beta <- para_est[1:(length(para_est)-3)]
para_copula <- tail(para_est,3)
c(rho, psi, kappa) %<-% para_copula

# Given the model, calculate marginal F1 = 1 - P(y_kj|x_j) = 1/(1+exp(X_j beta)
result <- F1 <- 1/( 1 + exp(c(X_event%*%para_beta)))


if (r_max == 0) {
cat("r_max equals 0, just use marginal probability!\n")
return (1- F1)
}

if (n == 1 || n_claim == 0) {
result <- 1 - F1 # 3 cases here: only one data point or no claim or (one claim + y = 1)
} else {

# another case with multiple records in event





# loop here for all records in the event
for (i in 1:n) {


if (n_claim  > 1) {

	if (y_event[i] == 0) {
		# all pair about claim = 1 and the index within the distance
		pairwise_distance <- data.frame(expand.grid(ind1 = i, ind2 = claim_ind)) %>%
			mutate(distance=fun(longlat[ind1,], longlat[ind2,])/1e3) %>% filter(distance <= within_dist) %>% select(ind1,ind2)
		
		if (NROW(pairwise_distance) > 0) {
		n_claim_within <- sum(pairwise_distance == i)  # find claim within the distance
		
		tmp_ind <- unique(c(pairwise_distance[,1], pairwise_distance[,2]))
		
		pairwise_distance <- data.frame(expand.grid(ind1 = tmp_ind, ind2 = tmp_ind)) %>% filter(ind1 != ind2) %>%
			mutate(distance=fun(longlat[ind1,], longlat[ind2,])/1e3)
		} else { n_claim_within <- 0}
		
		# Case 0: When y is 0, can use all 1 neighbor
		r0 <- min(r_max, n_claim_within)
		all_comb_y0 <- binary.combination(r0 + 1)
		sign_C_y0 <- (-1)^apply(all_comb_y0,1,sum)
	
	
		# Case 1: When y = 1
		r1 <- min(r_max, max(n_claim_within - 1 ,0))

		# Case 1a: at least two 1
		if (r1 > 0) {
		all_comb_y1 <- binary.combination(r1 + 1)
		sign_C_y1 <- (-1)^apply(all_comb_y1,1,sum)
		}
	
		if( r0 > 0) {
		ind_tmp <- c(i, unname(unlist(
					pairwise_distance %>% filter(ind1 == i) %>% slice_min(distance, n=r0, with_ties=FALSE) %>% select(ind2)))
					)
	
	

	
		F_tmp <- F1[ind_tmp]
		F_matrix <- F_combination(F_tmp, all_comb_y0)
	
		R_matrix <- R_calculation(ind_tmp, para_copula, pairwise_distance)
		numerator <- sum(Cop.cdf.general(F_matrix, R_matrix=R_matrix, 
						type=type, abseps=abseps, maxpts=maxpts, steps=steps)*sign_C_y0)
					
		if (r0 > 1) {
		denominator <- sum(Cop.cdf.general(F_matrix[-c(1:2^r0),-1], R_matrix=R_matrix[-1,-1], 
							type=type, abseps=abseps, maxpts=maxpts, steps=steps)*sign_C_y0[-c(1:2^r0)])
		} else { denominator <- 1 - F_tmp[2]}
		result[i] <- numerator/denominator
			
		} else { result[i] <- 1 - F1[i] }
	
		
	
	} else {
		# all pair about claim = 1 and the index within the distance
		pairwise_distance <- data.frame(expand.grid(ind1 = i, ind2 = claim_ind)) %>% filter(ind1 != ind2) %>%
			mutate(distance=fun(longlat[ind1,], longlat[ind2,])/1e3) %>% filter(distance <= within_dist) %>% select(ind1,ind2)
		
		if (NROW(pairwise_distance) > 0) {
			n_claim_within <- sum(pairwise_distance == i)  # find claim within the distance
		
		tmp_ind <- unique(c(pairwise_distance[,1], pairwise_distance[,2]))
		
		pairwise_distance <- data.frame(expand.grid(ind1 = tmp_ind, ind2 = tmp_ind)) %>% filter(ind1 != ind2) %>%
			mutate(distance=fun(longlat[ind1,], longlat[ind2,])/1e3)
		} else { n_claim_within <- 0}
	
		# Case 0: When y is 0, can use all  neighbor
		r0 <- min(r_max, n_claim_within)
		all_comb_y0 <- binary.combination(r0 + 1)
		sign_C_y0 <- (-1)^apply(all_comb_y0,1,sum)
	
	
		# Case 1: When y = 1
		r1 <- min(r_max, max(n_claim_within ,0))

		# Case 1a: at least two 1
		if (r1 > 0) {
		all_comb_y1 <- binary.combination(r1 + 1)
		sign_C_y1 <- (-1)^apply(all_comb_y1,1,sum)
		}
	
	

	
		if( r0 > 0) {
	
		ind_tmp <- c(i, unname(unlist(
					pairwise_distance %>% filter(ind1 == i) %>% slice_min(distance, n=r1, with_ties=FALSE) %>% select(ind2)))
					)
				
		F_tmp <- F1[ind_tmp]
		F_matrix <- F_combination(F_tmp, all_comb_y1)


		R_matrix <- R_calculation(ind_tmp, para_copula, pairwise_distance)
		numerator <- sum(Cop.cdf.general(F_matrix, R_matrix=R_matrix, 
					type=type, abseps=abseps, maxpts=maxpts, steps=steps)*sign_C_y1)
	
		if (r1 > 1) {
		denominator <- sum(Cop.cdf.general(F_matrix[-c(1:2^r1),-1], R_matrix=R_matrix[-1,-1], 
					type=type, abseps=abseps, maxpts=maxpts, steps=steps)*sign_C_y1[-c(1:2^r1)])
		} else {denominator <- 1 - F_tmp[2] }
				
		result[i] <- numerator/denominator
	
		} else {result[i] <- 1 - F1[i] }
		}

} else { # i.e. only one claim
	if (y_event[i] == 0) {
		# all pair about claim = 1 and the index
		pairwise_distance <- data.frame(expand.grid(ind1 = c(i,claim_ind), ind2 = claim_ind)) %>% filter (ind1 != ind2) %>%
			mutate(distance=fun(longlat[ind1,], longlat[ind2,])/1e3) %>% filter(distance <= within_dist | (ind1 != i & ind2 != i) ) 
	
		if (NROW(pairwise_distance) > 0) {
		n_claim_within <- sum(pairwise_distance == i)  # find claim within the distance
		
		tmp_ind <- unique(c(pairwise_distance[,1], pairwise_distance[,2]))
		
		pairwise_distance <- data.frame(expand.grid(ind1 = tmp_ind, ind2 = tmp_ind)) %>% filter(ind1 != ind2) %>%
			mutate(distance=fun(longlat[ind1,], longlat[ind2,])/1e3)
		} else { n_claim_within <- 0}
	
		# Case 0: When y is 0, can use all  neighbor
		r0 <- min(r_max, n_claim_within)
		all_comb_y0 <- binary.combination(r0 + 1)
		sign_C_y0 <- (-1)^apply(all_comb_y0,1,sum)
	
	
		# Case 1: When y = 1
		r1 <- min(r_max, max(n_claim_within - 1 ,0))

		# Case 1a: at least two 1
		if (r1 > 0) {
		all_comb_y1 <- binary.combination(r1 + 1)
		sign_C_y1 <- (-1)^apply(all_comb_y1,1,sum)
		}
	
		if (r0 > 0) {
	
		ind_tmp <- c(i, unname(unlist(
					pairwise_distance %>% filter(ind1 == i) %>% slice_min(distance, n=r0, with_ties=FALSE) %>% select(ind2)))
					)
	
		F_tmp <- F1[ind_tmp]
		F_matrix <- F_combination(F_tmp, all_comb_y0)
	
		R_matrix <- R_calculation(ind_tmp, para_copula, pairwise_distance)
		numerator <- sum(Cop.cdf.general(F_matrix, R_matrix=R_matrix, 
					type=type, abseps=abseps, maxpts=maxpts, steps=steps)*sign_C_y0)
				
		denominator <- 1 - F_tmp[2] # onyly one claim left, (-1)^0*C(F(1)) + (-1)*C(F(0))
		result[i] <- numerator/denominator
	
		} else { result[i] <- 1 - F1[i] }



	} else { result[i] <- 1 - F1[i] 	} # coz only claim + y_event[i] = 1 => no more claim case left => 1 - F1 as prediction

}
if (progress) cat("r0, r1, Pr(Y=1):", r0, r1, round(result[i],4), i,"out of", n, "\n")
}
}

result	
}



# Two stage estimation for k = 2
# dat		: raw dataset contains all variables (WITH index of event group and row index)
# pair.dat	: matrix of incies of pair and distance
# form		: formula type for the first stage GLM fitting
# grad		: logical indicator if  require APPROXIMATE gradient for two stages
# hessian	: locgical indicator if require hessian matrix of two stages
# ini		: initial value of stage 2
# tol.bound	: Small bound for the parameter space
# maxeval	: maximum number of function evaluation for optimization (default: 2e3)
# factr		: relative improment of objective function for declaring convergence for optimization (default: 1e8)
# progress	: logical indicator if require printing progress

est2.twostage <- function(dat, pair.dat, form,  grad=TRUE, hessian=TRUE,
	ini=c(0.05, 5, 0.3), tol.bound=1e-5, maxeval=2e3, factr=1e8, progress=TRUE) {

# Stage 1 by Weibull AFT model
m.wei <- survival::survreg(form, data=dat, dist="weibull",
	score=TRUE,
	x=TRUE, y=TRUE)

X <- model.matrix(m.wei) ; y <- m.wei$y[,1]

beta.wei <- m.wei$coef
scale.wei <- m.wei$scale

scale.inv <- 1/scale.wei


if (progress) cat("stage 1 Weibull AFT done, estimate are:", c(beta.wei, scale.wei), "\n")

# Survival at observed t: exp(-t^scale.inv*exp(-scale.inv* X%*%beta ))
# Hazard at t: scale.inv*t^(scale.inv - 1)*exp(-scale.inv* X%*%beta )

# F2 =  P(Y <= y) = 1 - Survival
# f2 = hazard*(1- F)

stage1 <- data.frame(y=y) %>% mutate(e.scale.X.beta = exp(-scale.inv * X%*%beta.wei)) %>%
			mutate( F2 = 1 - exp(-y^scale.inv*e.scale.X.beta), 
					f2 = scale.inv*y^(scale.inv - 1)*e.scale.X.beta*(1-F2)) %>% 
			mutate(f2 = ifelse(f2 <= 0, .Machine$double.xmin,f2))


# prepare the pair matrix for lookup
pair.dat <- pair.dat %>% mutate(
				f21 = stage1$f2[ind1], f22 = stage1$f2[ind2],
				F21 = stage1$F2[ind1], F22 = stage1$F2[ind2]
		)


logCL <- function(para, pair.dat=pair.dat, progress) {


ans <- sum(pair.dat %>% mutate(C2 = c.fun(F21, F22, distance, para=para)) %>%
		mutate(C2 = ifelse(C2 <= 0, .Machine$double.xmin,C2)) %>%
		transmute(logf = log(f21) + log(f22) + log(C2))) 


if (progress) cat("logCL:", ans, "; para: ", para,"\n")

ans
}

nlogCL <- function(para, pair.dat=pair.dat, progress) {
-logCL(para, pair.dat=pair.dat, progress)
}

# Stage 2

#obj.bfgs <- nloptr::lbfgs(ini, fn=nlogCL, gr=NULL, lower=rep(tol.bound, 3), pair.dat=pair.dat,
#		upper=c(1 - tol.bound, Inf, 1 - tol.bound))


obj <- nloptr::sbplx(ini, fn=nlogCL, lower=rep(tol.bound, 3), pair.dat=pair.dat, progress=progress,
	upper=c(1 - tol.bound, Inf, 1 - tol.bound),
	control=list(xtol_rel = 1e-5, # stop on small optimization step
				 maxeval = maxeval, # stop on this many function evaluations
				 ftol_rel = 1/factr # stop on change times function value
			))
			
est <- c(beta.wei, log(scale.wei), obj$par)
names(est) <- c(names(beta.wei), "log.scale", "rho", "psi", "kappa")

result <- list(obj=obj, est=est)

if (grad) {

if (progress) cat("Estiamtion done, calculating numerical gradient\n")

result$score <- c(m.wei$score, numDeriv::grad(logCL,  c(obj$par), pair.dat=pair.dat,
	progress=progress, method="Richardson"))
}


if (hessian) {

if (progress) cat("Estiamtion done, calculating numerical hessian\n")

result$hessian <- list(stage1= -solve(vcov(m.wei)), 
					stage2 = numDeriv::hessian(logCL, c(obj$par), pair.dat=pair.dat,
					progress=progress, method="Richardson"))
}
	
result
}


# Full estimation log-likelihood when k = 2
# dat		: raw dataset contains all variables (WITH index of event group and row index)
# pair.dat	: matrix of incies of pair and distance
# form		: formula type for the first stage GLM fitting
# grad		: logical indicator if  require APPROXIMATE gradient for two stages
# hessian	: locgical indicator if require hessian matrix of two stages
# ini.full	: initial value probably from two stage approach
# tol.bound	: Small bound for the parameter space
# maxeval	: maximum number of function evaluation for optimization (default: 5e3)
# factr		: relative improment of objective function for declaring convergence for optimization (default: 1e8)
# progress	: logical indicator if require printing progress

est2.full <- function(dat, pair.dat, form,  grad=TRUE, hessian=TRUE,
	ini.full, tol.bound=1e-5, maxeval=5e3, factr=1e8, progress=TRUE) {

# function that evaluate Weibull model at given para.vec (without copula parameters)
wei.noitr <- function(para.vec, form, dat) {
# force by no iteration
obj <- survival::survreg(form, data=dat, dist="weibull",
	score = TRUE, init = para.vec, X=X.out, y=TRUE,
	control=survreg.control(maxiter=0))

list(obj=obj, score = obj$score, loglik=obj$loglik[2])
}



logCL <- function(para.copula, pair.dat=pair.dat) {

ans <- sum(pair.dat %>% mutate(C2 = c.fun(F21, F22, distance, para=para.copula)) %>%
		mutate(C2 = ifelse(C2 <= 0, .Machine$double.xmin,C2)) %>%
		transmute(logf = log(f21) + log(f22) + log(C2))) 

ans
}



obj.tmp <- survival::survreg(form, data=dat, dist="weibull",
	score = TRUE, init = head(ini.full, length(ini.full) - 3), x=TRUE, y=TRUE,
	control=survreg.control(maxiter=0))

X <- model.matrix(obj.tmp) ; y <- obj.tmp$y[,1]


loglik.full <- function(para.full, X, y, form, dat, pair.dat, progress) {

para.vec <- head(para.full,length(para.full) - 3) # scale in log scale!
para.copula <- tail(para.full,3)

stage1.result <- wei.noitr(head(para.full,length(para.full) - 3), form, dat)


beta.wei <- stage1.result$obj$coef
scale.wei <- stage1.result$obj$scale

scale.inv <- 1/scale.wei


# Survival at observed t: exp(-t^scale.inv*exp(-scale.inv* X%*%beta ))
# Hazard at t: scale.inv*t^(scale.inv - 1)*exp(-scale.inv* X%*%beta )

# F2 =  P(Y <= y) = 1 - Survival
# f2 = hazard*(1- F)

stage1 <- data.frame(y=y) %>% mutate(e.scale.X.beta = exp(-scale.inv * X%*%beta.wei)) %>%
			mutate( F2 = 1 - exp(-y^scale.inv*e.scale.X.beta), 
					f2 = scale.inv*y^(scale.inv - 1)*e.scale.X.beta*(1-F2)) %>% 
			mutate(f2 = ifelse(f2 <= 0, .Machine$double.xmin,f2))


# prepare the pair matrix for lookup
pair.dat <- pair.dat %>% mutate(
				f21 = stage1$f2[ind1], f22 = stage1$f2[ind2],
				F21 = stage1$F2[ind1], F22 = stage1$F2[ind2]
		)


f.value <- stage1.result$loglik + logCL(para.copula, pair.dat=pair.dat)
if (progress) cat("logLik:", f.value , "; para:", round(para.full,3),"\n")

f.value 
}

nloglik.full <- function(para.full, X, y, form, dat, pair.dat, progress) {
	-loglik.full(para.full, X, y, form, dat, pair.dat, progress)
}


CI <- confint(obj.tmp, level= 1 - 1e-8)



lb <- c(c(obj.tmp$coef, log(obj.tmp$scale)) - qnorm(1 - 1e-10)*sqrt(diag(vcov(obj.tmp)))
		, rep(tol.bound, 3))
ub <- c(c(obj.tmp$coef, log(obj.tmp$scale)) + qnorm(1 - 1e-10)*sqrt(diag(vcov(obj.tmp)))
		, c(1 - tol.bound, Inf, 1 - tol.bound))

# obj.bfgs <- #optim(ini.full, fn=nloglik.full, gr=NULL, lower=lb, upper=ub,method="L-BFGS-B",
		#X=X, y=y, form=form, dat=dat, pair.dat=pair.dat,
		#control = list(parscale = c(rep(1e6,length(ini.full)-3),rep(1,3))))
		
obj <- nloptr::sbplx(ini.full, fn=nloglik.full, lower=lb, upper=ub,
	X=X, y=y, form=form, dat=dat, pair.dat=pair.dat, progress=progress,
	control=list(xtol_rel = 1e-5, # stop on small optimization step
				 maxeval = maxeval, # stop on this many function evaluations
				 ftol_rel = 1/factr # stop on change times function value
			))



est <- obj$par
names(est) <- c(colnames(X), "log.scale", "rho", "psi", "kappa")

result <- list(obj=obj, est=est)


if (grad) {

if (progress) cat("Estimation done, calculating numerical gradient\n")

result$score <- numDeriv::grad(loglik.full, obj$par, X=X, y=y, form=form, dat=dat, pair.dat=pair.dat, 
	progress=FALSE, method="Richardson")
}

if (hessian) {

if (progress) cat("Estimation done, calculating numerical hessian\n")

result$hessian <- numDeriv::hessian(loglik.full, obj$par, X=X, y=y, form=form, dat=dat, pair.dat=pair.dat, 
	progress=FALSE, method="Richardson")
}

result
}



# Full estimation weighted log-likelihood when k = 2
# dat		: raw dataset contains all variables (WITH index of event group and row index)
# pair.dat	: matrix of incies of pair and distance
# form		: formula type for the first stage GLM fitting
# grad		: logical indicator if  require APPROXIMATE gradient for two stages
# hessian	: locgical indicator if require hessian matrix of two stages
# ini.full	: initial value probably from two stage approach
# tol.bound	: Small bound for the parameter space
# maxeval	: maximum number of function evaluation for optimization (default: 5e3)
# factr		: relative improment of objective function for declaring convergence for optimization (default: 1e8)
# progress	: logical indicator if require printing progress

est2.full.weight <- function(dat, pair.dat, form, grad=TRUE, hessian=TRUE,
	ini.full, tol.bound=1e-5, maxeval=5e3, factr=1e8, progress=TRUE) {

W1 <- NROW(dat)
W2 <- NROW(pair.dat)

# function that evaluate Weibull model at given para.vec (without copula parameters)
wei.noitr <- function(para.vec, form, dat) {
# force by no iteration
obj <- survival::survreg(form, data=dat, dist="weibull",
	score = TRUE, init = para.vec, X=X.out, y=TRUE,
	control=survreg.control(maxiter=0))

list(obj=obj, score = obj$score, loglik=obj$loglik[2])
}



logCL <- function(para.copula, pair.dat=pair.dat) {

ans <- sum(pair.dat %>% mutate(C2 = c.fun(F21, F22, distance, para=para.copula)) %>%
		mutate(C2 = ifelse(C2 <= 0, .Machine$double.xmin,C2)) %>%
		transmute(logf = log(f21) + log(f22) + log(C2))) 

ans
}



obj.tmp <- survival::survreg(form, data=dat, dist="weibull",
	score = TRUE, init = head(ini.full, length(ini.full) - 3), x=TRUE, y=TRUE,
	control=survreg.control(maxiter=0))

X <- model.matrix(obj.tmp) ; y <- obj.tmp$y[,1]


loglik.full.weight <- function(para.full, X, y, form, dat, pair.dat, W1,W2, progress) {

para.vec <- head(para.full,length(para.full) - 3) # scale in log scale!
para.copula <- tail(para.full,3)

stage1.result <- wei.noitr(head(para.full,length(para.full) - 3), form, dat)


beta.wei <- stage1.result$obj$coef
scale.wei <- stage1.result$obj$scale

scale.inv <- 1/scale.wei


# Survival at observed t: exp(-t^scale.inv*exp(-scale.inv* X%*%beta ))
# Hazard at t: scale.inv*t^(scale.inv - 1)*exp(-scale.inv* X%*%beta )

# F2 =  P(Y <= y) = 1 - Survival
# f2 = hazard*(1- F)

stage1 <- data.frame(y=y) %>% mutate(e.scale.X.beta = exp(-scale.inv * X%*%beta.wei)) %>%
			mutate( F2 = 1 - exp(-y^scale.inv*e.scale.X.beta), 
					f2 = scale.inv*y^(scale.inv - 1)*e.scale.X.beta*(1-F2)) %>% 
			mutate(f2 = ifelse(f2 <= 0, .Machine$double.xmin,f2))


# prepare the pair matrix for lookup
pair.dat <- pair.dat %>% mutate(
				f21 = stage1$f2[ind1], f22 = stage1$f2[ind2],
				F21 = stage1$F2[ind1], F22 = stage1$F2[ind2]
		)


f.value <- stage1.result$loglik/W2 + logCL(para.copula, pair.dat=pair.dat)/W2
if (progress) cat("logLik:", f.value , "; para:", round(para.full,3),"\n")

f.value 
}

nloglik.full.weight <- function(para.full, X, y, form, dat, pair.dat, W1, W2, progress) {
	-loglik.full.weight(para.full, X, y, form, dat, pair.dat, W1, W2, progress)
}


CI <- confint(obj.tmp, level= 1 - 1e-8)



lb <- c(c(obj.tmp$coef, log(obj.tmp$scale)) - qnorm(1 - 1e-10)*sqrt(diag(vcov(obj.tmp)))
		, rep(tol.bound, 3))
ub <- c(c(obj.tmp$coef, log(obj.tmp$scale)) + qnorm(1 - 1e-10)*sqrt(diag(vcov(obj.tmp)))
		, c(1 - tol.bound, Inf, 1 - tol.bound))

# obj.bfgs <- #optim(ini.full, fn=nloglik.full, gr=NULL, lower=lb, upper=ub,method="L-BFGS-B",
		#X=X, y=y, form=form, dat=dat, pair.dat=pair.dat,
		#control = list(parscale = c(rep(1e6,length(ini.full)-3),rep(1,3))))
		
obj <- nloptr::sbplx(ini.full, fn=nloglik.full.weight, lower=lb, upper=ub,
	X=X, y=y, form=form, dat=dat, pair.dat=pair.dat, W1=W1, W2=W2, progress=progress,
	control=list(xtol_rel = 1e-5, # stop on small optimization step
				 maxeval = maxeval, # stop on this many function evaluations
				 ftol_rel = 1/factr # stop on change times function value
			))



est <- obj$par
names(est) <- c(colnames(X), "log.scale", "rho", "psi", "kappa")

result <- list(obj=obj, est=est)


if (grad) {

if (progress) cat("Estimation done, calculating numerical gradient\n")

result$score <- numDeriv::grad(loglik.full.weight, obj$par, X=X, y=y, form=form, dat=dat, pair.dat=pair.dat, 
		W1=W1, W2=W2, progress=progress, method="Richardson")
}


if (hessian) {

if (progress) cat("Estimation done, calculating numerical hessian\n")

result$hessian <- numDeriv::hessian(loglik.full.weight, obj$par, X=X, y=y, form=form, dat=dat, pair.dat=pair.dat, 
		W1=W1, W2=W2, progress=progress, method="Richardson")
}

result
}


# Predicted mean for k = 2 based on (at most) r nearest claim neighbor within a small distance
# for SINGLE hail event
# n_node: number of nodes in Gaussian-Laguerre quadrature (larger the better, 50 seems large enough)
# alpha: choice of weight function Gaussian-Laguerre quadrature (result is robust, control the "decay")
# Complexity: O(n_node*event*r_max^2) ~ O(n_node*event^3)
pred2.nnwithin.single.event <- function(X_event, evt_ind, y_event, within_dist=0.5,  
		para_est, longlat, r_max=5,
		n_node=50, alpha=2,
		fun=geosphere::distGeo,  pred_density=FALSE, pred_seq=c(1e-2,1e3, by=1e-2), progress=TRUE) {


pair_to_matrix <- function(pairwise_distance, n) {
pairwise_matrix <- matrix(0, nrow=n, ncol=n)
pairwise_matrix[(pairwise_distance$ind2 - 1)*n + pairwise_distance$ind1] <- pairwise_distance$distance
pairwise_matrix + t(pairwise_matrix) + diag(Inf,n)
}



n <- NROW(X_event)

if (length(unique(evt_ind)) != 1) { 
	cat("Not single event data is supported!")
	return(rep(NA, n))
}

para_beta <- para_est[1:(length(para_est)-3)]

beta.wei <- matrix(para_beta[1:(length(para_beta)-1)],ncol=1)
scale.wei <- exp(tail(para_beta,1))
scale.inv <- 1/scale.wei

para_copula <- tail(para_est,3)
c(rho, psi, kappa) %<-% para_copula


# use survreg notation here and note that
# survreg's scale = 1/(rweibull shape)
# survreg's intercept = log(rweibull scale) !!!

# CDF: F2 =  P(Y <= t) = 1 - Survival = 1 - exp(-(t/b)^a
# PDF: f2 = d P(Y <= t)/dt = a/b*(t/b)^(a-1)*exp(-(t/b)^a
# E(T): b*Gamma(1+1/a),
# where a = 1/survreg$scale ; b = exp(X%*%survreg$coef)

a <- 1/scale.wei ; b <- c(exp(X_event%*%beta.wei))

yFf <- data.frame(y=y_event) %>%
	mutate(F2= 1 - exp(-(y/b)^a), f2=a/b*(y/b)^(a-1)*exp(-(y/b)^a)) %>%
	mutate(f2 = ifelse(f2 <= 0, .Machine$double.xmin,f2))

result_marginal_mean 	<- b*gamma(1 + 1/a)
result_marginal_median 	<- b*log(2)^(1/a)

result <- result_marginal_mean

out <- gauss.quad(n_node ,"laguerre",alpha=alpha)

if (r_max == 0 || n == 1) {
cat("r_max equals 0 or single obs event, just use marginal distribution!\n")


} else {

# loop here

if (pred_density) {
	len_seq <- length(pred_seq)
	pred_density_matrix <- matrix(NA, nrow=n, ncol=len_seq)
}

for (i in 1:n) {


pairwise_distance <- data.frame(ind1=i, ind2=(1:n)[-i]) %>%
			mutate(distance=fun(longlat[ind1,], longlat[ind2,])/1e3) %>% 
				filter( distance <= within_dist)

n_claim_within <- NROW(pairwise_distance) # how many record within with_dist

r0 <- min(r_max, n_claim_within) + 1 # number of terms in numerator
r1 <- max(r0 - 1, 0) # number of term in denominator
			
if (NROW(pairwise_distance) > 0) {
	neighbor_within <- unname(unlist(pairwise_distance %>% 
				slice_min(distance, n=r1, with_ties=FALSE) %>% 
				select(ind2)))
				
	ind_tmp <- c(i, neighbor_within)

	if (length(neighbor_within) >=2)  {
		add_ind <- combn(neighbor_within, 2)
	
	# add more rows to pairwise_distance
	pairwise_distance <- bind_rows(
		pairwise_distance,
		data.frame(ind1=add_ind[1,], ind2=add_ind[2,]) %>%
			mutate(distance=fun(longlat[ind1,], longlat[ind2,])/1e3)
	)
	}
	
	R_matrix <- R_calculation(ind_tmp, para_copula, pairwise_distance)
	
	F_matrix <-  matrix(yFf$F2[ind_tmp], nrow=n_node, ncol=r0, byrow=TRUE)
	
	X_event_tmp <- X_event[rep(i,n_node),]
	
	a <- 1/scale.wei ; b_tmp <- c(exp(X_event_tmp%*%beta.wei))

	tmp <- data.frame(y=out$nodes) %>%
		mutate(F2= 1 - exp(-(y/b_tmp)^a), f2=a/b_tmp*(y/b_tmp)^(a-1)*exp(-(y/b_tmp)^a)) %>%
		mutate(f2 = ifelse(f2 <= 0, .Machine$double.xmin,f2))
	
	
	F_matrix[,1] <- tmp$F2
	
	numerator <- Cop.pdf.general(F_matrix, R_matrix)
	
	if (r1 > 0) { # at least one neighbor in denominator
		denominator <- Cop.pdf.general(yFf$F2[ind_tmp[-1]], R_matrix[-1,-1])
	} else { # only one neighbor in denominator c = dC(u)/du = du/du = 1
		denominator <- 1
	}
	denominator_vec <- rep(denominator, n_node)
	
	# E(Y) = \int_0^inf y*f(y|...)dy ~ \sum_y* w(y*) g(y*) by Gaussian-Laguerre quadrature
	# where g(y) = y*f(y|...)/w(y) = y*f(y)*numerator/denominator/w(y) and w(y) = y^alpha*exp(-y)
	result[i] <- tmp %>% select(y, f2) %>% mutate(
			numerator=numerator, denominator=denominator) %>%
			transmute(finally = exp((1-alpha)*log(y) + log(f2) + log(numerator)-log(denominator_vec) + y + log(out$weights))) %>%
			sum()
			
	if (pred_density) {
		R_matrix <- R_calculation(ind_tmp, para_copula, pairwise_distance)
	
		F_matrix <-  matrix(yFf$F2[ind_tmp], nrow=len_seq, ncol=r0, byrow=TRUE)
	
		X_event_tmp <- X_event[rep(i,len_seq),]
		
		a <- 1/scale.wei ; b_tmp <- c(exp(X_event_tmp%*%beta.wei))
		
		tmp <- data.frame(y=pred_seq) %>%
			mutate(F2= 1 - exp(-(y/b_tmp)^a), f2=a/b_tmp*(y/b_tmp)^(a-1)*exp(-(y/b_tmp)^a)) %>%
			mutate(f2 = ifelse(f2 <= 0, .Machine$double.xmin,f2))
			
		F_matrix[,1] <- tmp$F2
		numerator <- Cop.pdf.general(F_matrix, R_matrix)
		
		pred_density_matrix[i,] <- tmp$f2*numerator/denominator
		
	}
} else {
	if (pred_density) {
		# use marginal for density
		X_event_tmp <- X_event[rep(i,len_seq),]
		
		a <- 1/scale.wei ; b_tmp <- c(exp(X_event_tmp%*%beta.wei))
		
		tmp <- data.frame(y=pred_seq) %>%
			mutate(F2= 1 - exp(-(y/b_tmp)^a), f2=a/b_tmp*(y/b_tmp)^(a-1)*exp(-(y/b_tmp)^a)) %>%
			mutate(f2 = ifelse(f2 <= 0, .Machine$double.xmin,f2))
		
		pred_density_matrix[i,] <- tmp$f2
	}
}

if (progress) {
	cat("r0, r1,  Y, E(Y), E(Y|...):", r0, r1, round(c(y_event[i], result_marginal_mean[i], 
		result[i]),2), i,"out of", n, "\n")
}

} 


}

obj_pred <- data.frame(pred_mean_marginal=result_marginal_mean, 
	pred_median_marginal=result_marginal_median, 
	pred_mean_copula=result)
	
if (pred_density) {
	obj_pred <- list(predicted_value=obj_pred, 
		pred_density_list=list(pred_density_matrix=pred_density_matrix, pred_seq=pred_seq))
}
	
obj_pred
}






# Two stage estimation for k = 3
# dat		: raw dataset contains all variables (WITH index of event group and row index)
# pair.dat	: matrix of incies of pair and distance
# form_T: formula for k = 2 (F_T)
# form_Z: formula for k = 3 (F_Z) 
# grad		: logical indicator if  require APPROXIMATE gradient for two stages
# hessian	: locgical indicator if require hessian matrix of two stages
# ini		: initial value of stage 2
# tol.bound	: Small bound for the parameter space
# maxeval	: maximum number of function evaluation for optimization (default: 1e3)
# factr		: relative improment of objective function for declaring convergence for optimization (default: 1e8)
# progress	: logical indicator if require printing progress

est3.twostage <- function(dat, pair.dat, form_T, form_Z, grad=TRUE, hessian=TRUE,
	ini=c(0.05, 5, 0.3), tol.bound=1e-5, maxeval=1e3, factr=1e8, progress=TRUE) {

##### Stage 1 #####

# Weibull AFT model
m.wei <- survival::survreg(form_T, data=dat, dist="weibull",
	score=TRUE,
	x=TRUE, y=TRUE)

X_T <- model.matrix(m.wei) ; y_T <- m.wei$y[,1]

beta.wei <- m.wei$coef
names(beta.wei) <- paste0(names(beta.wei), "_T")
scale.wei <- m.wei$scale
names(scale.wei) <- "scale_T"

scale.inv <- 1/scale.wei ;  names(scale.inv) <- "scale_inv_T"

if (progress) cat("Stage 1 Weibull done, estimate are:", c(beta.wei, scale.wei), "\n")

F_T <- 1 - exp(-y_T^scale.inv*exp(-scale.inv * X_T%*%beta.wei))

# Gamma regression
m.gamma <- glm(form_Z, dat=dat, family = Gamma(link="log"))

X_Z <- model.matrix(m.gamma); y_Z <-  m.gamma$y

beta.gamma <- coef(m.gamma); names(beta.gamma) <- paste0(names(beta.gamma),"_Z")
disp.gamma <- summary(m.gamma)$dispersion

nu <- 1/disp.gamma 
log_mu <- X_Z%*%beta.gamma
mu <- exp(log_mu)

if (progress) cat("Stage 1 Gamma done, estimate are:", round(c(beta.wei, disp.gamma),3), "\n")


# Note dispersion parameter of Gamma regression is not estimated by MLE but method of moment
# nu <- (NROW(dat) - length(beta.gamma))/sum((y_Z/exp(log_mu) - 1)^2)  # MOM here
# logf_Z <- -nu*(y_Z/exp(log_mu) + log_mu) + (nu -1)*log(y_Z) + nu*log(nu) - lgamma(nu) # log-likelihood 

f_Z <- dgamma(y_Z, shape = nu, scale= mu/nu)
logf_Z <- dgamma(y_Z, shape = nu, scale= mu/nu, log=TRUE) # smaller rounding error than log(f_Z)?

F_Z <- pgamma(y_Z, shape = nu, scale= mu/nu)


# Bivariate copula (use Gaussian (1) here), seems Frank (5) matched tau slightly better
# general structure (probably slower)
cop <- VineCopula::BiCopEst(F_T, F_Z, family=1)

# small function for log-likelihood of bivariate copula at given parameter value (assume no par2)
logLik.cop <- function(para, cop) {

cop.tmp <- VineCopula::BiCop(cop$family, para)

sum(log(VineCopula::BiCopPDF(F_T, F_Z, cop.tmp)))
}

if (progress) cat("Stage 1 Bivariate copula done, estimate:", round(cop$par,3), "\n")

h_Z <- VineCopula::BiCopHfunc1(F_T, F_Z, cop) # d H(u,v)/dv (or du?) (i.e. conditional H function)
h <- VineCopula::BiCopPDF(F_T, F_Z, cop) # d^2 H(u,v)/ dudv (i.e. Copula density)

F3 <- h_Z
# f3 <- exp(logf_Z + log(h)) # f_Z*h(F_T, F_Z) # use exp(log) for numerical stability
log_f3 <- logf_Z + log(h)

stage1 <- data.frame(F3 = F3, log_f3 = log_f3) %>% 
			mutate( F3 = ifelse(F3 <= 0 | is.nan(F3), .Machine$double.xmin,F3))


# prepare the pair matrix for lookup
pair.dat <- pair.dat %>% mutate(
				log_f31 = stage1$log_f3[ind1], log_f32 = stage1$log_f3[ind2],
				F31 = stage1$F3[ind1], F32 = stage1$F3[ind2]
		)


logCL <- function(para, pair.dat=pair.dat, progress) {


ans <- sum(pair.dat %>% mutate(C3 = c.fun(F31, F32, distance, para=para)) %>%
		mutate(C3 = ifelse(C3 <= 0, .Machine$double.xmin, C3)) %>%
		transmute(logf = log_f31 + log_f32 + log(C3))) 


if (progress) cat("logCL:", ans, "; para: ", para,"\n")

ans
}

nlogCL <- function(para, pair.dat=pair.dat, progress) {
-logCL(para, pair.dat=pair.dat, progress)
}

# Stage 2

#obj.bfgs <- nloptr::lbfgs(ini, fn=nlogCL, gr=NULL, lower=rep(tol.bound, 3), pair.dat=pair.dat,
#		upper=c(1 - tol.bound, Inf, 1 - tol.bound))


obj <- nloptr::sbplx(ini, fn=nlogCL, lower=rep(tol.bound, 3), pair.dat=pair.dat, progress=progress,
	upper=c(1 - tol.bound, Inf, 1 - tol.bound),
	control=list(xtol_rel = 1e-5, # stop on small optimization step
				 maxeval = maxeval, # stop on this many function evaluations
				 ftol_rel = 1/factr # stop on change times function value
			))
			
est <- c(beta.wei, log(scale.wei), beta.gamma, disp.gamma, cop$par, obj$par)
names(est) <- c(names(beta.wei), "log_scale_T", 
				names(beta.gamma), "disp_Z",
				"copula_TZ",
				"rho", "psi", "kappa")

result <- list(obj=obj, est=est)

if (grad) {

if (progress) cat("Estimation done, calculating numerical gradient\n")

# Score for Weibull regression
score.wei <- m.wei$score; names(score.wei) <- c(names(beta.wei),"scale_T")

# "Approximate" score of Gamma regression
z_tmp <- y_Z/c(mu) - 1
score.gamma <- nu*c(sum(z_tmp),c(t(X_Z[,-1]) %*%z_tmp))
names(score.gamma) <- paste0(names(beta.gamma), "_Z")


score.cop <- numDeriv::grad(logLik.cop, cop$par, cop=cop)
names(score.cop) <- "copula_TZ"


score.spatial <- numDeriv::grad(logCL,  c(obj$par), pair.dat=pair.dat,
	progress=progress, method="Richardson")

names(score.spatial ) <- c("rho", "psi", "kappa")
	
result$score <- c(score.wei, score.gamma, score.cop, score.spatial)
}


if (hessian) {

if (progress) cat("Estimation done, calculating numerical hessian\n")
	
result$hessian <- list(
	stage1=list(weibull = -solve(vcov(m.wei)), gamma= -solve(vcov(m.gamma)),
	copula=numDeriv::hessian(logLik.cop, cop$par, cop=cop)),
			stage2 = numDeriv::hessian(logCL,  c(obj$par), pair.dat=pair.dat,
	progress=progress, method="Richardson"))
}
	
result
}




# Full estimation log-likelihood when k = 3
# dat		: raw dataset contains all variables (WITH index of event group and row index)
# pair.dat	: matrix of incies of pair and distance
# form_T: formula for k = 2 (F_T)
# form_Z: formula for k = 3 (F_Z) 
# grad		: logical indicator if  require APPROXIMATE gradient for two stages
# hessian	: locgical indicator if require hessian matrix of two stages
# ini		: initial value of stage 2
# tol.bound	: Small bound for the parameter space
# maxeval	: maximum number of function evaluation for optimization (default: 1e4)
# factr		: relative improment of objective function for declaring convergence for optimization (default: 1e8)
# progress	: logical indicator if require printing progress

est3.full <- function(dat, pair.dat, form_T, form_Z, grad=TRUE, hessian=FALSE,
	cop.family=1,
	ini.full, tol.bound=1e-5, factr=1e10, maxeval=1e4, progress=TRUE) {

############### Stage 1 ###############

# function that evaluate Weibull model at given para.vec (without copula parameters)
logLik.Weibull <- function(para.vec, form, dat, X_T) {
# force by no iteration
obj <- survival::survreg(form, data=dat, dist="weibull",
		init = para.vec, X=X.out, y=TRUE,
		control=survreg.control(maxiter=0))
	
beta.wei <- obj$coef
scale.wei <- obj$scale
scale.inv <- 1/scale.wei

F_T <- 1 - exp(-y_T^scale.inv*exp(-scale.inv * X_T%*%beta.wei))

list(obj=obj, loglik=obj$loglik[2], F_T=c(F_T))
}

# function that evaluate Gamma model at given para.vec (beta.gamma, disp.gamma)
logLik.Gamma <- function(para.vec, X_Z, y_Z) {

beta.gamma <- matrix(para.vec[1:NCOL(X_Z)],ncol=1)

disp.gamma <- tail(para.vec,1)

nu <- 1/disp.gamma 
log_mu <- c(X_Z%*%beta.gamma)
mu <- exp(log_mu)

f_Z <- dgamma(y_Z, shape = nu, scale= mu/nu)

logf_Z <- dgamma(y_Z, shape = nu, scale= mu/nu, log=TRUE) # smaller rounding error than log(f_Z)?

F_Z <- pgamma(y_Z, shape = nu, scale= mu/nu)

list(f_Z = f_Z, logf_Z = logf_Z, F_Z = F_Z, loglik=sum(logf_Z))
}


# function for log-likelihood of bivariate copula at given parameter value (assume no par2)
logLik.TZ <- function(para, F_T, F_Z, cop.family=1) {
cop.tmp <- VineCopula::BiCop(cop.family, para)
list(cop = cop.tmp, loglik=sum(log(VineCopula::BiCopPDF(F_T, F_Z, cop.tmp))))
}


# function for log CL for stage 2
logCL <- function(para, pair.dat=pair.dat, progress) {


ans <- sum(pair.dat %>% mutate(C3 = c.fun(F31, F32, distance, para=para)) %>%
		mutate(C3 = ifelse(C3 <= 0, .Machine$double.xmin, C3)) %>%
		transmute(logf = log_f31 + log_f32 + log(C3))) 
ans
}


logLik.full <- function(para.full, dat, pair.dat, X_T, y_T, X_Z, y_Z, cop.family=1, progress=TRUE) {

# Break parameter vector
para.tmp <- para.full

para.wei 	<- para.tmp[1:(NCOL(X_T)+1)] ; para.tmp <- para.tmp[-(1:(NCOL(X_T)+1))]

para.gamma  <- para.tmp[1:(NCOL(X_Z)+1)] ; para.tmp <- para.tmp[-(1:(NCOL(X_Z)+1))]

para.TZ <- para.tmp[1] ; para.tmp <- para.tmp[-1]

para.spatial <- para.tmp

###### Stage 1 #####
ans.wei 	<- logLik.Weibull(para.wei, form_T, dat, X_T)
ans.gamma 	<- logLik.Gamma(para.gamma, X_Z, y_Z)

ans.gamma.global <<- ans.gamma

F_T 	<- ans.wei$F_T

F_Z  	<- ans.gamma$F_Z
f_Z  	<- ans.gamma$f_Z
logf_Z 	<- ans.gamma$logf_Z


ans.TZ <- logLik.TZ(para.TZ, F_T, F_Z, cop.family)
cop <- ans.TZ$cop

h_Z <- VineCopula::BiCopHfunc1(F_T, F_Z, cop) # d H(u,v)/dv (or du?) (i.e. conditional H function)
h <- VineCopula::BiCopPDF(F_T, F_Z, cop) # d^2 H(u,v)/ dudv (i.e. Copula density)

loglik.stage1 <- ans.wei$loglik + ans.gamma$loglik + ans.TZ$loglik

############### Stage 2 ###############
# prepare pair.dat from stage 1 values
F3 <- h_Z
log_f3 <- logf_Z + log(h)

stage1 <- data.frame(F3 = F3, log_f3 = log_f3) %>% 
			mutate( F3 = ifelse(F3 <= 0 | is.nan(F3), .Machine$double.xmin,F3))


# prepare the pair matrix for lookup
pair.dat <- pair.dat %>% mutate(
				log_f31 = stage1$log_f3[ind1], log_f32 = stage1$log_f3[ind2],
				F31 = stage1$F3[ind1], F32 = stage1$F3[ind2]
		)
		

value <- loglik.stage1 + logCL(para.spatial, pair.dat=pair.dat, FALSE)
if(progress) {
cat("logLik:", value, "para:", round(para.full,3),"\n")
}

value
}





nlogLik.full <- function(para.full, dat, pair.dat, X_T, y_T, X_Z, y_Z, cop.family=1, progress=TRUE) {
-logLik.full(para.full, dat, pair.dat, X_T, y_T, X_Z, y_Z, cop.family, progress=progress)
}


# Prepare model matrices
m.wei <- survival::survreg(form_T, data=dat, dist="weibull", x=TRUE, y=TRUE)

X_T <- model.matrix(m.wei) ; y_T <- m.wei$y[,1] ; no_para_T <- NCOL(X_T)

m.gamma <- glm(form_Z, dat=dat, family = Gamma(link="log"))

X_Z <- model.matrix(m.gamma); y_Z <-  m.gamma$y


no.se <- 2

lb <- c(
summary(m.wei)$table[,1]  - no.se*summary(m.wei)$table[,2],		# beta_T, log_scale_T
summary(m.gamma)$coef[,1] - no.se*summary(m.gamma)$coef[,2], 	# beta_z
max(summary(m.gamma)$dispersion - 0.2, 0),						# dispersion_Z
-1 + tol.bound,													# bivariate copula
rep(tol.bound, 3)												# spatial copula
)

ub <- c(
summary(m.wei)$table[,1]  + no.se*summary(m.wei)$table[,2],		# beta_T, log_scale_T
summary(m.gamma)$coef[,1] + no.se*summary(m.gamma)$coef[,2], 	# beta_z
max(summary(m.gamma)$dispersion + 0.2, 0),						# dispersion_Z
1 - tol.bound,													# bivariate copula
c(1 - tol.bound, Inf, 1 - tol.bound)							# spatial copula 
)



obj <- nloptr::sbplx(ini.full, fn=nlogLik.full, 
	lower=lb, upper=ub,
	dat=dat, pair.dat=pair.dat, X_T=X_T, y_T=y_T, X_Z=X_Z, y_Z=y_Z, cop.family=cop.family, progress=progress,
	control=list(xtol_rel = 1e-5, # stop on small optimization step
				 maxeval = maxeval, # stop on this many function evaluations
				 ftol_rel = 1/factr # stop on change times function value
			))

est <- obj$par			
names(est) <- c(paste0(names(coef(m.wei)),"_T"), "log_scale_T", 
				paste0(names(coef(m.gamma)),"_Z"), "disp_Z",
				"copula_TZ",
				"rho", "psi", "kappa")

result <- list(obj=obj, est=est)

if (grad) {

if (progress) cat("Estimation done, calculating numerical gradient\n")

result$score <- numDeriv::grad(logLik.full, est, 
	dat=dat, pair.dat=pair.dat, X_T=X_T, y_T=y_T, X_Z=X_Z, y_Z=y_Z, cop.family=cop.family, progress=progress) 
}

if (hessian) {

if (progress) cat("Estimation done, calculating numerical hessian\n")

result$hessian <- numDeriv::hessian(logLik.full, est, 
	dat=dat, pair.dat=pair.dat, X_T=X_T, y_T=y_T, X_Z=X_Z, y_Z=y_Z, cop.family=cop.family, progress=progress) 
}


	
result
}



# Full estimation weighted log-likelihood when k = 3
# dat		: raw dataset contains all variables (WITH index of event group and row index)
# pair.dat	: matrix of incies of pair and distance
# form_T: formula for k = 2 (F_T)
# form_Z: formula for k = 3 (F_Z) 
# grad		: logical indicator if  require APPROXIMATE gradient for two stages
# hessian	: locgical indicator if require hessian matrix of two stages
# ini		: initial value of stage 2
# tol.bound	: Small bound for the parameter space
# maxeval	: maximum number of function evaluation for optimization (default: 1e4)
# factr		: relative improment of objective function for declaring convergence for optimization (default: 1e8)
# progress	: logical indicator if require printing progress

est3.full.weight <- function(dat, pair.dat, form_T, form_Z, grad=TRUE, hessian=FALSE,
	cop.family=1,
	ini.full, tol.bound=1e-5, factr=1e10, maxeval=1e4, progress=TRUE) {

W1 <- NROW(dat) ; W2 <- NROW(pair.dat)

############### Stage 1 ###############

# function that evaluate Weibull model at given para.vec (without copula parameters)
logLik.Weibull <- function(para.vec, form, dat, X_T) {
# force by no iteration
obj <- survival::survreg(form, data=dat, dist="weibull",
		init = para.vec, X=X.out, y=TRUE,
		control=survreg.control(maxiter=0))
	
beta.wei <- obj$coef
scale.wei <- obj$scale
scale.inv <- 1/scale.wei

F_T <- 1 - exp(-y_T^scale.inv*exp(-scale.inv * X_T%*%beta.wei))

list(obj=obj, loglik=obj$loglik[2], F_T=c(F_T))
}

# function that evaluate Gamma model at given para.vec (beta.gamma, disp.gamma)
logLik.Gamma <- function(para.vec, X_Z, y_Z) {

beta.gamma <- matrix(para.vec[1:NCOL(X_Z)],ncol=1)

disp.gamma <- tail(para.vec,1)

nu <- 1/disp.gamma 
log_mu <- c(X_Z%*%beta.gamma)
mu <- exp(log_mu)

f_Z <- dgamma(y_Z, shape = nu, scale= mu/nu)

logf_Z <- dgamma(y_Z, shape = nu, scale= mu/nu, log=TRUE) # smaller rounding error than log(f_Z)?

F_Z <- pgamma(y_Z, shape = nu, scale= mu/nu)

list(f_Z = f_Z, logf_Z = logf_Z, F_Z = F_Z, loglik=sum(logf_Z))
}


# function for log-likelihood of bivariate copula at given parameter value (assume no par2)
logLik.TZ <- function(para, F_T, F_Z, cop.family=1) {
cop.tmp <- VineCopula::BiCop(cop.family, para)
list(cop = cop.tmp, loglik=sum(log(VineCopula::BiCopPDF(F_T, F_Z, cop.tmp))))
}


# function for log CL for stage 2
logCL <- function(para, pair.dat=pair.dat, progress) {


ans <- sum(pair.dat %>% mutate(C3 = c.fun(F31, F32, distance, para=para)) %>%
		mutate(C3 = ifelse(C3 <= 0, .Machine$double.xmin, C3)) %>%
		transmute(logf = log_f31 + log_f32 + log(C3))) 
ans
}


logLik.full.weight <- function(para.full, dat, pair.dat, X_T, y_T, X_Z, y_Z, cop.family=1, W1=W1, W2=W2, progress=TRUE) {

# Break parameter vector
para.tmp <- para.full

para.wei 	<- para.tmp[1:(NCOL(X_T)+1)] ; para.tmp <- para.tmp[-(1:(NCOL(X_T)+1))]

para.gamma  <- para.tmp[1:(NCOL(X_Z)+1)] ; para.tmp <- para.tmp[-(1:(NCOL(X_Z)+1))]

para.TZ <- para.tmp[1] ; para.tmp <- para.tmp[-1]

para.spatial <- para.tmp

###### Stage 1 #####
ans.wei 	<- logLik.Weibull(para.wei, form_T, dat, X_T)
ans.gamma 	<- logLik.Gamma(para.gamma, X_Z, y_Z)

ans.gamma.global <<- ans.gamma

F_T 	<- ans.wei$F_T

F_Z  	<- ans.gamma$F_Z
f_Z  	<- ans.gamma$f_Z
logf_Z 	<- ans.gamma$logf_Z


ans.TZ <- logLik.TZ(para.TZ, F_T, F_Z, cop.family)
cop <- ans.TZ$cop

h_Z <- VineCopula::BiCopHfunc1(F_T, F_Z, cop) # d H(u,v)/dv (or du?) (i.e. conditional H function)
h <- VineCopula::BiCopPDF(F_T, F_Z, cop) # d^2 H(u,v)/ dudv (i.e. Copula density)

loglik.stage1 <- ans.wei$loglik + ans.gamma$loglik + ans.TZ$loglik

############### Stage 2 ###############
# prepare pair.dat from stage 1 values
F3 <- h_Z
log_f3 <- logf_Z + log(h)

stage1 <- data.frame(F3 = F3, log_f3 = log_f3) %>% 
			mutate( F3 = ifelse(F3 <= 0 | is.nan(F3), .Machine$double.xmin,F3))


# prepare the pair matrix for lookup
pair.dat <- pair.dat %>% mutate(
				log_f31 = stage1$log_f3[ind1], log_f32 = stage1$log_f3[ind2],
				F31 = stage1$F3[ind1], F32 = stage1$F3[ind2]
		)
		

value <- loglik.stage1/W1 + logCL(para.spatial, pair.dat=pair.dat, FALSE)/W2
if(progress) {
cat("Weighted logLik:", value, "para:", round(para.full,3),"\n")
}

value
}





nlogLik.full <- function(para.full, dat, pair.dat, X_T, y_T, X_Z, y_Z, cop.family=1, W1=W1, W2=W2, progress=TRUE) {
-logLik.full.weight(para.full, dat, pair.dat, X_T, y_T, X_Z, y_Z, cop.family=cop.family, W1=W1, W2=W2, progress=TRUE) 
}


# Prepare model matrices
m.wei <- survival::survreg(form_T, data=dat, dist="weibull", x=TRUE, y=TRUE)

X_T <- model.matrix(m.wei) ; y_T <- m.wei$y[,1] ; no_para_T <- NCOL(X_T)

m.gamma <- glm(form_Z, dat=dat, family = Gamma(link="log"))

X_Z <- model.matrix(m.gamma); y_Z <-  m.gamma$y


no.se <- 2

lb <- c(
summary(m.wei)$table[,1]  - no.se*summary(m.wei)$table[,2],		# beta_T, log_scale_T
summary(m.gamma)$coef[,1] - no.se*summary(m.gamma)$coef[,2], 	# beta_z
max(summary(m.gamma)$dispersion - 0.2, 0),						# dispersion_Z
-1 + tol.bound,													# bivariate copula
rep(tol.bound, 3)												# spatial copula
)

ub <- c(
summary(m.wei)$table[,1]  + no.se*summary(m.wei)$table[,2],		# beta_T, log_scale_T
summary(m.gamma)$coef[,1] + no.se*summary(m.gamma)$coef[,2], 	# beta_z
max(summary(m.gamma)$dispersion + 0.2, 0),						# dispersion_Z
1 - tol.bound,													# bivariate copula
c(1 - tol.bound, Inf, 1 - tol.bound)							# spatial copula 
)



obj <- nloptr::sbplx(ini.full, fn=nlogLik.full, 
	lower=lb, upper=ub,
	dat=dat, pair.dat=pair.dat, X_T=X_T, y_T=y_T, X_Z=X_Z, y_Z=y_Z, cop.family=cop.family,
	W1=W1, W2=W2, progress=progress,
	control=list(xtol_rel = 1e-5, # stop on small optimization step
				 maxeval = maxeval, # stop on this many function evaluations
				 ftol_rel = 1/factr # stop on change times function value
			))

est <- obj$par			
names(est) <- c(paste0(names(coef(m.wei)),"_T"), "log_scale_T", 
				paste0(names(coef(m.gamma)),"_Z"), "disp_Z",
				"copula_TZ",
				"rho", "psi", "kappa")

result <- list(obj=obj, est=est)

if (grad) {

if (progress) cat("Estimation done, calculating numerical gradient\n")

result$score <- numDeriv::grad(logLik.full.weight, est, 
	dat=dat, pair.dat=pair.dat, X_T=X_T, y_T=y_T, X_Z=X_Z, y_Z=y_Z, cop.family=cop.family, W1=W1, W2=W2, progress=progress) 
}

if (hessian) {

if (progress) cat("Estimation done, calculating numerical hessian\n")

result$hessian <- numDeriv::hessian(logLik.full.weight, est, 
	dat=dat, pair.dat=pair.dat, X_T=X_T, y_T=y_T, X_Z=X_Z, y_Z=y_Z, cop.family=cop.family, W1=W1, W2=W2, progress=progress) 
}


	
result
}


# Predicted mean for k = 3 based on (at most) r nearest claim neighbor within a small distance
# for SINGLE hail event
# n_node: number of nodes in Gaussian-Laguerre quadrature (larger the better, 50 seems large enough)
# alpha: choice of weight function Gaussian-Laguerre quadrature (result is robust, control the "decay")
# Complexity: O(n_node*event*r_max^2) ~ O(n_node*event^3)
pred3.nnwithin.single.event <- function(
		X_event_T, X_event_Z, T_event, Z_event,
		evt_ind, within_dist=0.5,  
		para_est, longlat, r_max=5,
		in_lower=0, in_upper=1e7, no_subdivision=1e3,
		fun=geosphere::distGeo, pred_density=FALSE, pred_seq=seq(1,1e5, by=1), progress=TRUE) {


pair_to_matrix <- function(pairwise_distance, n) {
pairwise_matrix <- matrix(0, nrow=n, ncol=n)
pairwise_matrix[(pairwise_distance$ind2 - 1)*n + pairwise_distance$ind1] <- pairwise_distance$distance
pairwise_matrix + t(pairwise_matrix) + diag(Inf,n)
}



n <- NROW(X_event_T)

if (length(unique(evt_ind)) != 1) { 
	cat("Not single event data is supported!")
	return(rep(NA, n))
}


NCOL_T <- NCOL(X_event_T)
NCOL_Z <- NCOL(X_event_Z)


beta.wei <- para_est[1:NCOL_T]
scale.wei <- exp(para_est[NCOL_T+1])
scale.inv <- 1/scale.wei

beta.gamma <- para_est[(NCOL_T+1+1):(NCOL_T+1+NCOL_Z)]

disp.gamma <- para_est[NCOL_T+NCOL_Z+2]

# Copula parameters
para_copula <- tail(para_est,4)
c(copula_TZ, rho, psi, kappa) %<-% para_copula

para_copula <- tail(para_copula,3)

# bivariate copula
cop <- VineCopula::BiCop(1, copula_TZ)

TZFf <- data.frame(T_event=T_event, Z_event=Z_event) %>%
			mutate(	F_T = 1 - exp(-T_event^scale.inv*exp(-scale.inv * c(X_event_T%*%beta.wei))),
					log_mu = c(X_event_Z%*%beta.gamma),
					mu = exp(log_mu),
					nu = rep(1/disp.gamma,n)) %>%
			mutate(	f_Z = dgamma(Z_event, shape = nu, scale= mu/nu),
					logf_Z = dgamma(Z_event, shape = nu, scale= mu/nu, log=TRUE),
					F_Z = pgamma(Z_event, shape = nu, scale= mu/nu),
				) %>%
			mutate( h_Z = VineCopula::BiCopHfunc1(F_T, F_Z, cop), # d H(u,v)/dv (or du?) (i.e. conditional H function, see help(VineCopula::BiCopHfunc)
					h = VineCopula::BiCopPDF(F_T, F_Z, cop) # d^2 H(u,v)/ dudv (i.e. Copula density)
				) %>%
			mutate(	F3 = h_Z, log_f3 = logf_Z + log(h)) %>%
			mutate(	f3 = exp(log_f3)) %>%
			mutate(	F3 = ifelse(F3 <= 0 | is.nan(F3), .Machine$double.xmin, F3),
					f3 = ifelse(f3 < 0 | is.nan(f3), .Machine$double.xmin, f3)
				)

result_marginal_mean <-  result_bivariate_cop <- result <- unname(unlist(TZFf$mu))


if (r_max == 0 || n == 1) {
cat("r_max equals 0 or single obs event, just use marginal distribution!\n")

} else {

if (pred_density) {
	len_seq <- length(pred_seq)
	pred_density_matrix <- matrix(NA, nrow=n, ncol=len_seq)
}

# loop here
for (i in 1:n) {

X_event_T_tmp <- X_event_T[i,] ; X_event_Z_tmp <- X_event_Z[i,]
T_event_tmp <- T_event[i]
F_T <- 1 - exp(-T_event_tmp^scale.inv*exp(-scale.inv * c(X_event_T_tmp%*%beta.wei)))
log_mu <- c(X_event_Z_tmp%*%beta.gamma) ; mu <- exp(log_mu)
nu <- 1/disp.gamma
	
result_bivariate_cop[i] <- integrate(Z_times_new_f3, 
	lower=in_lower, upper=in_upper, subdivisions=no_subdivision,
	F_T=F_T, nu=nu, mu=mu, cop=cop)$value	

pairwise_distance <- data.frame(ind1=i, ind2=(1:n)[-i]) %>%
			mutate(distance=fun(longlat[ind1,], longlat[ind2,])/1e3) %>% 
				filter( distance <= within_dist)

n_claim_within <- NROW(pairwise_distance) # how many record within with_dist

r0 <- min(r_max, n_claim_within) + 1 # number of terms in numerator
r1 <- max(r0 - 1, 0) # number of term in denominator

	
if (NROW(pairwise_distance) > 0) {
	neighbor_within <- unname(unlist(pairwise_distance %>% 
				slice_min(distance, n=r1, with_ties=FALSE) %>% 
				select(ind2)))
				
	ind_tmp <- c(i, neighbor_within)

	if (length(neighbor_within) >=2)  {
		add_ind <- combn(neighbor_within, 2)
	
	# add more rows to pairwise_distance
	pairwise_distance <- bind_rows(
		pairwise_distance,
		data.frame(ind1=add_ind[1,], ind2=add_ind[2,]) %>%
			mutate(distance=fun(longlat[ind1,], longlat[ind2,])/1e3)
	)
	}
	
	R_matrix <- R_calculation(ind_tmp, para_copula, pairwise_distance)
	
	F_matrix <-  TZFf$F3[ind_tmp]
	

	# E(Z) = \int_0^inf z*f(z|...)dz ~ \sum_z* w(z*) g(z*) by Gaussian-Laguerre quadrature, but FAIL numerically!!!!
	# Use truncated domain and hard integration here!
	result[i] <- integrate(Z_times_new_f_Z, lower=in_lower, upper=in_upper, subdivisions=no_subdivision,
		F_T=F_T, nu=nu, mu=mu, cop=cop, F_matrix=F_matrix, R_matrix=R_matrix,r1=r1)$value	
	
	
	
	if (pred_density)  {
		pred_density_matrix[i,] <- new_f_Z(pred_seq, F_T, nu, mu, cop, F_matrix, R_matrix, r1) 
	}
} else {
	if (pred_density)  {
		# use marginal for density (only bivaiate copula with T, no Gaussian copula)
		mu <- exp(X_event_Z[i,] %*%beta.gamma) ; nu <- 1/disp.gamma
		F_T <- rep(1 - exp(-T_event[i]^scale.inv*exp(-scale.inv * c(X_event_T[i,]%*%beta.wei))), len_seq)
		F_Z = pgamma(pred_seq, shape = nu, scale= mu/nu)
		
		
		logf_Z = dgamma(Z_event, shape = nu, scale= mu/nu, log=TRUE)
		h <- VineCopula::BiCopPDF(F_T, F_Z, cop)
		
		f3 <- exp(logf_Z + log(h))
		pred_density_matrix[i,] <- ifelse(f3 < 0 | is.nan(f3), .Machine$double.xmin, f3)
	}
}

if (progress) cat("r0, r1, Z, E(Z), E(Z|T), E(Z|T...):", r0, r1, 
		round(c(Z_event[i], result_marginal_mean[i], result_bivariate_cop[i],
			result[i]),4), i,"out of", n, "\n")
} 


}

obj_pred <- data.frame(pred_mean_marginal=result_marginal_mean, pred_mean_bivariate_copula=result_bivariate_cop,
	pred_mean_copula=result)
	
if (pred_density) {
	obj_pred <- list(predicted_value=obj_pred, 
		pred_density_list=list(pred_density_matrix=pred_density_matrix, pred_seq=pred_seq))
}
	
obj_pred
}

# assess average of continuous ranked probability score (CRPS) distance
# reduce to Cramér–von Mises criterion(?)
predictive_assess <- function(predicted, observed) {

mspe <- mean((predicted- observed)^2)
mape <- mean(abs((predicted- observed)/observed))

# Average CRPS by CVM distance
D_obs <- distr::DiscreteDistribution(supp=observed)
CRPS <- CvMDist(predicted, D_obs)

list(MSPE=mspe, MAPE=mape,CRPS=CRPS)
}

