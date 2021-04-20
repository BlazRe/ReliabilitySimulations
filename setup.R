# General options
options(scipen = 999, digits = 4, max.print=99999)
Sys.setenv(TMPDIR="C:/Rtmp")
Sys.setlocale("LC_ALL","Croatian")

# Load required packages
required_packages<-c("psy", "psych", "Lambda4", "purrr", "GPArotation")
lapply(required_packages, require, character.only = TRUE)

# Calculate the different reliability coefficients - safe omega total
mu2 <- function(x) {mu2 = psych::tenberge(x)$mu2}
calc_relaibilities_safe <- function(x, decimals=4) {
  alpha <- psy::cronbach(x)$alpha
  lambda2 <- Lambda4::lambda2(x)
  lambda4<-Lambda4::cov.lambda4(x, method = "Osburn")$l4
  mu2<-mu2(x)
  glb<-psych::glb.fa(x)$glb
  glb_a<-psych::glb.algebraic(x)$glb
  omegasafe<-possibly(Lambda4::omega.tot, otherwise = NA) 
  omega<-omegasafe(x)
  rel<-c(alpha, lambda2, lambda4, mu2, glb, glb_a, omega)
  names(rel)<-c("alpha", "lambda2", "lambda4", "mu2", "glb", "glb_a", "omega")
  rel<-round(rel, decimals)
  return(rel)
}

# Main function to run the simulation
run_sim<-function(conditions, iqr_cutoffs, replications){
  # Create lists for data storage
  sim<-vector(mode = "list", length = nrow(conditions))
  # For each condition
  for (i in 1:nrow(conditions)){
    mat<-matrix(rep(NA, replications*7), ncol=7)
    colnames(mat)<-coeffs
    # Extract info for convenience
    k<-conditions[i,"k"]
    alpha_pop<-conditions[i,"alpha_pop"]
    dist<-conditions[i,"dist"]
    categ<-conditions[i,"categ"]
    N<-conditions[i,"N"]
    # For each replication
    for (j in 1:replications){
      # Make truescores
      xt<-rnorm(N,0,1)
      pop_t<-truescoredf(xt,k)
      # Add errors
      pop_e<-errorsdf_ort(k,alpha_pop,N)
      pop_x<-pop_t+pop_e
      # apply distribution shape
      pop_x<-reshape(pop_x,dist)
      # Replace outliars
      # pop_x<-scale(pop_x)
      iqr_cutoff<-iqr_cutoffs[dist]
      ro<-rep_outliars(pop_x, iqr_cutoff)
      pop_x<-ro[[1]]
      # divide continuous variables into categories
      cat_x<-categorize(pop_x, categ)
      # calculate reliabilities
      mat[j,]<-calc_relaibilities_safe(cat_x)
    }
    # save in a list
    sim[[i]]<-mat
    print(i)
  }
  l<-list(sim=sim,conditions=conditions)
  return(l)
}

# Helper functions
# Headrick polynomial transform to introduce skew and kurtosis
transform_6params <- function(x,p0,p1,p2,p3,p4,p5){
  x<-scale(x)
  y <- p0 + p1*x + p2*x^2 + p3*x^3 + p4*x^4 + p5*x^5
  y<-scale(y)
  return(y)
}

# Reshape the distribution using specified parameters
reshape<-function(df,dist){
  dft<-df
  dft[]<-NA
  if (dist=="norm"){
    params<-list(p0=0,p1=1,p2=0,p3=0,p4=0,p5=0)
  } else if (dist=="asim1"){
    params<-list(p0=0,p1=-1,p2=0.15,p3=0,p4=0,p5=0)
  } else if (dist=="asim2"){
    params<-list(p0=0,p1=-1,p2=0.4,p3=0.02,p4=0,p5=0)
  } else if (dist=="plat"){
    params<-list(p0=0,p1=1,p2=0,p3=-0.055,p4=0,p5=0)
  } else if (dist=="lept"){
    params<-list(p0=0,p1=1,p2=0,p3=0.5,p4=0,p5=0)
  }
  
  for (i in 1:ncol(df)){
    
    dft[,i]<-do.call(transform_6params, c(list(x=df[,i]),params))
  }
  return(dft)
}

# Create data frame of true scores for different number of items
truescoredf<-function(xt, k){
  v<-rep(xt, k)
  df<-data.frame(matrix(v,ncol=k))
  return(df)
}

### Spearman-Brown solve for reliability of one item based on desired alpha_pop and k
rel_single<-function(k, alpha_pop){
  rs<--alpha_pop/((alpha_pop*k)-alpha_pop-k)
  return(rs)
}

### Calculate standard deviation of the error to be added to truescore
errorize<-function(rel_s){
  var_err<-(1/rel_s)-1
  sd_err<-sqrt(var_err)
  return(sd_err)
}
### Generate errors for one item
errors_one<-function(k,alpha_pop,pop.size){
  rs<-rel_single(k,alpha_pop)
  sde<-errorize(rs)
  err<-rnorm(pop.size,0,sd=sde)
  return(err)
}

# Generate data frame of errors for different number of items
errorsdf<-function(k,alpha_pop,pop.size){
  edf<-as.data.frame(replicate(k,errors_one(k,alpha_pop,pop.size)))
  return(edf)
}

# Generate data frame of orthogonal errors for different number of items
errorsdf_ort<-function(k,alpha_pop,pop.size){
  rs<-rel_single(k,alpha_pop)
  sde<-errorize(rs)
  X<-replicate(k,rnorm(pop.size))
  X.ort<-principal(X, nfactors=k, scores=TRUE)
  X.scores<-X.ort$scores
  edf<-scale(X.scores)*sde
  return(edf)
}

# Replace outliars with border values based on IQR
rep_outliars<-function(df, iqr_cutoff=1.5){
  lcount<-rep(NA,ncol(df))
  hcount<-rep(NA,ncol(df))
  for (i in 1:ncol(df)){
    q3<-quantile(df[,i], 3/4)
    q1<-quantile(df[,i], 1/4)
    iqr<-q3-q1
    lw<-q1-(iqr_cutoff*iqr)
    hw<-q3+(iqr_cutoff*iqr)
    lcount[i]<-sum(df[,i]<lw)
    hcount[i]<-sum(df[,i]>hw)
    df[df[,i]<lw,i]<-lw
    df[df[,i]>hw,i]<-hw
  }
  lr<-list(df=df,lcount=lcount,hcount=hcount)
  return(lr)
}

# Cut continuous variables in a number of discreet categories
categorize<-function(df,categ){
  dfc<-df
  dfc[]<-NA
  for (i in 1:ncol(df)){
    dfc[,i]<-as.numeric(cut(df[,i],breaks=categ))
  }
  return(dfc)
}
