# Simulation parameters
replications<-1000
coeffs<-c("alpha", "lambda2", "lambda4", "mu2", "glb", "glb_a", "omega")
iqr_cutoffs<-c("norm" = 2, "asim1" = 2, "asim2" = 2.5, "plat" = 0.7, "lept" = 5)

# Conditions
k<-c(4,8,12,20)
alpha_pop<-c(0.5,0.7,0.9)
dist<-c("norm","asim1","asim2","plat","lept")
categ<-c(3,5,7,9)
N<-c(50,100,200,500,1000)

# Combine conditions into the control object
cond_list<-list(k,categ,N,dist,alpha_pop)
conditions<-expand.grid(cond_list)
names(conditions)<-c("k","categ","N","dist","alpha_pop")

# Add single description for every condition
desc_cond<-apply(expand.grid(cond_list), 1, paste, collapse = "_")
desc_cond<-gsub(' ','',desc_cond)
desc_cond<-factor(desc_cond, levels = desc_cond)
conditions<-data.frame(conditions,desc_cond)

# Set seed
set.seed(100)

# Initiate simulation and estimate the elapsed time
start.time <- Sys.time()
sim100<-run_sim(conditions, iqr_cutoffs, replications)
end.time <- Sys.time()
elapsed.time <- end.time - start.time
elapsed.time

# Save result
save(sim100, file="sim100.Rdata")
