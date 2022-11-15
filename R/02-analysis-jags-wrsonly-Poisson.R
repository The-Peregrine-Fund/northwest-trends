library (jagsUI)
load("/scratch/brolek/northwest-wrs/data/data.rdata")
# load("./data/data.rdata")

sink("/scratch/brolek/northwest-wrs/JAGSmodel-wrsonly.txt")
# sink("./R/JAGSmodel.txt")
cat ("  
     model {
     ################################
     # indices: k=survey, i=route, 
     #          t=survey year e.g. 04-05
     # params: lmu.N=log(abundance)
     #         r= population growth rate
     ################################         
     #### Priors ###### 
psd <- 2
r.mu ~ dnorm( 0, 1/(psd*psd) )
lam.mu <- exp(r.mu)
sig.r.mu ~ dnorm(0, 1/(psd*psd))T(0,)
sig.proc.strat ~ dnorm(0, 1/(psd*psd))T(0,)
beta ~ dnorm(0, 1/(10*10) )
sig.obs[1] ~ dnorm(0, 1/(psd*psd) )T(0,)
sig.obs[2] ~ dnorm(0, 1/(psd*psd) )T(0,)

for (s in 1:nstrata){
  for( t3 in 1:(ntime2-1) ) {
     r.strat[t3, s] ~ dnorm(r.mu, 1/(sig.r.mu*sig.r.mu))
     lam.strat[t3, s] <- exp(r.strat[t3, s])
}} # t3 s

     #### data model ###### 
     # WRS
     for(k in 1:ncounts) {
      lN.est[k] <- lmu.N[time[k],rt[k]] + log(dist[k]) + beta*sp[k]
      logtheta[k] ~ dnorm (lN.est[k], 1/(sig.obs[1]*sig.obs[1]))
      log(theta[k]) <- logtheta[k]
      count[k] ~ dpois(theta[k])
     }  # k
     
     # # CBC
     # for(l in 1:ncounts2) {
     #  lN.est2[l] <- lmu.N2[time2[l],rt2[l]] + log(dist2[l])
     #  logtheta2[l] ~ dnorm (lN.est2[l], 1/(sig.obs[2]*sig.obs[2]))
     #  log(theta2[l]) <- logtheta2[l]
     #  count2[l] ~ dpois(theta2[l])
     # }  # l

     ### 1st year and dynamics ######
     for(i in 1:nroutes) {   #### priors for first year and strata ######
      lmn1st[i] <- log( mn1st[i]+0.0001 ) # add a small constant to prevent log(zero) which will result in failure to run
      lmu.N[ frst[i], i ] ~ dnorm( lmn1st[i], 1/(10*10) ) # priors for 1st yr abundance are roughly near observed values
      for( t in frst[i]:(ntime-1) ) { # dynamics
        lmu.N[t+1,i] <- lmu.N[t,i] + r[ t, i ] # WRS
        r[t,i] ~ dnorm( r.strat[t, strat1[i] ], 1/(sig.proc.strat*sig.proc.strat) )
     } } #t
     
    # for(j in 1:nroutes2) {  #### priors for first year and strata ######
    #   lmn1st2[j] <- log( mn1st2[j]+0.0001 ) # add a small constant to prevent log(zero) which will result in failure to run
    #   lmu.N2[ frst2[j], j ] ~ dnorm( lmn1st2[j], 1/(10*10) ) # priors for 1st yr abundance are roughly near observed values
    #  for( t2 in frst2[j]:(ntime2-1) ) { # dynamics
    #   lmu.N2[t2+1,j] <- lmu.N2[t2,j] + r2[ t2, j] # CBC
    #   r2[t2,j] ~ dnorm( r.strat[t2, strat2[j] ], 1/(sig.proc.strat*sig.proc.strat) )
    # } } # t2,j
   
     #####################
     # derived parameters
     ##################### 
     sig.ratio[1] <- sig.obs[1]/sig.proc.strat
     sig.ratio[2] <- sig.obs[2]/sig.proc.strat
     
}
     ", fill=TRUE)
sink()

params <- c("r.mu", "lam.mu","sig.r.mu", # highest level params  
            "r.strat", "sig.proc.strat", # shared params
            "beta", "sig.obs", "sig.ratio", 
            "r", "r2", 
            "lmu.N", "lmu.N2" 
)

nc <- 4; na <- 10000; ni <- 500000;  nb <- 400000; nt <- 400  
# nc <- 4; na <- 10000; ni <- 100000;  nb <- 50000; nt <- 50 
#nc <- 3; na <- 100; ni <- 200;  nb <- 100; nt <- 1 
for (i in 1:11){
  inits <- function(){
    list(
      mn1st=dall[[i]]$mn1st+0.0001,
      mn1st2=dall[[i]]$mn1st2+0.0001
    )}
  
out <- jags(model.file="/scratch/brolek/northwest-wrs/JAGSmodel-wrsonly.txt", 
            #model.file= "./R/JAGSmodel-cbconly.txt",          
            #inits = inits,
            dall[[i]], parameters.to.save = params,
               n.chains = nc, n.iter = ni,
               n.burnin=nb, n.thin=nt,
               n.adapt = na, parallel=T,
            codaOnly = c("r", "lmu.N", "lmu.N2"))

save(out, file=paste("/scratch/brolek/northwest-wrs/outputs/",
                     names(dall)[i],"-wrsonly",
                     ".rdata", sep="" ))
}