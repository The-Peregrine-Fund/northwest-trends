## ---- poisson ----------
library (jagsUI)
#load("/scratch/brolek/northwest-wrs/data/data.rdata")
load("./data/data.rdata")

#sink("/scratch/brolek/northwest-wrs/JAGSmodel.txt")
sink("./R/JAGSmodel.txt")
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
     
     # CBC
     for(l in 1:ncounts2) {
     lN.est2[l] <- lmu.N2[time2[l],rt2[l]] + log(dist2[l])
     logtheta2[l] ~ dnorm(lN.est2[l], 1/(sig.obs[2]*sig.obs[2]))
     log(theta2[l]) <- logtheta2[l]
     count2[l] ~ dpois(theta2[l])
     }  # l
     
     #### 1st year and dynamics ######
     for(i in 1:nroutes) {   #### priors for first year and strata ######
     lmn1st[i] <- log( mn1st[i]+0.0001 ) # add a small constant to prevent log(zero) which will result in failure to run
     lmu.N[ frst[i], i ] ~ dnorm( lmn1st[i], 1/(10*10) ) # priors for 1st yr abundance are roughly near observed values
     for( t in frst[i]:(ntime-1) ) { # dynamics
     lmu.N[t+1,i] <- lmu.N[t,i] + r[ t, i ] # WRS
     r[t,i] ~ dnorm( r.strat[t, strat1[i] ], 1/(sig.proc.strat*sig.proc.strat) ) 
     } } #t     
     
     for(j in 1:nroutes2) {  #### priors for first year and strata ######
     lmn1st2[j] <- log( mn1st2[j]+0.0001 ) # add a small constant to prevent log(zero) which will result in failure to run
     lmu.N2[ frst2[j], j ] ~ dnorm( lmn1st2[j], 1/(10*10) ) # priors for 1st yr abundance are roughly near observed values
     for( t2 in frst2[j]:(ntime2-1) ) { # dynamics
     lmu.N2[t2+1,j] <- lmu.N2[t2,j] + r2[ t2, j] # CBC
     r2[t2,j] ~ dnorm( r.strat[t2, strat2[j] ], 1/(sig.proc.strat*sig.proc.strat) ) 
     } } # t2,j
     
     #####################
     # Derived parameters
     # for model diagnostics
     ##################### 
     sig.ratio[1] <- sig.obs[1]/sig.proc.strat
     sig.ratio[2] <- sig.obs[2]/sig.proc.strat

   ###################
    # Assess GOF of the state-space models for counts
    # Step 1: Compute statistic for observed data
    # Step 2: Use discrepancy measure: mean absolute error
    # Step 3: Use test statistic: number of turns
    ###################
    # GOF for WRS model
    for(k in 1:ncounts) {
    c.exp.wrs[k] <- theta[k] # expected counts adult breeder
    c.obs.wrs[k] <- count[k] # observed counts
    c.rep.wrs[k] ~ dpois(theta[k]) # expected counts
    # Compute fit statistics, Mean absolute error
    dssm.obs.wrs[k] <- abs( ( (c.obs.wrs[k]) - (c.exp.wrs[k]) ) / (c.obs.wrs[k]+0.001)  )
    dssm.rep.wrs[k] <- abs( ( (c.rep.wrs[k]) - (c.exp.wrs[k]) ) / (c.rep.wrs[k]+0.001) )
    } # k
    dmape.obs.wrs <- sum(dssm.obs.wrs)
    dmape.rep.wrs <- sum(dssm.rep.wrs)
    # variance-mean ratio wrs
    tvm.rep.wrs <- sd(c.rep.wrs)^2/mean(c.rep.wrs)
    tvm.obs.wrs <- sd(count)^2/mean(count)
    
    # GOF for CBC model
    for(l in 1:ncounts2) {
    c.exp.cbc[l] <- theta2[l] # expected counts adult breeder
    c.obs.cbc[l] <- count2[l]
    c.rep.cbc[l] ~ dpois(theta2[l]) # expected counts
    # Compute fit statistics, Mean absolute error
    dssm.obs.cbc[l] <- abs( ( (c.obs.cbc[l]) - (c.exp.cbc[l]) ) / (c.obs.cbc[l]+0.001)  )
    dssm.rep.cbc[l] <- abs( ( (c.rep.cbc[l]) - (c.exp.cbc[l]) ) / (c.rep.cbc[l]+0.001) )
    } # l
    dmape.obs.cbc <- sum(dssm.obs.cbc)
    dmape.rep.cbc <- sum(dssm.rep.cbc)
    # variance-mean ratio cbc
    tvm.rep.cbc <- sd(c.rep.cbc)^2/mean(c.rep.cbc)
    tvm.obs.cbc <- sd(count2)^2/mean(count2)
        
     } # model end
     ", fill=TRUE)
sink()

params <- c("r.mu", "lam.mu","sig.r.mu",  # highest level params  
            "r.strat", "sig.proc.strat", # shared params
            "beta", "sig.obs", "sig.ratio", 
            "r", "r2", 
            "lmu.N", "lmu.N2",
            "dmape.obs.wrs", "dmape.rep.wrs", "dmape.obs.cbc", "dmape.rep.cbc",
            "tvm.rep.wrs", "tvm.obs.wrs", "tvm.rep.cbc", "tvm.obs.cbc"
)

nc <- 4; na <- 10000; ni <- 500000;  nb <- 400000; nt <- 400  
#nc <- 4; na <- 10000; ni <- 100000;  nb <- 50000; nt <- 50 
#nc <- 3; na <- 100; ni <- 200;  nb <- 100; nt <- 1 

for (i in 1:11){ # loop through species
  inits <- function(){
    list(
      theta=dall[[i]]$count,
      theta2=dall[[i]]$count2
    )}
  
out <- jags(#model.file="/scratch/brolek/northwest-wrs/JAGSmodel.txt", 
            model.file= "./R/JAGSmodel.txt",          
            dall[[i]], parameters.to.save = params,
               n.chains = nc, n.iter = ni,
               n.burnin=nb, n.thin=nt,
               n.adapt = na, parallel=T,
            codaOnly = c("r", "lmu.N", "lmu.N2"))

# save(out, file=paste("/scratch/brolek/northwest-wrs/outputs/",
#                      names(dall)[i],"-wrscbc",
#                      ".rdata", sep="" ))
}
