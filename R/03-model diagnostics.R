#*********************
# Checks GOF (top) and convergence (bottom)
#********************

#*********************
# Check GOF ----
#********************
source("./R/plotdiag-function.R")
# create file names
floc <- "G:\\Projects_ARCHIVED\\nw-wrs\\noCali\\"
fldr <- c("poisson\\", "negbin\\", "zip\\")
sp <- c("AMKE", "BAEA", "COHA", "FEHA", "GOEA", "NOHA", "PRFA", "RLHA", "RSHA", "RTHA", "WTKI")
sffx <- c("-wrscbc.rdata", "-wrscbc-nb.rdata", "-wrscbc-zip.rdata")
flnms <- list()
for (i in 1:length(sp)) {
  flnms[[i]] <- paste0(floc, fldr, sp[i], sffx )
}

i <- 1 # 
sp[i]
load(flnms[[i]][1])
plot.diag(out, ratio=T, lab=paste0(sp[i], "poisson") )
load(flnms[[i]][2])
plot.diag(out, ratio=F, lab=paste0(sp[i], "negbin") )
load(flnms[[i]][3])
plot.diag(out, ratio=F, lab=paste0(sp[i], "zip") )

i <- 2 # 
sp[i]
load(flnms[[i]][1])
plot.diag(out, ratio=T, lab=paste0(sp[i], "poisson") )
load(flnms[[i]][2])
plot.diag(out, ratio=F, lab=paste0(sp[i], "negbin") )
load(flnms[[i]][3])
plot.diag(out, ratio=F, lab=paste0(sp[i], "zip") )

i <- 3 # 
sp[i]
load(flnms[[i]][1])
plot.diag(out, ratio=T, lab=paste0(sp[i], "poisson") )
load(flnms[[i]][2])
plot.diag(out, ratio=F, lab=paste0(sp[i], "negbin") )
load(flnms[[i]][3])
plot.diag(out, ratio=F, lab=paste0(sp[i], "zip") )

i <- 4 # 
sp[i]
load(flnms[[i]][1])
plot.diag(out, ratio=T, lab=paste0(sp[i], "poisson") )
load(flnms[[i]][2])
plot.diag(out, ratio=F, lab=paste0(sp[i], "negbin") )
load(flnms[[i]][3])
plot.diag(out, ratio=F, lab=paste0(sp[i], "zip") )

i <- 5 # 
sp[i]
load(flnms[[i]][1])
plot.diag(out, ratio=T, lab=paste0(sp[i], "poisson") )
load(flnms[[i]][2])
plot.diag(out, ratio=F, lab=paste0(sp[i], "negbin") )
load(flnms[[i]][3])
plot.diag(out, ratio=F, lab=paste0(sp[i], "zip") )

i <- 6 # 
sp[i]
load(flnms[[i]][1])
plot.diag(out, ratio=T, lab=paste0(sp[i], "poisson") )
load(flnms[[i]][2])
plot.diag(out, ratio=F, lab=paste0(sp[i], "negbin") )
load(flnms[[i]][3])
plot.diag(out, ratio=F, lab=paste0(sp[i], "zip") )

i <- 7 # 
sp[i]
load(flnms[[i]][1])
plot.diag(out, ratio=T, lab=paste0(sp[i], "poisson") )
load(flnms[[i]][2])
plot.diag(out, ratio=F, lab=paste0(sp[i], "negbin") )
load(flnms[[i]][3])
plot.diag(out, ratio=F, lab=paste0(sp[i], "zip") )

i <- 8 # 
sp[i]
load(flnms[[i]][1])
plot.diag(out, ratio=T, lab=paste0(sp[i], "poisson") )
load(flnms[[i]][2])
plot.diag(out, ratio=F, lab=paste0(sp[i], "negbin") )
load(flnms[[i]][3])
plot.diag(out, ratio=F, lab=paste0(sp[i], "zip") )

i <- 9 # 
sp[i]
load(flnms[[i]][1])
plot.diag(out, ratio=T, lab=paste0(sp[i], "poisson") )
load(flnms[[i]][2])
plot.diag(out, ratio=F, lab=paste0(sp[i], "negbin") )
load(flnms[[i]][3])
plot.diag(out, ratio=F, lab=paste0(sp[i], "zip") )

i <- 10 # 
sp[i]
load(flnms[[i]][1])
plot.diag(out, ratio=T, lab=paste0(sp[i], "poisson") )
load(flnms[[i]][2])
plot.diag(out, ratio=F, lab=paste0(sp[i], "negbin") )
load(flnms[[i]][3])
plot.diag(out, ratio=F, lab=paste0(sp[i], "zip") )

i <- 11 # 
sp[i]
load(flnms[[i]][1])
plot.diag(out, ratio=T, lab=paste0(sp[i], "poisson") )
load(flnms[[i]][2])
plot.diag(out, ratio=F, lab=paste0(sp[i], "negbin") )
load(flnms[[i]][3])
plot.diag(out, ratio=F, lab=paste0(sp[i], "zip") )

#*********************
# Check convergence ----
#********************
library(jagsUI)
fls <- list.files("C:\\Users\\rolek.brian\\Documents\\Projects\\northwest-wrs\\outputs\\withGOF\\final",
           full.names = TRUE)
paramszi <- c("r.mu", "lam.mu","sig.r.mu",  # highest level params  
            #"r.strat", 
            "sig.proc.strat", # shared params
            "beta", "sig.obs",  
            "psi", "psi2")
paramsnb <- c("r.mu", "lam.mu","sig.r.mu",  # highest level params  
            #"r.strat", 
            "sig.proc.strat", # shared params
            "beta", "sig.obs",  
            "rr", "rr2")
fls[[1]]
load(fls[[1]])
traceplot(out, paramszi)

fls[[2]]
load(fls[[2]])
traceplot(out, paramszi)

fls[[3]]
load(fls[[3]])
traceplot(out, paramsnb)

fls[[4]]
load(fls[[4]])
traceplot(out, paramszi)

fls[[5]]
load(fls[[5]])
traceplot(out, paramsnb)

fls[[6]]
load(fls[[6]])
traceplot(out, paramszi)

fls[[7]]
load(fls[[7]])
traceplot(out, paramsnb)

fls[[8]]
load(fls[[8]])
traceplot(out, paramsnb)

fls[[9]]
load(fls[[9]])
traceplot(out, paramsnb)

fls[[10]]
load(fls[[10]])
traceplot(out, paramszi)

fls[[11]]
load(fls[[11]])
traceplot(out, paramsnb)

# Make a table for GOF bayesian p values
bp.tab <- data.frame("Species"=NA, "Data set"=rep(c("WRS", "CBC"), 11), 
                     "Poisson"=NA, "negative.binomial"=NA, "zero.inflated.Poisson"=NA)
j <- 1
for (i in 1:11){
bp.tab[j:(j+1), "Species"] <- sp[i]
load(flnms[[i]][1])
bp.tab[j:(j+1), "Poisson"] <- bp.func(out)
load(flnms[[i]][2])
bp.tab[j:(j+1), "negative.binomial"] <- bp.func(out)
load(flnms[[i]][3])
bp.tab[j:(j+1), "zero.inflated.Poisson"] <- bp.func(out)
j <- j+2
}

write.csv()