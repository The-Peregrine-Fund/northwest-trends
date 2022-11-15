library(jagsUI)
library (viridis)
library(MCMCvis)
library(tidybayes)
library(HDInterval)
library(reshape2)
library(ggplot2)
library(dplyr)
load("C:\\Users\\rolek.brian\\Documents\\GitHub\\northwest-wrs-trends\\data\\data.rdata")

#######################
# Convergence checks
#######################
flpath <- list.files("G:\\nw-wrs\\outputs7\\both\\", full.names=T)
flnm <- list.files("G:\\nw-wrs\\outputs7\\both\\")
conv <- data.frame(spp=rep(NA,length(flpath)), r.mu=NA, sig.r.mu=NA, 
                   proc.mu=NA, sig.proc.mu=NA, 
                   sig.obs1=NA, sig.obs2=NA)

for (i in 1:length(flpath)){
  load(flpath[i])
  conv$spp[i] <- flnm[i]
  conv$r.mu[i] <- out$summary["r.mu", "Rhat"]
  conv$sig.r.mu[i] <- out$summary["sig.r.mu", "Rhat"]
  conv$sig.proc.mu[i] <- out$summary["sig.proc.strat", "Rhat"]
  # conv$proc.mu[i] <- out$summary["r.strat", "Rhat"]
  # conv$sig.proc.mu[i] <- out$summary["sig.proc.strat", "Rhat"]
  conv$sig.obs1[i] <- out$summary["sig.obs[1]", "Rhat"]
  conv$sig.obs2[i] <- out$summary["sig.obs[2]", "Rhat"]
  par(mfrow=c(2,3), oma=c(0,0,1,0))
  print(flpath[i])
  traceplot(out, "r.mu")
  title(main=substr( flnm[i], 1, 4), outer=TRUE)
  traceplot(out, "sig.r.mu")
  traceplot(out, "sig.proc.strat")
  traceplot(out, "sig.obs[1]")
  traceplot(out, "sig.obs[2]")
}

###############
# Data summaries
################
# number of surveys in CBC and WRS surveys per year
nu.wrs <- table(all.dat$surv_yr[!is.na(all.dat$total)])
wrs.yr <- paste(2004:2019, "-", 2005:2020, sep="")
wrsdf <- data.frame(yr=wrs.yr, wrs=nu.wrs)[,-2]
nu.cbc <- table(csp.sp@data$Count_yr[!is.na(csp.sp@data$RTHA)])
cbc.yr <- paste(2004:2019, "-", 2005:2020, sep="")
cbcdf <- data.frame(yr=cbc.yr, cbc=nu.cbc)[,-2]
both <- merge(cbcdf, wrsdf, by="yr", all.x=T)
colnames(both) <- c("Year", "No. CBC", "No. WRS")
write.csv(both, "C:\\Users\\rolek.brian\\Documents\\Projects\\northwest-wrs\\docs\\surveys-per-yr.csv" )

################
# NEW Plots
################
flnms <- list.files("G:\\nw-wrs\\outputs7\\both\\", full.names=T)
fl <- list.files("G:\\nw-wrs\\outputs7\\both\\")
sppnms <- substr( fl, 1, 4)

# calculate model estimates for each
# species in each strata, and 
# summarize data
strata_tab <- strata_tab[order(strata_tab$strata_num), ]
colnames(strata_tab)[6] <- "BAEA"

cts.wrs <- cts.cbc <- cd <- wd <- list()
for (i in 1:length(flnms)){
  cd[[i]] <- wd[[i]] <- list()
  # data manipulation
  spdat <- dall[[ sppnms[i] ]]
  str_tab <- strata_tab[order(strata_tab[, sppnms[i]]), ] # sort strata for each sp
  sub.wrs <- !is.na(spdat$count) # subset to remove NAs
  wd[[i]]$count <- spdat$count[sub.wrs] 
  wd[[i]]$dist <- spdat$dist[sub.wrs]
  wd[[i]]$rt <- spdat$rt[sub.wrs]
  tstrat1 <- spdat$strat1
  tstrat1 <- tstrat1[wd[[i]]$rt]
  wd[[i]]$strat <- str_tab$strata_num[ tstrat1 ] # reassign strata for some sp
  wd[[i]]$time <- spdat$time[sub.wrs]
  sub.cbc <- !is.na(spdat$count2) # subset to remove NAs
  cd[[i]]$count <- spdat$count2[sub.cbc] 
  cd[[i]]$dist <- spdat$dist2[sub.cbc]
  cd[[i]]$rt <- spdat$rt2[sub.cbc]
  tstrat2 <- spdat$strat2
  tstrat2 <- tstrat2[cd[[i]]$rt]
  cd[[i]]$strat <- str_tab$strata_num[ tstrat2 ] # reassign strata for some sp
  cd[[i]]$time <- spdat$time2[sub.cbc]
  # calculate mean count/100km in each strata from raw data 
  temp.wrs <- with(wd[[i]], tapply( count/dist, 
          list(time, strat), mean) )
  arw <- array(NA, dim=c(16, 10))# add NAs for yrs not surveyed by WRS
  for (j in 1:ncol(temp.wrs)){
    cnum <- as.numeric(colnames(temp.wrs)[j])
    arw[, cnum] <- temp.wrs[,j]
  }
  cts.wrs[[i]] <- arw
  
  temp.cbc <- with(cd[[i]], tapply( count/dist, 
                                    list(time, strat), 
                                    mean, na.rm=T) )
  arc <- array(NA, dim=c(16, 10))
  for (j in 1:ncol(temp.wrs)){
    cnum <- as.numeric(colnames(temp.cbc)[j])
    arc[, cnum] <- temp.cbc[,j]
  }
  cts.cbc[[i]] <- arc
}
# extract abundance model estimates for each strata  
N.cbc.all <- N.wrs.all <- N.wrs <- N.wrs.hdis <- N.cbc <- N.cbc.hdis <- list()
w.rts <- all.dat[!duplicated(all.dat$route_num), c("route_num", "strata_num")]
c.rts <- csp.sp[!duplicated(csp.sp$route_num2 ), c("route_num2", "strata_num2")]
for (i in 1:length(flnms)){  
  load(flnms[[i]])
  # average abundance over all routes 
  # for each data set
  #N.cbc.all <- N.wrs.all <- N.wrs.hdis <- N.cbc.hdis <- list()
  wrs.sims <- apply(out$sims.list$lmu.N, c(1,2) , mean, na.rm=T )
  N.wrs.all[[i]] <- exp(apply(wrs.sims, 2, mean, na.rm=T ))
  N.wrs.hdis[[i]] <- exp(apply(wrs.sims, 2, hdi, na.rm=T, credMass=0.95 ))
  cbc.sims <- apply(out$sims.list$lmu.N2, c(1,2) , mean, na.rm=T )
  N.cbc.all[[i]] <- exp(apply(cbc.sims, 2, mean, na.rm=T ))
  N.cbc.hdis[[i]] <- exp(apply(cbc.sims, 2, hdi, na.rm=T, credMass=0.95 ))
  # WRS average abundance for each strata 
  mdf <- melt(out$sims.list$lmu.N)
  colnames(mdf) <- c("iter", "yr", "route_num", "value") 
  mdf$strata_num <- dall[[ sppnms[i] ]]$strat1[ mdf$route_num ]
  str_tab <- strata_tab[order(strata_tab[, sppnms[i]]), ]
  mdf$strata_num2 <- str_tab$strata_num[ mdf$strata_num ]
  temp.wrs <- with(mdf, exp(tapply(value, list(yr,strata_num2), mean, na.rm=T )))
  arw2 <- array(NA, dim=c(16, 10))
  for (j in 1:ncol(temp.wrs)){
    cnum <- as.numeric(colnames(temp.wrs)[j])
    arw2[1:16, cnum] <- temp.wrs[,j]
  }
  N.wrs[[i]] <- arw2
  # CBC average abundance for each strata
  mdf2 <- melt(out$sims.list$lmu.N2)
  colnames(mdf2) <- c("iter", "yr", "route_num", "value") 
  mdf2$strata_num <- dall[[ sppnms[i] ]]$strat2[ mdf2$route_num ]
  mdf2$strata_num2 <- str_tab$strata_num[ mdf2$strata_num ]
  temp.cbc <- with(mdf2, exp(tapply(value, list(yr,strata_num2), mean, na.rm=T )))
  arc2 <- array(NA, dim=c(16, 10))
  for (j in 1:ncol(temp.cbc)){
    cnum <- as.numeric(colnames(temp.cbc)[j])
    arc2[1:16, cnum] <- temp.cbc[,j]
  }
  N.cbc[[i]] <- arc2
}
# z score abundance in each strata for each species 
cts.wrsz <- N.wrsz <- cts.cbcz <- N.cbcz <- list()
for (i in 1:length(flnms)){
  cts.wrsz[[i]] <- N.wrsz[[i]] <- cts.cbcz[[i]] <- N.cbcz[[i]] <- array(NA, dim=dim(N.cbc[[i]]))
  for (j in 1:10){
    cts.wrsz[[i]][,j] <- ( cts.wrs[[i]][,j]-mean(cts.wrs[[i]][,j], na.rm=T) )/ sd(cts.wrs[[i]][,j], na.rm=T)
    N.wrsz[[i]][,j] <- ( N.wrs[[i]][,j]-mean(N.wrs[[i]][,j], na.rm=T) )/ sd(N.wrs[[i]][,j], na.rm=T)
    cts.cbcz[[i]][,j] <- ( cts.cbc[[i]][,j]-mean(cts.cbc[[i]][,j], na.rm=T) )/ sd(cts.cbc[[i]][,j], na.rm=T)
    N.cbcz[[i]][,j] <- ( N.cbc[[i]][,j]-mean(N.cbc[[i]][,j], na.rm=T) )/ sd(N.cbc[[i]][,j], na.rm=T)
  }}
# plot relative abundance for each strata on the natural scale
# plot CBC first
spp <- c("American \nKestrel", "Bald \nEagle", "Cooper's \nHawk", 
         "Ferruginous \nHawk", "Golden \nEagle", "Northern \nHarrier",
         "Prairie \nFalcon", "Rough-legged \nHawk", "Red-shouldered \nHawk",
         "Red-tailed \nHawk", "White-tailed \nKite")
tiff(height=6, width=10, res=300, units="in",
     "C:\\Users\\rolek.brian\\Documents\\Projects\\northwest-wrs\\docs\\figs\\abund_CBC_strata.tiff")
par(mfrow=c(3,4), mar=c(0,3,0,0), oma=c(6,5,5,5))
for (i in 1:length(flnms)){
  mn.Ns <- rbind(N.cbcz[[i]], N.wrsz[[i]])
  mn <- min(mn.Ns, na.rm=T)
  mx <- max(mn.Ns, na.rm=T)
  md <- (mx-mn)/2
  preddates <- 2004:2019
  plot(preddates, rep(NA, length(preddates)), 
       ylim=c(-2, mx*1.2),
       xlab="", ylab="",
       xaxt="n", yaxt="n")
  xdates.at <- seq(2004, 2019, by=5)
  xdates <- c(2005, rep(NA, length(xdates.at)-2), 2019)
  axis(1, at=xdates.at, labels=NA, cex.axis=2)
  axis(2, at=c(-2, 0, 2), labels=c(-2, 0, 2), cex.axis=2)
  for (j in 1:10){ lines(preddates, N.cbcz[[i]][,j], lwd=2, col=viridis(10)[j]) }
  #for (j in 1:10){ lines(preddates, N.wrsz[[i]][,j], lwd=2, col=viridis(10)[j], lty=2) } 
  title(spp[i], font.main=1, line=-4, cex.main=2)
  if (i >=8){  axis(1, at=xdates.at, labels=xdates, cex.axis=2)   }
}
plot.new()
legend(0.2, 0.85, legend=1:10, lty=1, col=viridis(10), 
       title="Strata", ncol=2, cex=1.5, bty="n",
       lwd=2, xpd=T, seg.len=0.5)
mtext(side=2, "Abundance index", outer=T, cex=2, line=1.5)
mtext(side=1, "Year", outer=T, cex=2, line=4.5)
dev.off()

# plot WRS abundance for each strata
tiff(height=6, width=10, res=300, units="in",
     "C:\\Users\\rolek.brian\\Documents\\Projects\\northwest-wrs\\docs\\figs\\abund_WRS_strata.tiff")
par(mfrow=c(3,4), mar=c(0,3,0,0), oma=c(6,5,5,5))
for (i in 1:length(flnms)){
  mn.Ns <- rbind(N.cbcz[[i]], N.wrsz[[i]])
  mn <- min(mn.Ns, na.rm=T)
  mx <- max(mn.Ns, na.rm=T)
  md <- (mx-mn)/2
  preddates <- 2004:2019
  plot(preddates, rep(NA, length(preddates)), 
       ylim=c(-2, mx*1.2),
       xlab="", ylab="",
       xaxt="n", yaxt="n")
  xdates.at <- seq(2004, 2019, by=5)
  xdates <- c(2005, rep(NA, length(xdates.at)-2), 2019)
  axis(1, at=xdates.at, labels=NA, cex.axis=2)
  axis(2, at=c(-2, 0, 2), labels=c(-2, 0, 2), cex.axis=2)
  #for (j in 1:10){ lines(preddates, N.cbcz[[i]][,j], lwd=2, col=viridis(10)[j]) }
  for (j in 1:10){ lines(preddates, N.wrsz[[i]][,j], lwd=2, col=viridis(10)[j], lty=1) } 
  title(spp[i], font.main=1, line=-4, cex.main=2)
  if (i >=8){  axis(1, at=xdates.at, labels=xdates, cex.axis=2)   }
}
plot.new()
legend(0.2, 0.85, legend=1:10, lty=1, col=viridis(10), 
       title="Strata", ncol=2, cex=1.5, bty="n",
       lwd=2, xpd=T, seg.len=0.5)
mtext(side=2, "Abundance index", outer=T, cex=2, line=1.5)
mtext(side=1, "Year", outer=T, cex=2, line=4.5)
dev.off()

#################
# Plot abundance for entire region
#################
tiff(height=6, width=10, res=300, units="in",
     "C:\\Users\\rolek.brian\\Documents\\Projects\\northwest-wrs\\docs\\figs\\abund_wrs-and-cbc.tiff")
par(mfrow=c(3,4), mar=c(0,3,0,0), oma=c(6,5,5,5))
Nc <- Nw <- Nc.hdis <- Nw.hdis <- list()
for (i in 1:length(flnms)){
  Nc[[i]] <- (N.cbc.all[[i]]-mean(N.cbc.all[[i]], na.rm=T))/sd(N.cbc.all[[i]], na.rm=T)
  Nc.hdis[[i]] <- (N.cbc.hdis[[i]]-mean(N.cbc.all[[i]], na.rm=T))/sd(N.cbc.all[[i]], na.rm=T)
  Nw[[i]] <- (N.wrs.all[[i]]-mean(N.wrs.all[[i]], na.rm=T))/sd(N.wrs.all[[i]], na.rm=T)   
  Nw.hdis[[i]] <- (N.wrs.hdis[[i]]-mean(N.wrs.all[[i]], na.rm=T))/sd(N.wrs.all[[i]], na.rm=T)       
  mn.Ns <- c(Nc[[i]], Nw[[i]])
  mn <- floor(min(mn.Ns, na.rm=T))
  mx <- ceiling(max(mn.Ns, na.rm=T))
  md <- (mx-mn)/2
  preddates <- 2004:2019
  plot(preddates, rep(NA, length(preddates)), 
       ylim=c(-2.5, 4),
       xlab="", ylab="",
       xaxt="n", yaxt="n")
  xdates.at <- seq(2004, 2019, by=5)
  xdates <- c(2005, rep(NA, length(xdates.at)-2), 2019)
  axis(1, at=xdates.at, labels=NA, cex.axis=2)
  if (i %in% c(1,5,9)){
  axis(2, at=c(-4, -2, 0, 2, 4), labels=c(NA, -2, 0, 2, NA), cex.axis=2) 
  } else{
    axis(2, at=c(-4, -2, 0, 2, 4), labels=c(NA, NA, NA, NA, NA), cex.axis=2) 
    }
  lines(preddates, Nc[[i]], lwd=2, col=viridis(2)[1]) 
  lines(preddates, Nw[[i]], lwd=2, col=viridis(2)[2])  
  title(spp[i], font.main=1, line=-4, cex.main=2)
  if (i >=8){  axis(1, at=xdates.at, labels=xdates, cex.axis=2)   }
}
plot.new()
legend(0.2, 0.6, legend=c("CBC", "WRS"), lty=1, col=viridis(2), 
       title="Survey", cex=2, bty="n",
       lwd=2, xpd=T, seg.len=0.5)
mtext(side=2, "Abundance index", outer=T, cex=2, line=1.5)
mtext(side=1, "Year", outer=T, cex=2, line=4.5)
dev.off()

#################
# Plot pop growth for entire region each year
#################
# derive r params for plotting
ryr <- ryr.hdis <- list()
for (i in 1:length(flnms)){  
  load(flnms[[i]])
  # mean r each year
  mnr.sims <- apply(out$sims.list$r.strat, c(1,2), mean, na.rm=T)
  ryr[[i]] <- apply(mnr.sims, 2, mean, na.rm=T)
  ryr.hdis[[i]] <- apply(mnr.sims, 2, hdi, na.rm=T)
}

tiff(height=6, width=10, res=300, units="in",
     "C:\\Users\\rolek.brian\\Documents\\Projects\\northwest-wrs\\docs\\figs\\pop-growth.tiff")
par(mfrow=c(3,4), mar=c(0,3,0,0), oma=c(6,5,5,5))
for (i in 1:length(flnms)){
  mn <- round(min(ryr.hdis[[i]], na.rm=T),2)
  mx <- round(max(ryr.hdis[[i]], na.rm=T),2)
  md <- (mx-mn)/2
  if (i %in% c(2,3,4,5)){  
    int <- round(min(c(abs(mn), abs(mx))),2) } else{
    int <- round(min(c(abs(mn), abs(mx))),1) }
  preddates <- 2004:2019
  plot(preddates, rep(NA, length(preddates)), 
       ylim=c(mn*1.1, mx*1.4),
       xlab="", ylab="",
       xaxt="n", yaxt="n")
  xdates.at <- seq(2005, 2019, by=5)
  xdates <- c(2005, rep(NA, length(xdates.at)-2), 2019)
  axis(1, at=xdates.at, labels=NA, cex.axis=2)
  if(i==2){
    axis(2, at=c(-2*int, 0, int*2, int*4, int*6), labels=c( NA, 0, NA, int*4, NA), cex.axis=1.5)
  } else{
  axis(2, at=c(-2*int, -int, 0, int, 2*int, 3*int, 4*int), labels=c(NA, -int, NA, int, NA, 3*int, NA), cex.axis=1.5)
  }
  polygon(c(preddates, rev(preddates)), c(ryr.hdis[[i]][1,], rev(ryr.hdis[[i]][2,])),
          col="gray60", border=NA)
  lines(preddates, ryr[[i]], lwd=2)
  abline(h=0, lty=2, lwd=2)
  title(spp[i], font.main=1, line=-4, cex.main=2)
  if (i >=8){  axis(1, at=xdates.at, labels=xdates.at, cex.axis=2)   }
}
mtext(side=2, "Population growth rate (r)", outer=T, cex=2, line=1.5)
mtext(side=1, "Year", outer=T, cex=2, line=4.5)
dev.off()

######################
# Create a table of mean overall growth rates
#####################
spp2 <- c("American Kestrel", "Bald Eagle", "Cooper's Hawk", 
         "Ferruginous Hawk", "Golden Eagle", "Northern Harrier",
         "Prairie Falcon", "Rough-legged Hawk", "Red-shouldered Hawk",
         "Red-tailed Hawk", "White-tailed Kite")
tab <- data.frame(Species=spp2, Mean=NA, SE=NA, LHDI=NA, UHDI=NA)
for (i in 1:length(flnms)){  
  load(flnms[[i]])
  tab$Median[i] <- round(median(out$sims.list$r.mu, na.rm=T),3)
  tab$UHDI[i] <- round(hdi(out$sims.list$r.mu, na.rm=T)[1],3)
  tab$LHDI[i] <- round(hdi(out$sims.list$r.mu, na.rm=T)[2],3)
  tab$p_dir[i] <- ifelse(tab$Median[i]<0, mean(out$sims.list$r.mu<0), 
         mean(out$sims.list$r.mu>0))
}
# which are sig
tab$sig <- tab$p_dir>0.975
write.csv(tab, "C:\\Users\\rolek.brian\\Documents\\Projects\\northwest-wrs\\docs\\pop-growth-rates.csv" )

# plot the average growth rates
spp3 <- c("American Kestrel*", "Bald Eagle*", "Cooper's Hawk*", 
          "Ferruginous Hawk*", "Golden Eagle", "Northern Harrier",
          "Prairie Falcon*", "Rough-legged Hawk", "Red-shouldered Hawk*",
          "Red-tailed Hawk*", "White-tailed Kite*")
ldf <- list()
for (i in 1:length(flnms)){  
  load(flnms[[i]])
  ldf[[i]] <-  gather_draws( out, r.mu)
  ldf[[i]]$sp <- spp3[i]
}
df <- do.call(rbind, ldf)
# sort names by median values
mds <- tapply(df$.value, df$sp, FUN=median)
ord <- rev(rownames(mds)[order(mds)])
df$sp <- factor(df$sp, levels = ord )
 
pgr <- df %>% # survival for 6 months
  ggplot(aes(x = .value, y = sp)) +
  stat_gradientinterval(fill_type = "gradient",
                        point_interval = median_hdi,
                        .width = c(0.67, 0.95)) +
  # stat_halfeye(point_interval = median_hdi,
  #              .width = c(0.67, 0.95) )+
  geom_vline(xintercept=0, linetype="dashed", size=1) +
  theme_minimal() +
  ylab("Species common name") + xlab("Population\ngrowth rate (r)") +
  xlim(-0.15,0.15) +
  theme(text = element_text(size = 20))  

ggsave("C:\\Users\\rolek.brian\\Documents\\Projects\\northwest-wrs\\docs\\figs\\average_growth_rts.tiff",
       pgr, dpi=400, height=5, width=7, units="in",
       device="tiff")


########################
# Compare r over time
# among data integrating models
# and those without data integration
########################
flnms <- list() 
flnms[[1]] <- list.files("G:\\nw-wrs\\outputs6\\both\\", full.names=T)
flnms[[2]] <- list.files("G:\\nw-wrs\\outputs6\\wrs\\", full.names=T)
flnms[[3]] <- list.files("G:\\nw-wrs\\outputs6\\cbc\\", full.names=T)
# extract mean r's over all sites
r.wrs <- r.cbc <- list()
for ( sp in 1:length(flnms[[2]]) ){
  load(flnms[[2]][sp])
  r.wrs[[sp]] <- apply(out$sims.list$r, c(1,2), mean, na.rm=T)

  load(flnms[[3]][sp])
  r.cbc[[sp]] <- apply(out$sims.list$r2, c(1,2), mean, na.rm=T)
}
med <- function(x){ apply(x, 2, median, na.rm=T)}
hd <- function(x){ apply(x, 2, hdi, na.rm=T)}
wm <- lapply(r.wrs, med)
wh <- lapply(r.wrs, hd)
cm <- lapply(r.cbc, med)
ch <- lapply(r.cbc, hd)
# subsetting to remove NAs
fst <- function(x){ min(which(!is.na(x))) }
lapply(wm, fst)
sub <- function(x){ x[40:length(x)]}
wm <- lapply(wm, sub)
cm <- lapply(cm, sub)
sub2 <- function(x){ x[,40:ncol(x)]}
wh <- lapply(wh, sub2)
ch <- lapply(ch, sub2)

# calculate correlation coefficients
n.iter <- nrow(r.wrs[[1]])
nsp <- length(r.wrs)
cc <- array(NA, dim = c(n.iter,nsp))
for (sp in 1:nsp){
rw <- r.wrs[[sp]][, 40:54]  
cw <- r.cbc[[sp]][, 40:54]
for (i in 1:n.iter){
  cc[i,sp] <- cor(rw[i,], cw[i,])
}
}

corr <- round(apply(cc, 2, median),3)
corh <- round(apply(cc, 2, hdi),3)

# Compute the posterior modes of correlation coefficients
for (j in 1:6){
  m <- density(correl.h[,j], na.rm = TRUE)
  correl.est[j,1]<- m$x[which(m$y==max(m$y))]
}
correl.est <- round(correl.est,2)
n.iter <- dim(out$sims.list$l.mu.F)[[1]]
P <- c()
P[1] <- sum(correl.h[!is.na(correl.h[,1]),1]>0)/n.iter

tiff(height=6, width=10, res=300, units="in",
     "C:\\Users\\rolek.brian\\Documents\\Projects\\northwest-wrs\\docs\\figs\\r-cors-WRSvsCBC.tiff")
par(mfrow=c(3,4), mar=c(2,2,1,1), oma=c(5,5,1,1))
for (sp in 1:length(wm)){
  mn.x <- min(wh[[sp]], na.rm=T)
  mx.x <- max(wh[[sp]], na.rm=T)
  mn.y <- min(ch[[sp]], na.rm=T)
  mx.y <- max(ch[[sp]], na.rm=T)
  plot(wm[[sp]], cm[[sp]], 
       xlim=c(mn.x, mx.x)*1.1, ylim=c(mn.y, mx.y)*1.1,
       xaxt="n", yaxt="n", type="n")
  abline(h=0, col="gray")
  abline(v=0, col="gray")
  points(wm[[sp]], cm[[sp]], pch=16)
  xlim <- if( par("usr")[2]-par("usr")[1]>0.1){
                  c(round(par("usr")[1:2]*0.85, 1)) } else{
                  c(round(par("usr")[1:2]*0.85, 2))}
  
  ylim <- if( par("usr")[4]-par("usr")[3]>0.1){
                  c(round(par("usr")[3:4]*0.85, 1)) } else{
                  c(round(par("usr")[3:4]*0.85, 2))}
  axis(1, xaxs="s", xaxp=c(xlim, 2))
  axis(2, xaxs="s", yaxp=c(ylim, 2))
  maint <- paste(spp2[sp],"\n", corr[sp], " (", 
        corh[,sp][1], " - ", corh[,sp][2], ")", sep="")
  title(main=maint , line=-3, cex.main=1.5)
for (t in 1:length(wm[[1]])){
  lines(x=c( wh[[sp]][1,t], wh[[sp]][2,t]), y=c(cm[[sp]][t], cm[[sp]][t]))
  lines(x=c( wm[[sp]][t], wm[[sp]][t]), y=c( ch[[sp]][1,t], ch[[sp]][2,t]) )
} #t
  } #sp
mtext("Population growth rate (r), WRS data",
      outer=T, side=1, line=2, cex=2)
mtext("Population growth rate (r), CBC data",
      outer=T, side=2, line=2, cex=2)
dev.off()

##################
# plot overall means
###################
# extract mean r's over all sites
# r.comp <- r.wrs <- r.cbc <- list()
# for ( sp in 1:length(flnms[[2]]) ){
#   load(flnms[[1]][sp])
#   r.comp[[sp]] <- apply(out$sims.list$r, c(1), mean, na.rm=T)
#   
#   load(flnms[[2]][sp])
#   r.wrs[[sp]] <- apply(out$sims.list$r, c(1), mean, na.rm=T)
#   
#   load(flnms[[3]][sp])
#   r.cbc[[sp]] <- apply(out$sims.list$r2, c(1), mean, na.rm=T)
# }
# wm <- lapply(r.wrs, median)
# wh <- lapply(r.wrs, hdi)
# cm <- lapply(r.cbc, median)
# ch <- lapply(r.cbc, hdi)

r.comp <- r.wrs <- r.cbc <- list()
for ( sp in 1:length(flnms[[2]]) ){
  load(flnms[[1]][sp])
  r.comp[[sp]] <- out$sims.list$r.mu
  
  load(flnms[[2]][sp])
  r.wrs[[sp]] <- out$sims.list$r.mu
  
  load(flnms[[3]][sp])
  r.cbc[[sp]] <- out$sims.list$r.mu
}

wm <- lapply(r.wrs, median)
wh <- lapply(r.wrs, hdi)
cm <- lapply(r.cbc, median)
ch <- lapply(r.cbc, hdi)
bm <- lapply(r.comp, median)
bh <- lapply(r.comp, hdi)
cols <- viridis(4)
spp4 <- c("AMKE", "BAEA", "COHA", 
          "FEHA", "GOEA", "NOHA",
          "PRFA", "RLHA", "RSHA",
          "RTHA", "WTKI")
plot(1:length(spp4), cm, type="n", 
     xlim=c(0, length(spp4)+1), ylim=c(-0.6, 0.2),
     xlab="Species", 
     ylab="Average population growth rate",
     xaxt="n", yaxt="n")
abline(h=seq(-0.6, 0.2, by=0.1), col="gray80")
abline(v=seq(0.5, 11.5, by=1), col="gray80")

points((1:length(spp4))-0.3, wm, pch="-", cex=2, col=cols[1])
points(1:length(spp4), cm, pch="-", cex=2, col=cols[2])
points((1:length(spp4))+0.3, bm, pch="-", cex=2, col=cols[3])
for (i in 1:length(wm)){
  segments(x0=(1:11-0.3)[i], 
           x1=(1:11-0.3)[i], 
           y0=wh[[i]][1], 
           y1=wh[[i]][2], 
           lwd=2, col=cols[1] )
  segments(x0=c(1:11)[i], 
           x1=c(1:11)[i], 
           y0=ch[[i]][1], 
           y1=ch[[i]][2] , 
            lwd=2, col=cols[2] )
  segments(x0=(1:11+0.3)[i], 
           x1=(1:11+0.3)[i], 
           y0=bh[[i]][1], 
           y1=bh[[i]][2], 
           lwd=2, col=cols[3] )
}
abline(h=0, lty=2, lwd=2)
axis(1, at=1:11, labels=spp4, las=2)
axis(2, at=c(-0.6, -0.4, -0.2, 0, 0.2), 
     labels=c(-0.6, -0.4, -0.2, 0, 0.2))
legend(1,-0.2, legend=c("WRS only", "CBC only", "Composite"), 
       lwd=2, pch="-", col=cols[1:3])
#######################
# OLD Plots
#####################
flnms <- list.files("C:\\Users\\rolek.brian\\Documents\\Projects\\northwest-wrs\\outputs\\", full.names=T)
sppnms <- list.files("C:\\Users\\rolek.brian\\Documents\\Projects\\northwest-wrs\\outputs\\")

for (i in 1:length(flnms)){
  load(flnms[[i]])
  par(mfrow=c(2,3))
  traceplot(out, "mean.r")
  mtext(side=3, sppnms[[i]])
  traceplot(out, "sig.proc")
  mtext(side=3, sppnms[[i]])
  traceplot(out, "sig.obs")
  mtext(side=3, sppnms[[i]])
}

uniq.time$d <- ifelse(uniq.time$m==2, 14, 16)
preddates <- as.Date( with(uniq.time, paste(y,m,d, sep="/")), format="%Y/%m/%d")
mycol <- rgb(t(col2rgb("gray50")), max = 255, alpha = 150)
survs <- t(tapply(surv.dat$route_num, list(surv.dat$route_num, surv.dat$year_month) , length, default=0))
surveyed <- survs>0
surveyed[surveyed==F] <- NA

spplabels <- c("American\nKestrel", "Bald\nEagle", "Barn\nOwl", "Cooper's\nHawk",
               "Ferruginous\nHawk", "Great-horned\nOwl", "Golden\nEagle", "Merlin",
               "Northern\nHarrier", "Peregrine\nFalcon", "Prairie\nFalcon", "Rough-legged\nHawk",
               "Red-shouldered\nHawk", "Red-tailed\nHawk", "Sharp-shinned\nHawk", "White-tailed\nKite")
# plot on the natural scale
par(mfrow=c(4,4), mar=c(0,3,0,0), oma=c(5,5,1,1))
for (i in 1:length(flnms)){
  load(flnms[[i]])
  mn.N <- exp(apply(out$sims.list$lmu.N, 2, mean ))
  mn.Ns<- mn.N*surveyed
  uci.Ns <- exp(apply(mn.Ns, 2, quantile, 0.975 ))
  lci.Ns <- exp(apply(mn.Ns, 2, quantile, 0.025 ))
  mn <- floor(min(mn.Ns, na.rm=T))
  mx <- ceiling(max(mn.Ns, na.rm=T))
  md <- (mx-mn)/2
plot(preddates, rep(NA, length(preddates)), 
     ylim=c(mn, mx*1.2),
     xlab="", ylab="",
     xaxt="n", yaxt="n")
xdates <- as.Date( paste(c(2005, 2010, 2015, 2020), "/1/1", sep=""), format="%Y/%m/%d")
axis(2, at=c(mn, mn+md, mx), labels=c(mn, md, mx), cex.axis=2)
if(i %in% c(13:16)){
axis(1, at=xdates, labels=c(2005, 2010, 2015, 2020), cex.axis=2)}
if(!i %in% c(13:16)){
  axis(1, at=xdates, labels=c(NA, NA, NA, NA))}
#for (r in 1:389){ lines(preddates, mn.Ns[,r], col=mycol ) }
polygon(c(preddates, rev(preddates)), c(uci.N, rev(lci.N)), col=mycol, border=NA )
lines(preddates, mn.N, lwd=4)
title(spplabels[i], font=1, line=-4, cex.main=2)
}


# plot on the log scale
par(mfrow=c(4,4), mar=c(0,0,0,0), oma=c(10,10,1,1))
for (i in 1:length(flnms)){
  load(flnms[[i]])
  mn.N <- exp(apply(out$sims.list$lmu.N, 2, mean ))
  mn.Ns<- mn.N*surveyed
  mn <- min(mn.Ns, na.rm=T)
  mx <- ceiling(max(mn.Ns, na.rm=T))
  md <- (mx-mn)/2
  uci.N <- exp(apply(out$sims.list$lmu.N, 2, quantile, 0.975 ))
  lci.N <- exp(apply(out$sims.list$lmu.N, 2, quantile, 0.025 ))
  plot(preddates, rep(NA, length(preddates)), 
       ylim=c(log(0.000001), log(100)*1.2),
       xlab="", ylab="",
       xaxt="n", yaxt="n")
  xdates <- as.Date( paste(c(2005, 2010, 2015, 2020), "/1/1", sep=""), format="%Y/%m/%d")
  if(i %in% c(1,5,9,13)){
    axis(2, at=log(c(0.000001, 0.00001, 0.0001, 0.001, 0.01, 0.1, 1, 10, 100)), 
         labels=c(0.000001, NA, NA, NA, NA, NA, NA, NA, 100), 
         cex.axis=2, las=1)}
  if(i %in% c(13:16)){
    axis(1, at=xdates, labels=c(NA, 2010, NA, 2020), cex.axis=2)}
  # if(!i %in% c(13:16)){
  #   axis(1, at=xdates, labels=c(NA, NA, NA, NA))}
  polygon(c(preddates, rev(preddates)), c(uci.N, rev(lci.N)), col=mycol, border=NA )
  lines(preddates, log(apply(mn.Ns, 1, mean, na.rm=T)), lwd=4)
  title(spplabels[i], font=1, line=-4, cex.main=2)
}
mtext(side=2, "Abundance (spaced on log scale)", outer=T, line=7, cex=2)
mtext(side=1, "Year", outer=T, line=5, cex=2)

#############
# plot with CIs insted of routes
##############
par(mfrow=c(4,4), mar=c(0,0,0,0), oma=c(10,10,1,1))
for (i in 1:length(flnms)){
  load(flnms[[i]])
  l.Nr <- apply(out$sims.list$lmu.N, c(1,2), mean )
  l.N <- exp(apply(l.Nr, 2, mean ))
  uci.N <- exp(apply(l.Nr, 2, quantile, 0.975 ))
  lci.N <- exp(apply(l.Nr, 2, quantile, 0.025 ))
  mn <- min(lci.N)
  mx <- max(uci.N)
  md <- (mx-mn)/2
  plot(preddates, rep(NA, length(preddates)), 
       ylim=c(mn, mx*1.2),
       xlab="", ylab="",
       xaxt="n", yaxt="n")
  xdates <- as.Date( paste(c(2005, 2010, 2015, 2020), "/1/1", sep=""), format="%Y/%m/%d")
  if(i %in% c(1,5,9,13)){
    axis(2, at=c(mn, md, mx), 
         labels=c(mn, md, mx), 
         cex.axis=2, las=1)}
  if(i %in% c(13:16)){
    axis(1, at=xdates, labels=c(NA, 2010, NA, 2020), cex.axis=2)}
  polygon(c(preddates, rev(preddates)), c(uci.N, rev(lci.N)), col=mycol, border=NA )
  lines(preddates, l.N, lwd=4)
  title(spplabels[i], font=1, line=-4, cex.main=2)
}
mtext(side=2, "Abundance (spaced on log scale)", outer=T, line=7, cex=2)
mtext(side=1, "Year", outer=T, line=5, cex=2)
