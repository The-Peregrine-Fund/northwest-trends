plot.diag <- function(out, ratio=TRUE, lab=""){
  par(mfrow=c(2,2))
  # plot WRS mean absolute percentage error
  mx <- max(c(out$sims.list$dmape.rep.wrs,out$sims.list$dmape.obs.wrs))
  mn <- min(c(out$sims.list$dmape.rep.wrs,out$sims.list$dmape.obs.wrs))
  plot(jitter(out$sims.list$dmape.obs.wrs, amount=300), 
       jitter(out$sims.list$dmape.rep.wrs, amount=300),
       main=paste0("Mean absolute percentage error\nWRS model\n",lab),
       ylab="Discrepancy replicate values",
       xlab="Discrepancy observed values", 
       xlim=c(mn,mx), ylim=c(mn,mx), 
       pch=16, cex=0.5, col="gray10")
  curve(1*x, from=mn, to=mx, add=T, lty=2, lwd=2, col="blue")
  bp1 <- round(mean(out$sims.list$dmape.rep.wrs > out$sims.list$dmape.obs.wrs),2)
  loc <- ifelse(bp1 < 0.5, "topleft", "bottomright")
  legend(loc, legend=bquote(p[B]==.(bp1)), bty="n", cex=2)
  
  # plot CBC mean absolute percentage error
  mx <- max(c(out$sims.list$dmape.rep.cbc,out$sims.list$dmape.obs.cbc))
  mn <- min(c(out$sims.list$dmape.rep.cbc,out$sims.list$dmape.obs.cbc))
  plot(jitter(out$sims.list$dmape.obs.cbc), 
       jitter(out$sims.list$dmape.rep.cbc),
       main=paste0("Mean absolute percentage error\nCBC model\n",lab),
       ylab="Discrepancy replicate values",
       xlab="Discrepancy observed values", 
       xlim=c(mn,mx), ylim=c(mn,mx), 
       pch=16, cex=0.5, col="gray10")
  curve(1*x, from=mn, to=mx, add=T, lty=2, lwd=2, col="blue")
  bp2 <- round(mean(out$sims.list$dmape.rep.cbc > out$sims.list$dmape.obs.cbc),2)
  loc <- ifelse(bp2 < 0.5, "topleft", "bottomright")
  legend(loc, legend=bquote(p[B]==.(bp2)), bty="n", cex=2)
  
  if (ratio==TRUE){
    # plot WRS variance/mean ratio
    hist(out$sims.list$tvm.rep.wrs, nclass=50,
         xlab="variance/mean ", main=NA, axes=FALSE)
    abline(v=out$mean$tvm.obs.wrs, col="red")
    axis(1); axis(2)
    # plot CBC variance/mean ratio
    hist(out$sims.list$tvm.rep.cbc, nclass=50,
         xlab="variance/mean ", main=NA, axes=FALSE)
    abline(v=out$mean$tvm.obs.cbc, col="red")
    axis(1); axis(2)
  }
  
  print(paste0("WRS ratio of observation and process error=", round(out$mean$sig.ratio[1],1) ))
  print(paste0("CBC ratio of observation and process error=", round(out$mean$sig.ratio[2],1) ))
  return(c(WRS=bp1,CBC=bp2))
}

bp.func <- function(out){
  bp1 <- round(mean(out$sims.list$dmape.rep.wrs > out$sims.list$dmape.obs.wrs),2)
  bp2 <- round(mean(out$sims.list$dmape.rep.cbc > out$sims.list$dmape.obs.cbc),2)
return(c(WRS=bp1,CBC=bp2))
}