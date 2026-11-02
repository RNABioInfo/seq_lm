suppressPackageStartupMessages(library(edgeR))
source('bin/de_analysis/profile_intervals.R')
set.seed(941)
y <- matrix(rnbinom(600*6,mu=40,size=10),600,6);rownames(y)<-paste0('g',1:600)
y[1,4:6]<-0;y[2,1:3]<-0;y[3,4:6]<-c(1e-6,0,0)
group <- factor(rep(c('T0','T1'),each=3));design<-model.matrix(~0+group);colnames(design)<-levels(group)
d<-DGEList(y,group=group);d<-normLibSizes(d);d<-estimateDisp(d,design,robust=TRUE)
f<-glmQLFit(d,design,robust=TRUE,keep.unit.mat=TRUE)
r<-profile_intervals(f,c(0,0,0,1,1,1))
stopifnot(is.infinite(r$lower[1]),r$lower[1]<0,is.finite(r$upper[1]),
 is.infinite(r$upper[2]),is.finite(r$lower[2]),is.infinite(r$R[1]),is.infinite(r$R[2]),
 is.finite(r$R[3]),r$R[3]>10000,
 max(c(r$lower_root_relative_error,r$upper_root_relative_error),na.rm=TRUE)<1e-4)
off<-as.matrix(f$offset)+matrix(0,nrow(y),ncol(y));phi<-r$effective_phi
for(b in c(-3,0,2)) {
 z<-glmFit(f$counts,design=matrix(1,6,1),offset=off+outer(rep(b,nrow(y)),c(0,0,0,1,1,1)),dispersion=phi,prior.count=0)
 stable<-stable_null_deviance(f$counts,off,rep(b,nrow(y)),c(0,0,0,1,1,1),phi)
 stopifnot(max(abs(z$deviance-stable))<1e-5)
}
q<-glmQLFTest(f,contrast=c(-1,1));null<-stable_null_deviance(f$counts,off,rep(0,nrow(y)),c(0,0,0,1,1,1),phi)
stopifnot(max(abs(pmax(0,null-f$deviance)/f$s2.post-q$table$F))<1e-4)
# Pointwise profile CI membership equals the installed QL test at zero.
contains<-r$lower<=0&r$upper>=0
stopifnot(identical(unname(contains),unname(q$table$PValue>=.05)))
# Poisson limiting model; no spurious NB cancellation.
z<-glmFit(f$counts,design=matrix(1,6,1),offset=off,dispersion=0,prior.count=0)
stable<-stable_null_deviance(f$counts,off,rep(0,nrow(y)),c(0,0,0,1,1,1),rep(0,nrow(y)))
stopifnot(max(abs(z$deviance-stable))<1e-5)
# Larger sample exposure at fixed QL variance/df improves interior precision.
f2<-f;f2$counts<-f$fitted.values;f2$deviance<-rep(0,nrow(y))
r2<-profile_intervals(f2,c(0,0,0,1,1,1))
f3<-f2;f3$counts[,4:6]<-f3$counts[,4:6]*2;f3$offset<-off;f3$offset[,4:6]<-off[,4:6]+log(2)
r3<-profile_intervals(f3,c(0,0,0,1,1,1))
stopifnot(all(r3$R[4:600]<=r2$R[4:600]*(1+1e-6),na.rm=TRUE))
cat('PASS: boundary, extreme fractional count, two-sided root accuracy, NB profile agreement, QLF statistic, CI inversion, Poisson, depth sensitivity\n')
