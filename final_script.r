@@ -0,0 +1,207 @@
##Script for running CV_protocol 100 times on datasets

library(phyloseq)
library(ROCR)
library(randomForest)
library(kernlab)
library(Hmisc)

load('/ifs/home/leeb01/otus.rdata')
source('/ifs/home/leeb01/pesame.r')

##FUNCTION
rf.kfoldAUC = function(x, y, k=8, ...){
  n=length(y)
  split.idx = split(unlist(tapply(1:n, y, sample)), rep(1:k, length=n))
  aucs = rep(NA, k)
  split.sample = list()
  x[!is.finite(x)] = 0
  y = factor(y)
  selected = matrix(NA, nrow=dim(x)[2], ncol=k)
  for(i in 1:k){
    sps = table(1:n %in% split.idx[[i]], y)
    if(!any(sps==0)){
      otus.x = table.wilcox.auc(x[-split.idx[[i]],], y[-split.idx[[i]]]); 
      sig_otu <- which(t(otus.x)[,"p.value"]<0.05)
      ## run the model only on OTUs w/p-value<0.05
      selected[sig_otu, i] = 1
      if(length(sig_otu)==0){
        aucs[i] = 0.5
        next
      }
      if(length(sig_otu)==1){
        aucs[i]=wilcox.auc(c(x[split.idx[[i]], sig_otu]), y[split.idx[[i]]])["auc"]   
        next
      }

      rfm <- randomForest(x[-split.idx[[i]], sig_otu, drop=F], y[-split.idx[[i]]], ...)
      aucs[i] = performance(prediction(predict(rfm, x[split.idx[[i]], sig_otu, drop=F], type="prob")[,2], 
                                       y[split.idx[[i]]]), measure='auc')@y.values[[1]]
    }
    
    split.sample = append(split.sample, sps)
  }
  list(aucs=aucs, splits = split.idx, splits.sample.size = split.sample, selected.var = selected)  
}

##MALES

   ##OTU + HOSP
d.merg.hosp.results <- replicate(100, rf.kfoldAUC(dd1, sample_data(otu.subs)$CASECTL, k=8))
d.merg.hosp.aucs <- d.merg.hosp.results[seq(1, 400, 4)]
##Evaluation
d.merg.hosp.aucs.q <- quantile(unlist(d.merg.hosp.aucs), probs=c(0, 0.5, 0.95, 1), na.rm = TRUE)
d.merg.hosp.aucs.mean <- mean(unlist(d.merg.hosp.aucs), na.rm = TRUE)


   ##OTU + AGE
d.merg.age.results <- replicate(100, rf.kfoldAUC(dd2, sample_data(otu.subs)$CASECTL, k=8))
d.merg.age.aucs <- d.merg.age.results[seq(1, 400, 4)]
##Evaluation
d.merg.age.aucs.q <- quantile(unlist(d.merg.age.aucs), probs=c(0, 0.5, 0.95, 1), na.rm = TRUE)
d.merg.age.aucs.mean <- mean(unlist(d.merg.age.aucs), na.rm = TRUE)


   ##OTU + RACE
d.merg.race.results <-replicate(100, rf.kfoldAUC(dd3, sample_data(otu.subs)$CASECTL, k=8))
d.merg.race.aucs <- d.merg.race.results[seq(1, 400, 4)]
##Evaluation
d.merg.race.aucs.q <- quantile(unlist(d.merg.race.aucs), probs=c(0, 0.5, 0.95, 1), na.rm = TRUE)
d.merg.race.aucs.mean <- mean(unlist(d.merg.race.aucs), na.rm = TRUE)


   ##ALL LINKED
d.merg.all.results <- replicate(100, rf.kfoldAUC(dd.merg.all, sample_data(otu.subs)$CASECTL, k=8))
d.merg.all.aucs <- d.merg.all.results[seq(1, 400, 4)]
##Evaluation
d.merg.all.aucs.q <- quantile(unlist(d.merg.all.aucs), probs=c(0, 0.5, 0.95, 1), na.rm = TRUE)
d.merg.all.aucs.mean <- mean(unlist(d.merg.all.aucs), na.rm = TRUE)


   ##UNLINKED
d.merg.unlinked.results <- replicate(100, rf.kfoldAUC(t(otu_table(otu.subs)), sample_data(otu.subs)$CASECTL, k=8))
d.merg.unlinked.aucs <- d.merg.unlinked.results[seq(1, 400, 4)]
##Evaluation
d.merg.unlinked.aucs.q <- quantile(unlist(d.merg.unlinked.aucs), probs=c(0, 0.5, 0.95, 1), na.rm = TRUE)
d.merg.unlinked.aucs.mean <- mean(unlist(d.merg.unlinked.aucs), na.rm = TRUE)


  ##PHYLUM
d.phylum.mat = matrix(t(otu_table(phylum.subs)), ncol=nspecies(phylum.subs), nrow=nsamples(phylum.subs))
d.phylum.results <- replicate(100, rf.kfoldAUC(d.phylum.mat, sample_data(phylum.subs)$CASECTL, k=8))
d.phylum.aucs <- d.phylum.results[seq(1, 400, 4)]
##Evaluation
d.phylum.aucs.q <- quantile(unlist(d.phylum.aucs), probs=c(0, 0.5, 0.95, 1), na.rm = TRUE)
d.phylum.aucs.mean <- mean(unlist(d.phylum.aucs), na.rm = TRUE)

   
   ##CLASS
d.class.mat = matrix(t(otu_table(class.subs)), ncol=nspecies(class.subs), nrow=nsamples(class.subs))
d.class.results <- replicate(100, rf.kfoldAUC(d.class.mat, sample_data(class.subs)$CASECTL, k=8))
d.class.aucs <- d.class.results[seq(1, 400, 4)]
##Evaluation
d.class.aucs.q <- quantile(unlist(d.class.aucs), probs=c(0, 0.5, 0.95, 1), na.rm = TRUE)
d.class.aucs.mean <- mean(unlist(d.class.aucs), na.rm = TRUE)


   ##ORDER
d.order.mat = matrix(t(otu_table(order.subs)), ncol=nspecies(order.subs), nrow=nsamples(order.subs))
d.order.results <- replicate(100, rf.kfoldAUC(d.order.mat, sample_data(order.subs)$CASECTL, k=8))
d.order.aucs <- d.order.results[seq(1, 400, 4)]
##Evaluation
d.order.aucs.q <- quantile(unlist(d.order.aucs), probs=c(0, 0.5, 0.95, 1), na.rm = TRUE)
d.order.aucs.mean <- mean(unlist(d.order.aucs), na.rm = TRUE)


   ##FAMILY
d.family.mat = matrix(t(otu_table(family.subs)), ncol=nspecies(family.subs), nrow=nsamples(family.subs))
d.family.results <- replicate(100, rf.kfoldAUC(d.family.mat, sample_data(family.subs)$CASECTL, k=8))
d.family.aucs <- d.family.results[seq(1, 400, 4)]
##Evaluation
d.family.aucs.q <- quantile(unlist(d.family.aucs), probs=c(0, 0.5, 0.95, 1), na.rm = TRUE)
d.family.aucs.mean <- mean(unlist(d.family.aucs), na.rm = TRUE)


##FEMALES


   ##OTU + HOSP
f.merg.hosp.results <- replicate(100, rf.kfoldAUC(ff1, sample_data(otu.subs.f)$CASECTL, k=8))
f.merg.hosp.aucs <- f.merg.hosp.results[seq(1, 400, 4)]
##Evaluation
f.merg.hosp.aucs.q <- quantile(unlist(f.merg.hosp.aucs), probs=c(0, 0.5, 0.95, 1), na.rm = TRUE)
f.merg.hosp.aucs.mean <- mean(unlist(f.merg.hosp.aucs), na.rm = TRUE)


   ##OTU + AGE
f.merg.age.results <- replicate(100, rf.kfoldAUC(ff2, sample_data(otu.subs.f)$CASECTL, k=8))
f.merg.age.aucs <- f.merg.age.results[seq(1, 400, 4)]
##Evaluation
f.merg.age.aucs.q <- quantile(unlist(f.merg.age.aucs), probs=c(0, 0.5, 0.95, 1), na.rm = TRUE)
f.merg.age.aucs.mean <- mean(unlist(f.merg.age.aucs), na.rm = TRUE)


   #OTU + RACE
f.merg.race.results <- replicate(100, rf.kfoldAUC(ff3, sample_data(otu.subs.f)$CASECTL, k=8))
f.merg.race.aucs <- f.merg.race.results[seq(1, 400, 4)]
##Evaluation
f.merg.race.aucs.q <- quantile(unlist(f.merg.race.aucs), probs=c(0, 0.5, 0.95, 1), na.rm = TRUE)
f.merg.race.aucs.mean <- mean(unlist(f.merg.race.aucs), na.rm = TRUE)


   ##ALL LINKED
f.merg.all.results <- replicate(100, rf.kfoldAUC(ff.merg.all, sample_data(otu.subs.f)$CASECTL, k=8))
f.merg.all.aucs <- f.merg.all.results[seq(1, 400, 4)]
##Evaluation
f.merg.all.aucs.q <- quantile(unlist(f.merg.all.aucs), probs=c(0, 0.5, 0.95, 1), na.rm = TRUE)
f.merg.all.aucs.mean <- mean(unlist(f.merg.all.aucs), na.rm = TRUE)


   ##UNLINKED
f.merg.unlinked.results <- replicate(100, rf.kfoldAUC(t(otu_table(otu.subs.f)), sample_data(otu.subs.f)$CASECTL, k=8))
f.merg.unlinked.aucs <- f.merg.unlinked.results[seq(1, 400, 4)]
##Evaluation
f.merg.unlinked.aucs.q <- quantile(unlist(f.merg.unlinked.aucs), probs=c(0, 0.5, 0.95, 1), na.rm = TRUE)
f.merg.unlinked.aucs.mean <- mean(unlist(f.merg.unlinked.aucs), na.rm = TRUE)


   ##PHYLUM
f.phylum.mat = matrix(t(otu_table(phylum.subs.f)), ncol=nspecies(phylum.subs.f), nrow=nsamples(phylum.subs.f))
f.phylum.results <- replicate(100, rf.kfoldAUC(f.phylum.mat, sample_data(phylum.subs.f)$CASECTL, k=8))
f.phylum.aucs <- f.phylum.results[seq(1, 400, 4)]
##Evaluation
f.phylum.aucs.q <- quantile(unlist(f.phylum.aucs), probs=c(0, 0.5, 0.95, 1), na.rm = TRUE)
f.phylum.aucs.mean <- mean(unlist(f.phylum.aucs), na.rm = TRUE)


   ##CLASS
f.class.mat = matrix(t(otu_table(class.subs.f)), ncol=nspecies(class.subs.f), nrow=nsamples(class.subs.f))
f.class.results <- replicate(100, rf.kfoldAUC(f.class.mat, sample_data(class.subs.f)$CASECTL, k=8))
f.class.aucs <- f.class.results[seq(1, 400, 4)]
##Evaluation
f.class.aucs.q <- quantile(unlist(f.class.aucs), probs=c(0, 0.5, 0.95, 1), na.rm = TRUE)
f.class.aucs.mean <- mean(unlist(f.class.aucs), na.rm = TRUE)


   ##ORDER
f.order.mat = matrix(t(otu_table(order.subs.f)), ncol=nspecies(order.subs.f), nrow=nsamples(order.subs.f))
f.order.results <- replicate(100, rf.kfoldAUC(f.order.mat, sample_data(order.subs.f)$CASECTL, k=8))
f.order.aucs <- f.order.results[seq(1, 400, 4)]
##Evaluation
f.order.aucs.q <- quantile(unlist(f.order.aucs), probs=c(0, 0.5, 0.95, 1), na.rm = TRUE)
f.order.aucs.mean <- mean(unlist(f.order.aucs), na.rm = TRUE)


   ##FAMILY
f.family.mat = matrix(t(otu_table(family.subs.f)), ncol=nspecies(family.subs.f), nrow=nsamples(family.subs.f))
f.family.results <- replicate(100, rf.kfoldAUC(f.family.mat, sample_data(family.subs.f)$CASECTL, k=8))
f.family.aucs <- f.family.results[seq(1, 400, 4)]
##Evaluation
f.family.aucs.q <- quantile(unlist(f.family.aucs), probs=c(0, 0.5, 0.95, 1), na.rm = TRUE)
f.family.aucs.mean <- mean(unlist(f.family.aucs), na.rm = TRUE)


   ##AUC EVALUATION AT DIFFERENT PERCENTILES
quantile(unlist(d.merg.hosp.aucs))
save.image('/Users/brianlee/Documents/CRC/CV_protocol_w_f_s_output/CV analyses results/merg.analysis.all.rdata')
