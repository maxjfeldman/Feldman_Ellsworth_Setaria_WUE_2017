##############################################
# Load in dependencies
##############################################
library(lme4)
library(ggplot2)
library(lattice)
library(lmomco)

## The first thing is to set your R session to the base directory you just downloaded from github
## insert path below...
setwd()

## Tester
#setwd("~/Dropbox/Feldman_Ellsworth_Setaria_WUE_2017/")

##### CREATE DIRECTORY PATHS ##### 

## Make the directory of the folder you downloaded the current working directory
home.dir<-getwd()


##############################################
# Lets define a function to get a loess fit of timeseries data
##############################################

get.loess.fit<-function(fit, times, geno, cond) {
  return_df<-c()
  predict.vals<-predict(fit, times, se=T)
  genotype<-rep(geno, length(times))
  condition<-rep(as.character(cond), length(times))
  M<-predict.vals$fit
  M.lo<-M - predict.vals$se.fit
  M.hi<-M + predict.vals$se.fit
  slope<-c(0)
  for(s in 2:length(times)) {
    s.temp<-(M[s] - M[s-1]) / 2
    slope<-c(slope, s.temp)
  }
  slope<-slope*10
  return_df<-cbind(genotype, condition,times, M, M.lo, M.hi, slope)
  return(return_df)
}


##############################################
# Make a function to do loess fit, make plots and output data and report
##############################################

report.loess.values<-function(rawdata, trait, genos, treatments, days, from, to, plotname){
  genos<-as.character(sort(unique(genos)))
  treatments<-as.character(treatments)
  treatments<-unique(treatments)
  treatments<-sort(treatments)
  t1<-treatments[1]
  t2<-treatments[2]
  print(t1)
  print(t2)
  dap_i<-unique(days)
  dap_i<-sort(as.numeric(as.character(dap_i)))
  times = seq(from = from, to = to, by=0.1)
  colnumber <- which(colnames(rawdata) %in% trait)

  # Get loess fits for each genotype using get.loess.fit() fxn
  ril_loess_model_fit<-c()
  for (i in 1:length(genos)) {
    r<-genos[i]
    temp1<-rawdata[rawdata$genotype == as.character(r),]
    per.ril<-c()
    for (j in 1:length(treatments)) {
      t<-treatments[j]
      per.t<-c()
      temp2<-temp1[temp1$treatment == as.character(t),]
      if (nrow(temp2) < 1) {next;}
      colnumber2 <- which(colnames(temp2) %in% trait)
      # Log of 0 is INF need to replace with another small #
      temp2[,colnumber2]<-replace(temp2[,colnumber2], temp2[,colnumber2] <= 0, 1)
      out.loess<-loess(get(trait)~dap_i, data=temp2)
      output<-get.loess.fit(out.loess, times, r, t)
      ril_loess_model_fit<-rbind(ril_loess_model_fit, output)
    }
  }
  
  # Now make sure the resulting dataframe has 
  ril_loess_model_fit<-ril_loess_model_fit[complete.cases(ril_loess_model_fit),]
  colnames(ril_loess_model_fit)<-c('ril', 'treatment', 'dap_i', 'M', 'M.lo', 'M.hi', 'AGR')
  ril_loess_model_fit<-as.data.frame(ril_loess_model_fit)
  ril_loess_model_fit$ril<-as.character(ril_loess_model_fit$ril)
  ril_loess_model_fit$treatment<-as.character(ril_loess_model_fit$treatment)
  ril_loess_model_fit$dap_i<-as.numeric(as.character(ril_loess_model_fit$dap_i))
  ril_loess_model_fit$M<-as.numeric(as.character(ril_loess_model_fit$M))
  ril_loess_model_fit$M.lo<-as.numeric(as.character(ril_loess_model_fit$M.lo))
  ril_loess_model_fit$M.hi<-as.numeric(as.character(ril_loess_model_fit$M.hi))
  ril_loess_model_fit$AGR<-as.numeric(as.character(ril_loess_model_fit$AGR))
  ril_loess_model_fit<-ril_loess_model_fit[ril_loess_model_fit$dap_i %in% dap_i,]
  # Lets remove M.lo and M.hi 
  # Sometimes these are difficult to estimate and end up NA
  # Creates plotting problems below
  ril_loess_model_fit<-ril_loess_model_fit[,c(1:4,7)]
  
  rate_id<-sort(unique(as.character(ril_loess_model_fit$ril)))
  growth_rate_report<-c()
  pdf(plotname, paper = "a4", width = 21/2.54, height = (7*5.3)/2.54)
  layout(matrix(c(1:4), ncol = 1, byrow = T), heights = c(1, 1, 1, 1))
  for(i in 1:length(rate_id)) {
    r<-rate_id[i]
    ril_set<-rawdata[(rawdata$genotype == r),]
    ril_rates<-ril_loess_model_fit[(ril_loess_model_fit$ril == r),]
    plant_ids<-unique(ril_set$plantbarcode)
    colnumber3 <- which(colnames(ril_set) %in% trait)
    max.b<-max(max(ril_set[,colnumber3],na.rm=T), max(as.numeric(as.character(ril_rates$M)), na.rm=T),na.rm=T)
    min.b<-min(min(ril_set[,colnumber3],na.rm=T), min(as.numeric(as.character(ril_rates$M)), na.rm=T),na.rm=T)
    set<-ril_set[ril_set$plantbarcode == plant_ids[1],]
    # Start making plots
    if (set[1,'treatment'] == t2) {l.color<-c("light blue")}
    if (set[1,'treatment'] == t1) {l.color<-c("gold")}
    
    plot(set[,colnumber3]~set$dap_i, type='p',  pch=19, xlim=c(from,to), ylim=c(min.b, max.b), col=l.color, xlab="Days after planting", ylab=trait, main=r)
    if(length(plant_ids) >1) {
      for (j in 2:length(plant_ids)) {
        set<-ril_set[ril_set$plantbarcode == plant_ids[j],]
        if (set[1,'treatment'] == t2) {l.color<-c("light blue")}
        if (set[1,'treatment'] == t1) {l.color<-c("gold")}
        points(set[,colnumber3]~set$dap_i, type='p',  pch=19, xlim=c(from,to), ylim=c(min.b, max.b), col=l.color, xlab="Day", ylab=trait)
      }
    }
    
    rate.t2<-ril_rates[ril_rates$treatment == t2, ]
    if (nrow(rate.t2) > 0) {
      l.color<-c("blue")
      p.color<-c("dark blue")
      max.rate.t2<-max(as.numeric(as.character(rate.t2$AGR)),na.rm=T)
      day.t2<-rate.t2[rate.t2$AGR == max.rate.t2, 'dap_i']
      max.val.t2<-as.numeric(as.character(rate.t2[rate.t2$AGR == max.rate.t2, 'M']))
      lines(rate.t2$M~rate.t2$dap_i, lwd=2, col=l.color)
      points(day.t2, max.val.t2, cex=1.5, col=p.color, pch=18)
      # Can optionally plot confidence interval
      #lines(rate.t2$M.hi~rate.t2$dap_i, lty=2, col=c('navy'))
      #lines(rate.t2$M.lo~rate.t2$dap_i, lty=2, col=c('navy'))
    }
    
    rate.t1<-ril_rates[ril_rates$treatment == t1, ]  
    if (nrow(rate.t1) > 0) {
      l.color<-c("orange")
      p.color<-c("dark orange")
      max.rate.t1<-max(rate.t1$AGR,na.rm=T)
      day.t1<-rate.t1[rate.t1$AGR == max.rate.t1, 'dap_i']
      max.val.t1<-rate.t1[rate.t1$AGR == max.rate.t1, 'M'] 
      lines(rate.t1$M~rate.t1$dap_i, lwd=2, col=l.color)
      points(day.t1, max.val.t1, cex=1.5, col=p.color, pch=18)
      #lines(rate.t1$M.hi~rate.t1$dap_i, lty=2, col=c('brown'))
      #lines(rate.t1$M.lo~rate.t1$dap_i, lty=2, col=c('brown'))
    }
    # treatment drought
    rate.t1<-ril_rates[ril_rates$treatment == t1, ]  
    
    if (nrow(rate.t1) > 0) {
      max.day_trait.t1<-max(rate.t1$M, na.rm=T)
      max.agr.t1<-max(rate.t1$AGR,na.rm=T)
      rate.t1<-rate.t1[complete.cases(rate.t1),]
      max.agr.day.t1<-rate.t1[rate.t1$AGR == max.agr.t1, 'dap_i']
      max.agr.day_trait.t1<-rate.t1[rate.t1$AGR == max.agr.t1, 'M']
    }
    
    # treatment well watered
    rate.t2<-ril_rates[ril_rates$treatment == t2, ]
    if (nrow(rate.t2) > 0) {
      max.day_trait.t2<-max(rate.t2$M,na.rm=T)
      max.agr.t2<-max(rate.t2$AGR,na.rm=T)
      rate.t2<-rate.t2[complete.cases(rate.t2),]
      max.agr.day.t2<-rate.t2[rate.t2$AGR == max.agr.t2, 'dap_i']
      max.agr.day_trait.t2<-rate.t2[rate.t2$AGR == max.agr.t2, 'M']
    }
    
    # Generate the report on a per/ril basis
    if (length(unique(ril_rates$treatment)) > 1) {
      ril_entry<-c(r, max.day_trait.t2, max.day_trait.t1, max.day_trait.t2 - max.day_trait.t1, max.agr.t2, max.agr.t1, max.agr.t2 - max.agr.t1, max.agr.day.t2, max.agr.day.t1, max.agr.day.t2 - max.agr.day.t1, max.agr.day_trait.t2, max.agr.day_trait.t1, max.agr.day_trait.t2 - max.agr.day_trait.t1)
      growth_rate_report<-rbind(growth_rate_report, ril_entry)
    }
    # Plot rates
    if (length(unique(ril_rates$treatment)) > 1) {
      max.r<-max(max(rate.t2$AGR,na.rm=T), max(rate.t1$AGR,na.rm=T),na.rm=T)
      min.r<-min(min(rate.t2$AGR,na.rm=T), min(rate.t1$AGR,na.rm=T),na.rm=T)
      rate.t2<<-rate.t2
      rate.t1<<-rate.t1
      max.r<<-max.r
      min.r<<-min.r
      max.agr.day.t2<<-max.agr.day.t2
      max.agr.t2<<-max.agr.t2
      max.agr.day.t1<<-max.agr.day.t1
      max.agr.t1<<-max.agr.t1
      plot(rate.t2$AGR~rate.t2$dap_i, type="l", col="blue", xlab='Days after planting', ylab='Rate', ylim=c(min.r,max.r))
      lines(rate.t1$AGR~rate.t1$dap_i, col="orange")
      points(max.agr.day.t2, max.agr.t2, pch=18, col="dark blue")
      points(max.agr.day.t1, max.agr.t1, pch=18, col="dark orange")
      
    }
  }
  dev.off()
  
 
  #ril_loess_model_fit<-ril_loess_model_fit[,c(1:4,7)]
  colnames(ril_loess_model_fit)<-c("genotype", "treatment", "dap_i", trait, paste(trait, '_rate', sep=""))
  # Give column names to growth report
  rownames(growth_rate_report)<-c(1:nrow(growth_rate_report))
  growth_rate_report<-as.data.frame(growth_rate_report)
  colnames(growth_rate_report)<-c("genotype", "max_value.t2", "max_value.t1","max_value.diff","max_rate.t2","max_rate.t1","max_rate.diff","max_day.t2","max_day.t1","max_day.diff","value_max_rate_day.t2","value_max_rate_day.t1", "value_max_rate_day.diff")
  ril_loess_model_fit<<-ril_loess_model_fit
  growth_rate_report<<-growth_rate_report
}


##############################################
# Lets define a function to calculate heritability
##############################################

# Broad sense heritability
get_h2<-function(data){
  i.treat<-unique(data$treatment)
  pheno<-data
  year<-unique(pheno[,3])
  exp<-unique(pheno[,2])
  pheno<-pheno[,c(7,4,9:length(colnames(pheno)))]
  colnames(pheno)[c(1,2)]<-c("id", "treatment")
  pheno[pheno == "."] <- NA
  colnames(pheno)[3:ncol(pheno)]<-paste(colnames(pheno)[3:ncol(pheno)] , exp, sep="_")
  H2<-c()
  for (i in 3:length(colnames(pheno))){
    # Get complete cases
    cc.pheno<-pheno[complete.cases(pheno[,i]),c(1:2,i)]
    # Build linear model each cofactor is a random effect
    model<-lmer(cc.pheno[,3]~(1|id)+(1|treatment)+(1|id:treatment), data=cc.pheno)
    
    # Extract variance from model object, save individual components in vector 're' and residual variance as a scalar named 'res'
    re<-as.numeric(VarCorr(model))
    res<-attr(VarCorr(model), "sc")^2
    
    # Extract individual components (order will remain the same)
    gxt.var<-re[1]
    geno.var<-re[2]
    treat.var<-re[3]
    # Total variance is sum of all variances
    tot.var<-sum(re, res)
    
    reps.t1<-table(pheno[pheno$treatment == i.treat[1], 'id'])
    reps.t2<-table(pheno[pheno$treatment == i.treat[2], 'id'])
    reps.treatment<-c(reps.t1, reps.t2)
    
    reps.t1<-as.character(unique(pheno[pheno$treatment == i.treat[1], 'id']))
    reps.t2<-as.character(unique(pheno[pheno$treatment == i.treat[2], 'id']))
    unique.combined <- c(as.character(reps.t1), as.character(reps.t2))
    
    freq.unique.combined <- table(unique.combined)
    
    # Calculate the harmonic mean replication within treatment blocks
    hm_treatment<-harmonic.mean(freq.unique.combined)$harmean
    
    # Now get a count of total genotypic replication
    reps.total<-table(pheno[,'id'])
    # Get the harmonic mean of this quantity
    hm_total<-harmonic.mean(reps.total)$harmean
    
    # Calculate heritability as described by AEL 
    # H2 = geno.var/(geno.var + (gxt.var/harmonic mean of treatment block replication) + (residual.var/harmonic mean of total genotype replication) )
    h2<-((geno.var)/(geno.var + (gxt.var/hm_treatment) + (res/hm_total)))
    
    # This is the heritability
    H2<-c(H2,h2)
    
  }
  names(H2)<-colnames(pheno)[3:ncol(pheno)]
  return(H2)
}

# Heritability within treatment
get_h2_in_treatment<-function(data){
  i.treat<-unique(data$treatment)
  pheno<-data
  year<-unique(pheno[,3])
  exp<-unique(pheno[,2])
  pheno<-pheno[,c(7,4,9:length(colnames(pheno)))]
  colnames(pheno)[c(1,2)]<-c("id", "treatment")
  pheno[pheno == "."] <- NA
  colnames(pheno)[3:ncol(pheno)]<-paste(colnames(pheno)[3:ncol(pheno)] , exp, sep="_")
  variance.out<-c()
  for (t in 1:length(i.treat)) {
    # Create variables to store values
    treatment.pheno<-pheno[pheno$treatment == i.treat[t],]
    H2<-c()
    e2<-c()
    
    # For each treatment.phenotype calculate variance
    for(i in 3:length(colnames(treatment.pheno))){
      # Use only RILs with all measurements for each treatment.phenotype
      cc.treatment.pheno<-treatment.pheno[complete.cases(treatment.pheno[,i]),c(1:2,i)]
      # Build linear model each cofactor is a random effect
      model<-lmer(cc.treatment.pheno[,3]~(1|id), data=cc.treatment.pheno, control=lmerControl(check.nobs.vs.nlev = "ignore", check.nobs.vs.rankZ = "ignore",check.nobs.vs.nRE="ignore"))
      # Extract variance from model object, save individual components in vector 're' and residual variance as a scalar named 'res'
      re<-as.numeric(VarCorr(model))
      res<-attr(VarCorr(model), "sc")^2
      # Extract individual components (order will remain the same)
      geno.var<-re[1]
      # Total variance is sum of all variances
      tot.var<-sum(re, res)
      # Get proportion of variance
      h<-geno.var/tot.var
      e<-res/tot.var
      # Append variables to a vector of variables
      H2<-c(H2,h)
      e2<-c(e2,e)
    }
    
    variance<-rbind(H2, e2)
    colnames(treatment.pheno)[3:length(treatment.pheno)]<-paste(i.treat[t], colnames(treatment.pheno)[3:length(treatment.pheno)], sep="_")
    colnames(variance)<-colnames(treatment.pheno)[3:length(treatment.pheno)]
    rownames(variance)<-c('Genotype', 'Error')
    assign(paste('variance', i.treat[t], sep="_"), variance)
    variance.out<-cbind(variance.out, variance)
  }
  return(variance.out)
}  

# Total variance partition
get_total_var<-function(data){
  i.treat<-unique(data$treatment)
  pheno<-data
  year<-unique(pheno[,3])
  exp<-unique(pheno[,2])
  pheno<-pheno[,c(7,4,9:length(colnames(pheno)))]
  colnames(pheno)[c(1,2)]<-c("id", "treatment")
  pheno[pheno == "."] <- NA
  colnames(pheno)[3:ncol(pheno)]<-paste(colnames(pheno)[3:ncol(pheno)] , exp, sep="_")
  
  H2<-c()
  t2<-c()
  e2<-c()
  gxt2<-c()
  
  for (i in 3:length(colnames(pheno))){
    # Get complete cases
    cc.pheno<-pheno[complete.cases(pheno[,i]),c(1:2,i)]
    
    # Build linear model each cofactor is a random effect
    model<-lmer(cc.pheno[,3]~(1|id)+(1|treatment)+(1|id:treatment), data=cc.pheno)
    
    # Extract variance from model object, save individual components in vector 're' and residual variance as a scalar named 'res'
    re<-as.numeric(VarCorr(model))
    res<-attr(VarCorr(model), "sc")^2
    
    # Extract individual components (order will remain the same)
    gxt.var<-re[1]
    geno.var<-re[2]
    treat.var<-re[3]
    # Total variance is sum of all variances
    tot.var<-sum(re, res)
    
    # Get proportion of variance
    h<-geno.var/tot.var
    t<-treat.var/tot.var
    e<-res/tot.var
    gxt<-gxt.var/tot.var
    # Append variables to a vector of variables
    H2<-c(H2,h)
    t2<-c(t2,t)
    e2<-c(e2,e)
    gxt2<-c(gxt2, gxt)
    
  }
  variance<-rbind(H2, t2, gxt2, e2)
  colnames(variance)<-colnames(pheno)[3:length(pheno)]
  rownames(variance)<-c('Genotype', 'Treatment', 'G x Treatment', 'Error')
  return(variance)
}


# Total variance partition field (includes plot)
get_total_var_field<-function(data){
  i.treat<-unique(data$treatment)
  pheno<-data
  year<-unique(pheno[,3])
  exp<-unique(pheno[,2])
  pheno<-pheno[,c(7,4,9:length(colnames(pheno)))]
  colnames(pheno)[c(1,2)]<-c("id", "treatment")
  pheno[pheno == "."] <- NA
  colnames(pheno)[3:ncol(pheno)]<-paste(colnames(pheno)[3:ncol(pheno)] , exp, sep="_")
  
  H2<-c()
  t2<-c()
  e2<-c()
  p2<-c()
  gxt2<-c()
  
  for (i in 5:length(colnames(pheno))){
    # Get complete cases
    cc.pheno<-pheno[complete.cases(pheno[,i]),c(1:2,i)]
    
    # Build linear model each cofactor is a random effect
    model<-lmer(cc.pheno[,3]~(1|id)+(1|treatment)+(1|plot)+(1|id:treatment), data=cc.pheno)
    
    # Extract variance from model object, save individual components in vector 're' and residual variance as a scalar named 'res'
    re<-as.numeric(VarCorr(model))
    res<-attr(VarCorr(model), "sc")^2
    
    # Extract individual components (order will remain the same)
    gxt.var<-re[1]
    geno.var<-re[2]
    plot.var<-re[3]
    treat.var<-re[4]
    # Total variance is sum of all variances
    tot.var<-sum(re, res)
    
    # Get proportion of variance for all factors
    h<-geno.var/tot.var
    t<-treat.var/tot.var
    e<-res/tot.var
    gxt<-gxt.var/tot.var
    p<-plot.var/tot.var
    # Append variables to a vector of variables
    H2<-c(H2,h)
    t2<-c(t2,t)
    e2<-c(e2,e)
    gxt2<-c(gxt2, gxt)
    p2<-c(p2, p)
    
  }
  variance<-rbind(H2, t2, p2, gxt2, e2)
  colnames(variance)<-colnames(pheno)[3:length(pheno)]
  rownames(variance)<-c('Genotype', 'Treatment','Plot', 'G x Treatment', 'Error')
  return(variance)
}


loess.fit.for.h2<-function(data, trait){
  barcodes<-unique(data$plantbarcode)
  dap_i<-sort(unique(data$dap_i))
  plant_id_loess_model_fit<-c()
  for(b in barcodes){
    print(b)
    temp<-data[data$plantbarcode == b,]
    genotype<-unique(temp$genotype)
    treatment<-unique(temp$treatment)
    #colnumber <- which(colnames(temp) %in% trait)
    out.loess<-loess(get(trait)~dap_i, data=temp)
    times = seq(from = min(temp$dap_i), to = max(temp$dap_i), by=0.1)
    output<-get.loess.estimates(out.loess, times, b)
    output<-as.data.frame(output)
    output<-output[output$times %in% dap_i,]
    output$genotype<-rep(genotype, nrow(output))
    output$treatment<-rep(treatment, nrow(output))
    colnames(output)<-c("plantbarcode", "dap_i", trait,paste(trait, ".lo" ,sep=""),paste(trait, ".hi", sep=""), paste(trait, ".slope", sep=""), 'genotype', 'treatment') 
    output<-output[,c(1,7,8,2:6)]
    plant_id_loess_model_fit<-rbind(plant_id_loess_model_fit, output)
  }
  return(plant_id_loess_model_fit)
}

get.loess.estimates<-function(fit, times, id) {
  return_df<-c()
  predict.vals<-predict(fit, times, se=T)
  ids<-rep(as.character(id), length(times))
  #condition<-rep(as.character(cond), length(times))
  M<-predict.vals$fit
  M.lo<-M - predict.vals$se.fit
  M.hi<-M + predict.vals$se.fit
  slope<-c(0)
  for(s in 2:length(times)) {
    s.temp<-(M[s] - M[s-1]) / 2
    slope<-c(slope, s.temp)
  }
  slope<-slope*10
  return_df<-cbind(ids,times, M, M.lo, M.hi, slope)
  return(return_df)
}


##############################################
# Lets define a function to merge the QTL results based upon a single marker
##############################################

merged_table<-c()
unify_chr<-function(temp){
  all_qtl<-sort(table(temp$marker), decreasing=T)
  if(length(all_qtl) > 1){
    #m.name<-names(all_qtl)[1]
    m.name<-temp[temp$lod == max(temp$lod),'marker']
    m.name<-m.name[1]
    ave.pos<-mean(temp[temp$marker == m.name, 'pos'])
    #cr<-unique(temp[temp$marker == names(all_qtl)[1],'chr'])
    #po<-unique(temp[temp$marker == names(all_qtl)[1],'pos'])
    cr<-unique(temp[temp$marker == m.name,'chr'])
    po<-unique(temp[temp$marker == m.name,'pos'])
    max.pos<-ave.pos+10
    min.pos<-ave.pos-10
    subset<-temp[temp$pos > min.pos & temp$pos < max.pos,]
    subset$marker<-rep(m.name, length(subset$marker))
    subset$chr<-rep(cr, length(subset$chr))
    subset$pos<-rep(po, length(subset$pos))
    temp<-temp[temp$pos < min.pos | temp$pos > max.pos,]
    merged_table<<-rbind(merged_table, subset)
    unify_chr(temp)
  } 
  if(length(all_qtl) == 1){
    merged_table<<-rbind(merged_table, temp)
    #unify_chr(temp)
  }
}


unify_marker<-function(input){
  
  chrs<-sort(unique(input$chr))
  merged_table<<-c()
  for(ch in chrs) {
    temp<-input[input$chr == ch,]
    temp$marker<-as.character(temp$marker)
    unify_chr(temp)
  }
  return(merged_table)
}


###### This is a function to get unique qtl from a qtl summary table
# Basically collpase redundant QTL into 10 cM intervals

remove_dup_qtl<-function(temp){
  all_qtl<-sort(table(temp$marker), decreasing=T)
  if (length(all_qtl) == 1) {
    treatments<-as.character(unique(temp$treatment))
    if(length(treatments) == 1) { 
      m.names<<-c(m.names, names(all_qtl)[1])
      # <<- means change the global variable (chr<<-max) changes the global variable chr to local variable max
      chr<<-c(chr,unique(temp[temp$marker == names(all_qtl)[1],'chr']))
      pos<<-c(pos,unique(temp[temp$marker == names(all_qtl)[1],'pos']))
      t<-as.character(unique(temp$treatment))
      condition<<-c(condition, t)
      qtl_count<<-c(qtl_count, 1)
      
      max.lod<<-c(max.lod, max(temp[temp$marker == names(all_qtl)[1],'lod'], na.rm = T))
      max.prop.var<<-c(max.prop.var, max(temp[temp$marker == names(all_qtl)[1],'prop.var'], na.rm = T))
      max.fx<<-c(max.fx, max(temp[temp$marker == names(all_qtl)[1],'additive.fx'], na.rm = T))
      max.fx_se<<-c(max.fx_se, max(temp[temp$marker == names(all_qtl)[1],'additive.fx_se'], na.rm = T))
      max.L.CI_pos<<-c(max.L.CI_pos, max(temp[temp$marker == names(all_qtl)[1],'L.CI_pos'], na.rm = T))
      max.R.CI_pos<<-c(max.R.CI_pos, max(temp[temp$marker == names(all_qtl)[1],'R.CI_pos'], na.rm = T))
      
      med.lod<<-c(med.lod, median(temp[temp$marker == names(all_qtl)[1],'lod'], na.rm = T))
      med.prop.var<<-c(med.prop.var, median(temp[temp$marker == names(all_qtl)[1],'prop.var'], na.rm = T))
      med.fx<<-c(med.fx, median(temp[temp$marker == names(all_qtl)[1],'additive.fx'], na.rm = T))
      med.fx_se<<-c(med.fx_se, median(temp[temp$marker == names(all_qtl)[1],'additive.fx_se'], na.rm = T))
      med.L.CI_pos<<-c(med.L.CI_pos, median(temp[temp$marker == names(all_qtl)[1],'L.CI_pos'], na.rm = T))
      med.R.CI_pos<<-c(med.R.CI_pos, median(temp[temp$marker == names(all_qtl)[1],'R.CI_pos'], na.rm = T))
      
      min.lod<<-c(min.lod, min(temp[temp$marker == names(all_qtl)[1],'lod'], na.rm = T))
      min.prop.var<<-c(min.prop.var, min(temp[temp$marker == names(all_qtl)[1],'prop.var'], na.rm = T))
      min.fx<<-c(min.fx, min(temp[temp$marker == names(all_qtl)[1],'additive.fx'], na.rm = T))
      min.fx_se<<-c(min.fx_se, min(temp[temp$marker == names(all_qtl)[1],'additive.fx_se'], na.rm = T))
      min.L.CI_pos<<-c(min.L.CI_pos, min(temp[temp$marker == names(all_qtl)[1],'L.CI_pos'], na.rm = T))
      min.R.CI_pos<<-c(min.R.CI_pos, min(temp[temp$marker == names(all_qtl)[1],'R.CI_pos'], na.rm = T))
      
      print(chr) 
      print(pos)
      print(qtl_count)
      print(med.lod)
    } 
    if(length(treatments) > 1){
      for(t in treatments) {
        t.temp<-temp[temp$treatment == t,]
        m.names<<-c(m.names, names(all_qtl)[1])
        # <<- means change the global variable (chr<<-max) changes the global variable chr to local variable max
        chr<<-c(chr,unique(t.temp[t.temp$marker == names(all_qtl)[1],'chr']))
        pos<<-c(pos,unique(t.temp[t.temp$marker == names(all_qtl)[1],'pos']))
        condition<<-c(condition, t)
        qtl_count<<-c(qtl_count, 1)
        
        max.lod<<-c(max.lod, max(t.temp[t.temp$marker == names(all_qtl)[1],'lod'], na.rm=T))
        max.prop.var<<-c(max.prop.var, max(t.temp[t.temp$marker == names(all_qtl)[1],'prop.var'], na.rm = T))
        max.fx<<-c(max.fx, max(t.temp[t.temp$marker == names(all_qtl)[1],'additive.fx'], na.rm = T))
        max.fx_se<<-c(max.fx_se, max(t.temp[t.temp$marker == names(all_qtl)[1],'additive.fx_se'], na.rm = T))
        max.L.CI_pos<<-c(max.L.CI_pos, max(t.temp[t.temp$marker == names(all_qtl)[1],'L.CI_pos'], na.rm = T))
        max.R.CI_pos<<-c(max.R.CI_pos, max(t.temp[t.temp$marker == names(all_qtl)[1],'R.CI_pos'], na.rm = T))
        
        med.lod<<-c(med.lod, median(t.temp[t.temp$marker == names(all_qtl)[1],'lod'], na.rm=T))
        med.prop.var<<-c(med.prop.var, median(t.temp[t.temp$marker == names(all_qtl)[1],'prop.var'], na.rm = T))
        med.fx<<-c(med.fx, median(t.temp[t.temp$marker == names(all_qtl)[1],'additive.fx'], na.rm = T))
        med.fx_se<<-c(med.fx_se, median(t.temp[t.temp$marker == names(all_qtl)[1],'additive.fx_se'], na.rm = T))
        med.L.CI_pos<<-c(med.L.CI_pos, median(t.temp[t.temp$marker == names(all_qtl)[1],'L.CI_pos'], na.rm = T))
        med.R.CI_pos<<-c(med.R.CI_pos, median(t.temp[t.temp$marker == names(all_qtl)[1],'R.CI_pos'], na.rm = T))
        
        min.lod<<-c(min.lod, min(t.temp[t.temp$marker == names(all_qtl)[1],'lod'], na.rm=T))
        min.prop.var<<-c(min.prop.var, min(t.temp[t.temp$marker == names(all_qtl)[1],'prop.var'], na.rm = T))
        min.fx<<-c(min.fx, min(t.temp[t.temp$marker == names(all_qtl)[1],'additive.fx'], na.rm = T))
        min.fx_se<<-c(min.fx_se, min(t.temp[t.temp$marker == names(all_qtl)[1],'additive.fx_se'], na.rm = T))
        min.L.CI_pos<<-c(min.L.CI_pos, min(t.temp[t.temp$marker == names(all_qtl)[1],'L.CI_pos'], na.rm = T))
        min.R.CI_pos<<-c(min.R.CI_pos, min(t.temp[t.temp$marker == names(all_qtl)[1],'R.CI_pos'], na.rm = T))
        
        print(chr) 
        print(pos)
        print(qtl_count)
        print(med.lod)
      }
    }
  }
  if (length(all_qtl) > 1) {
    
    #name<-names(all_qtl)[1]
    name<-temp[temp$lod == max(temp$lod),'marker']
    name<-name[1]
    tester<-temp[temp$marker == name,]
    treatments<-as.character(unique(tester$treatment))
    
    if(length(treatments) == 1) { 
      ave.pos<-mean(temp[temp$marker == name, 'pos'])
      #m.name<-names(all_qtl)[1]
      m.name<-name
      #cr<-unique(temp[temp$marker == names(all_qtl)[1],'chr'])
      #po<-unique(temp[temp$marker == names(all_qtl)[1],'pos'])
      cr<-unique(temp[temp$marker == m.name,'chr'])
      po<-unique(temp[temp$marker == m.name,'pos'])
      max.pos<-ave.pos+10
      min.pos<-ave.pos-10
      subset<-temp[temp$pos > min.pos & temp$pos < max.pos,]
      if(length(unique(subset$treatment)) == 1){
        m.names<<-c(m.names, m.name)
        chr<<-c(chr, cr)
        pos<<-c(pos, po)
        qtl_c<-nrow(subset)
        
        x.lod<-max(subset$lod, na.rm = T)
        x.prop.var<-max(subset$prop.var, na.rm = T)
        x.fx<-max(subset$additive.fx, na.rm = T)
        x.fx_se<-max(subset$additive.fx_se, na.rm = T)
        x.L.CI_pos<-max(subset$L.CI_pos, na.rm = T)
        x.R.CI_pos<-max(subset$R.CI_pos, na.rm = T)
        
        m.lod<-median(subset$lod, na.rm = T)
        m.prop.var<-median(subset$prop.var, na.rm = T)
        m.fx<-median(subset$additive.fx, na.rm = T)
        m.fx_se<-median(subset$additive.fx_se, na.rm = T)
        m.L.CI_pos<-median(subset$L.CI_pos, na.rm = T)
        m.R.CI_pos<-median(subset$R.CI_pos, na.rm = T)
        
        n.lod<-min(subset$lod, na.rm = T)
        n.prop.var<-min(subset$prop.var, na.rm = T)
        n.fx<-min(subset$additive.fx, na.rm = T)
        n.fx_se<-min(subset$additive.fx_se, na.rm = T)
        n.L.CI_pos<-min(subset$L.CI_pos, na.rm = T)
        n.R.CI_pos<-min(subset$R.CI_pos, na.rm = T)
        
        temp<-temp[temp$pos < min.pos | temp$pos > max.pos,]
        print(ave.pos) 
        print(chr) 
        print(pos)
        #print(collapsed_qtl)
        print(med.lod)
        condition<<-c(condition, treatments[1])
        qtl_count<<-c(qtl_count, qtl_c)
        
        max.lod<<-c(max.lod, x.lod)
        max.prop.var<<-c(max.prop.var, x.prop.var)
        max.fx<<-c(max.fx, x.fx)
        max.fx_se<<-c(max.fx_se, x.fx_se)
        max.L.CI_pos<<-c(max.L.CI_pos, x.L.CI_pos)
        max.R.CI_pos<<-c(max.R.CI_pos, x.R.CI_pos)
        
        med.lod<<-c(med.lod, m.lod)
        med.prop.var<<-c(med.prop.var, m.prop.var)
        med.fx<<-c(med.fx, m.fx)
        med.fx_se<<-c(med.fx_se, m.fx_se)
        med.L.CI_pos<<-c(med.L.CI_pos, m.L.CI_pos)
        med.R.CI_pos<<-c(med.R.CI_pos, m.R.CI_pos)
        
        min.lod<<-c(min.lod, n.lod)
        min.prop.var<<-c(min.prop.var, n.prop.var)
        min.fx<<-c(min.fx, n.fx)
        min.fx_se<<-c(min.fx_se, n.fx_se)
        min.L.CI_pos<<-c(min.L.CI_pos, n.L.CI_pos)
        min.R.CI_pos<<-c(min.R.CI_pos, n.R.CI_pos)
        
        remove_dup_qtl(temp)
      }
      if(length(unique(subset$treatment)) > 1){
        subset.ts<-unique(subset$treatment)
        for (t in subset.ts) {
          m.names<<-c(m.names, m.name)
          chr<<-c(chr, cr)
          pos<<-c(pos, po)
          t.subset<-subset[subset$treatment == t,]
          qtl_c<-nrow(t.subset)
          
          x.lod<-max(t.subset$lod, na.rm = T)
          x.prop.var<-max(t.subset$prop.var, na.rm = T)
          x.fx<-max(t.subset$additive.fx, na.rm = T)
          x.fx_se<-max(t.subset$additive.fx_se, na.rm = T)
          x.L.CI_pos<-max(t.subset$L.CI_pos, na.rm = T)
          x.R.CI_pos<-max(t.subset$R.CI_pos, na.rm = T)
          
          m.lod<-median(t.subset$lod, na.rm = T)
          m.prop.var<-median(t.subset$prop.var, na.rm = T)
          m.fx<-median(t.subset$additive.fx, na.rm = T)
          m.fx_se<-median(t.subset$additive.fx_se, na.rm = T)
          m.L.CI_pos<-median(t.subset$L.CI_pos, na.rm = T)
          m.R.CI_pos<-median(t.subset$R.CI_pos, na.rm = T)
          
          n.lod<-min(t.subset$lod, na.rm = T)
          n.prop.var<-min(t.subset$prop.var, na.rm = T)
          n.fx<-min(t.subset$additive.fx, na.rm = T)
          n.fx_se<-min(t.subset$additive.fx_se, na.rm = T)
          n.L.CI_pos<-min(t.subset$L.CI_pos, na.rm = T)
          n.R.CI_pos<-min(t.subset$R.CI_pos, na.rm = T)
          
          print(ave.pos) 
          print(chr) 
          print(pos)
          #print(collapsed_qtl)
          print(med.lod)
          
          condition<<-c(condition, t)
          qtl_count<<-c(qtl_count, qtl_c)
          
          max.lod<<-c(max.lod, x.lod)
          max.prop.var<<-c(max.prop.var, x.prop.var)
          max.fx<<-c(max.fx, x.fx)
          max.fx_se<<-c(max.fx_se, x.fx_se)
          max.L.CI_pos<<-c(max.L.CI_pos, x.L.CI_pos)
          max.R.CI_pos<<-c(max.R.CI_pos, x.R.CI_pos)
          
          med.lod<<-c(med.lod, m.lod)
          med.prop.var<<-c(med.prop.var, m.prop.var)
          med.fx<<-c(med.fx, m.fx)
          med.fx_se<<-c(med.fx_se, m.fx_se)
          med.L.CI_pos<<-c(med.L.CI_pos, m.L.CI_pos)
          med.R.CI_pos<<-c(med.R.CI_pos, m.R.CI_pos)
          
          min.lod<<-c(min.lod, n.lod)
          min.prop.var<<-c(min.prop.var, n.prop.var)
          min.fx<<-c(min.fx, n.fx)
          min.fx_se<<-c(min.fx_se, n.fx_se)
          min.L.CI_pos<<-c(min.L.CI_pos, n.L.CI_pos)
          min.R.CI_pos<<-c(min.R.CI_pos, n.R.CI_pos)
          
        }
        temp<-temp[temp$pos < min.pos | temp$pos > max.pos,]
        remove_dup_qtl(temp)
      }
    }
    if(length(treatments) > 1) {
      #for (t in treatments) {
      ave.pos<-mean(temp[temp$marker == name, 'pos'])
      #m.name<-names(all_qtl)[1]
      #cr<-unique(temp[temp$marker == names(all_qtl)[1],'chr'])
      #po<-unique(temp[temp$marker == names(all_qtl)[1],'pos'])
      m.name<-name
      cr<-unique(temp[temp$marker == m.name,'chr'])
      po<-unique(temp[temp$marker == m.name,'pos'])
      max.pos<-ave.pos+10
      min.pos<-ave.pos-10
      subset<-temp[temp$pos > min.pos & temp$pos < max.pos,]
      subset.ts<-unique(subset$treatment)
      for (t in subset.ts) {
        m.names<<-c(m.names, m.name)
        chr<<-c(chr, cr)
        pos<<-c(pos, po)
        t.subset<-subset[subset$treatment == t,]
        qtl_c<-nrow(t.subset)
        
        x.lod<-max(t.subset$lod, na.rm = T)
        x.prop.var<-max(t.subset$prop.var, na.rm = T)
        x.fx<-max(t.subset$additive.fx, na.rm = T)
        x.fx_se<-max(t.subset$additive.fx_se, na.rm = T)
        x.L.CI_pos<-max(t.subset$L.CI_pos, na.rm = T)
        x.R.CI_pos<-max(t.subset$R.CI_pos, na.rm = T)
        
        m.lod<-median(t.subset$lod, na.rm = T)
        m.prop.var<-median(t.subset$prop.var, na.rm = T)
        m.fx<-median(t.subset$additive.fx, na.rm = T)
        m.fx_se<-median(t.subset$additive.fx_se, na.rm = T)
        m.L.CI_pos<-median(t.subset$L.CI_pos, na.rm = T)
        m.R.CI_pos<-median(t.subset$R.CI_pos, na.rm = T)
        
        n.lod<-min(t.subset$lod, na.rm = T)
        n.prop.var<-min(t.subset$prop.var, na.rm = T)
        n.fx<-min(t.subset$additive.fx, na.rm = T)
        n.fx_se<-min(t.subset$additive.fx_se, na.rm = T)
        n.L.CI_pos<-min(t.subset$L.CI_pos, na.rm = T)
        n.R.CI_pos<-min(t.subset$R.CI_pos, na.rm = T)
        
        print(ave.pos) 
        print(chr) 
        print(pos)
        #print(collapsed_qtl)
        print(med.lod)
        
        condition<<-c(condition, t)
        qtl_count<<-c(qtl_count, qtl_c)
        
        max.lod<<-c(max.lod, x.lod)
        max.prop.var<<-c(max.prop.var, x.prop.var)
        max.fx<<-c(max.fx, x.fx)
        max.fx_se<<-c(max.fx_se, x.fx_se)
        max.L.CI_pos<<-c(max.L.CI_pos, x.L.CI_pos)
        max.R.CI_pos<<-c(max.R.CI_pos, x.R.CI_pos)
        
        med.lod<<-c(med.lod, m.lod)
        med.prop.var<<-c(med.prop.var, m.prop.var)
        med.fx<<-c(med.fx, m.fx)
        med.fx_se<<-c(med.fx_se, m.fx_se)
        med.L.CI_pos<<-c(med.L.CI_pos, m.L.CI_pos)
        med.R.CI_pos<<-c(med.R.CI_pos, m.R.CI_pos)
        
        min.lod<<-c(min.lod, n.lod)
        min.prop.var<<-c(min.prop.var, n.prop.var)
        min.fx<<-c(min.fx, n.fx)
        min.fx_se<<-c(min.fx_se, n.fx_se)
        min.L.CI_pos<<-c(min.L.CI_pos, n.L.CI_pos)
        min.R.CI_pos<<-c(min.R.CI_pos, n.R.CI_pos)
        
      }
      temp<-temp[temp$pos < min.pos | temp$pos > max.pos,]
      remove_dup_qtl(temp)
    }
  }
}


condense_qtl<-function(input){
  
  chrs<-sort(unique(input$chr))
  m.names<<-c()
  chr<<-c()
  pos<<-c()
  condition<<-c()
  qtl_count<<-c()
  med.lod<<-c()
  med.prop.var<<-c()
  med.fx<<-c()
  med.fx_se<<-c()
  med.L.CI_pos<<-c()
  med.R.CI_pos<<-c()
  
  max.lod<<-c()
  max.prop.var<<-c()
  max.fx<<-c()
  max.fx_se<<-c()
  max.L.CI_pos<<-c()
  max.R.CI_pos<<-c()
  
  min.lod<<-c()
  min.prop.var<<-c()
  min.fx<<-c()
  min.fx_se<<-c()
  min.L.CI_pos<<-c()
  min.R.CI_pos<<-c()
  
  for(ch in chrs) {
    temp<-input[input$chr == ch,]
    temp$marker<-as.character(temp$marker)
    remove_dup_qtl(temp)
  }
  
  input.collapsed<-as.data.frame(cbind(m.names, chr, pos, condition, qtl_count, max.lod, max.prop.var, max.fx, max.fx_se, max.L.CI_pos, max.R.CI_pos,med.lod, med.prop.var, med.fx, med.fx_se, med.L.CI_pos, med.R.CI_pos,min.lod, min.prop.var, min.fx, min.fx_se, min.L.CI_pos, min.R.CI_pos))        
  return(input.collapsed)
}




##############################################
# Lets define a function to make common QTL plots
##############################################


make_qtl_common_plot<-function(all.qtl, plotname) {

  all.qtl$chr<-factor(all.qtl$chr, levels=c(1,2,3,4,5,6,7,8,9))

  fx.size<-all.qtl$additive.fx
  fx.size<-as.numeric(as.character(fx.size))

  plot.char<-c()
  for(i in 1:length(fx.size)){
    if (fx.size[i] > 0) {plot.char<-c(plot.char, '24')}
    if (fx.size[i] < 0) {plot.char<-c(plot.char, '25')}
  }

  all.qtl$plot.char<-plot.char
  all.qtl$plot.char<-as.factor(all.qtl$plot.char)

  all.qtl$group<-paste(all.qtl$exp, all.qtl$year, all.qtl$treatment, sep="_")

  treatments<-as.character(all.qtl$treatment)
  treatment.name<-unique(treatments)
  plot.col<-c()
  for(i in 1:length(treatments)){
    logical<-treatments[i] == treatment.name
    col<-which(logical, arr.ind=TRUE)
    plot.col<-c(plot.col, col)
  }

  all.qtl$plot.col<-plot.col
  
  pdf(plotname)
  p<-ggplot() + geom_point(data = all.qtl, aes(x = pos, y = prop.var, shape=plot.char, colour=as.character(plot.col), fill=as.character(plot.col)),size=3, alpha=0.5) + geom_blank(data = blank_data, aes(x = x, y = y)) + facet_wrap(~chr, scales = "free_x") + expand_limits(x = 0) + scale_x_continuous(expand = c(0, 0)) + theme_bw() +   scale_shape_manual(values=c(24,25)) 
  print(p + scale_color_manual(values=c("1" = "orange", "2" = "blue")) + scale_fill_manual(values=c("1" = "orange", "2" = "blue")) + ylab("% Variance")  + xlab("Genome Position") + theme(legend.position = "none"))
  dev.off()

}




make_qtl_common_plot_diff<-function(all.qtl, plotname) {
  
  all.qtl$chr<-factor(all.qtl$chr, levels=c(1,2,3,4,5,6,7,8,9))
  
  fx.size<-all.qtl$additive.fx
  fx.size<-as.numeric(as.character(fx.size))
  
  plot.char<-c()
  for(i in 1:length(fx.size)){
    if (fx.size[i] > 0) {plot.char<-c(plot.char, '24')}
    if (fx.size[i] < 0) {plot.char<-c(plot.char, '25')}
  }
  
  all.qtl$plot.char<-plot.char
  all.qtl$plot.char<-as.factor(all.qtl$plot.char)
  
  all.qtl$group<-paste(all.qtl$exp, all.qtl$year, all.qtl$treatment, sep="_")
  
  treatments<-as.character(all.qtl$treatment)
  treatment.name<-unique(treatments)
  plot.col<-c()
  for(i in 1:length(treatments)){
    logical<-treatments[i] == treatment.name
    col<-which(logical, arr.ind=TRUE)
    plot.col<-c(plot.col, col)
  }
  
  all.qtl$plot.col<-plot.col
  
  pdf(plotname)
  p<-ggplot() + geom_point(data = all.qtl, aes(x = pos, y = prop.var, shape=plot.char, colour=as.character(plot.col), fill=as.character(plot.col)),size=3, alpha=0.5) + geom_blank(data = blank_data, aes(x = x, y = y)) + facet_wrap(~chr, scales = "free_x") + expand_limits(x = 0) + scale_x_continuous(expand = c(0, 0)) + theme_bw() +   scale_shape_manual(values=c(24,25)) 
  print(p + scale_color_manual(values=c("1" = "orange", "2" = "blue", "3" = "grey")) + scale_fill_manual(values=c("1" = "orange", "2" = "blue", "3" = "grey")) + ylab("% Variance")  + xlab("Genome Position") + theme(legend.position = "none"))
  dev.off()
  
}




##############################################
# Lets define a function to calculated predicted values and residuals from a major axis model
##############################################
get.lmodel2.values<-function(model, rma){
  if(rma == 'N') {
    # Get model intercepts
    ols.int<-model$regression.results$Intercept[1]
    ma.int<-model$regression.results$Intercept[2]
    sma.int<-model$regression.results$Intercept[3]
    # Get model slope
    ols.slope<-model$regression.results$Slope[1]
    ma.slope<-model$regression.results$Slope[2]
    sma.slope<-model$regression.results$Slope[3]
    
    # Get values you specified as X
    x<-model$x
    y<-model$y
    
    # Get predicted values
    y.ols.pred<-(ols.slope * x) + ols.int
    y.ma.pred<-(ma.slope * x) + ma.int
    y.sma.pred<-(sma.slope * x) + sma.int
    
    # Get residuals from the fit
    y.ols.res<-(y-y.ols.pred)
    y.ma.res<-(y-y.ma.pred)
    y.sma.res<-(y-y.sma.pred)
    
    # Format results
    out<-cbind(x,y,y.ols.pred,y.ma.pred,y.sma.pred,y.ols.res,y.ma.res,y.sma.res)
    colnames(out)<-c('x','y','ols.pred','ma.pred','sma.pred','ols.res','ma.res','sma.res')
  }
  
  if(rma == 'N') {
    # Get model intercepts
    ols.int<-model$regression.results$Intercept[1]
    ma.int<-model$regression.results$Intercept[2]
    sma.int<-model$regression.results$Intercept[3]
    rma.int<-model$regression.results$Intercept[4]
    # Get model slope
    ols.slope<-model$regression.results$Slope[1]
    ma.slope<-model$regression.results$Slope[2]
    sma.slope<-model$regression.results$Slope[3]
    rma.slope<-model$regression.results$Slope[4]
    
    # Get values you specified as X
    x<-model$x
    y<-model$y
    
    # Get predicted values
    y.ols.pred<-(ols.slope * x) + ols.int
    y.ma.pred<-(ma.slope * x) + ma.int
    y.sma.pred<-(sma.slope * x) + sma.int
    y.rma.pred<-(rma.slope * x) + rma.int
    
    # Get residuals from the fit
    y.ols.res<-(y.ols.pred-y)
    y.ma.res<-(y.ma.pred-y)
    y.sma.res<-(y.sma.pred-y)
    y.rma.res<-(y.rma.pred-y)
    
    # Format results
    out<-cbind(x,y,y.ols.pred,y.ma.pred,y.sma.pred,y.rma.pred,y.ols.res,y.ma.res,y.sma.res,y.rma.res)
    colnames(out)<-c('x','y','ols.pred','ma.pred','sma.pred','rma.pred','ols.res','ma.res','sma.res','rma.res')
  }
return(out)
}

setwd(home.dir)
save.image('analysis_fxns.Rdata')

