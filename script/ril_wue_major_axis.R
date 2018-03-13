
library(ggplot2)
library(lmodel2)

## The first thing is to set your R session to the base directory you just downloaded from github
setwd()


## Tester
#setwd("~/Dropbox/Feldman_Ellsworth_Setaria_WUE_2017/")

##### CREATE DIRECTORIES ##### 

## Make the directory of the folder you downloaded the current working directory
home.dir<-getwd()
setwd(home.dir)
load('analysis_fxns.Rdata')

wue_data.dir<-paste(home.dir, '/data/', sep="")

## If this directory doesn't exist the download is not complete
if (file.exists(wue_data.dir)){
  print("Directory exists!")
} else {
  dir.create(file.path(wue_data.dir))
}

## Make the results directory
wue_results.dir<-paste(home.dir, '/results/', sep="")

## If this directory doesn't exist the download is not complete
if (file.exists(wue_results.dir)){
  print("Directory exists!")
} else {
  dir.create(file.path(wue_results.dir))
}

## Make a directory for plant size results
wue_results.plant_size.dir<-paste(wue_results.dir, 'plant_size/', sep="")

## If this directory doesn't exist the download is not complete
if (file.exists(wue_results.plant_size.dir)){
  print("Directory exists!")
} else {
  dir.create(file.path(wue_results.plant_size.dir))
}

## Make a directory for water use results
wue_results.water_use.dir<-paste(wue_results.dir, 'water_use/', sep="")

## If this directory doesn't exist the download is not complete
if (file.exists(wue_results.water_use.dir)){
  print("Directory exists!")
} else {
  dir.create(file.path(wue_results.water_use.dir))
}

## Make a directory for water use results
wue_results.te.dir<-paste(wue_results.dir, 'te/', sep="")

## If this directory doesn't exist the download is not complete
if (file.exists(wue_results.te.dir)){
  print("Directory exists!")
} else {
  dir.create(file.path(wue_results.te.dir))
}

##### END ##### 


setwd(wue_results.te.dir)
load("ril_transpiration_efficiency.Rdata")

# Define fxn to get values from lmodel2 value
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
  if(rma == 'Y') {
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
    y.ols.res<-(y-y.ols.pred)
    y.ma.res<-(y-y.ma.pred)
    y.sma.res<-(y-y.sma.pred)
    y.rma.res<-(y-y.rma.pred)
    
    # Format results
    out<-cbind(x,y,y.ols.pred,y.ma.pred,y.sma.pred,y.rma.pred,y.ols.res,y.ma.res,y.sma.res,y.rma.res)
    colnames(out)<-c('x','y','ols.pred','ma.pred','sma.pred','rma.pred','ols.res','ma.res','sma.res','rma.res')
  }
  return(out)
}

# Lets run a loop and fit model sv_area ~ water used 
days<-sort(unique(day_size$dap_i))
te_residual.ma<-c()
for(d in days){
  temp<-day_size[day_size$dap_i == d,]
  treats<-sort(unique(temp$treatment))
  temp.t<-c()
  for(t in treats){
    temp2<-temp[temp$treatment == t,]
    # Get fit for cumulative
    te_total.mdl<-lmodel2(sv_area~water_lost_total, data=temp2, range.y="relative", range.x="relative", nperm=99)
    fit_te_total.val<-get.lmodel2.values(te_total.mdl, 'Y')
    fit_te_total.val<-as.data.frame(fit_te_total.val)
    ma_te_fit_total<-fit_te_total.val$ma.pred
    ma_te_residual_total<-fit_te_total.val$ma.res
    # Get fit for daily/rate
    te_day.mdl<-lmodel2(sv_area_day~water_lost_day, data=temp2, range.y="relative", range.x="relative", nperm=99)
    fit_te_day.val<-get.lmodel2.values(te_day.mdl, 'Y')
    fit_te_day.val<-as.data.frame(fit_te_day.val)
    ma_te_fit_day<-fit_te_day.val$ma.pred
    ma_te_residual_day<-fit_te_day.val$ma.res
    temp3<-cbind(temp2,ma_te_fit_total,ma_te_residual_total,ma_te_fit_day,ma_te_residual_day)
    temp.t<-rbind(temp.t, temp3)
  }
  te_residual.ma<-rbind(te_residual.ma, temp.t)
}





# Format for te_fit_total QTL pipeline
wue.l<-te_residual.ma
wue.l<-wue.l[,c('genotype', 'treatment', 'dap_i', 'ma_te_fit_total')]

days<-sort(unique(wue.l$dap_i))

ril_te_fit_total_qtl<-c()
for(d in 1:length(days)) {
  day<-days[d]
  temp.data<-wue.l[wue.l$dap_i == day, ]
  colnames(temp.data)[4]<-paste('ma_te_fit_total', day, sep="_")
  temp.data<-temp.data[,c(1,2,4)]
  if (d == 1) {
    ril_te_fit_total_qtl<-temp.data
  }
  
  if (d > 1) {
    ril_te_fit_total_qtl<-merge(ril_te_fit_total_qtl, temp.data, by=c('genotype', 'treatment'), all=T)
  }
}

# Add misc columns
ril_te_fit_total_qtl$Obs<-c(1:nrow(ril_te_fit_total_qtl))
ril_te_fit_total_qtl$experiment<-rep("BP14", nrow(ril_te_fit_total_qtl))
ril_te_fit_total_qtl$year<-rep("2014", nrow(ril_te_fit_total_qtl))
ril_te_fit_total_qtl$plot<-rep("bellweater", nrow(ril_te_fit_total_qtl))
ril_te_fit_total_qtl$plot_id<-rep("bellweater", nrow(ril_te_fit_total_qtl))
ril_te_fit_total_qtl$sampling<-rep("bellweater", nrow(ril_te_fit_total_qtl))



# Reorder columns
ril_te_fit_total_qtl<-ril_te_fit_total_qtl[,c(20,21,22,2,23,24,1,25,3:19)]

colnames(ril_te_fit_total_qtl)[7]<-c("id")
write.csv(ril_te_fit_total_qtl, file="ril_ma_te_fit_total_qtl.csv", quote=F, row.names=F)


# Format for te_residual_total QTL pipeline
wue.l<-te_residual.ma
wue.l<-wue.l[,c('genotype', 'treatment', 'dap_i', 'ma_te_residual_total')]

days<-sort(unique(wue.l$dap_i))

ril_te_residual_total_qtl<-c()
for(d in 1:length(days)) {
  day<-days[d]
  temp.data<-wue.l[wue.l$dap_i == day, ]
  colnames(temp.data)[4]<-paste('ma_te_residual_total', day, sep="_")
  temp.data<-temp.data[,c(1,2,4)]
  if (d == 1) {
    ril_te_residual_total_qtl<-temp.data
  }
  
  if (d > 1) {
    ril_te_residual_total_qtl<-merge(ril_te_residual_total_qtl, temp.data, by=c('genotype', 'treatment'), all=T)
  }
}

# Add misc columns
ril_te_residual_total_qtl$Obs<-c(1:nrow(ril_te_residual_total_qtl))
ril_te_residual_total_qtl$experiment<-rep("BP14", nrow(ril_te_residual_total_qtl))
ril_te_residual_total_qtl$year<-rep("2014", nrow(ril_te_residual_total_qtl))
ril_te_residual_total_qtl$plot<-rep("bellweater", nrow(ril_te_residual_total_qtl))
ril_te_residual_total_qtl$plot_id<-rep("bellweater", nrow(ril_te_residual_total_qtl))
ril_te_residual_total_qtl$sampling<-rep("bellweater", nrow(ril_te_residual_total_qtl))

# Reorder columns
ril_te_residual_total_qtl<-ril_te_residual_total_qtl[,c(20,21,22,2,23,24,1,25,3:19)]

colnames(ril_te_residual_total_qtl)[7]<-c("id")
write.csv(ril_te_residual_total_qtl, file="ril_ma_te_residual_total_qtl.csv", quote=F, row.names=F)

#################################################



# Format te_fit_day for QTL pipeline
wue.l<-te_residual.ma
wue.l<-wue.l[,c('genotype', 'treatment', 'dap_i', 'ma_te_fit_day')]

days<-sort(unique(wue.l$dap_i))

ril_te_fit_day_qtl<-c()
for(d in 1:length(days)) {
  day<-days[d]
  temp.data<-wue.l[wue.l$dap_i == day, ]
  colnames(temp.data)[4]<-paste('ma_te_fit_day', day, sep="_")
  temp.data<-temp.data[,c(1,2,4)]
  if (d == 1) {
    ril_te_fit_day_qtl<-temp.data
  }
  
  if (d > 1) {
    ril_te_fit_day_qtl<-merge(ril_te_fit_day_qtl, temp.data, by=c('genotype', 'treatment'), all=T)
  }
}

# Add misc columns
ril_te_fit_day_qtl$Obs<-c(1:nrow(ril_te_fit_day_qtl))
ril_te_fit_day_qtl$experiment<-rep("BP14", nrow(ril_te_fit_day_qtl))
ril_te_fit_day_qtl$year<-rep("2014", nrow(ril_te_fit_day_qtl))
ril_te_fit_day_qtl$plot<-rep("bellweater", nrow(ril_te_fit_day_qtl))
ril_te_fit_day_qtl$plot_id<-rep("bellweater", nrow(ril_te_fit_day_qtl))
ril_te_fit_day_qtl$sampling<-rep("bellweater", nrow(ril_te_fit_day_qtl))

# Reorder columns
ril_te_fit_day_qtl<-ril_te_fit_day_qtl[,c(20,21,22,2,23,24,1,25,3:19)]

colnames(ril_te_fit_day_qtl)[7]<-c("id")
write.csv(ril_te_fit_day_qtl, file="ril_ma_te_fit_day_qtl.csv", quote=F, row.names=F)

# Format te_residual_day for QTL pipeline
wue.l<-te_residual.ma
wue.l<-wue.l[,c('genotype', 'treatment', 'dap_i', 'ma_te_residual_day')]

days<-sort(unique(wue.l$dap_i))

ril_te_residual_day_qtl<-c()
for(d in 1:length(days)) {
  day<-days[d]
  temp.data<-wue.l[wue.l$dap_i == day, ]
  colnames(temp.data)[4]<-paste('ma_te_residual_day', day, sep="_")
  temp.data<-temp.data[,c(1,2,4)]
  if (d == 1) {
    ril_te_residual_day_qtl<-temp.data
  }
  
  if (d > 1) {
    ril_te_residual_day_qtl<-merge(ril_te_residual_day_qtl, temp.data, by=c('genotype', 'treatment'), all=T)
  }
}

# Add misc columns
ril_te_residual_day_qtl$Obs<-c(1:nrow(ril_te_residual_day_qtl))
ril_te_residual_day_qtl$experiment<-rep("BP14", nrow(ril_te_residual_day_qtl))
ril_te_residual_day_qtl$year<-rep("2014", nrow(ril_te_residual_day_qtl))
ril_te_residual_day_qtl$plot<-rep("bellweater", nrow(ril_te_residual_day_qtl))
ril_te_residual_day_qtl$plot_id<-rep("bellweater", nrow(ril_te_residual_day_qtl))
ril_te_residual_day_qtl$sampling<-rep("bellweater", nrow(ril_te_residual_day_qtl))

# Reorder columns
ril_te_residual_day_qtl<-ril_te_residual_day_qtl[,c(20,21,22,2,23,24,1,25,3:19)]

colnames(ril_te_residual_day_qtl)[7]<-c("id")
write.csv(ril_te_residual_day_qtl, file="ril_ma_te_residual_day_qtl.csv", quote=F, row.names=F)


library(ggplot2)
library(gplots)




# Combine biomass with water data
#setwd(wue_results.plant_size.dir)
#load('ril_biomass.Rdata')
#setwd(wue_results.water_use.dir)
#load('ril_water_use.Rdata')

setwd(wue_results.te.dir)



# Calculate the relative magnitude as % of maximum value
#ids<-unique(te_residual.ma$genotype)
#ma.te_residual.r<-c()
#for (i in ids){
#  temp<-te_residual.ma[te_residual.ma$genotype == i,]
#  fw_max<-max(temp$fw_biomass_day, na.rm=T)
#  temp$fw_biomass_day.r<-temp$fw_biomass_day/fw_max
#  sv_max<-max(temp$sv_area_day, na.rm=T)
#  temp$sv_area_day.r<-temp$sv_area_day/sv_max
#  fw_water_max<-max(temp$fw_wue_ratio_day, na.rm=T)
#  temp$fw_water_day.r<-temp$fw_wue_ratio_day/fw_water_max
#  # WUE RATIO
#  wue_total_max<-max(temp$wue_total, na.rm=T)
#  temp$wue_total.r<-temp$wue_total/wue_total_max
#  wue_day_max<-max(temp$wue_day, na.rm=T)
#  temp$wue_day.r<-temp$wue_day/wue_day_max
#  # TE RESIDUAL VALUES
#  ma_te_residual_total_max<-max(temp$ma_te_residual_total, na.rm=T)
#  temp$ma_te_residual_total.r<-temp$ma_te_residual_total/ma_te_residual_total_max
#  ma_te_residual_day_max<-max(temp$ma_te_residual_day, na.rm=T)
#  temp$ma_te_residual_day.r<-temp$ma_te_residual_day/ma_te_residual_day_max
#  # TE FIT VALUES
#  ma_te_fit_total_max<-max(temp$ma_te_fit_total, na.rm=T)
#  temp$ma_te_fit_total.r<-temp$ma_te_fit_total/ma_te_fit_total_max
#  ma_te_fit_day_max<-max(temp$ma_te_fit_day, na.rm=T)
#  temp$ma_te_fit_day.r<-temp$ma_te_fit_day/ma_te_fit_day_max
#  ma.te_residual.r<-rbind(ma.te_residual.r, temp)
  
#}

# Now lets combine loess summaries

#ril_loess_trait.r<-merge(ril_loess_trait.r, ma.te_residual.r[,c(1:3,16:ncol(ma.te_residual.r))], 
#                  by = c("genotype", "treatment", "dap_i"))


# Lets do some correlation analysis
#days<-sort(unique(ril_loess_trait.r$dap_i))
#corr.table.all<-c()
#corr.table.dry<-c()
#corr.table.wet<-c()
#for(d in days){
#  temp<-ril_loess_trait.r[ril_loess_trait.r$dap_i == d,]
#  temp2<-merge(temp, isotope, by=c("genotype", "treatment"))
#  temp2.d<-temp2[temp2$treatment=='dry',]
#  temp2.w<-temp2[temp2$treatment=='wet',]
#  cor.a<-cor(temp2[,4:38])
#  cor.a<-round(cor.a,2)
#  rownames(cor.a)<-paste(rownames(cor.a), d, sep="+")
#  cor.d<-cor(temp2.d[,4:38])
#  cor.d<-round(cor.d,2)
#  rownames(cor.d)<-paste(rownames(cor.d), d, sep="+")
#  cor.w<-cor(temp2.w[,4:38])
#  cor.w<-round(cor.w,2)
#  rownames(cor.w)<-paste(rownames(cor.w), d, sep="+")
#  corr.table.all<-rbind(corr.table.all, cor.a)
#  corr.table.dry<-rbind(corr.table.dry, cor.d)
#  corr.table.wet<-rbind(corr.table.wet, cor.w)
#}

## Not a surprise STDEV is zero. Think about what the model fit is doing

# Lets add columns with condition, trait and day for plotting
#trait<-c()
#days<-c()
#treatment<-rep('all', nrow(corr.table.all))
#for(r in 1:nrow(corr.table.all)){
#  temp<-rownames(corr.table.all)[r]
#  field<-strsplit(temp, "\\+")
#  t<-field[[1]][1]
#  trait<-c(trait,t)
#  d<-field[[1]][2]
#  days<-c(days,d)
#}
#corr.table.all<-cbind(trait, days, treatment, corr.table.all)
#corr.table.all<-as.data.frame(corr.table.all)
#corr.table.all$days<-as.numeric(as.character(corr.table.all$days))
#for(c in 4:ncol(corr.table.all)){
#  corr.table.all[,c]<-as.numeric(as.character(corr.table.all[,c]))
#}



# Now dry
#trait<-c()
#days<-c()
#treatment<-rep('dry', nrow(corr.table.dry))
#for(r in 1:nrow(corr.table.dry)){
#  temp<-rownames(corr.table.dry)[r]
#  field<-strsplit(temp, "\\+")
#  t<-field[[1]][1]
#  trait<-c(trait,t)
#  d<-field[[1]][2]
#  days<-c(days,d)
#}
#corr.table.dry<-cbind(trait, days, treatment, corr.table.dry)
#corr.table.dry<-as.data.frame(corr.table.dry)
#corr.table.dry$days<-as.numeric(as.character(corr.table.dry$days))
#for(c in 4:ncol(corr.table.dry)){
#  corr.table.dry[,c]<-as.numeric(as.character(corr.table.dry[,c]))
#}

# Lets add columns with condition, trait and day for plotting
#trait<-c()
#days<-c()
#treatment<-rep('wet', nrow(corr.table.wet))
#for(r in 1:nrow(corr.table.wet)){
#  temp<-rownames(corr.table.wet)[r]
#  field<-strsplit(temp, "\\+")
#  t<-field[[1]][1]
#  trait<-c(trait,t)
#  d<-field[[1]][2]
#  days<-c(days,d)
#}
#corr.table.wet<-cbind(trait, days, treatment, corr.table.wet)
#corr.table.wet<-as.data.frame(corr.table.wet)
#corr.table.wet$days<-as.numeric(as.character(corr.table.wet$days))
#for(c in 4:ncol(corr.table.wet)){
#  corr.table.wet[,c]<-as.numeric(as.character(corr.table.wet[,c]))
#}

#corr.table<-rbind(corr.table.all, corr.table.dry, corr.table.wet)

###########################################################
# Lets look at ratio derivation of TE 
###########################################################

#temp<-corr.table[corr.table$trait == 'wue_day' | corr.table$trait == 'wue_total',]

# Make a plot of how d13C correlates with TE calculated on a per day basis and TE calculated across all timepoints
# Could make a plot of all by changing temp to corr.table but the plot is too large and hard to examine

#pdf("Correlation_of_d13C_with_wue_ratio_per_day_and_cumulative_wue_ratio.pdf")
#p<-ggplot(temp, aes(x=days, y=d13C, colour=treatment)) + geom_line() + facet_wrap(~trait) + theme_bw() + ylab("Correlation with d13C") + xlab("Day after planting") + scale_colour_manual(values = c("red","orange","navy"))
#p
#dev.off()

# It is interesting that there is such a difference between wet and dry for total
# Lets examine those values more closely with a scatter plot

#days<-sort(unique(ril_loess_trait.r$dap_i))
#pdf("d13C.v.wue_ratio_day.pdf")
#for(d in days){
#  temp<-ril_loess_trait.r[ril_loess_trait.r$dap_i == d,]
#  temp<-merge(temp, isotope, by=c("genotype", "treatment"))
#  p<-ggplot(temp, aes(x=wue_day, y=d13C, colour=treatment)) + geom_point() + theme_bw() + labs(title=d) + scale_colour_manual(values = c("orange","navy"))
#  print(p)
#}
#dev.off()

#days<-sort(unique(ril_loess_trait.r$dap_i))
#pdf("d13C.v.wue_ratio_total.pdf")
#for(d in days){
#  temp<-ril_loess_trait.r[ril_loess_trait.r$dap_i == d,]
#  temp<-merge(temp, isotope, by=c("genotype", "treatment"))
#  p<-ggplot(temp, aes(x=wue_total, y=d13C, colour=treatment)) + geom_point() + theme_bw() + labs(title=d) + scale_colour_manual(values = c("orange","navy"))
#  print(p)
#}
#dev.off()



###########################################################
# Lets look at alternative derivation of TE based upon residual
###########################################################

#temp<-corr.table[corr.table$trait == 'ma_te_residual_day' | corr.table$trait == 'ma_te_residual_total',]

# Make a plot of how d13C correlates with TE calculated on a per day basis and TE calculated across all timepoints
# Could make a plot of all by changing temp to corr.table but the plot is too large and hard to examine

#pdf("Correlation_of_d13C_with_residual_calculaetd_te_per_day_and_cumulative_ET.pdf")
#p<-ggplot(temp, aes(x=days, y=d13C, colour=treatment)) + geom_line() + facet_wrap(~trait) + theme_bw() + ylab("Correlation with d13C") + xlab("Day after planting") + scale_colour_manual(values = c("red","orange","navy"))
#p
#dev.off()

# It is interesting that there is such a difference between wet and dry for total
# Lets examine those values more closely with a scatter plot

#days<-sort(unique(ril_loess_trait.r$dap_i))
#pdf("d13C.v.te_residual_day.pdf")
#for(d in days){
#  temp<-ril_loess_trait.r[ril_loess_trait.r$dap_i == d,]
#  temp<-merge(temp, isotope, by=c("genotype", "treatment"))
#  p<-ggplot(temp, aes(x=ma_te_residual_day, y=d13C, colour=treatment)) + geom_point() + theme_bw() + labs(title=d) + scale_colour_manual(values = c("orange","navy"))
#  print(p)
#}
#dev.off()

#days<-sort(unique(ril_loess_trait.r$dap_i))
#pdf("d13C.v.te_residual_total.pdf")
#for(d in days){
#  temp<-ril_loess_trait.r[ril_loess_trait.r$dap_i == d,]
#  temp<-merge(temp, isotope, by=c("genotype", "treatment"))
#  p<-ggplot(temp, aes(x=ma_te_residual_total, y=d13C, colour=treatment)) + geom_point() + theme_bw() + labs(title=d) + scale_colour_manual(values = c("orange","navy"))
#  print(p)
#}
#dev.off()







###########################################################
# Lets look at alternative derivation of TE based upon residual (Model fit this time)
###########################################################

#temp<-corr.table[corr.table$trait == 'ma_te_fit_day' | corr.table$trait == 'ma_te_fit_total',]

# Make a plot of how d13C correlates with TE calculated on a per day basis and TE calculated across all timepoints
# Could make a plot of all by changing temp to corr.table but the plot is too large and hard to examine

#pdf("Correlation_of_d13C_with_residual_calculaetd_te_per_day_and_cumulative_ET.pdf")
#p<-ggplot(temp, aes(x=days, y=d13C, colour=treatment)) + geom_line() + facet_wrap(~trait) + theme_bw() + ylab("Correlation with d13C") + xlab("Day after planting") + scale_colour_manual(values = c("red","orange","navy"))
#p
#dev.off()

# It is interesting that there is such a difference between wet and dry for total
# Lets examine those values more closely with a scatter plot

#days<-sort(unique(ril_loess_trait.r$dap_i))
#pdf("d13C.v.te_residual_day.pdf")
#for(d in days){
#  temp<-ril_loess_trait.r[ril_loess_trait.r$dap_i == d,]
#  temp<-merge(temp, isotope, by=c("genotype", "treatment"))
#  p<-ggplot(temp, aes(x=te_residual_day, y=d13C, colour=treatment)) + geom_point() + theme_bw() + labs(title=d) + scale_colour_manual(values = c("orange","navy"))
#  print(p)
#}
#dev.off()

#days<-sort(unique(ril_loess_trait.r$dap_i))
#pdf("d13C.v.te_residual_total.pdf")
#for(d in days){
#  temp<-ril_loess_trait.r[ril_loess_trait.r$dap_i == d,]
#  temp<-merge(temp, isotope, by=c("genotype", "treatment"))
#  p<-ggplot(temp, aes(x=te_residual_total, y=d13C, colour=treatment)) + geom_point() + theme_bw() + labs(title=d) + scale_colour_manual(values = c("orange","navy"))
#  print(p)
#}
#dev.off()

#save.image("te_major_axis_supplemental.Rdata")







#temp<-corr.table[corr.table$trait == 'ma_te_fit_day' | corr.table$trait == 'ma_te_fit_total',]

# Make a plot of how d13C correlates with TE calculated on a per day basis and TE calculated across all timepoints
# Could make a plot of all by changing temp to corr.table but the plot is too large and hard to examine

#pdf("Correlation_of_d13C_with_residual_calculaetd_te_per_day_and_cumulative_ET.pdf")
#p<-ggplot(temp, aes(x=days, y=d13C, colour=treatment)) + geom_line() + facet_wrap(~trait) + theme_bw() + ylab("Correlation with d13C") + xlab("Day after planting") + scale_colour_manual(values = c("red","orange","navy"))
#p
#dev.off()

#temp<-corr.table[corr.table$trait == 'ma_te_fit_day' | corr.table$trait == 'ma_te_fit_total' | corr.table$trait == 'ma_te_residual_day' | corr.table$trait == 'ma_te_residual_total',]

# Make a plot of how d13C correlates with TE calculated on a per day basis and TE calculated across all timepoints
# Could make a plot of all by changing temp to corr.table but the plot is too large and hard to examine

#pdf("Correlation_of_d13C_with_MA_te_per_day_and_cumulative_ET.pdf")
#p<-ggplot(temp, aes(x=days, y=d13C, colour=treatment)) + geom_line() + facet_wrap(~trait) + theme_bw() + ylab("Correlation with d13C") + xlab("Day after planting") + scale_colour_manual(values = c("red","orange","navy"))
#p
#dev.off()


#### SV_AREA + WATER LOST

#temp<-corr.table[corr.table$trait == 'sv_area_day' | corr.table$trait == 'sv_area' | corr.table$trait == 'water_lost_day' | corr.table$trait == 'water_lost_total',]

# Make a plot of how d13C correlates with TE calculated on a per day basis and TE calculated across all timepoints
# Could make a plot of all by changing temp to corr.table but the plot is too large and hard to examine

#pdf("Correlation_of_d13C_with_sv_area_and_water_lost.pdf")
#p<-ggplot(temp, aes(x=days, y=d13C, colour=treatment)) + geom_line() + facet_wrap(~trait) + theme_bw() + ylab("Correlation with d13C") + xlab("Day after sewing") + scale_colour_manual(values = c("red","orange","navy"))
#p
#dev.off()



##############################################################################################
# Lets do the fit for all of these model types
##############################################################################################
setwd(wue_results.te.dir)

# Define fxn to get values from lmodel2 value
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
  if(rma == 'Y') {
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
    y.ols.res<-(y-y.ols.pred)
    y.ma.res<-(y-y.ma.pred)
    y.sma.res<-(y-y.sma.pred)
    y.rma.res<-(y-y.rma.pred)
    
    # Format results
    out<-cbind(x,y,y.ols.pred,y.ma.pred,y.sma.pred,y.rma.pred,y.ols.res,y.ma.res,y.sma.res,y.rma.res)
    colnames(out)<-c('x','y','ols.pred','ma.pred','sma.pred','rma.pred','ols.res','ma.res','sma.res','rma.res')
  }
  return(out)
}


# Lets run a loop and fit model sv_area ~ water used 
# TE will be the residual of the model fit
days<-sort(unique(day_size$dap_i))
te_residual.ma<-c()
for(d in days){
  temp<-day_size[day_size$dap_i == d,]
  treats<-sort(unique(temp$treatment))
  temp.t<-c()
  for(t in treats){
    temp2<-temp[temp$treatment == t,]
    # Get fit for cumulative
    te_total.mdl<-lmodel2(sv_area~water_lost_total, data=temp2, range.y="relative", range.x="relative", nperm=99)
    fit_te_total.val<-get.lmodel2.values(te_total.mdl, 'Y')
    fit_te_total.val<-as.data.frame(fit_te_total.val)
    ma_te_fit_total<-fit_te_total.val$ma.pred
    ma_te_residual_total<-fit_te_total.val$ma.res
    sma_te_fit_total<-fit_te_total.val$sma.pred
    sma_te_residual_total<-fit_te_total.val$sma.res
    rma_te_fit_total<-fit_te_total.val$sma.pred
    rma_te_residual_total<-fit_te_total.val$sma.res
    
    # Get fit for daily/rate
    te_day.mdl<-lmodel2(sv_area_day~water_lost_day, data=temp2, range.y="relative", range.x="relative", nperm=99)
    fit_te_day.val<-get.lmodel2.values(te_day.mdl, 'Y')
    fit_te_day.val<-as.data.frame(fit_te_day.val)
    ma_te_fit_day<-fit_te_day.val$ma.pred
    ma_te_residual_day<-fit_te_day.val$ma.res
    sma_te_fit_day<-fit_te_day.val$sma.pred
    sma_te_residual_day<-fit_te_day.val$sma.res
    rma_te_fit_day<-fit_te_day.val$sma.pred
    rma_te_residual_day<-fit_te_day.val$sma.res
    
    
    temp3<-cbind(temp2,ma_te_fit_total,ma_te_residual_total,ma_te_fit_day,ma_te_residual_day,sma_te_fit_total,sma_te_residual_total,sma_te_fit_day,sma_te_residual_day,rma_te_fit_total,rma_te_residual_total,rma_te_fit_day,rma_te_residual_day)
    temp.t<-rbind(temp.t, temp3)
  }
  te_residual.ma<-rbind(te_residual.ma, temp.t)
}


###########################################################################
# MA regression
###########################################################################
# Format for te_fit_total QTL pipeline
wue.l<-te_residual.ma
wue.l<-wue.l[,c('genotype', 'treatment', 'dap_i', 'ma_te_fit_total')]

days<-sort(unique(wue.l$dap_i))

ril_te_fit_total_qtl<-c()
for(d in 1:length(days)) {
  day<-days[d]
  temp.data<-wue.l[wue.l$dap_i == day, ]
  colnames(temp.data)[4]<-paste('ma_te_fit_total', day, sep="_")
  temp.data<-temp.data[,c(1,2,4)]
  if (d == 1) {
    ril_te_fit_total_qtl<-temp.data
  }
  
  if (d > 1) {
    ril_te_fit_total_qtl<-merge(ril_te_fit_total_qtl, temp.data, by=c('genotype', 'treatment'), all=T)
  }
}

# Add misc columns
ril_te_fit_total_qtl$Obs<-c(1:nrow(ril_te_fit_total_qtl))
ril_te_fit_total_qtl$experiment<-rep("BP14", nrow(ril_te_fit_total_qtl))
ril_te_fit_total_qtl$year<-rep("2014", nrow(ril_te_fit_total_qtl))
ril_te_fit_total_qtl$plot<-rep("bellwether", nrow(ril_te_fit_total_qtl))
ril_te_fit_total_qtl$plot_id<-rep("bellwether", nrow(ril_te_fit_total_qtl))
ril_te_fit_total_qtl$sampling<-rep("bellwether", nrow(ril_te_fit_total_qtl))



# Reorder columns
ril_te_fit_total_qtl<-ril_te_fit_total_qtl[,c(20,21,22,2,23,24,1,25,3:19)]

colnames(ril_te_fit_total_qtl)[7]<-c("id")
write.csv(ril_te_fit_total_qtl, file="ril_ma_te_fit_total_qtl.csv", quote=F, row.names=F)


# Format for te_residual_total QTL pipeline
wue.l<-te_residual.ma
wue.l<-wue.l[,c('genotype', 'treatment', 'dap_i', 'ma_te_residual_total')]

days<-sort(unique(wue.l$dap_i))

ril_te_residual_total_qtl<-c()
for(d in 1:length(days)) {
  day<-days[d]
  temp.data<-wue.l[wue.l$dap_i == day, ]
  colnames(temp.data)[4]<-paste('ma_te_residual_total', day, sep="_")
  temp.data<-temp.data[,c(1,2,4)]
  if (d == 1) {
    ril_te_residual_total_qtl<-temp.data
  }
  
  if (d > 1) {
    ril_te_residual_total_qtl<-merge(ril_te_residual_total_qtl, temp.data, by=c('genotype', 'treatment'), all=T)
  }
}

# Add misc columns
ril_te_residual_total_qtl$Obs<-c(1:nrow(ril_te_residual_total_qtl))
ril_te_residual_total_qtl$experiment<-rep("BP14", nrow(ril_te_residual_total_qtl))
ril_te_residual_total_qtl$year<-rep("2014", nrow(ril_te_residual_total_qtl))
ril_te_residual_total_qtl$plot<-rep("bellweater", nrow(ril_te_residual_total_qtl))
ril_te_residual_total_qtl$plot_id<-rep("bellweater", nrow(ril_te_residual_total_qtl))
ril_te_residual_total_qtl$sampling<-rep("bellweater", nrow(ril_te_residual_total_qtl))

# Reorder columns
ril_te_residual_total_qtl<-ril_te_residual_total_qtl[,c(20,21,22,2,23,24,1,25,3:19)]

colnames(ril_te_residual_total_qtl)[7]<-c("id")
write.csv(ril_te_residual_total_qtl, file="ril_ma_te_residual_total_qtl.csv", quote=F, row.names=F)

#################################################




###########################################################################
# SMA regression
###########################################################################
# Format for te_fit_total QTL pipeline
wue.l<-te_residual.ma
wue.l<-wue.l[,c('genotype', 'treatment', 'dap_i', 'sma_te_fit_total')]

days<-sort(unique(wue.l$dap_i))

ril_te_fit_total_qtl<-c()
for(d in 1:length(days)) {
  day<-days[d]
  temp.data<-wue.l[wue.l$dap_i == day, ]
  colnames(temp.data)[4]<-paste('sma_te_fit_total', day, sep="_")
  temp.data<-temp.data[,c(1,2,4)]
  if (d == 1) {
    ril_te_fit_total_qtl<-temp.data
  }
  
  if (d > 1) {
    ril_te_fit_total_qtl<-merge(ril_te_fit_total_qtl, temp.data, by=c('genotype', 'treatment'), all=T)
  }
}

# Add misc columns
ril_te_fit_total_qtl$Obs<-c(1:nrow(ril_te_fit_total_qtl))
ril_te_fit_total_qtl$experiment<-rep("BP14", nrow(ril_te_fit_total_qtl))
ril_te_fit_total_qtl$year<-rep("2014", nrow(ril_te_fit_total_qtl))
ril_te_fit_total_qtl$plot<-rep("bellweater", nrow(ril_te_fit_total_qtl))
ril_te_fit_total_qtl$plot_id<-rep("bellweater", nrow(ril_te_fit_total_qtl))
ril_te_fit_total_qtl$sampling<-rep("bellweater", nrow(ril_te_fit_total_qtl))



# Reorder columns
ril_te_fit_total_qtl<-ril_te_fit_total_qtl[,c(20,21,22,2,23,24,1,25,3:19)]

colnames(ril_te_fit_total_qtl)[7]<-c("id")
write.csv(ril_te_fit_total_qtl, file="ril_sma_te_fit_total_qtl.csv", quote=F, row.names=F)


# Format for te_residual_total QTL pipeline
wue.l<-te_residual.ma
wue.l<-wue.l[,c('genotype', 'treatment', 'dap_i', 'sma_te_residual_total')]

days<-sort(unique(wue.l$dap_i))

ril_te_residual_total_qtl<-c()
for(d in 1:length(days)) {
  day<-days[d]
  temp.data<-wue.l[wue.l$dap_i == day, ]
  colnames(temp.data)[4]<-paste('sma_te_residual_total', day, sep="_")
  temp.data<-temp.data[,c(1,2,4)]
  if (d == 1) {
    ril_te_residual_total_qtl<-temp.data
  }
  
  if (d > 1) {
    ril_te_residual_total_qtl<-merge(ril_te_residual_total_qtl, temp.data, by=c('genotype', 'treatment'), all=T)
  }
}

# Add misc columns
ril_te_residual_total_qtl$Obs<-c(1:nrow(ril_te_residual_total_qtl))
ril_te_residual_total_qtl$experiment<-rep("BP14", nrow(ril_te_residual_total_qtl))
ril_te_residual_total_qtl$year<-rep("2014", nrow(ril_te_residual_total_qtl))
ril_te_residual_total_qtl$plot<-rep("bellweater", nrow(ril_te_residual_total_qtl))
ril_te_residual_total_qtl$plot_id<-rep("bellweater", nrow(ril_te_residual_total_qtl))
ril_te_residual_total_qtl$sampling<-rep("bellweater", nrow(ril_te_residual_total_qtl))

# Reorder columns
ril_te_residual_total_qtl<-ril_te_residual_total_qtl[,c(20,21,22,2,23,24,1,25,3:19)]

colnames(ril_te_residual_total_qtl)[7]<-c("id")
write.csv(ril_te_residual_total_qtl, file="ril_sma_te_residual_total_qtl.csv", quote=F, row.names=F)

#################################################







###########################################################################
# RMA regression
###########################################################################
# Format for te_fit_total QTL pipeline
wue.l<-te_residual.ma
wue.l<-wue.l[,c('genotype', 'treatment', 'dap_i', 'rma_te_fit_total')]

days<-sort(unique(wue.l$dap_i))

ril_te_fit_total_qtl<-c()
for(d in 1:length(days)) {
  day<-days[d]
  temp.data<-wue.l[wue.l$dap_i == day, ]
  colnames(temp.data)[4]<-paste('rma_te_fit_total', day, sep="_")
  temp.data<-temp.data[,c(1,2,4)]
  if (d == 1) {
    ril_te_fit_total_qtl<-temp.data
  }
  
  if (d > 1) {
    ril_te_fit_total_qtl<-merge(ril_te_fit_total_qtl, temp.data, by=c('genotype', 'treatment'), all=T)
  }
}

# Add misc columns
ril_te_fit_total_qtl$Obs<-c(1:nrow(ril_te_fit_total_qtl))
ril_te_fit_total_qtl$experiment<-rep("BP14", nrow(ril_te_fit_total_qtl))
ril_te_fit_total_qtl$year<-rep("2014", nrow(ril_te_fit_total_qtl))
ril_te_fit_total_qtl$plot<-rep("bellweater", nrow(ril_te_fit_total_qtl))
ril_te_fit_total_qtl$plot_id<-rep("bellweater", nrow(ril_te_fit_total_qtl))
ril_te_fit_total_qtl$sampling<-rep("bellweater", nrow(ril_te_fit_total_qtl))



# Reorder columns
ril_te_fit_total_qtl<-ril_te_fit_total_qtl[,c(20,21,22,2,23,24,1,25,3:19)]

colnames(ril_te_fit_total_qtl)[7]<-c("id")
write.csv(ril_te_fit_total_qtl, file="ril_rma_te_fit_total_qtl.csv", quote=F, row.names=F)


# Format for te_residual_total QTL pipeline
wue.l<-te_residual.ma
wue.l<-wue.l[,c('genotype', 'treatment', 'dap_i', 'rma_te_residual_total')]

days<-sort(unique(wue.l$dap_i))

ril_te_residual_total_qtl<-c()
for(d in 1:length(days)) {
  day<-days[d]
  temp.data<-wue.l[wue.l$dap_i == day, ]
  colnames(temp.data)[4]<-paste('rma_te_residual_total', day, sep="_")
  temp.data<-temp.data[,c(1,2,4)]
  if (d == 1) {
    ril_te_residual_total_qtl<-temp.data
  }
  
  if (d > 1) {
    ril_te_residual_total_qtl<-merge(ril_te_residual_total_qtl, temp.data, by=c('genotype', 'treatment'), all=T)
  }
}

# Add misc columns
ril_te_residual_total_qtl$Obs<-c(1:nrow(ril_te_residual_total_qtl))
ril_te_residual_total_qtl$experiment<-rep("BP14", nrow(ril_te_residual_total_qtl))
ril_te_residual_total_qtl$year<-rep("2014", nrow(ril_te_residual_total_qtl))
ril_te_residual_total_qtl$plot<-rep("bellweater", nrow(ril_te_residual_total_qtl))
ril_te_residual_total_qtl$plot_id<-rep("bellweater", nrow(ril_te_residual_total_qtl))
ril_te_residual_total_qtl$sampling<-rep("bellweater", nrow(ril_te_residual_total_qtl))

# Reorder columns
ril_te_residual_total_qtl<-ril_te_residual_total_qtl[,c(20,21,22,2,23,24,1,25,3:19)]

colnames(ril_te_residual_total_qtl)[7]<-c("id")
write.csv(ril_te_residual_total_qtl, file="ril_rma_te_residual_total_qtl.csv", quote=F, row.names=F)

#################################################

save.image("te_major_axis_supplemental.Rdata")

