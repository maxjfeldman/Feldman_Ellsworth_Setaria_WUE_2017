
library(ggplot2)
library(gplots)

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

## Make a directory for te model results
wue_results.te.dir<-paste(wue_results.dir, 'te/', sep="")

## If this directory doesn't exist the download is not complete
if (file.exists(wue_results.te.dir)){
  print("Directory exists!")
} else {
  dir.create(file.path(wue_results.te.dir))
}

wue_results.te.correlation.dir<-paste(wue_results.te.dir, 'correlation/', sep="")

## If this directory doesn't exist the download is not complete
if (file.exists(wue_results.te.correlation.dir)){
  print("Directory exists!")
} else {
  dir.create(file.path(wue_results.te.correlation.dir))
}

## Make a directory for isotope results
wue_results.isotope.dir<-paste(wue_results.dir, 'isotope/', sep="")

## If this directory doesn't exist the download is not complete
if (file.exists(wue_results.isotope.dir)){
  print("Directory exists!")
} else {
  dir.create(file.path(wue_results.isotope.dir))
}



# Combine biomass with water data
setwd(wue_results.plant_size.dir)
load('ril_biomass.Rdata')

setwd(wue_results.water_use.dir)
load('ril_water_use.Rdata')

setwd(wue_results.te.dir)
load('ril_transpiration_efficiency.Rdata')

setwd(wue_data.dir)
isotope_qtl<-read.csv("ril_stable_isotope_qtl.csv")
isotope<-isotope_qtl[,c(4,7,9:16)]

setwd(wue_results.te.dir)

isotope.ag<-aggregate(isotope[,3:ncol(isotope)],by=list(isotope$genotype, isotope$treatment), mean, na.action = na.pass, na.rm=TRUE)
colnames(isotope.ag)[1:2]<-c("genotype", "treatment")

## Calculate the relative magnitude as % of maximum value
ids<-unique(te_model$genotype)
te_model.r<-c()
for (i in ids){
  temp<-te_model[te_model$genotype == i,]
  fw_max<-max(temp$fw_biomass_day, na.rm=T)
  temp$fw_biomass_day.r<-temp$fw_biomass_day/fw_max
  sv_max<-max(temp$sv_area_day, na.rm=T)
  temp$sv_area_day.r<-temp$sv_area_day/sv_max
  fw_wue_ratio_max<-max(temp$fw_wue_ratio_day, na.rm=T)
  temp$fw_wue_ratio_day.r<-temp$fw_wue_ratio_day/fw_wue_ratio_max
  ## WUE RATIO
  wue_total_max<-max(temp$wue_total, na.rm=T)
  temp$wue_total.r<-temp$wue_total/wue_total_max
  wue_day_max<-max(temp$wue_day, na.rm=T)
  temp$wue_day.r<-temp$wue_day/wue_day_max
  ## TE RESIDUAL VALUES
  te_residual_total_max<-max(temp$te_residual_total, na.rm=T)
  temp$te_residual_total.r<-temp$te_residual_total/te_residual_total_max
  te_residual_day_max<-max(temp$te_residual_day, na.rm=T)
  temp$te_residual_day.r<-temp$te_residual_day/te_residual_day_max
  ## TE FIT VALUES
  te_fit_total_max<-max(temp$te_fit_total, na.rm=T)
  temp$te_fit_total.r<-temp$te_fit_total/te_fit_total_max
  te_fit_day_max<-max(temp$te_fit_day, na.rm=T)
  temp$te_fit_day.r<-temp$te_fit_day/te_fit_day_max
  te_model.r<-rbind(te_model.r, temp)
  
}

## Now lets combine loess summaries
ril_loess_trait.r<-merge(ril_loess_trait.r, te_model.r[,c(1:3,16:ncol(te_model.r))], by=c("genotype", "treatment", "dap_i"))

p<-ggplot(ril_loess_trait.r, aes(x=dap_i, y=sv_area_day, colour=treatment)) + geom_smooth()
p<-ggplot(ril_loess_trait.r, aes(x=dap_i, y=fw_wue_ratio_day, colour=treatment)) + geom_smooth()
p<-ggplot(ril_loess_trait.r, aes(x=dap_i, y=te_residual_total, colour=treatment)) + geom_smooth()
p<-ggplot(ril_loess_trait.r, aes(x=dap_i, y=te_residual_day, colour=treatment)) + geom_smooth()


## Lets do some correlation analysis
days<-sort(unique(ril_loess_trait.r$dap_i))
corr.table.all<-c()
corr.table.dry<-c()
corr.table.wet<-c()
for(d in days){
  temp<-ril_loess_trait.r[ril_loess_trait.r$dap_i == d,]
  temp2<-merge(temp, isotope.ag, by=c("genotype", "treatment"))
  temp2.d<-temp2[temp2$treatment=='dry',]
  temp2.w<-temp2[temp2$treatment=='wet',]
  cor.a<-cor(temp2[,4:ncol(temp2)])
  cor.a<-round(cor.a,2)
  rownames(cor.a)<-paste(rownames(cor.a), d, sep="+")
  cor.d<-cor(temp2.d[,4:ncol(temp2)])
  cor.d<-round(cor.d,2)
  rownames(cor.d)<-paste(rownames(cor.d), d, sep="+")
  cor.w<-cor(temp2.w[,4:ncol(temp2)])
  cor.w<-round(cor.w,2)
  rownames(cor.w)<-paste(rownames(cor.w), d, sep="+")
  corr.table.all<-rbind(corr.table.all, cor.a)
  corr.table.dry<-rbind(corr.table.dry, cor.d)
  corr.table.wet<-rbind(corr.table.wet, cor.w)
}

## Lets add columns with condition, trait and day for plotting
trait<-c()
days<-c()
treatment<-rep('all', nrow(corr.table.all))
for(r in 1:nrow(corr.table.all)){
  temp<-rownames(corr.table.all)[r]
  field<-strsplit(temp, "\\+")
  t<-field[[1]][1]
  trait<-c(trait,t)
  d<-field[[1]][2]
  days<-c(days,d)
}
corr.table.all<-cbind(trait, days, treatment, corr.table.all)
corr.table.all<-as.data.frame(corr.table.all)
corr.table.all$days<-as.numeric(as.character(corr.table.all$days))
for(c in 4:ncol(corr.table.all)){
  corr.table.all[,c]<-as.numeric(as.character(corr.table.all[,c]))
}



## Now dry
trait<-c()
days<-c()
treatment<-rep('dry', nrow(corr.table.dry))
for(r in 1:nrow(corr.table.dry)){
  temp<-rownames(corr.table.dry)[r]
  field<-strsplit(temp, "\\+")
  t<-field[[1]][1]
  trait<-c(trait,t)
  d<-field[[1]][2]
  days<-c(days,d)
}
corr.table.dry<-cbind(trait, days, treatment, corr.table.dry)
corr.table.dry<-as.data.frame(corr.table.dry)
corr.table.dry$days<-as.numeric(as.character(corr.table.dry$days))
for(c in 4:ncol(corr.table.dry)){
  corr.table.dry[,c]<-as.numeric(as.character(corr.table.dry[,c]))
}

## Lets add columns with condition, trait and day for plotting
trait<-c()
days<-c()
treatment<-rep('wet', nrow(corr.table.wet))
for(r in 1:nrow(corr.table.wet)){
  temp<-rownames(corr.table.wet)[r]
  field<-strsplit(temp, "\\+")
  t<-field[[1]][1]
  trait<-c(trait,t)
  d<-field[[1]][2]
  days<-c(days,d)
}
corr.table.wet<-cbind(trait, days, treatment, corr.table.wet)
corr.table.wet<-as.data.frame(corr.table.wet)
corr.table.wet$days<-as.numeric(as.character(corr.table.wet$days))
for(c in 4:ncol(corr.table.wet)){
  corr.table.wet[,c]<-as.numeric(as.character(corr.table.wet[,c]))
}

corr.table<-rbind(corr.table.all, corr.table.dry, corr.table.wet)


###########################################################
## Lets look at WUE 
###########################################################

temp<-corr.table[corr.table$trait == 'wue_day' | corr.table$trait == 'wue_total',]

## Make a plot of how d13C correlates with WUE calculated on a per day basis and WUE calculated across all timepoints
## Could make a plot of all by changing temp to corr.table but the plot is too large and hard to examine

setwd(wue_results.te.correlation.dir)

pdf("Correlation_of_d13C_with_wue_ratio_per_day_and_cumulative_wue_ratio.pdf")
p<-ggplot(temp, aes(x=days, y=d13C, colour=treatment)) + geom_line() + facet_wrap(~trait) + theme_bw() + ylab("Correlation with d13C") + xlab("Day after planting") + scale_colour_manual(values = c("red","orange","navy"))
p
dev.off()

## It is interesting that there is such a difference between wet and dry for total
## Lets examine those values more closely with a scatter plot

days<-sort(unique(ril_loess_trait.r$dap_i))
pdf("d13C.v.wue_ratio_day.pdf")
for(d in days){
  temp<-ril_loess_trait.r[ril_loess_trait.r$dap_i == d,]
  temp<-merge(temp, isotope, by=c("genotype", "treatment"))
  p<-ggplot(temp, aes(x=wue_day, y=d13C, colour=treatment)) + geom_point() + theme_bw() + labs(title=d) + scale_colour_manual(values = c("orange","navy"))
  print(p)
}
dev.off()

days<-sort(unique(ril_loess_trait.r$dap_i))
pdf("d13C.v.wue_ratio_total.pdf")
for(d in days){
  temp<-ril_loess_trait.r[ril_loess_trait.r$dap_i == d,]
  temp<-merge(temp, isotope, by=c("genotype", "treatment"))
  p<-ggplot(temp, aes(x=wue_total, y=d13C, colour=treatment)) + geom_point() + theme_bw() + labs(title=d) + scale_colour_manual(values = c("orange","navy"))
  print(p)
}
dev.off()

###########################################################
## Lets look at TE residual
###########################################################

temp<-corr.table[corr.table$trait == 'te_residual_day' | corr.table$trait == 'te_residual_total',]

setwd(wue_results.te.correlation.dir)


## Make a plot of how d13C correlates with TE residual calculated on a per day basis and TE residual calculated across all timepoints
## Could make a plot of all by changing temp to corr.table but the plot is too large and hard to examine

pdf("Correlation_of_d13C_with_TE_residual_per_day_and_cumulative_TE_residual.pdf")
p<-ggplot(temp, aes(x=days, y=d13C, colour=treatment)) + geom_line() + facet_wrap(~trait) + theme_bw() + ylab("Correlation with d13C") + xlab("Day after planting") + scale_colour_manual(values = c("red","orange","navy"))
p
dev.off()

## It is interesting that there is such a difference between wet and dry for total
## Lets examine those values more closely with a scatter plot

days<-sort(unique(ril_loess_trait.r$dap_i))
pdf("d13C.v.te_residual_day.pdf")
for(d in days){
  temp<-ril_loess_trait.r[ril_loess_trait.r$dap_i == d,]
  temp<-merge(temp, isotope, by=c("genotype", "treatment"))
  p<-ggplot(temp, aes(x=te_residual_day, y=d13C, colour=treatment)) + geom_point() + theme_bw() + labs(title=d) + scale_colour_manual(values = c("orange","navy"))
  print(p)
}
dev.off()

days<-sort(unique(ril_loess_trait.r$dap_i))
pdf("d13C.v.te_residual_total.pdf")
for(d in days){
  temp<-ril_loess_trait.r[ril_loess_trait.r$dap_i == d,]
  temp<-merge(temp, isotope, by=c("genotype", "treatment"))
  p<-ggplot(temp, aes(x=te_residual_total, y=d13C, colour=treatment)) + geom_point() + theme_bw() + labs(title=d) + scale_colour_manual(values = c("orange","navy"))
  print(p)
}
dev.off()



###########################################################
## Lets look at fit of TE model
###########################################################

temp<-corr.table[corr.table$trait == 'te_fit_day' | corr.table$trait == 'te_fit_total',]

## Make a plot of how d13C correlates with TE calculated on a per day basis and TE calculated across all timepoints
## Could make a plot of all by changing temp to corr.table but the plot is too large and hard to examine

pdf("Correlation_of_d13C_with_TE_fit_per_day_and_cumulative_TE_fit.pdf")
p<-ggplot(temp, aes(x=days, y=d13C, colour=treatment)) + geom_line() + facet_wrap(~trait) + theme_bw() + ylab("Correlation with d13C") + xlab("Day after planting") + scale_colour_manual(values = c("red","orange","navy"))
p
dev.off()

days<-sort(unique(ril_loess_trait.r$dap_i))
pdf("d13C.v.te_fit_day.pdf")
for(d in days){
  temp<-ril_loess_trait.r[ril_loess_trait.r$dap_i == d,]
  temp<-merge(temp, isotope, by=c("genotype", "treatment"))
  p<-ggplot(temp, aes(x=te_fit_day, y=d13C, colour=treatment)) + geom_point() + theme_bw() + labs(title=d) + scale_colour_manual(values = c("orange","navy"))
  print(p)
}
dev.off()

days<-sort(unique(ril_loess_trait.r$dap_i))
pdf("d13C.v.te_fit_total.pdf")
for(d in days){
  temp<-ril_loess_trait.r[ril_loess_trait.r$dap_i == d,]
  temp<-merge(temp, isotope, by=c("genotype", "treatment"))
  p<-ggplot(temp, aes(x=te_fit_total, y=d13C, colour=treatment)) + geom_point() + theme_bw() + labs(title=d) + scale_colour_manual(values = c("orange","navy"))
  print(p)
}
dev.off()


###########################################################
## Lets look at correlation with sv_area
###########################################################

temp<-corr.table[corr.table$trait == 'sv_area_day' | corr.table$trait == 'sv_area',]

## Make a plot of how d13C correlates with plant size on a per day basis (rate) and plant size calculated across all timepoints
## Could make a plot of all by changing temp to corr.table but the plot is too large and hard to examine

pdf("Correlation_of_d13C_with_sv_area_per_day_and_cumulative_sv_area.pdf")
p<-ggplot(temp, aes(x=days, y=d13C, colour=treatment)) + geom_line() + facet_wrap(~trait) + theme_bw() + ylab("Correlation with d13C") + xlab("Day after sewing") + scale_colour_manual(values = c("red","orange","navy"))
p
dev.off()

## It is interesting that there is such a difference between wet and dry for total
## Lets examine those values more closely with a scatter plot

days<-sort(unique(ril_loess_trait.r$dap_i))
pdf("d13C.v.sv_area.pdf")
for(d in days){
  temp<-ril_loess_trait.r[ril_loess_trait.r$dap_i == d,]
  temp<-merge(temp, isotope, by=c("genotype", "treatment"))
  p<-ggplot(temp, aes(x=sv_area_day, y=d13C, colour=treatment)) + geom_point() + theme_bw() + labs(title=d) + scale_colour_manual(values = c("orange","navy"))
  print(p)
}
dev.off()

days<-sort(unique(ril_loess_trait.r$dap_i))
pdf("d13C.v.sv_area.pdf")
for(d in days){
  temp<-ril_loess_trait.r[ril_loess_trait.r$dap_i == d,]
  temp<-merge(temp, isotope, by=c("genotype", "treatment"))
  p<-ggplot(temp, aes(x=sv_area, y=d13C, colour=treatment)) + geom_point() + theme_bw() + labs(title=d) + scale_colour_manual(values = c("orange","navy"))
  print(p)
}
dev.off()



###########################################################
## Lets look at correlation with water_lost
###########################################################

temp<-corr.table[corr.table$trait == 'water_lost_day' | corr.table$trait == 'water_lost_total',]

# Make a plot of how d13C correlates with water use calculated on a per day basis and water use calculated across all timepoints
# Could make a plot of all by changing temp to corr.table but the plot is too large and hard to examine

pdf("Correlation_of_d13C_with_water_lost_per_day_and_cumulative_water_lost.pdf")
p<-ggplot(temp, aes(x=days, y=d13C, colour=treatment)) + geom_line() + facet_wrap(~trait) + theme_bw() + ylab("Correlation with d13C") + xlab("Day after planting") + scale_colour_manual(values = c("red","orange","navy"))
p
dev.off()


pdf("water_lost_day_v.d13C_scatter_plot.pdf")
p<-ggplot(temp, aes(x=water_lost_day, y=d13C, color=treatment)) + geom_point(size=2) + theme_bw()
p2<- p + theme(legend.position="none", axis.text=element_text(size=20), axis.title=element_text(size=24,face="bold"), plot.title = element_text(size = 24, face = "bold")) + scale_color_manual(values=c("red", "orange", "navy"))+ xlab("Water lost (mL)")
print(p2)
dev.off()


pdf("water_lost_total_v.d13C_scatter_plot.pdf")
p<-ggplot(temp, aes(x=water_lost_total, y=d13C, color=treatment)) + geom_point(size=2) + theme_bw()
p2<- p + theme(legend.position="none", axis.text=element_text(size=20), axis.title=element_text(size=24,face="bold"), plot.title = element_text(size = 24, face = "bold")) + scale_color_manual(values=c("red", "orange", "navy"))+ xlab("Water lost (mL)")
print(p2)
dev.off()


setwd(wue_results.te.dir)
save.image("ril_transpiration_analysis.Rdata")







