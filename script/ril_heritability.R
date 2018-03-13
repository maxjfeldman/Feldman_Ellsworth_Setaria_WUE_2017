library(ggplot2)
library(gplots)
library(lme4)
library(lmomco)

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

wue_results.h2.dir<-paste(wue_results.dir, 'heritability/', sep="")

if (file.exists(wue_results.h2.dir)){
  print("Directory exists!")
} else {
  dir.create(file.path(wue_results.h2.dir))
}

##### END ##### 

## Combine biomass with water data
setwd(wue_results.plant_size.dir)
load('ril_biomass.Rdata')
setwd(wue_results.water_use.dir)
load('ril_water_use.Rdata')
setwd(wue_results.te.dir)
load('ril_transpiration_efficiency.Rdata')

## Write heritability results to the correct directory
setwd(wue_results.h2.dir)


## The first is to model each plant_id to get the value of the trait using loess regression
daily.sv_area<-loess.fit.for.h2(vis.wue.17, 'sv_area')

days<-sort(unique(daily.sv_area$dap_i))
days<-sort(as.numeric(as.character(days)))

## Calculate how much each plant id changes in size at each day based on projected values
## This is basically so we can calculate heritability for the rate statistic
day_size_id<-c()
for (d in 1:length(days)){
  day<-days[d]
  today<-daily.sv_area[daily.sv_area$dap_i == day,]
  print("day is:")
  print(day)
  if (d == 1) {
    print("In d == 1 loop!")
    plantbarcode<-today$plantbarcode
    genotype<-today$genotype
    treatment<-today$treatment
    dap_i<-as.character(today$dap_i)
    sv_area_day<-today$sv_area
    day_size_id<-cbind(as.character(plantbarcode), genotype, treatment, dap_i, sv_area_day)
    colnames(day_size_id)[1]<-c("plantbarcode")
  }
  
  if (d > 1) {
    print('now in other loop')
    yesterday<-daily.sv_area[daily.sv_area$dap_i == (as.numeric(as.character(day)) - 1),]
    yesterday<-yesterday[,-c(4)]
    colnames(yesterday)[c(4)]<-c('sv_area.y')
    today<-merge(today[,c(1:5)], yesterday[,c(1:4)], by=c("plantbarcode", "genotype", "treatment"))
    plantbarcode<-today$plantbarcode
    genotype<-today$genotype 
    treatment<-today$treatment
    dap_i<-as.character(today$dap_i)
    today$sv_area<-as.numeric(as.character(today$sv_area))
    today$sv_area.y<-as.numeric(as.character(today$sv_area.y))
    sv_area_day<-today$sv_area - today$sv_area.y
    temp<-cbind(as.character(plantbarcode), genotype, treatment, dap_i, sv_area_day)
    colnames(temp)[1]<-c("plantbarcode")
    temp<-as.data.frame(temp)
    day_size_id<-rbind(day_size_id, temp)
  }
}

day_size_id<-as.data.frame(day_size_id)
day_size_id$dap_i<-as.numeric(as.character(day_size_id$dap_i))
day_size_id$sv_area_day<-as.numeric(as.character(day_size_id$sv_area_day))
day_size_id<-day_size_id[day_size_id$dap_i > 16,]


## Merge the cumulative and daily values
te_based_on_id<-merge(daily.sv_area, day_size_id, by=c("plantbarcode", "genotype", "treatment", "dap_i"))
water.m$plantbarcode<-as.character(water.m$plantbarcode)
water.m$genotype<-as.character(water.m$genotype)
water.m$treatment<-as.character(water.m$treatment)
water.m$dap_i<-as.numeric(as.character(water.m$dap_i))
water.m<-water.m[water.m$dap_i > 16, ]

## Keep a copy before doing merge
te_based_on_id.old<-te_based_on_id

## Merge in the data for water used
te_based_on_id<-merge(te_based_on_id, water.m, by=c("plantbarcode", "genotype", "treatment", "dap_i"), all=T)


te_based_on_id<-te_based_on_id[,c(1:5,9,11:12)]
colnames(te_based_on_id)[c(5,7)]<-c("sv_area_total", "water_lost_day")

## Change factors to numeric
for(c in 4:ncol(te_based_on_id)){
  te_based_on_id[,c]<-as.numeric(as.character(te_based_on_id[,c]))
}



## Lets get the te_model fit
days<-sort(unique(te_based_on_id$dap_i))
te_model<-c()
for(d in days){
  temp<-te_based_on_id[te_based_on_id$dap_i == d,]
  temp<-temp[complete.cases(temp),]
  treats<-sort(unique(temp$treatment))
  temp.t<-c()
  for(t in treats){
    temp2<-temp[temp$treatment == t,]
    te_total.mdl<-lm(sv_area_total~water_lost_total, data=temp2)
    te_fit_total<-predict(te_total.mdl)
    te_residual_total<-residuals(te_total.mdl)
    te_day.mdl<-lm(sv_area_day~water_lost_day, data=temp2)
    te_fit_day<-predict(te_day.mdl)
    te_residual_day<-residuals(te_day.mdl)
    temp3<-cbind(temp2,te_fit_total,te_residual_total,te_fit_day,te_residual_day)
    temp.t<-rbind(temp.t, temp3)
  }
  te_model<-rbind(te_model, temp.t)
}

## Lets add a column for wue_ratio
te_model$wue_ratio_total<-te_model$sv_area_total/te_model$water_lost_total
te_model$wue_ratio_day<-te_model$sv_area_day/te_model$water_lost_day



## Need to convert to QTL format 

## Here we need to reformat for get_h2() fxn
te_model$experiment<-rep('DP14', nrow(te_model))
te_model$year<-rep('2014', nrow(te_model))
te_model$plot<-rep('bellweather', nrow(te_model))
te_model$plot_id<-rep('bellweather', nrow(te_model))
te_model<-te_model[,c(15:16,3,17,18,2,1,5:14,4)]



## First lets work on sv_water_total (TE calculated cumulatively throuhout experiment)
te_residual_total_only<-te_model[,c(1:7,13,18)]
days<-sort(unique(te_residual_total_only$dap_i))
te_residual_total_4_h2<-c()
for (d in 1:length(days)){
  day<-days[d]
  temp<-te_residual_total_only[te_residual_total_only$dap_i == day,]
  temp<-temp[,c(1:8)]
  colnames(temp)[8]<-paste(colnames(temp)[8], day, sep="_")
  if(d == 1) {
    te_residual_total_4_h2<-temp
  }
  if (d > 1){
    te_residual_total_4_h2<-merge(te_residual_total_4_h2, temp, by=c("experiment", "year", "treatment", "plot", "plot_id", "genotype", "plantbarcode"), all=T)  
  }
}

# Add the final column and adjust
te_residual_total_4_h2$Obs<-c(1:nrow(te_residual_total_4_h2))
te_residual_total_4_h2<-te_residual_total_4_h2[,c(25,1:24)]

# Calculate H2
te_residual_total_h2<-get_h2(te_residual_total_4_h2)
te_residual_total_h2_in_treatment<-get_h2_in_treatment(te_residual_total_4_h2)

pdf('te_residual_total_h2.pdf')

plot(te_residual_total_h2~c(17:33), type="l", ylim=c(0,1), col="red", ylab="H2 of cumulative biomass~water residual", xlab="Days after sewing")
lines(te_residual_total_h2_in_treatment[1,1:17]~c(17:33), ylim=c(0,1), col="orange")
lines(te_residual_total_h2_in_treatment[1,18:34]~c(17:33), ylim=c(0,1), col="navy")
#legend("topright",legend=c("all", "dry", "wet"), col=c("red","orange","navy"), lty=c(1,1,1), bty=c('n'), cex=0.5)

dev.off()

te_fit_total_only<-te_model[,c(1:7,12,18)]
days<-sort(unique(te_fit_total_only$dap_i))
te_fit_total_4_h2<-c()
for (d in 1:length(days)){
  day<-days[d]
  temp<-te_fit_total_only[te_fit_total_only$dap_i == day,]
  temp<-temp[,c(1:8)]
  colnames(temp)[8]<-paste(colnames(temp)[8], day, sep="_")
  if(d == 1) {
    te_fit_total_4_h2<-temp
  }
  if (d > 1){
    te_fit_total_4_h2<-merge(te_fit_total_4_h2, temp, by=c("experiment", "year", "treatment", "plot", "plot_id", "genotype", "plantbarcode"), all=T)  
  }
}

# Add the final column and adjust
te_fit_total_4_h2$Obs<-c(1:nrow(te_fit_total_4_h2))
te_fit_total_4_h2<-te_fit_total_4_h2[,c(25,1:24)]

# Calculate H2
te_fit_total_h2<-get_h2(te_fit_total_4_h2)
te_fit_total_h2_in_treatment<-get_h2_in_treatment(te_fit_total_4_h2)

pdf('te_fit_total_h2.pdf')

plot(te_fit_total_h2~c(17:33), type="l", ylim=c(0,1), col="red", ylab="H2 of cumulative biomass~water fit", xlab="Days after sewing")
lines(te_fit_total_h2_in_treatment[1,1:17]~c(17:33), ylim=c(0,1), col="orange")
lines(te_fit_total_h2_in_treatment[1,18:34]~c(17:33), ylim=c(0,1), col="navy")
#legend("topright",legend=c("all", "dry", "wet"), col=c("red","orange","navy"), lty=c(1,1,1), bty=c('n'), cex=0.5)

dev.off()


## Calculate H2 for WUE ratio total

wue_ratio_total_only<-te_model[,c(1:7,16,18)]
days<-sort(unique(wue_ratio_total_only$dap_i))
wue_ratio_total_4_h2<-c()
for (d in 1:length(days)){
  day<-days[d]
  temp<-wue_ratio_total_only[wue_ratio_total_only$dap_i == day,]
  temp<-temp[,c(1:8)]
  colnames(temp)[8]<-paste(colnames(temp)[8], day, sep="_")
  if(d == 1) {
    wue_ratio_total_4_h2<-temp
  }
  if (d > 1){
    wue_ratio_total_4_h2<-merge(wue_ratio_total_4_h2, temp, by=c("experiment", "year", "treatment", "plot", "plot_id", "genotype", "plantbarcode"), all=T)  
  }
}

# Add the final column and adjust
wue_ratio_total_4_h2$Obs<-c(1:nrow(wue_ratio_total_4_h2))
wue_ratio_total_4_h2<-wue_ratio_total_4_h2[,c(25,1:24)]

# Calculate H2
wue_ratio_total_h2<-get_h2(wue_ratio_total_4_h2)
wue_ratio_total_h2_in_treatment<-get_h2_in_treatment(wue_ratio_total_4_h2)

pdf('wue_ratio_total_h2.pdf')

plot(wue_ratio_total_h2~c(17:33), type="l", ylim=c(0,1), col="red", ylab="H2 of WUE ratio total", xlab="Days after sewing")
lines(wue_ratio_total_h2_in_treatment[1,1:17]~c(17:33), ylim=c(0,1), col="orange")
lines(wue_ratio_total_h2_in_treatment[1,18:34]~c(17:33), ylim=c(0,1), col="navy")
#legend("topright",legend=c("all", "dry", "wet"), col=c("red","orange","navy"), lty=c(1,1,1), bty=c('n'), cex=0.5)

dev.off()

## Calculate H2 for plant size 
sv_area_only<-te_model[,c(1:7,8,18)]
days<-sort(unique(sv_area_only$dap_i))
sv_area_4_h2<-c()
for (d in 1:length(days)){
  day<-days[d]
  temp<-sv_area_only[sv_area_only$dap_i == day,]
  temp<-temp[,c(1:8)]
  colnames(temp)[8]<-paste(colnames(temp)[8], day, sep="_")
  if(d == 1) {
    sv_area_4_h2<-temp
  }
  if (d > 1){
    sv_area_4_h2<-merge(sv_area_4_h2, temp, by=c("experiment", "year", "treatment", "plot", "plot_id", "genotype", "plantbarcode"), all=T)  
  }
}

# Add the final column and adjust
sv_area_4_h2$Obs<-c(1:nrow(sv_area_4_h2))
sv_area_4_h2<-sv_area_4_h2[,c(25,1:24)]

# Calculate H2
sv_area_h2<-get_h2(sv_area_4_h2)
sv_area_h2_in_treatment<-get_h2_in_treatment(sv_area_4_h2)

pdf('sv_area_h2.pdf')

plot(sv_area_h2~c(17:33), type="l", ylim=c(0,1), col="red", ylab="H2 of sv_area (pixels)", xlab="Days after planting")
lines(sv_area_h2_in_treatment[1,1:17]~c(17:33), ylim=c(0,1), col="orange")
lines(sv_area_h2_in_treatment[1,18:34]~c(17:33), ylim=c(0,1), col="navy")
#legend("topright",legend=c("all", "dry", "wet"), col=c("red","orange","navy"), lty=c(1,1,1), bty=c('n'), cex=0.5)

dev.off()

##  Now water used total
water_lost_total_only<-te_model[,c(1:7,11,18)]
days<-sort(unique(water_lost_total_only$dap_i))
water_lost_total_4_h2<-c()
for (d in 1:length(days)){
  day<-days[d]
  temp<-water_lost_total_only[water_lost_total_only$dap_i == day,]
  temp<-temp[,c(1:8)]
  colnames(temp)[8]<-paste(colnames(temp)[8], day, sep="_")
  if(d == 1) {
    water_lost_total_4_h2<-temp
  }
  if (d > 1){
    water_lost_total_4_h2<-merge(water_lost_total_4_h2, temp, by=c("experiment", "year", "treatment", "plot", "plot_id", "genotype", "plantbarcode"), all=T)  
  }
}

# Add the final column and adjust
water_lost_total_4_h2$Obs<-c(1:nrow(water_lost_total_4_h2))
water_lost_total_4_h2<-water_lost_total_4_h2[,c(25,1:24)]

# Calculate H2
water_lost_total_h2<-get_h2(water_lost_total_4_h2)
water_lost_total_h2_in_treatment<-get_h2_in_treatment(water_lost_total_4_h2)

pdf('water_lost_total_h2.pdf')

plot(water_lost_total_h2~c(17:33), type="l", ylim=c(0,1), col="red", ylab="H2 of water_lost_total (pixels)", xlab="Days after planting")
lines(water_lost_total_h2_in_treatment[1,1:17]~c(17:33), ylim=c(0,1), col="orange")
lines(water_lost_total_h2_in_treatment[1,18:34]~c(17:33), ylim=c(0,1), col="navy")
#legend("topright",legend=c("all", "dry", "wet"), col=c("red","orange","navy"), lty=c(1,1,1), bty=c('n'), cex=0.5)

dev.off()

############################
## Cool, lets make some more plots using ggplot2
############################


## Make a data frame for each subset
dap_i<-c(17:33)
comp.h2<-cbind(dap_i,rep('all',length(dap_i)), sv_area_h2, water_lost_total_h2,te_residual_total_h2, te_fit_total_h2, wue_ratio_total_h2)
colnames(comp.h2)<-c("dap_i", "treatment", "sv_area_total", "water_lost_total", "te_residual_total", "te_fit_total", "wue_ratio_total")
## Dry
dry.comp.h2<-cbind(dap_i,rep('dry', length(dap_i)),sv_area_h2_in_treatment[1,1:17],water_lost_total_h2_in_treatment[1,1:17],te_residual_total_h2_in_treatment[1,1:17], te_fit_total_h2_in_treatment[1,1:17], wue_ratio_total_h2_in_treatment[1,1:17])
colnames(dry.comp.h2)<-c("dap_i", "treatment", "sv_area_total", "water_lost_total", "te_residual_total", "te_fit_total", "wue_ratio_total")
## Wet
wet.comp.h2<-cbind(dap_i,rep('wet', length(dap_i)),sv_area_h2_in_treatment[1,18:34],water_lost_total_h2_in_treatment[1,18:34],te_residual_total_h2_in_treatment[1,18:34], te_fit_total_h2_in_treatment[1,18:34], wue_ratio_total_h2_in_treatment[1,18:34])
colnames(wet.comp.h2)<-c("dap_i", "treatment", "sv_area_total", "water_lost_total", "te_residual_total", "te_fit_total", "wue_ratio_total")

## Combine them
all.comp.h2<-rbind(comp.h2, dry.comp.h2, wet.comp.h2)

rownames(all.comp.h2)<-c(1:nrow(all.comp.h2))

gg.comp.h2<-c()
for(i in 3:ncol(all.comp.h2)){
  temp<-as.data.frame(all.comp.h2[,c(i)])
  temp$trait<-rep(colnames(all.comp.h2)[i], nrow(all.comp.h2))
  temp<-cbind(temp, all.comp.h2[,1:2])
  colnames(temp)<-c("value", "trait", "dap_i", 'treatment')
  gg.comp.h2<-rbind(gg.comp.h2, temp)
}

gg.comp.h2$value<-as.numeric(as.character(gg.comp.h2$value))
gg.comp.h2$dap_i<-as.numeric(as.character(gg.comp.h2$dap_i))


pdf("ril_heritability_ggplot2.pdf")
p<-ggplot(gg.comp.h2, aes(x=dap_i, y=value, color=treatment)) + geom_line(size=2)
g<-p+ scale_color_manual(values=c("red", "orange", "navy")) + ylab("Broad sense H2") + xlab("Days after planting") + facet_wrap(~trait)
d<-g +theme_bw() + ylim(0,1) + theme(legend.position="none", axis.text=element_text(size=20), axis.title=element_text(size=24,face="bold"), plot.title = element_text(size = 24, face = "bold"))
d
dev.off()

save.image('ril_heritability.Rdata')


write.csv(all.comp.h2, "cumulative_broadsense_H2.csv", quote=F, row.names=F)



###################################################
## Now lets work with the daily/rate values
###################################################

## TE RESIDUAL DAILY/RATE

te_residual_daily_only <- te_model[, c("experiment", "year", "plot", "plot_id", "plantbarcode", "genotype", "treatment", "dap_i", 
                                          "water_lost_day", "sv_area_day", "wue_ratio_day", "te_fit_day", "te_residual_day")]
days <- sort(unique(te_residual_daily_only$dap_i))
te_residual_daily_4_h2 <- c()
for (d in 1:length(days)){
  day <- days[d]
  temp <- te_residual_daily_only[te_residual_daily_only$dap_i == day,]
  temp <- temp[,c(1:7, 13)]
  colnames(temp)[8] <- paste(colnames(temp)[8], day, sep="_")
  if(d == 1) {
    te_residual_daily_4_h2<-temp
  }
  if (d > 1){
    te_residual_daily_4_h2 <- merge(te_residual_daily_4_h2, temp, by = c("experiment", "year", "treatment", "plot", 
                                                                         "plot_id", "genotype", "plantbarcode"), all = T)  
  }
}

## Add the final column and adjust
te_residual_daily_4_h2$Obs<-c(1:nrow(te_residual_daily_4_h2))
te_residual_daily_4_h2<-te_residual_daily_4_h2[,c(25,1:24)]

## Calculate H2
te_residual_daily_h2<-get_h2(te_residual_daily_4_h2)
te_residual_daily_h2_in_treatment<-get_h2_in_treatment(te_residual_daily_4_h2)


pdf('te_residual_daily_h2.pdf')

plot(te_residual_daily_h2~c(17:33), type="l", ylim=c(0,1), col="red", ylab="H2 of cumulative biomass~water residual", xlab="Days after sewing")
lines(te_residual_daily_h2_in_treatment[1,1:17]~c(17:33), ylim=c(0,1), col="orange")
lines(te_residual_daily_h2_in_treatment[1,18:34]~c(17:33), ylim=c(0,1), col="navy")
#legend("topright",legend=c("all", "dry", "wet"), col=c("red","orange","navy"), lty=c(1,1,1), bty=c('n'), cex=0.5)

dev.off()

## TE FIT DAILY/RATE

te_fit_daily_only <- te_model[, c("experiment", "year", "plot", "plot_id", "plantbarcode", "genotype", "treatment", "dap_i", 
                                     "water_lost_day", "sv_area_day", "wue_ratio_day", "te_fit_day", "te_residual_day")]
days<-sort(unique(te_fit_daily_only$dap_i))
te_fit_daily_4_h2<-c()
for (d in 1:length(days)){
  day<-days[d]
  temp<-te_fit_daily_only[te_fit_daily_only$dap_i == day,]
  temp<-temp[,c(1:7, 12)]
  colnames(temp)[8]<-paste(colnames(temp)[8], day, sep="_")
  if(d == 1) {
    te_fit_daily_4_h2<-temp
  }
  if (d > 1){
    te_fit_daily_4_h2<-merge(te_fit_daily_4_h2, temp, by=c("experiment", "year", "treatment", "plot", "plot_id", "genotype", "plantbarcode"), all=T)  
  }
}

## Add the final column and adjust
te_fit_daily_4_h2$Obs<-c(1:nrow(te_fit_daily_4_h2))
te_fit_daily_4_h2<-te_fit_daily_4_h2[,c(25,1:24)]

## Calculate H2
te_fit_daily_h2<-get_h2(te_fit_daily_4_h2)
te_fit_daily_h2_in_treatment<-get_h2_in_treatment(te_fit_daily_4_h2)

pdf('te_fit_daily_h2.pdf')

plot(te_fit_daily_h2~c(17:33), type="l", ylim=c(0,1), col="red", ylab="H2 of cumulative biomass~water fit", xlab="Days after sewing")
lines(te_fit_daily_h2_in_treatment[1,1:17]~c(17:33), ylim=c(0,1), col="orange")
lines(te_fit_daily_h2_in_treatment[1,18:34]~c(17:33), ylim=c(0,1), col="navy")
#legend("topright",legend=c("all", "dry", "wet"), col=c("red","orange","navy"), lty=c(1,1,1), bty=c('n'), cex=0.5)

dev.off()


## Calculate H2 for WUE ratio daily

wue_ratio_daily_only<-te_model[, c("experiment", "year", "plot", "plot_id", "plantbarcode", "genotype", "treatment", "dap_i", 
                                      "water_lost_day", "sv_area_day", "wue_ratio_day", "te_fit_day", "te_residual_day")]
days<-sort(unique(wue_ratio_daily_only$dap_i))
wue_ratio_daily_4_h2<-c()
for (d in 1:length(days)){
  day<-days[d]
  temp<-wue_ratio_daily_only[wue_ratio_daily_only$dap_i == day,]
  temp<-temp[,c(1:7, 11)]
  colnames(temp)[8]<-paste(colnames(temp)[8], day, sep="_")
  if(d == 1) {
    wue_ratio_daily_4_h2<-temp
  }
  if (d > 1){
    wue_ratio_daily_4_h2<-merge(wue_ratio_daily_4_h2, temp, by=c("experiment", "year", "treatment", "plot", "plot_id", "genotype", "plantbarcode"), all=T)  
  }
}

## Add the final column and adjust
wue_ratio_daily_4_h2$Obs<-c(1:nrow(wue_ratio_daily_4_h2))
wue_ratio_daily_4_h2<-wue_ratio_daily_4_h2[,c(25,1:24)]

## Calculate H2
wue_ratio_daily_h2<-get_h2(wue_ratio_daily_4_h2)
wue_ratio_daily_h2_in_treatment<-get_h2_in_treatment(wue_ratio_daily_4_h2)

pdf('wue_ratio_daily_h2.pdf')

plot(wue_ratio_daily_h2~c(17:33), type="l", ylim=c(0,1), col="red", ylab="H2 of WUE ratio daily", xlab="Days after sewing")
lines(wue_ratio_daily_h2_in_treatment[1,1:17]~c(17:33), ylim=c(0,1), col="orange")
lines(wue_ratio_daily_h2_in_treatment[1,18:34]~c(17:33), ylim=c(0,1), col="navy")
#legend("topright",legend=c("all", "dry", "wet"), col=c("red","orange","navy"), lty=c(1,1,1), bty=c('n'), cex=0.5)

dev.off()

## SV AREA PER DAY

sv_area_daily_only <- te_model[, c("experiment", "year", "plot", "plot_id", "plantbarcode", "genotype", "treatment", "dap_i", 
                                      "water_lost_day", "sv_area_day", "wue_ratio_day", "te_fit_day", "te_residual_day")]
days<-sort(unique(sv_area_daily_only$dap_i))
sv_area_daily_4_h2<-c()
for (d in 1:length(days)){
  day<-days[d]
  temp<-sv_area_daily_only[sv_area_daily_only$dap_i == day,]
  temp<-temp[,c(1:7, 10)]
  colnames(temp)[8]<-paste(colnames(temp)[8], day, sep="_")
  if(d == 1) {
    sv_area_daily_4_h2<-temp
  }
  if (d > 1){
    sv_area_daily_4_h2<-merge(sv_area_daily_4_h2, temp, by = c("experiment", "year", "treatment", "plot", "plot_id", 
                                                               "genotype", "plantbarcode"), all = T)  
  }
}

## Add the final column and adjust
sv_area_daily_4_h2$Obs<-c(1:nrow(sv_area_daily_4_h2))
sv_area_daily_4_h2<-sv_area_daily_4_h2[,c(25,1:24)]

## Calculate H2
sv_area_daily_h2<-get_h2(sv_area_daily_4_h2)
sv_area_daily_h2_in_treatment<-get_h2_in_treatment(sv_area_daily_4_h2)

pdf('sv_area_h2.pdf')

plot(sv_area_daily_h2~c(17:33), type="l", ylim=c(0,1), col="red", ylab="H2 of sv_area (pixels)", xlab="Days after planting")
lines(sv_area_daily_h2_in_treatment[1,1:17]~c(17:33), ylim=c(0,1), col="orange")
lines(sv_area_daily_h2_in_treatment[1,18:34]~c(17:33), ylim=c(0,1), col="navy")
#legend("topright",legend=c("all", "dry", "wet"), col=c("red","orange","navy"), lty=c(1,1,1), bty=c('n'), cex=0.5)

dev.off()

## WATER LOST DAILY

water_lost_daily_only <- te_model[, c("experiment", "year", "plot", "plot_id", "plantbarcode", "genotype", "treatment", "dap_i", 
                                         "water_lost_day", "sv_area_day", "wue_ratio_day", "te_fit_day", "te_residual_day")]
days<-sort(unique(water_lost_daily_only$dap_i))
water_lost_daily_4_h2<-c()
for (d in 1:length(days)){
  day<-days[d]
  temp<-water_lost_daily_only[water_lost_daily_only$dap_i == day,]
  temp<-temp[,c(1:7,9)]
  colnames(temp)[8]<-paste(colnames(temp)[8], day, sep="_")
  if(d == 1) {
    water_lost_daily_4_h2<-temp
  }
  if (d > 1){
    water_lost_daily_4_h2<-merge(water_lost_daily_4_h2, temp, by=c("experiment", "year", "treatment", "plot", "plot_id", "genotype", "plantbarcode"), all=T)  
  }
}

## Add the final column and adjust
water_lost_daily_4_h2$Obs<-c(1:nrow(water_lost_daily_4_h2))
water_lost_daily_4_h2<-water_lost_daily_4_h2[,c(25,1:24)]

# Calculate H2
water_lost_daily_h2<-get_h2(water_lost_daily_4_h2)
water_lost_daily_h2_in_treatment<-get_h2_in_treatment(water_lost_daily_4_h2)

pdf('water_lost_daily_h2.pdf')

plot(water_lost_daily_h2~c(17:33), type="l", ylim=c(0,1), col="red", ylab="H2 of water_lost_daily (mL)", xlab="Days after planting")
lines(water_lost_daily_h2_in_treatment[1,1:17]~c(17:33), ylim=c(0,1), col="orange")
lines(water_lost_daily_h2_in_treatment[1,18:34]~c(17:33), ylim=c(0,1), col="navy")
#legend("topright",legend=c("all", "dry", "wet"), col=c("red","orange","navy"), lty=c(1,1,1), bty=c('n'), cex=0.5)

dev.off()

####################
## Make some GGPLOT compatible data.frames for plotting
####################

## Make a data frame for each subset
dap_i<-c(17:33)
day_comp.h2<-cbind(dap_i,rep('all',length(dap_i)), sv_area_daily_h2, water_lost_daily_h2,te_residual_daily_h2, 
                   te_fit_daily_h2, wue_ratio_daily_h2)
colnames(day_comp.h2)<-c("dap_i", "treatment", "sv_area_daily", "water_lost_daily", "te_residual_daily", 
                         "te_fit_daily", "wue_ratio_daily")
## Dry
dry.day_comp.h2 <- cbind(dap_i, rep('dry', length(dap_i)), sv_area_h2_in_treatment[1, 1:17], 
                         water_lost_daily_h2_in_treatment[1, 1:17], te_residual_daily_h2_in_treatment[1, 1:17], 
                         te_fit_daily_h2_in_treatment[1, 1:17], wue_ratio_daily_h2_in_treatment[1, 1:17])
colnames(dry.day_comp.h2)<-c("dap_i", "treatment", "sv_area_daily", "water_lost_daily", "te_residual_daily", 
                             "te_fit_daily", "wue_ratio_daily")
## Wet
wet.day_comp.h2 <- cbind(dap_i,rep('wet', length(dap_i)), sv_area_daily_h2_in_treatment[1,18:34], 
                         water_lost_daily_h2_in_treatment[1, 18:34], te_residual_daily_h2_in_treatment[1, 18:34], 
                         te_fit_daily_h2_in_treatment[1, 18:34], wue_ratio_daily_h2_in_treatment[1, 18:34])
colnames(wet.day_comp.h2)<-c("dap_i", "treatment", "sv_area_daily", "water_lost_daily", "te_residual_daily", 
                             "te_fit_daily", "wue_ratio_daily")

## Combine them
all.day_comp.h2 <- rbind(day_comp.h2, dry.day_comp.h2, wet.day_comp.h2)
rownames(all.day_comp.h2)<-c(1:nrow(all.day_comp.h2))

gg.comp.daily.h2 <- c()
for(i in 3:ncol(all.day_comp.h2)){
  temp <- as.data.frame(all.day_comp.h2[,c(i)])
  temp$trait <- rep(colnames(all.day_comp.h2)[i], nrow(all.day_comp.h2))
  temp <- cbind(temp, all.day_comp.h2[,1:2])
  colnames(temp) <- c("value", "trait", "dap_i", 'treatment')
  gg.comp.daily.h2 <- rbind(gg.comp.daily.h2, temp)
}

gg.comp.daily.h2$value<-as.numeric(as.character(gg.comp.daily.h2$value))
gg.comp.daily.h2$dap_i<-as.numeric(as.character(gg.comp.daily.h2$dap_i))


pdf("ril_heritability_rate.pdf")
p<-ggplot(gg.comp.daily.h2, aes(x=dap_i, y=value, color=treatment)) + geom_line(size=2)
g<-p+ scale_color_manual(values=c("red", "orange", "navy")) + ylab("Broad sense H2") + xlab("Days after planting") + facet_wrap(~trait)
d<-g +theme_bw() + ylim(0,1) + theme(legend.position="none", axis.text=element_text(size=20), axis.title=element_text(size=24,face="bold"), plot.title = element_text(size = 24, face = "bold"))
d
dev.off()


setwd(wue_results.h2.dir)
save.image('ril_heritability.Rdata')


## Lets get average heritability for each trait across all treatments
gg.comp.h2$cat<-paste(gg.comp.h2$trait, gg.comp.h2$treatment, sep=".")


cats<-unique(gg.comp.h2$cat)
h2_aves<-c()
for (i in cats){
  temp<-gg.comp.h2[gg.comp.h2$cat == i,]
  val<-mean(temp$value)
  temp2<-c(i, val)
  h2_aves<-rbind(h2_aves, temp2)
}


## Now do daily or rate of change
gg.comp.daily.h2$cat<-paste(gg.comp.daily.h2$trait, gg.comp.daily.h2$treatment, sep=".")

cats<-unique(gg.comp.daily.h2$cat)
h2_rate_aves<-c()
for (i in cats){
  temp<-gg.comp.daily.h2[gg.comp.daily.h2$cat == i,]
  val<-mean(temp$value)
  temp2<-c(i, val)
  h2_rate_aves<-rbind(h2_rate_aves, temp2)
}

colnames(h2_aves)<-c("trait", "h2_ave")
colnames(h2_rate_aves)<-c("trait", "h2_rate_ave")

H2_aves<-cbind(h2_aves, h2_rate_aves[,2])
colnames(H2_aves)[3]<-c("h2_rate_ave")
rownames(H2_aves)<-c(1:nrow(H2_aves))
H2_aves<-as.data.frame(H2_aves)
H2_aves$h2_ave<-as.numeric(as.character(H2_aves$h2_ave))
H2_aves$h2_rate_ave<-as.numeric(as.character(H2_aves$h2_rate_ave))


H2_aves$diff<-H2_aves[,2] - H2_aves[,3]
mean(H2_aves$diff)

## Here you can see there is a 5% difference on average across all traits



###########################################################################################
# Lets do the same thing with % variance explained for all factors
###########################################################################################

sv_area_total_prop.var<-get_total_var(sv_area_4_h2)
water_lost_total_prop.var<-get_total_var(water_lost_total_4_h2)
te_residual_total_prop.var<-get_total_var(te_residual_total_4_h2)
te_fit_total_prop.var<-get_total_var(te_fit_total_4_h2)
wue_ratio_total_prop.var<-get_total_var(wue_ratio_total_4_h2)


#########################################
## Make a data frame for each subset

## sv_area_total
dap_i<-c(17:33)
gg.sv_area_total.pv<-c()
for(r in 1:nrow(sv_area_total_prop.var)){
  name<-rownames(sv_area_total_prop.var)[r]
  temp<-sv_area_total_prop.var[r,]
  temp2<-cbind(rep('sv_area_total',length(temp)), dap_i, rep(name, length(temp)), temp)
  gg.sv_area_total.pv<-rbind(gg.sv_area_total.pv, temp2)
}
colnames(gg.sv_area_total.pv)<-c("trait", "dap_i", "type", "value")

## water_lost_total
dap_i<-c(17:33)
gg.water_lost_total.pv<-c()
for(r in 1:nrow(water_lost_total_prop.var)){
  name<-rownames(water_lost_total_prop.var)[r]
  temp<-water_lost_total_prop.var[r,]
  temp2<-cbind(rep('water_lost_total',length(temp)), dap_i, rep(name, length(temp)), temp)
  gg.water_lost_total.pv<-rbind(gg.water_lost_total.pv, temp2)
}
colnames(gg.water_lost_total.pv)<-c("trait", "dap_i", "type", "value")

## te_residual_total
dap_i<-c(17:33)
gg.te_residual_total.pv<-c()
for(r in 1:nrow(te_residual_total_prop.var)){
  name<-rownames(te_residual_total_prop.var)[r]
  temp<-te_residual_total_prop.var[r,]
  temp2<-cbind(rep('te_residual_total',length(temp)), dap_i, rep(name, length(temp)), temp)
  gg.te_residual_total.pv<-rbind(gg.te_residual_total.pv, temp2)
}
colnames(gg.te_residual_total.pv)<-c("trait", "dap_i", "type", "value")

## te_fit_total
dap_i<-c(17:33)
gg.te_fit_total.pv<-c()
for(r in 1:nrow(te_fit_total_prop.var)){
  name<-rownames(te_fit_total_prop.var)[r]
  temp<-te_fit_total_prop.var[r,]
  temp2<-cbind(rep('te_fit_total',length(temp)), dap_i, rep(name, length(temp)), temp)
  gg.te_fit_total.pv<-rbind(gg.te_fit_total.pv, temp2)
}
colnames(gg.te_fit_total.pv)<-c("trait", "dap_i", "type", "value")

## wue_ratio_total
dap_i<-c(17:33)
gg.wue_ratio_total.pv<-c()
for(r in 1:nrow(wue_ratio_total_prop.var)){
  name<-rownames(wue_ratio_total_prop.var)[r]
  temp<-wue_ratio_total_prop.var[r,]
  temp2<-cbind(rep('wue_ratio_total',length(temp)), dap_i, rep(name, length(temp)), temp)
  gg.wue_ratio_total.pv<-rbind(gg.wue_ratio_total.pv, temp2)
}
colnames(gg.wue_ratio_total.pv)<-c("trait", "dap_i", "type", "value")


gg.prop.var<-rbind(gg.sv_area_total.pv, gg.water_lost_total.pv, gg.te_residual_total.pv, gg.te_fit_total.pv, gg.wue_ratio_total.pv)
rownames(gg.prop.var)<-c(1:nrow(gg.prop.var))
gg.prop.var<-as.data.frame(gg.prop.var)

gg.prop.var$value<-as.numeric(as.character(gg.prop.var$value))
gg.prop.var$dap_i<-as.numeric(as.character(gg.prop.var$dap_i))
#levels(gg.prop.var$type)<-c("Genotype", "Treatment", "G X Treatment", "Error")

pdf("ril_prop.var_total.pdf")
p<-ggplot(gg.prop.var, aes(x=dap_i, y=value, color=type)) + geom_line(size=2)
g<-p+ scale_color_manual(values=c("red", "blue", "purple", "black")) + ylab("% of variance") + xlab("Days after planting") + facet_wrap(~trait)
g<-p + ylab("% of variance") + xlab("Days after planting") + facet_wrap(~trait)
#d<-g +theme_bw() + ylim(0,1) + theme(legend.position="none", axis.text=element_text(size=20), axis.title=element_text(size=24,face="bold"), plot.title = element_text(size = 24, face = "bold"))
d<-g +theme_bw() + ylim(0,1) + theme(axis.text=element_text(size=20), axis.title=element_text(size=24,face="bold"), plot.title = element_text(size = 24, face = "bold"))
d
dev.off()

save.image('ril_heritability.Rdata')


###########################################################################################
# RATE/DAILY
###########################################################################################

sv_area_daily_prop.var<-get_total_var(sv_area_daily_4_h2)
water_lost_daily_prop.var<-get_total_var(water_lost_daily_4_h2)
te_residual_daily_prop.var<-get_total_var(te_residual_daily_4_h2)
te_fit_daily_prop.var<-get_total_var(te_fit_daily_4_h2)
wue_ratio_daily_prop.var<-get_total_var(wue_ratio_daily_4_h2)

#########################################
## Make a data frame for each subset

## sv_area_daily
dap_i<-c(17:33)
gg.sv_area_daily.pv<-c()
for(r in 1:nrow(sv_area_daily_prop.var)){
  name<-rownames(sv_area_daily_prop.var)[r]
  temp<-sv_area_daily_prop.var[r,]
  temp2<-cbind(rep('sv_area_daily',length(temp)), dap_i, rep(name, length(temp)), temp)
  gg.sv_area_daily.pv<-rbind(gg.sv_area_daily.pv, temp2)
}
colnames(gg.sv_area_daily.pv)<-c("trait", "dap_i", "type", "value")

## water_lost_daily
dap_i<-c(17:33)
gg.water_lost_daily.pv<-c()
for(r in 1:nrow(water_lost_daily_prop.var)){
  name<-rownames(water_lost_daily_prop.var)[r]
  temp<-water_lost_daily_prop.var[r,]
  temp2<-cbind(rep('water_lost_daily',length(temp)), dap_i, rep(name, length(temp)), temp)
  gg.water_lost_daily.pv<-rbind(gg.water_lost_daily.pv, temp2)
}
colnames(gg.water_lost_daily.pv)<-c("trait", "dap_i", "type", "value")

# te_residual_daily
dap_i<-c(17:33)
gg.te_residual_daily.pv<-c()
for(r in 1:nrow(te_residual_daily_prop.var)){
  name<-rownames(te_residual_daily_prop.var)[r]
  temp<-te_residual_daily_prop.var[r,]
  temp2<-cbind(rep('te_residual_daily',length(temp)), dap_i, rep(name, length(temp)), temp)
  gg.te_residual_daily.pv<-rbind(gg.te_residual_daily.pv, temp2)
}
colnames(gg.te_residual_daily.pv)<-c("trait", "dap_i", "type", "value")

## te_fit_daily
dap_i<-c(17:33)
gg.te_fit_daily.pv<-c()
for(r in 1:nrow(te_fit_daily_prop.var)){
  name<-rownames(te_fit_daily_prop.var)[r]
  temp<-te_fit_daily_prop.var[r,]
  temp2<-cbind(rep('te_fit_daily',length(temp)), dap_i, rep(name, length(temp)), temp)
  gg.te_fit_daily.pv<-rbind(gg.te_fit_daily.pv, temp2)
}
colnames(gg.te_fit_daily.pv)<-c("trait", "dap_i", "type", "value")

## wue_ratio_daily
dap_i<-c(17:33)
gg.wue_ratio_daily.pv<-c()
for(r in 1:nrow(wue_ratio_daily_prop.var)){
  name<-rownames(wue_ratio_daily_prop.var)[r]
  temp<-wue_ratio_daily_prop.var[r,]
  temp2<-cbind(rep('wue_ratio_daily',length(temp)), dap_i, rep(name, length(temp)), temp)
  gg.wue_ratio_daily.pv<-rbind(gg.wue_ratio_daily.pv, temp2)
}
colnames(gg.wue_ratio_daily.pv)<-c("trait", "dap_i", "type", "value")


gg.prop.var.daily<-rbind(gg.sv_area_daily.pv, gg.water_lost_daily.pv, gg.te_residual_daily.pv, gg.te_fit_daily.pv, gg.wue_ratio_daily.pv)
rownames(gg.prop.var.daily)<-c(1:nrow(gg.prop.var.daily))
gg.prop.var.daily<-as.data.frame(gg.prop.var.daily)

gg.prop.var.daily$value<-as.numeric(as.character(gg.prop.var.daily$value))
gg.prop.var.daily$dap_i<-as.numeric(as.character(gg.prop.var.daily$dap_i))
#levels(gg.prop.var$type)<-c("Genotype", "Treatment", "G X Treatment", "Error")

pdf("ril_prop.var_daily.pdf")
p<-ggplot(gg.prop.var.daily, aes(x=dap_i, y=value, color=type)) + geom_line(size=2)
g<-p+ scale_color_manual(values=c("red", "blue", "purple", "black")) + ylab("% of variance") + xlab("Days after planting") + facet_wrap(~trait)
g<-p + ylab("% of variance") + xlab("Days after planting") + facet_wrap(~trait)
#d<-g +theme_bw() + ylim(0,1) + theme(legend.position="none", axis.text=element_text(size=20), axis.title=element_text(size=24,face="bold"), plot.title = element_text(size = 24, face = "bold"))
d<-g +theme_bw() + ylim(0,1) + theme(axis.text=element_text(size=20), axis.title=element_text(size=24,face="bold"), plot.title = element_text(size = 24, face = "bold"))
d
dev.off()


setwd(wue_results.h2.dir)
save.image('ril_heritability.Rdata')














## Revisit after isotopes

##################
# Lets make a plot for the Maize poster
##################

load("~/Dropbox/ellsworth_feldman_water_use/results/isotope/ril_stable_isotope4.Rdata")

isotope_qtl<-isotope_qtl[,1:9]
hv<-get_h2(isotope_qtl)
hv_trt<-get_h2_in_treatment(isotope_qtl)
pdf("ril_heritability_3.pdf")

p<-ggplot(gg.comp.h2[gg.comp.h2$treatment == 'all',], aes(x=dap_i, y=value, color=trait)) + geom_line(size=1.3)
g<-p+ scale_color_manual(values=c("green", "red", "blue","black","orange")) + ylab("Heritability (H2)") + xlab("Days after planting") 
d<-g +theme_bw() + ylim(0,1) + theme(legend.position="none", axis.text=element_text(size=20), axis.title=element_text(size=24,face="bold"), plot.title = element_text(size = 24, face = "bold"))
d + geom_point(x = 32, y = hv, shape=8, color='purple', size=6)

print(d)
dev.off()


prop.var.30<-gg.prop.var[gg.prop.var$dap_i == 30,]
value<-get_total_var(isotope_qtl)
trait<-rep("d13C", 4)
dap_i<-rep('30', 4)
type<-c('Genotype', 'Treatment', 'G x Treatment', 'Error')
var.isotope<-cbind(trait, dap_i, type, value)
colnames(var.isotope)[4]<-c('value')

prop.var.30<-rbind(prop.var.30, var.isotope)
prop.var.30$trait<-as.character(prop.var.30$trait)
prop.var.30[prop.var.30$trait == 'sv_area_total', 'trait'] <- c("Biomass")
prop.var.30[prop.var.30$trait == 'water_lost_total', 'trait'] <- c("Water lost")
prop.var.30[prop.var.30$trait == 'te_fit_total', 'trait'] <- c("TE fit")
prop.var.30[prop.var.30$trait == 'te_residual_total', 'trait'] <- c("TE residual")
prop.var.30[prop.var.30$trait == 'wue_ratio_total', 'trait'] <- c("WUE ratio")

prop.var.30$trait<-factor(prop.var.30$trait, c('Biomass', "Water lost", "TE fit", "TE residual", "WUE ratio", "d13C"))
prop.var.30$type<-factor(prop.var.30$type, c("Genotype", "Treatment", "G x Treatment", "Error"))

prop.var.30$value<-as.numeric(as.character(prop.var.30$value))


p<-ggplot(data=prop.var.30, aes(x=factor(1), y=value, fill=factor(type))) + geom_bar(width=1, stat = 'identity') + theme_bw() + facet_grid(facets=.~trait)
q<-p + theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.text=element_text(size=20), axis.title=element_text(size=24,face="bold"), plot.title = element_text(size = 24, face = "bold"))
x<-q + ylab("% of total variance") + guides(fill = guide_legend(title = "Variance"))

pdf("prop.variance.day30.pdf", width=8, height=8)
print(x)
dev.off()
x



#### Load

load("~/Dropbox/ellsworth_feldman_water_use/results/heritability_2/ril_heritability_2.Rdata")

load("~/Dropbox/ellsworth_feldman_water_use/results/isotope/ril_stable_isotope2.Rdata")

isotope_qtl<-isotope_qtl[,1:9]
hv<-get_h2(isotope_qtl)

p<-ggplot(gg.comp.h2[gg.comp.h2$treatment == 'all',], aes(x=dap_i, y=value, color=trait)) + geom_line(size=1.3)
g<-p+ scale_color_manual(values=c("green", "red", "blue","black","orange")) + ylab("Heritability (H2)") + xlab("Days after planting") 
d<-g +theme_bw() + ylim(0,1) + theme(legend.position="none", axis.text=element_text(size=20), axis.title=element_text(size=24,face="bold"), plot.title = element_text(size = 24, face = "bold"))
d + geom_point(x = 32, y = hv, shape=8, color='purple', size=6)

print(d)
dev.off()


prop.var.30<-gg.prop.var[gg.prop.var$dap_i == 30,]
value<-get_total_var(isotope_qtl)





### For DDPSC Science Retreat
load("~/Dropbox/ellsworth_feldman_water_use/results/isotope/ril_stable_isotope4.Rdata")

isotope_qtl<-isotope_qtl[,1:9]
hv<-get_h2(isotope_qtl)
hv_trt<-get_h2_in_treatment(isotope_qtl)

h2.iso<-c(hv, hv_trt[1,])
cond<-c("all", "wet", "dry")
h2.iso<-rbind(h2.iso,cond)
h2.iso<-t(h2.iso)
colnames(h2.iso)<-c("h2", "treatment")
h2.iso<-as.data.frame(h2.iso)
h2.iso$h2<-as.numeric(as.character(h2.iso$h2))
h2.iso$trait<-rep("d13C", nrow(h2.iso))

pdf("d13C_h2.pdf", height=4, width=3)
p<-ggplot(h2.iso, aes(y=h2, x=trait, color=treatment)) + geom_point(size=3) + theme_bw() + theme(legend.position="none", axis.text=element_text(size=20), axis.title=element_text(size=24,face="bold"), plot.title = element_text(size = 24, face = "bold"))
p2<- p+ scale_color_manual(values=c("red", "orange", "navy")) + ylab("Broad sense H2") + xlab("") + facet_wrap(~trait)
print(p2)
dev.off()


prop.var_iso<-get_total_var(isotope_qtl)
type<-row.names(prop.var_iso)
trait<-rep("d13C", nrow(prop.var_iso))
prop.var_iso<-cbind(prop.var_iso, type, trait)
#prop.var_iso<-t(prop.var_iso)
prop.var_iso<-as.data.frame(prop.var_iso)
prop.var_iso$d13C_BP14<-as.numeric(as.character(prop.var_iso$d13C_BP14))

pdf("d13C_prop.var.pdf", height=4, width=3)
p<-ggplot(prop.var_iso, aes(y=d13C_BP14, x=trait, color=type)) + geom_point(size=3) + theme_bw() + theme(legend.position="none", axis.text=element_text(size=20), axis.title=element_text(size=24,face="bold"), plot.title = element_text(size = 24, face = "bold"))
p2<- p + ylab("% of variance") + xlab("") + facet_wrap(~trait)
print(p2)
dev.off()
