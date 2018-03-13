## Load packages one might use
library(ggplot2, warn.conflicts = FALSE)
library(lubridate, warn.conflicts = FALSE)
library(MASS, warn.conflicts = FALSE)
library(car, warn.conflicts = FALSE)
library(nlme, warn.conflicts = FALSE)
library(mvtnorm, warn.conflicts = FALSE)
library(grid, warn.conflicts = FALSE)

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

## Make a directory for water use results
wue_results.water_use.dir<-paste(wue_results.dir, 'water_use/', sep="")

## If this directory doesn't exist the download is not complete
if (file.exists(wue_results.water_use.dir)){
  print("Directory exists!")
} else {
  dir.create(file.path(wue_results.water_use.dir))
}

##### END ##### 


setwd(wue_data.dir)

## This is a standard output from LemnaTec instrument
w<-read.csv("SnapshotInfo_RIL_allday.csv")

w<-w[w$water.amount != -1,]

planting_date = as.POSIXct("2014-1-13")
w$date<-as.POSIXct(w$timestamp, origin = "1970-01-01")
w$dap <- as.numeric(w$date - planting_date)
w$dap_i<-as.integer(w$dap)

## Read a key that associates phenotyper ID (aka 'plant_id') to RIL id
key<-read.csv('phe2ril2.csv', header=F)
colnames(key)<-c('phe', 'ril')

## Make the RIL name same as in QTL map file
key$ril<-sprintf("%03d", key$ril)
key$ril<-paste('RIL_', key$ril, sep="")

## Add a new (empty column called genotype)
w$genotype<-NA
## Add genotypes 
for(i in 1:nrow(key)) {
  phe_name<-key[i,1]
  ril_name<-key[i,2]
  w$genotype[grep(paste('Dr', phe_name, 'A', sep=""), w$plant.barcode)]<-ril_name
}

## Add parental genotypes
w$genotype[grep('Dp1', w$plant.barcode)]<-c("A10")
w$genotype[grep('Dp2', w$plant.barcode)]<-c("B100")

## Add treatment
w$treatment<-c('none')
w$treatment[grep("AA", w$plant.barcode)]<-c('wet')
w$treatment[grep("AB", w$plant.barcode)]<-c('dry')

water_added_dap<-aggregate(w$water.amount, by=list(w$plant.barcode, w$car.tag, w$dap_i, w$genotype, w$treatment), sum)
colnames(water_added_dap)<-c('plantbarcode', 'cartag', 'dap_i', 'genotype', 'treatment','water_added')
weight_after_dap<-aggregate(w$weight.after, by=list(w$plant.barcode, w$car.tag, w$dap_i, w$genotype, w$treatment), min)
colnames(weight_after_dap)<-c('plantbarcode', 'cartag', 'dap_i', 'genotype', 'treatment','weight_total')

ril_water<-merge(water_added_dap, weight_after_dap, by=c('plantbarcode', 'cartag', 'dap_i', 'genotype', 'treatment'))

## Might think about how increases in set point for water increases throughout experiment 
p<-ggplot(data=weight_after_dap, aes(x=dap_i, y=weight_total, colour=treatment)) + geom_point(size=0.1)
p
## Looks like there are changes on days 19, and 25

dap<-sort(unique(ril_water$dap_i))
dap<-dap[dap <= 33 & dap > 14]

water_lost<-c()
for (d in 1:length(dap)){
  day<-dap[d]
  if (d == 1) {
    temp<-ril_water[ril_water$dap_i == day,]
    water_lost<-temp[,1:6]
    colnames(water_lost)[6]<-c("water_lost")
  }
  if (d > 1){
    temp<-ril_water[ril_water$dap_i == day,]
    if (day == 19 | day == 25){
      plant_ids<-temp$plantbarcode
      yesterday<-ril_water[ril_water$plantbarcode %in% plant_ids,]
      yesterday<-yesterday[yesterday$dap_i == (day-1),c(1:2,4:5,7)]
      colnames(yesterday)[5]<-c('weight_yesterday')
      temp3<-merge(temp, yesterday, by=c('plantbarcode', 'cartag', 'genotype', 'treatment'))
      temp3$weight_offset<-temp3$weight_total - temp3$weight_yesterday
      temp3$water_lost<-temp3$water_added - temp3$weight_offset
      temp3<-temp3[,c(1,2,5,3,4,10)]
      water_lost<-rbind(water_lost, temp3)
      next
    }
    temp2<-temp[temp$water_added > 0,]
    colnames(temp2)[6]<-c("water_lost")
    water_lost<-rbind(water_lost, temp2[,c(1:6)])
    dry_down<-temp[temp$water_added == 0,]
    if (nrow(dry_down) > 0){
      plant_ids<-dry_down$plantbarcode
      yesterday<-ril_water[ril_water$plantbarcode %in% plant_ids,]
      yesterday<-yesterday[yesterday$dap_i == (day-1),c(1:2,4:5,7)]
      colnames(yesterday)[5]<-c('weight_yesterday')
      temp3<-merge(dry_down, yesterday, by=c('plantbarcode', 'cartag', 'genotype', 'treatment'))
      temp3$water_lost<-temp3$weight_yesterday - temp3$weight_total
      temp3<-temp3[,c(1,2,5,3,4,9)]
      water_lost<-rbind(water_lost, temp3)
    }
  }
}
p<-ggplot(data=water_lost, aes(x=dap_i, y=water_lost, colour=treatment)) + geom_point(cex=0.3) 
p
p<-ggplot(data=water_lost, aes(x=dap_i, y=water_lost, colour=treatment)) + geom_point(cex=0.3) + ylim(0,50)
p

## Remove all values below zero and above 200 (day 17, 32, 33)
water_lost<-water_lost[!(water_lost$water_lost > 200 & water_lost$dap_i %in% c(17,32,33)),]
water_lost<-water_lost[water_lost$water_lost > 0,]
p<-ggplot(data=water_lost, aes(x=dap_i, y=water_lost, colour=treatment)) + geom_point(cex=0.3) 
p

p<-ggplot(data=water_lost[water_lost$treatment == 'dry',], aes(x=dap_i, y=water_lost, colour=treatment)) + geom_point(cex=0.3) 
p

## Remove outlier
water_lost<-water_lost[!(water_lost$treatment == 'dry' & water_lost$dap_i == 32 & water_lost$water_lost > 145),]
p<-ggplot(data=water_lost[water_lost$treatment == 'wet',], aes(x=dap_i, y=water_lost, colour=treatment)) + geom_point(cex=0.3) 
p

p<-ggplot(data=water_lost, aes(x=dap_i, y=water_lost, colour=treatment)) + geom_smooth() + ylim(0,60)
p



setwd(wue_results.water_use.dir)
save.image("ril_water_use.Rdata")

empty_pots<-c("Dr27AA001310", "Dr41AA001393", "Dr53AB001468", "Dr65AB001541", "Dr72AA001581")

## Lets break out the empty pots
## Rename them in the data.frame water total
water_lost$genotype<-as.character(water_lost$genotype)
water_lost[water_lost$plantbarcode %in% empty_pots, 'genotype']<-c("EMPTY")

## Make a subset of parental lines and empty pots call it parentals
parentals<-water_lost[water_lost$genotype == "A10" | water_lost$genotype == "B100" | water_lost$genotype == "EMPTY",]
## Plot the parental lines water use vs. empty
p<-ggplot(parentals, aes(dap_i, water_lost, color=genotype)) + geom_smooth() + facet_wrap(~treatment)
p

water_lost<-water_lost[water_lost$dap_i < 34 & water_lost$dap_i >14,]
water_days<-c(15:33)
water_total<-c()
for(d in 1:length(water_days)) {
  day<-water_days[d]
  temp<-water_lost[water_lost$dap_i <= day,]
  temp<-aggregate(temp[,6], by=list(temp$plantbarcode, temp$cartag, temp$genotype, temp$treatment), sum)
  temp$dap_i<-rep(day, nrow(temp))
  colnames(temp)<-c("plantbarcode", "cartag","genotype","treatment","water_lost_total","dap_i")
  water_total<-rbind(water_total, temp)
}

## Lets merge this total water lost with the daily water lost
water.m<-merge(water_lost, water_total, by=c("plantbarcode", "cartag", "genotype", "treatment", "dap_i"))

parentals<-water.m[water.m$genotype == "A10" | water.m$genotype == "B100" | water.m$genotype == "EMPTY",]
## Plot the parental lines water use vs. empty
p<-ggplot(parentals, aes(dap_i, water_lost_total, color=genotype)) + geom_smooth() + facet_wrap(~treatment)
p
p<-ggplot(parentals, aes(dap_i, water_lost, color=genotype)) + geom_smooth() + facet_wrap(~treatment)
p


## Need to subtract the water from empty from water added
days<-sort(unique(water.m$dap_i))
ril_water<-c()
for(d in 1:length(days)){
  ## Calculate the amount of water that evaporated from empty pots in both conditions
  temp<-water.m[water.m$dap_i == days[d],]
  evaporation_day.w<-min(temp[temp$genotype == 'EMPTY' & temp$treatment == 'wet', 'water_lost'], na.rm=T)
  evaporation_day.d<-min(temp[temp$genotype == 'EMPTY' & temp$treatment == 'dry', 'water_lost'], na.rm=T)
  evaporation_total.w<-min(temp[temp$genotype == 'EMPTY' & temp$treatment == 'wet', 'water_lost_total'], na.rm=T)
  evaporation_total.d<-min(temp[temp$genotype == 'EMPTY' & temp$treatment == 'dry', 'water_lost_total'], na.rm=T)
  
  ## If there is no value set equal to zero
  ## See warning message
  if(!is.finite(evaporation_day.w)) {evaporation_day.w <-0} 
  if(!is.finite(evaporation_day.d)) {evaporation_day.d <-0} 
  if(!is.finite(evaporation_total.w)) {evaporation_total.w <-0} 
  if(!is.finite(evaporation_total.d)) {evaporation_total.d <-0} 
  
  ## Get a value for transpiration (subptract the amount of water lost due to evaporation from empty pots)
  transpiration_day.w<-temp[temp$treatment == 'wet', 'water_lost'] - evaporation_day.w
  transpiration_day.d<-temp[temp$treatment == 'dry', 'water_lost'] - evaporation_day.d
  transpiration_total.w<-temp[temp$treatment == 'wet', 'water_lost_total'] - evaporation_total.w
  transpiration_total.d<-temp[temp$treatment == 'dry', 'water_lost_total'] - evaporation_total.d
  temp$transpiration_day<-c(NA)
  temp$transpiration_total<-c(NA)
  temp[temp$treatment == 'wet','transpiration_day']<-transpiration_day.w
  temp[temp$treatment == 'dry','transpiration_day']<-transpiration_day.d
  temp[temp$treatment == 'wet','transpiration_total']<-transpiration_total.w
  temp[temp$treatment == 'dry','transpiration_total']<-transpiration_total.d
  temp[temp$treatment == 'wet','evaporation_day']<-evaporation_day.w
  temp[temp$treatment == 'dry','evaporation_day']<-evaporation_day.d
  temp[temp$treatment == 'wet','evaporation_total']<-evaporation_total.w
  temp[temp$treatment == 'dry','evaporation_total']<-evaporation_total.d
  ## Now calculate the ratio of evaporation over total evapotranspiration
  temp$ev_to_evt_day<-temp$evaporation_day/temp$water_lost
  temp$ev_to_evt_total<-temp$evaporation_total/temp$water_lost_total
  ril_water<-rbind(ril_water, temp)
}


p<-ggplot(data=ril_water, aes(x=dap_i, y=water_lost_total, colour=treatment)) + geom_smooth()
p
p<-ggplot(data=ril_water, aes(x=dap_i, y=transpiration_total, colour=treatment)) + geom_smooth()
p
p<-ggplot(data=ril_water, aes(x=dap_i, y=evaporation_total, colour=treatment)) + geom_smooth()
p
p<-ggplot(data=ril_water, aes(x=dap_i, y=ev_to_evt_total, colour=treatment)) + geom_smooth()
p

## Save the data
save.image("ril_water_use.Rdata")


##############################################
## Now get the loess fit for various water types
##############################################

## WATER LOST TOTAL
## Now go for a loess fit
genos<-as.character(sort(unique(ril_water$genotype)))
treatments<-unique(ril_water$treatment)
dap_i<-unique(ril_water$dap_i)
dap_i<-sort(as.numeric(as.character(dap_i)))
times = seq(from = 15, to = 33, by=0.1)

report.loess.values(ril_water, 'water_lost_total', genos, treatments, dap_i, min(dap_i), max(dap_i), "ril_loess_estimates_water_lost_total.pdf")

## Give the output an appropriate name
ril_loess_water_lost_total<-ril_loess_model_fit
## Rename the growth report
water_lost_total_report<-growth_rate_report


#############################################
## Format for QTL pipeline
water.l<-ril_loess_water_lost_total

days<-sort(unique(water.l$dap_i))

ril_water_lost_total_qtl<-c()
for(d in 1:length(days)) {
  day<-days[d]
  temp.data<-water.l[water.l$dap_i == day, ]
  colnames(temp.data)[4]<-paste('water_lost_total', day, sep="_")
  temp.data<-temp.data[,c(1,2,4)]
  if (d == 1) {
    ril_water_lost_total_qtl<-temp.data
  }
  
  if (d > 1) {
    ril_water_lost_total_qtl<-merge(ril_water_lost_total_qtl, temp.data, by=c('genotype', 'treatment'), all=T)
  }
}

## Add misc columns
ril_water_lost_total_qtl$Obs<-c(1:nrow(ril_water_lost_total_qtl))
ril_water_lost_total_qtl$experiment<-rep("BP14", nrow(ril_water_lost_total_qtl))
ril_water_lost_total_qtl$year<-rep("2014", nrow(ril_water_lost_total_qtl))
ril_water_lost_total_qtl$plot<-rep("bellweater", nrow(ril_water_lost_total_qtl))
ril_water_lost_total_qtl$plot_id<-rep("bellweater", nrow(ril_water_lost_total_qtl))
ril_water_lost_total_qtl$sampling<-rep("bellweater", nrow(ril_water_lost_total_qtl))

ril_water_lost_total_qtl<-ril_water_lost_total_qtl[,c(22,23,24,2,25,26,1,27,5:21)]
colnames(ril_water_lost_total_qtl)[7]<-c("id")
write.csv(ril_water_lost_total_qtl, file="ril_loess_water_lost_total_qtl.csv", quote=F, row.names=F)


##############################################
## TRANSPIRATION TOTAL
## Now go for a loess fit
genos<-as.character(sort(unique(ril_water$genotype)))
treatments<-unique(ril_water$treatment)
dap_i<-unique(ril_water$dap_i)
dap_i<-sort(as.numeric(as.character(dap_i)))
times = seq(from = 15, to = 33, by=0.1)

report.loess.values(ril_water, 'transpiration_total', genos, treatments, dap_i, min(dap_i), max(dap_i), "ril_loess_estimates_transpiration_total.pdf")

# Give output and appropriate name
ril_loess_transpiration_total<-ril_loess_model_fit
# Rename the growth report
transpiration_total_report<-growth_rate_report

#############################################
## Format for QTL pipeline
transpiration.l<-ril_loess_transpiration_total

days<-sort(unique(transpiration.l$dap_i))

ril_transpiration_total_qtl<-c()
for(d in 1:length(days)) {
  day<-days[d]
  temp.data<-transpiration.l[transpiration.l$dap_i == day, ]
  colnames(temp.data)[4]<-paste('transpiration_total', day, sep="_")
  temp.data<-temp.data[,c(1,2,4)]
  if (d == 1) {
    ril_transpiration_total_qtl<-temp.data
  }
  
  if (d > 1) {
    ril_transpiration_total_qtl<-merge(ril_transpiration_total_qtl, temp.data, by=c('genotype', 'treatment'), all=T)
  }
}

# Add misc columns
ril_transpiration_total_qtl$Obs<-c(1:nrow(ril_transpiration_total_qtl))
ril_transpiration_total_qtl$experiment<-rep("BP14", nrow(ril_transpiration_total_qtl))
ril_transpiration_total_qtl$year<-rep("2014", nrow(ril_transpiration_total_qtl))
ril_transpiration_total_qtl$plot<-rep("bellweater", nrow(ril_transpiration_total_qtl))
ril_transpiration_total_qtl$plot_id<-rep("bellweater", nrow(ril_transpiration_total_qtl))
ril_transpiration_total_qtl$sampling<-rep("bellweater", nrow(ril_transpiration_total_qtl))

ril_transpiration_total_qtl<-ril_transpiration_total_qtl[,c(22,23,24,2,25,26,1,27,3:21)]
colnames(ril_transpiration_total_qtl)[7]<-c("id")
write.csv(ril_transpiration_total_qtl, file="ril_loess_transpiration_total_qtl.csv", quote=F, row.names=F)

save.image("ril_water_use.Rdata")


##############################################
## EVAPORATION TOTAL
## Now go for a loess fit
genos<-as.character(c("EMPTY"))
treatments<-unique(ril_water$treatment)
dap_i<-unique(ril_water$dap_i)
dap_i<-sort(as.numeric(as.character(dap_i)))
times = seq(from = 15, to = 33, by=0.1)
report.loess.values(ril_water, 'evaporation_total', genos, treatments, dap_i, min(dap_i), max(dap_i), "ril_loess_estimates_evaporation_total.pdf")

ril_loess_evaporation_total<-ril_loess_model_fit
# Rename the growth report
evaporation_total_report<-growth_rate_report

save.image("ril_water_use.Rdata")

#############################################################################################################
## Now lets summarize by day
#############################################################################################################

## WATER LOST DAY
## Now go for a loess fit
genos<-as.character(sort(unique(ril_water$genotype)))
treatments<-unique(ril_water$treatment)
dap_i<-unique(ril_water$dap_i)
dap_i<-sort(as.numeric(as.character(dap_i)))
times = seq(from = 15, to = 33, by=0.1)

report.loess.values(ril_water, 'water_lost', genos, treatments, dap_i, min(dap_i), max(dap_i), "ril_loess_estimates_water_lost_day.pdf")

## Give output an appropriate name
ril_loess_water_lost_day<-ril_loess_model_fit
## Rename the growth report
water_lost_day_report<-growth_rate_report



#############################################
## Format for QTL pipeline
water.l<-ril_loess_water_lost_day

days<-sort(unique(water.l$dap_i))

ril_water_lost_day_qtl<-c()
for(d in 1:length(days)) {
  day<-days[d]
  temp.data<-water.l[water.l$dap_i == day, ]
  colnames(temp.data)[4]<-paste('water_lost_day', day, sep="_")
  temp.data<-temp.data[,c(1,2,4)]
  if (d == 1) {
    ril_water_lost_day_qtl<-temp.data
  }
  
  if (d > 1) {
    ril_water_lost_day_qtl<-merge(ril_water_lost_day_qtl, temp.data, by=c('genotype', 'treatment'), all=T)
  }
}

# Add misc columns
ril_water_lost_day_qtl$Obs<-c(1:nrow(ril_water_lost_day_qtl))
ril_water_lost_day_qtl$experiment<-rep("BP14", nrow(ril_water_lost_day_qtl))
ril_water_lost_day_qtl$year<-rep("2014", nrow(ril_water_lost_day_qtl))
ril_water_lost_day_qtl$plot<-rep("bellweater", nrow(ril_water_lost_day_qtl))
ril_water_lost_day_qtl$plot_id<-rep("bellweater", nrow(ril_water_lost_day_qtl))
ril_water_lost_day_qtl$sampling<-rep("bellweater", nrow(ril_water_lost_day_qtl))

ril_water_lost_day_qtl<-ril_water_lost_day_qtl[,c(22,23,24,2,25,26,1,27,3:21)]
colnames(ril_water_lost_day_qtl)[7]<-c("id")
write.csv(ril_water_lost_day_qtl, file="ril_loess_water_lost_day_qtl.csv", quote=F, row.names=F)



##############################################
## TRANSPIRATION DAY
## Now go for a loess fit
genos<-as.character(sort(unique(ril_water$genotype)))
treatments<-unique(ril_water$treatment)
dap_i<-unique(ril_water$dap_i)
dap_i<-sort(as.numeric(as.character(dap_i)))
times = seq(from = 15, to = 33, by=0.1)

report.loess.values(ril_water, 'transpiration_day', genos, treatments, dap_i, min(dap_i), max(dap_i), "ril_loess_estimates_transpiration_day.pdf")

## Give output an appropriate name
ril_loess_transpiration_day<-ril_loess_model_fit
## Rename the growth report
transpiration_day_report<-growth_rate_report

#############################################
## Format for QTL pipeline
transpiration.l<-ril_loess_transpiration_day

days<-sort(unique(transpiration.l$dap_i))

ril_transpiration_day_qtl<-c()
for(d in 1:length(days)) {
  day<-days[d]
  temp.data<-transpiration.l[transpiration.l$dap_i == day, ]
  colnames(temp.data)[4]<-paste('transpiration_day', day, sep="_")
  temp.data<-temp.data[,c(1,2,4)]
  if (d == 1) {
    ril_transpiration_day_qtl<-temp.data
  }
  
  if (d > 1) {
    ril_transpiration_day_qtl<-merge(ril_transpiration_day_qtl, temp.data, by=c('genotype', 'treatment'), all=T)
  }
}

## Add misc columns
ril_transpiration_day_qtl$Obs<-c(1:nrow(ril_transpiration_day_qtl))
ril_transpiration_day_qtl$experiment<-rep("BP14", nrow(ril_transpiration_day_qtl))
ril_transpiration_day_qtl$year<-rep("2014", nrow(ril_transpiration_day_qtl))
ril_transpiration_day_qtl$plot<-rep("bellweater", nrow(ril_transpiration_day_qtl))
ril_transpiration_day_qtl$plot_id<-rep("bellweater", nrow(ril_transpiration_day_qtl))
ril_transpiration_day_qtl$sampling<-rep("bellweater", nrow(ril_transpiration_day_qtl))

ril_transpiration_day_qtl<-ril_transpiration_day_qtl[,c(22,23,24,2,25,26,1,27,3:21)]
colnames(ril_transpiration_day_qtl)[7]<-c("id")
write.csv(ril_transpiration_day_qtl, file="ril_loess_transpiration_day_qtl.csv", quote=F, row.names=F)

save.image("ril_water_use.Rdata")


##############################################
## EVAPORATION DAY
## Now go for a loess fit
genos<-as.character(c("EMPTY"))
treatments<-unique(ril_water$treatment)
dap_i<-unique(ril_water$dap_i)
dap_i<-sort(as.numeric(as.character(dap_i)))
times = seq(from = 8, to = 33, by=0.1)

report.loess.values(ril_water, 'evaporation_day', genos, treatments, dap_i, min(dap_i), max(dap_i), "ril_loess_estimates_evaporation_day.pdf")

## Give output an appropriate name
ril_loess_evaporation_day<-ril_loess_model_fit
## Rename the growth report
evaporation_day_report<-growth_rate_report

save.image("ril_water_use.Rdata")


#####################################################################
## Lets get some values for water loss throughout entire experiment (NOT JUST AFTER DAY 15)
#####################################################################

setwd(wue_data.dir)
w<-read.csv("SnapshotInfo_RIL_allday.csv")

w<-w[w$water.amount != -1,]

planting_date = as.POSIXct("2014-1-13")
w$date<-as.POSIXct(w$timestamp, origin = "1970-01-01")
w$dap <- as.numeric(w$date - planting_date)
w$dap_i<-as.integer(w$dap)

## Read a key that associates phenotyper ID (aka 'plant_id') to RIL id
key<-read.csv('phe2ril2.csv', header=F)
colnames(key)<-c('phe', 'ril')

## Make the RIL name same as in QTL map file
key$ril<-sprintf("%03d", key$ril)
key$ril<-paste('RIL_', key$ril, sep="")

## Add a new (empty column called genotype)
w$genotype<-NA
## Add genotypes 
for(i in 1:nrow(key)) {
  phe_name<-key[i,1]
  ril_name<-key[i,2]
  w$genotype[grep(paste('Dr', phe_name, 'A', sep=""), w$plant.barcode)]<-ril_name
}

# Add parental genotypes
w$genotype[grep('Dp1', w$plant.barcode)]<-c("A10")
w$genotype[grep('Dp2', w$plant.barcode)]<-c("B100")

## Add treatment
w$treatment<-c('none')
w$treatment[grep("AA", w$plant.barcode)]<-c('wet')
w$treatment[grep("AB", w$plant.barcode)]<-c('dry')

water_added_dap<-aggregate(w$water.amount, by=list(w$plant.barcode, w$car.tag, w$dap_i, w$genotype, w$treatment), sum)
colnames(water_added_dap)<-c('plantbarcode', 'cartag', 'dap_i', 'genotype', 'treatment','water_added')
weight_after_dap<-aggregate(w$weight.after, by=list(w$plant.barcode, w$car.tag, w$dap_i, w$genotype, w$treatment), min)
colnames(weight_after_dap)<-c('plantbarcode', 'cartag', 'dap_i', 'genotype', 'treatment','weight_total')

ril_water.all<-merge(water_added_dap, weight_after_dap, by=c('plantbarcode', 'cartag', 'dap_i', 'genotype', 'treatment'))

## Might think about how increases in set point for water increases throughout experiment 
p<-ggplot(data=weight_after_dap, aes(x=dap_i, y=weight_total, colour=treatment)) + geom_point(size=0.1)
p
## Looks like there are changes on days 19, and 25
dap<-sort(unique(ril_water.all$dap_i))

water_lost.all<-c()
for (d in 1:length(dap)){
  day<-dap[d]
  if (d == 1) {
    temp<-ril_water.all[ril_water.all$dap_i == day,]
    water_lost.all<-temp[,1:6]
    colnames(water_lost.all)[6]<-c("water_lost")
  }
  if (d > 1){
    temp<-ril_water.all[ril_water.all$dap_i == day,]
    if (day == 19 | day == 25){
      plant_ids<-temp$plantbarcode
      yesterday<-ril_water.all[ril_water.all$plantbarcode %in% plant_ids,]
      yesterday<-yesterday[yesterday$dap_i == (day-1),c(1:2,4:5,7)]
      colnames(yesterday)[5]<-c('weight_yesterday')
      temp3<-merge(temp, yesterday, by=c('plantbarcode', 'cartag', 'genotype', 'treatment'))
      temp3$weight_offset<-temp3$weight_total - temp3$weight_yesterday
      temp3$water_lost<-temp3$water_added - temp3$weight_offset
      temp3<-temp3[,c(1,2,5,3,4,10)]
      water_lost.all<-rbind(water_lost.all, temp3)
      next
    }
    temp2<-temp[temp$water_added > 0,]
    colnames(temp2)[6]<-c("water_lost")
    water_lost.all<-rbind(water_lost.all, temp2[,c(1:6)])
    dry_down<-temp[temp$water_added == 0,]
    if (nrow(dry_down) > 0){
      plant_ids<-dry_down$plantbarcode
      yesterday<-ril_water.all[ril_water.all$plantbarcode %in% plant_ids,]
      yesterday<-yesterday[yesterday$dap_i == (day-1),c(1:2,4:5,7)]
      colnames(yesterday)[5]<-c('weight_yesterday')
      temp3<-merge(dry_down, yesterday, by=c('plantbarcode', 'cartag', 'genotype', 'treatment'))
      temp3$water_lost<-temp3$weight_yesterday - temp3$weight_total
      temp3<-temp3[,c(1,2,5,3,4,9)]
      water_lost.all<-rbind(water_lost.all, temp3)
    }
  }
}
p<-ggplot(data=water_lost.all, aes(x=dap_i, y=water_lost, colour=treatment)) + geom_point(cex=0.3) 
p
p<-ggplot(data=water_lost.all, aes(x=dap_i, y=water_lost, colour=treatment)) + geom_point(cex=0.3) + ylim(0,50)
p
## Remove all values below zero and above 200 (day 17, 32, 33)
water_lost.all<-water_lost.all[!(water_lost.all$water_lost > 200 & water_lost.all$dap_i %in% c(17,32,33)),]
water_lost.all<-water_lost.all[water_lost.all$water_lost > 0,]
p<-ggplot(data=water_lost.all, aes(x=dap_i, y=water_lost, colour=treatment)) + geom_point(cex=0.3) 
p

p<-ggplot(data=water_lost.all[water_lost.all$treatment == 'dry',], aes(x=dap_i, y=water_lost, colour=treatment)) + geom_point(cex=0.3) 
p

## Remove outlier
water_lost.all<-water_lost.all[!(water_lost.all$treatment == 'dry' & water_lost.all$dap_i == 32 & water_lost.all$water_lost > 145),]
p<-ggplot(data=water_lost.all[water_lost.all$treatment == 'wet',], aes(x=dap_i, y=water_lost, colour=treatment)) + geom_point(cex=0.3) 
p

p<-ggplot(data=water_lost.all, aes(x=dap_i, y=water_lost, colour=treatment)) + geom_smooth() + ylim(0,60)
p

setwd(wue_results.water_use.dir)
save.image("ril_water_use.Rdata")

empty_pots<-c("Dr27AA001310", "Dr41AA001393", "Dr53AB001468", "Dr65AB001541", "Dr72AA001581")

## Lets break out the empty pots
## Rename them in the data.frame water total
water_lost.all$genotype<-as.character(water_lost.all$genotype)
water_lost.all[water_lost.all$plantbarcode %in% empty_pots, 'genotype']<-c("EMPTY")

parentals<-water_lost.all[water_lost.all$genotype == "A10" | water_lost.all$genotype == "B100" | water_lost.all$genotype == "EMPTY",]
## Plot the parental lines water use vs. empty
p<-ggplot(parentals, aes(dap_i, water_lost, color=genotype)) + geom_smooth() + facet_wrap(~treatment)
p

water_days<-sort(unique(water_lost.all$dap_i))
water_total.all<-c()
for(d in 1:length(water_days)) {
  day<-water_days[d]
  temp<-water_lost.all[water_lost.all$dap_i <= day,]
  temp<-aggregate(temp[,6], by=list(temp$plantbarcode, temp$cartag, temp$genotype, temp$treatment), sum)
  temp$dap_i<-rep(day, nrow(temp))
  colnames(temp)<-c("plantbarcode", "cartag","genotype","treatment","water_lost_total","dap_i")
  water_total.all<-rbind(water_total.all, temp)
}

## Lets merge this total water lost with the daily water lost
water.m.all<-merge(water_lost.all, water_total.all, by=c("plantbarcode", "cartag", "genotype", "treatment", "dap_i"))
p<-ggplot(water.m.all, aes(dap_i, water_lost_total, color=treatment)) + geom_smooth() 
p

## Build a loess model for each water lost daily and total water lost
## This is to summarize values by genotype. Aggregation may be a better choice.
## Purely for consistency in methods

## daily water lost caclulated here
genos<-as.character(sort(unique(water.m.all$genotype)))
treatments<-unique(water.m.all$treatment)
dap_i<-unique(water.m.all$dap_i)
dap_i<-sort(as.numeric(as.character(dap_i)))
times = seq(from = min(dap_i), to = max(dap_i), by=0.1)

report.loess.values(water.m.all, 'water_lost', genos, treatments, dap_i, min(dap_i), max(dap_i), "ril_loess_estimates_water_lost_day.all.pdf")

## Give output an appropriate name
ril_loess_water_lost_day.all<-ril_loess_model_fit
## Rename the growth report
water_lost_day.all_report<-growth_rate_report

############################################################
# Cumulative water loss
############################################################

genos<-as.character(sort(unique(water.m.all$genotype)))
treatments<-unique(water.m.all$treatment)
dap_i<-unique(water.m.all$dap_i)
dap_i<-sort(as.numeric(as.character(dap_i)))
times = seq(from = min(dap_i), to = max(dap_i), by=0.1)

report.loess.values(water.m.all, 'water_lost_total', genos, treatments, dap_i, min(dap_i), max(dap_i), "ril_loess_estimates_water_lost_total.all.pdf")

## Give output an appropriate name
ril_loess_water_lost_total.all<-ril_loess_model_fit
## Rename the growth report
water_lost_total.all_report<-growth_rate_report

all_day_water_lost<-merge(ril_loess_water_lost_day.all, ril_loess_water_lost_total.all, by=c("genotype", "treatment","dap_i"))

colnames(all_day_water_lost)[4:5]<-c("water_lost_day", "water_lost_day_rate")

save.image("ril_water_use.Rdata")
rm(list=ls())




