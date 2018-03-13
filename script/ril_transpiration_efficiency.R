library(ggplot2)
library(gridExtra)

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


## Combine biomass with water data
setwd(wue_results.plant_size.dir)
load('ril_biomass.Rdata')

setwd(wue_results.water_use.dir)
load('ril_water_use.Rdata')

## Make TE the working directory
setwd(wue_results.te.dir)

## Lets examine how much water was lost during the first two days
p<-ggplot(data=ril_water, aes(x=dap_i, y=water_lost_total,colour=treatment)) + geom_point() + xlim(15,17) +ylim(-5,20)
p

## Lets look at how big the plants are in milligrams during the first two dates
p<-ggplot(data=vis, aes(x=dap_i, y=fw_biomass,colour=treatment)) + geom_point() + xlim(15,17) +ylim(-5,1000)
p

## Lets look at how evaporation between empty pots and pots with plants in them differ
empties<-ril_water[ril_water$genotype == "EMPTY",]
p<-ggplot(data=ril_water, aes(x=dap_i, y=water_lost,colour=treatment)) + geom_smooth()
p + geom_point(data=empties, aes(x=dap_i, y=water_lost,colour=treatment))

#################################################################
# Load in
#################################################################
## Get days 15 and on
vis.all_days<-vis
vis<-vis[vis$dap_i > 14,]

## Lets merge the plan size and water use directories to assess WUE
vis.wue<-merge(vis, ril_water, by=c("plantbarcode", "cartag", "genotype", "treatment", "dap_i"))
vis.wue<-vis.wue[vis.wue$dap_i > 14, ]

## Calculate fresh weight biomass / evapotranspiration total 
vis.wue$fw_water_total<-vis.wue$fw_biomass/vis.wue$water_lost_total

## Calculate dry weight biomass / evapotranspiration total 
vis.wue$dw_water_total<-vis.wue$dw_biomass/vis.wue$water_lost_total

## Calculate sv_area / evapotranspiration total 
vis.wue$sv_water_total<-vis.wue$sv_area/vis.wue$water_lost_total

# Lets look at what evapotranspiration efficiency looks like on day 15
pdf("sv_area_over_water_lost_total_day15.pdf")

p<-ggplot(data=vis.wue, aes(x=dap_i, y=sv_water_total, colour = treatment)) + geom_point() + scale_colour_manual(values = c("orange","navy")) + xlab("Days after planting") + ylab("WUE ratio = SV_area/Total_water_added (ml) ")
p + theme_bw()

dev.off()

# Large artifacts on the first two days. Lets remove these and examine the distribution on day 17
pdf("sv_area_over_water_lost_total_day17.pdf")

p<-ggplot(data=vis.wue, aes(x=dap_i, y=sv_water_total, colour = treatment)) + geom_point() + scale_colour_manual(values = c("orange","navy")) + xlab("Days after planting") + ylab("WUE ratio = SV_area/Total_water_added (ml) ") +xlim(17, 33)
p + theme_bw() +  ylim(0,1000)

dev.off()

## Look at relationship between evapotranspiration and sv_area
p<-ggplot(data=vis.wue[vis.wue$dap_i == 15,], aes(x=sv_area, y=water_lost_total, colour = treatment)) + geom_point()
p

## A large % of the dry pots have not lost any water yet, try day 16
p<-ggplot(data=vis.wue[vis.wue$dap_i == 16,], aes(x=sv_area, y=water_lost_total, colour = treatment)) + geom_point()
p
## Day 16 looks better but still a fair number of 'dry' treatment pots with no water lost

p<-ggplot(data=vis.wue[vis.wue$dap_i == 17,], aes(x=sv_area, y=water_lost_total, colour = treatment)) + geom_point()
p
## Day 17 looks okay

p<-ggplot(data=vis.wue[vis.wue$dap_i == 18,], aes(x=sv_area, y=water_lost_total, colour = treatment)) + geom_point()
p
## Day 18 is for sure


## Plot WUE ratio at day 17
p<-ggplot(data=vis.wue[vis.wue$dap_i == 17,], aes(x=dap_i, y=sv_water_total, colour = treatment)) + geom_point()
p

## Lets remove data before day 15 for calculation of traits for QTL analysis

#vis.wue.all.days<-vis.wue
#vis.wue<-vis.wue[vis.wue$dap_i > 14,]



## Lets plot biomass vs water use over all days and on each day
p<-ggplot(data=vis.wue, aes(x=water_lost_total, y=sv_area, colour=treatment)) + geom_point() + scale_color_manual(values=c("orange","navy"))
p + theme_bw()

## Lets get another data frame for day 17 on
vis.wue.17<-vis.wue[vis.wue$dap_i > 16,]

## Lets facet the plot by days, WUE/TE is the deviation from the relationship between total sv_area and water_lost_total
pdf("sv_area_total_v.water_lost_total.pdf")
p<-ggplot(data=vis.wue.17, aes(x=water_lost_total, y=sv_area, colour=treatment)) + geom_point() + facet_wrap(~dap_i) + scale_color_manual(values=c("orange","navy"))
p + theme_bw()
dev.off()

# Lets examine the relationship between sv_area and wue_ratio
pdf("wue_ratio_v.sv_area.pdf")
p<-ggplot(data=vis.wue.17, aes(x=sv_area, y=sv_water_total, colour=treatment)) + geom_point() + facet_wrap(~dap_i) + scale_color_manual(values=c("orange","navy"))
p + theme_bw()
dev.off()

# Lets try making a series of histograms of each sv_water_total facetted by day
pdf("histogram_of_wue_ratio.pdf")
p<-ggplot(data=vis.wue.17, aes(sv_water_total, colour=treatment)) +  geom_density(alpha = 0.2) + facet_wrap(~dap_i) + scale_color_manual(values=c("orange","navy"))
p + theme_bw()
dev.off()


## Okay now that we see that day 17 is the best day 
## Lets get better estimates of water use from day 15 on...

setwd(wue_results.water_use.dir)

# WATER_DAY
# Now go for a loess fit
genos<-as.character(sort(unique(vis.wue$genotype)))
treatments<-unique(vis.wue$treatment)
dap_i<-unique(vis.wue$dap_i)
dap_i<-sort(as.numeric(as.character(dap_i)))
times = seq(from = 15, to = 33, by=0.1)


report.loess.values(ril_water, 'water_lost', genos, treatments, dap_i, min(times), max(times), "ril_loess_estimates_water_lost_day.15_dap_i.pdf")


# Give the output a more informative name
ril_loess_water_lost_day<-ril_loess_model_fit
colnames(ril_loess_water_lost_day)<-c("genotype", "treatment", "dap_i", "water_lost_day", "water_lost_day_rate")
# Rename the growth report
water_lost_day_report<-growth_rate_report

#############################################
## Format for QTL pipeline
water_day_reformat<-ril_loess_water_lost_day

days<-sort(unique(water_day_reformat$dap_i))

ril_water_lost_qtl<-c()
for(d in 1:length(days)) {
  day<-days[d]
  temp.data<-water_day_reformat[water_day_reformat$dap_i == day, ]
  colnames(temp.data)[4]<-paste('water_lost', day, sep="_")
  temp.data<-temp.data[,c(1,2,4)]
  if (d == 1) {
    ril_water_lost_qtl<-temp.data
  }
  
  if (d > 1) {
    ril_water_lost_qtl<-merge(ril_water_lost_qtl, temp.data, by=c('genotype', 'treatment'), all=T)
  }
}

# Add misc columns
ril_water_lost_qtl$Obs<-c(1:nrow(ril_water_lost_qtl))
ril_water_lost_qtl$experiment<-rep("BP14", nrow(ril_water_lost_qtl))
ril_water_lost_qtl$year<-rep("2014", nrow(ril_water_lost_qtl))
ril_water_lost_qtl$plot<-rep("bellweater", nrow(ril_water_lost_qtl))
ril_water_lost_qtl$plot_id<-rep("bellweater", nrow(ril_water_lost_qtl))
ril_water_lost_qtl$sampling<-rep("bellweater", nrow(ril_water_lost_qtl))

## Day 17
ril_water_lost_qtl<-ril_water_lost_qtl[,c(22,23,24,2,25,26,1,27,5:21)]

colnames(ril_water_lost_qtl)[7]<-c("id")
write.csv(ril_water_lost_qtl, file="ril_loess_water_lost_day_qtl.csv", quote=F, row.names=F)



##################################################################################
## Lets calculate biomass/size accumulated per day, and water use/transpiration efficiency
## Calculated on a relative basis
##################################################################################
setwd(wue_results.water_use.dir)

## Now lets combine loess summaries
colnames(ril_loess_fw_biomass)[c(1,4,5)]<-c('genotype', 'fw_biomass', 'fw_biomass_rate')
colnames(ril_loess_sv_area)[c(1,4,5)]<-c('genotype', 'sv_area', 'sv_area_rate')


ril_loess_trait<-merge(ril_loess_fw_biomass, ril_loess_sv_area, by=c('genotype', 'treatment', 'dap_i'))
ril_loess_trait<-merge(ril_loess_trait, ril_loess_water_lost_day, by=c('genotype', 'treatment', 'dap_i'))
ril_loess_trait<-merge(ril_loess_trait, ril_loess_water_lost_total, by=c('genotype', 'treatment', 'dap_i'))

## Lets generate a similar data frame for loess summaries
ids<-unique(ril_loess_trait$genotype)
ril_loess_trait.r<-c()
for (i in ids){
  temp<-ril_loess_trait[ril_loess_trait$genotype == i,]
  temp.w<-temp[temp$treatment == 'wet',]
  if(nrow(temp.w) > 0) {
    fw_max<-max(temp.w$fw_biomass)
    temp.w$fw_biomass.r<-temp.w$fw_biomass/fw_max
    sv_max<-max(temp.w$sv_area)
    temp.w$sv_area.r<-temp.w$sv_area/sv_max
    water_lost_total_max<-max(temp.w$water_lost_total)
    temp.w$water_lost_total.r<-temp.w$water_lost_total/water_lost_total_max
  
    water_lost_day_max<-max(temp.w$water_lost_day)
    temp.w$water_lost_day.r<-temp.w$water_lost_day/water_lost_day_max
   
  }
  # Now dry 
  temp.d<-temp[temp$treatment == 'dry',]
  if (nrow(temp.d) > 0) {
    fw_max<-max(temp.d$fw_biomass)
    temp.d$fw_biomass.r<-temp.d$fw_biomass/fw_max
    sv_max<-max(temp.d$sv_area)
    temp.d$sv_area.r<-temp.d$sv_area/sv_max
    water_lost_total_max<-max(temp.d$water_lost_total)
    temp.d$water_lost_total.r<-temp.d$water_lost_total/water_lost_total_max
    water_lost_day_max<-max(temp.d$water_lost_day)
    temp.d$water_lost_day.r<-temp.d$water_lost_day/water_lost_day_max
  }
  temp<-rbind(temp.d, temp.w)
  ril_loess_trait.r<-rbind(ril_loess_trait.r, temp)
  
}

################################################################
## Lets prepare some daily efficiency measurements
## *** NOTICE ***
## The way we calculate daily WUE ratio is different than previous water and biomass values. 
## In this case we fit before each trait before division rather than the opposite (divide (sv_area/water_lost)) then fit.
################################################################

days<-sort(unique(ril_loess_trait.r$dap_i))

day_size<-c()
for (d in 1:length(days)){
  day<-days[d]
  today<-ril_loess_trait.r[ril_loess_trait.r$dap_i == day,]
  if (d == 1) {
    genotype<-today$genotype
    treatment<-today$treatment
    dap_i<-today$dap_i
    fw_biomass_day<-today$fw_biomass
    sv_area_day<-today$sv_area
    day_size<-cbind(genotype, treatment, dap_i, fw_biomass_day, sv_area_day)
  }
  
  if (d > 1) {
    yesterday<-ril_loess_trait.r[ril_loess_trait.r$dap_i == (day - 1),]
    yesterday<-yesterday[,-c(3)]
    colnames(yesterday)[c(3,5)]<-c('fw_biomass.y', 'sv_area.y')
    today<-merge(today[,c(1:4,6)], yesterday[,c(1:3,5)], by=c("genotype", "treatment"))
    genotype<-today$genotype 
    treatment<-today$treatment
    dap_i<-today$dap_i
    fw_biomass_day<-today$fw_biomass - today$fw_biomass.y
    sv_area_day<-today$sv_area - today$sv_area.y
    temp<-cbind(genotype, treatment, dap_i, fw_biomass_day, sv_area_day)
    day_size<-rbind(day_size, temp)
  }
}

day_size<-as.data.frame(day_size)
day_size$dap_i<-as.numeric(as.character(day_size$dap_i))
day_size$sv_area_day<-as.numeric(as.character(day_size$sv_area_day))
day_size$fw_biomass_day<-as.numeric(as.character(day_size$fw_biomass_day))
day_size<-day_size[day_size$dap_i > 16,]

## Notice that a few of the values of sv_area_day and fw_biomass_day are negative
p<-ggplot(data=day_size, aes(x=dap_i, y=sv_area_day, colour=treatment)) + geom_point(size=0.3) 
p
p<-ggplot(data=day_size, aes(x=dap_i, y=fw_biomass_day, colour=treatment)) + geom_point(size=0.3) 
p

## This is impossible, lets set all negative values to 0
day_size[day_size$sv_area_day < 0, 'sv_area_day']<-c(0)
day_size[day_size$fw_biomass_day < 0, 'fw_biomass_day']<-c(0)

# Merge with other traits 
day_size<-merge(ril_loess_trait.r, day_size, by=c("genotype", "treatment", "dap_i"))

setwd(wue_results.te.dir)

# Calculate WUE ratio
day_size$wue_day<-day_size$sv_area_day/day_size$water_lost_day
day_size$wue_total<-day_size$sv_area/day_size$water_lost_total
day_size$fw_wue_ratio_total<-day_size$fw_biomass/day_size$water_lost_total
day_size$fw_wue_ratio_day<-day_size$fw_biomass_day/day_size$water_lost_day

## Look at plots of these
p<-ggplot(data=day_size, aes(x=dap_i, y=wue_day, colour=treatment)) + geom_smooth()
p
p<-ggplot(data=day_size, aes(x=dap_i, y=fw_wue_ratio_day, colour=treatment)) + geom_smooth()
p
p<-ggplot(data=day_size, aes(x=dap_i, y=wue_total, colour=treatment)) + geom_smooth()
p

#############################################
## Format for QTL pipeline
wue.day.l<-day_size
wue.day.l<-wue.day.l[,c('genotype', 'treatment', 'dap_i', 'wue_day')]

days<-sort(unique(wue.day.l$dap_i))

ril_wue_ratio_day_qtl<-c()
for(d in 1:length(days)) {
  day<-days[d]
  temp.data<-wue.day.l[wue.day.l$dap_i == day, ]
  colnames(temp.data)[4]<-paste('wue_day', day, sep="_")
  temp.data<-temp.data[,c(1,2,4)]
  if (d == 1) {
    ril_wue_ratio_day_qtl<-temp.data
  }
  
  if (d > 1) {
    ril_wue_ratio_day_qtl<-merge(ril_wue_ratio_day_qtl, temp.data, by=c('genotype', 'treatment'), all=T)
  }
}

# Add misc columns
ril_wue_ratio_day_qtl$Obs<-c(1:nrow(ril_wue_ratio_day_qtl))
ril_wue_ratio_day_qtl$experiment<-rep("BP14", nrow(ril_wue_ratio_day_qtl))
ril_wue_ratio_day_qtl$year<-rep("2014", nrow(ril_wue_ratio_day_qtl))
ril_wue_ratio_day_qtl$plot<-rep("bellweater", nrow(ril_wue_ratio_day_qtl))
ril_wue_ratio_day_qtl$plot_id<-rep("bellweater", nrow(ril_wue_ratio_day_qtl))
ril_wue_ratio_day_qtl$sampling<-rep("bellweater", nrow(ril_wue_ratio_day_qtl))

# Reorder columns
ril_wue_ratio_day_qtl<-ril_wue_ratio_day_qtl[,c(20,21,22,2,23,24,1,25,3:19)]

colnames(ril_wue_ratio_day_qtl)[7]<-c("id")
write.csv(ril_wue_ratio_day_qtl, file="ril_loess_wue_ratio_day_qtl.csv", quote=F, row.names=F)

####################################################
## Now fw_biomass_day/water_day
## Format for QTL pipeline

setwd(wue_results.te.dir)

wue.day.fw.l<-day_size
wue.day.fw.l<-wue.day.fw.l[,c('genotype', 'treatment', 'dap_i', 'fw_wue_ratio_day')]

days<-sort(unique(wue.day.fw.l$dap_i))

ril_fw_water_day_qtl<-c()
for(d in 1:length(days)) {
  day<-days[d]
  temp.data<-wue.day.fw.l[wue.day.fw.l$dap_i == day, ]
  colnames(temp.data)[4]<-paste('fw_water_day', day, sep="_")
  temp.data<-temp.data[,c(1,2,4)]
  if (d == 1) {
    ril_fw_water_day_qtl<-temp.data
  }
  
  if (d > 1) {
    ril_fw_water_day_qtl<-merge(ril_fw_water_day_qtl, temp.data, by=c('genotype', 'treatment'), all=T)
  }
}

# Add misc columns
ril_fw_water_day_qtl$Obs<-c(1:nrow(ril_fw_water_day_qtl))
ril_fw_water_day_qtl$experiment<-rep("BP14", nrow(ril_fw_water_day_qtl))
ril_fw_water_day_qtl$year<-rep("2014", nrow(ril_fw_water_day_qtl))
ril_fw_water_day_qtl$plot<-rep("bellweater", nrow(ril_fw_water_day_qtl))
ril_fw_water_day_qtl$plot_id<-rep("bellweater", nrow(ril_fw_water_day_qtl))
ril_fw_water_day_qtl$sampling<-rep("bellweater", nrow(ril_fw_water_day_qtl))

## Day 17
ril_fw_water_day_qtl<-ril_fw_water_day_qtl[,c(20,21,22,2,23,24,1,25,3:19)]

colnames(ril_fw_water_day_qtl)[7]<-c("id")
write.csv(ril_fw_water_day_qtl, file="ril_loess_fw_wue_day_qtl.csv", quote=F, row.names=F)


#########################################################################################
## Lets make some plots of the trends observed in each line
#########################################################################################

## ***** These take a long time to print because your doing a fit for each gentoype with geom_smooth()
setwd(wue_results.plant_size.dir)

## Size per day
genos<-(unique(day_size$genotype))
pdf("ril_loess_estimtates_sv_area_per_day.pdf", paper = "a4", width = 21/2.54, height = (7*5.3)/2.54)
layout(matrix(c(1:4), ncol = 1, byrow = T), heights = c(1, 1, 1, 1))
for (g in genos){
  temp<-day_size[day_size$genotype == g,]
  p<-ggplot(data=temp, aes(x=dap_i, y=sv_area_day, colour=treatment)) + geom_smooth() + scale_color_manual(values=c("orange","navy")) + labs(title=g)
  p<-p + theme_bw()
  print(p)
}
dev.off()

genos<-(unique(day_size$genotype))
pdf("ril_loess_estimates_fw_biomass_per_day.pdf", paper = "a4", width = 21/2.54, height = (7*5.3)/2.54)
layout(matrix(c(1:4), ncol = 1, byrow = T), heights = c(1, 1, 1, 1))
for (g in genos){
  temp<-day_size[day_size$genotype == g,]
  p<-ggplot(data=temp, aes(x=dap_i, y=fw_biomass_day, colour=treatment)) + geom_smooth() + scale_color_manual(values=c("orange","navy")) + labs(title=g)
  p<-p + theme_bw()
  print(p)
}
dev.off()

setwd(wue_results.te.dir)

# WUE ratio per day
genos<-(unique(day_size$genotype))
pdf("ril_loess_estimates_wue_day.pdf", paper = "a4", width = 21/2.54, height = (7*5.3)/2.54)
layout(matrix(c(1:4), ncol = 1, byrow = T), heights = c(1, 1, 1, 1))
for (g in genos){
  temp<-day_size[day_size$genotype == g,]
  p<-ggplot(data=temp, aes(x=dap_i, y=wue_day, colour=treatment)) + geom_smooth() + scale_color_manual(values=c("orange","navy")) + labs(title=g) + theme_bw()
  p<-p + theme_bw()
  print(p)
}
dev.off()

genos<-(unique(day_size$genotype))
pdf("ril_loess_estimates_fw_wue_ratio_per_day.pdf", paper = "a4", width = 21/2.54, height = (7*5.3)/2.54)
layout(matrix(c(1:4), ncol = 1, byrow = T), heights = c(1, 1, 1, 1))
for (g in genos){
  temp<-day_size[day_size$genotype == g,]
  p<-ggplot(data=temp, aes(x=dap_i, y=fw_wue_ratio_day, colour=treatment)) + geom_smooth() + scale_color_manual(values=c("orange","navy")) + labs(title=g)
  p<-p + theme_bw()
  print(p)
}
dev.off()


#############################################
# Need to write calculation of sv_area per day to .csv in QTL format
#############################################
# Format for QTL pipeline
setwd(wue_results.plant_size.dir)

# day_size
size.day.l<-day_size
size.day.l<-size.day.l[,c('genotype', 'treatment', 'dap_i', 'sv_area')]

days<-sort(unique(size.day.l$dap_i))

ril_sv_area_day_qtl<-c()
for(d in 1:length(days)) {
  day<-days[d]
  temp.data<-size.day.l[size.day.l$dap_i == day, ]
  colnames(temp.data)[4]<-paste('sv_area_day', day, sep="_")
  temp.data<-temp.data[,c(1,2,4)]
  if (d == 1) {
    ril_sv_area_day_qtl<-temp.data
  }
  
  if (d > 1) {
    ril_sv_area_day_qtl<-merge(ril_sv_area_day_qtl, temp.data, by=c('genotype', 'treatment'), all=T)
  }
}

# Add misc columns
ril_sv_area_day_qtl$Obs<-c(1:nrow(ril_sv_area_day_qtl))
ril_sv_area_day_qtl$experiment<-rep("BP14", nrow(ril_sv_area_day_qtl))
ril_sv_area_day_qtl$year<-rep("2014", nrow(ril_sv_area_day_qtl))
ril_sv_area_day_qtl$plot<-rep("bellweater", nrow(ril_sv_area_day_qtl))
ril_sv_area_day_qtl$plot_id<-rep("bellweater", nrow(ril_sv_area_day_qtl))
ril_sv_area_day_qtl$sampling<-rep("bellweater", nrow(ril_sv_area_day_qtl))

# Reorder columns
ril_sv_area_day_qtl<-ril_sv_area_day_qtl[,c(20,21,22,2,23,24,1,25,3:19)]

colnames(ril_sv_area_day_qtl)[7]<-c("id")
setwd(wue_results.plant_size.dir)
write.csv(ril_sv_area_day_qtl, file="ril_loess_sv_area_day_qtl.csv", quote=F, row.names=F)



#############################################################################################
## Lets use the ratio of the fit values to calculate the WUE ratio total

setwd(wue_results.te.dir)

days<-sort(unique(day_size$dap_i))
genos<-sort(unique(day_size$genotype))

## WUE ratio
genos<-(unique(day_size$genotype))
pdf("ril_loess_estimtates_wue_ratio_total.pdf", paper = "a4", width = 21/2.54, height = (7*5.3)/2.54)
layout(matrix(c(1:4), ncol = 1, byrow = T), heights = c(1, 1, 1, 1))
for (g in genos){
  temp<-day_size[day_size$genotype == g,]
  p<-ggplot(data=temp, aes(x=dap_i, y=wue_total, colour=treatment)) + geom_smooth() + scale_color_manual(values=c("orange","navy")) + labs(title=g) 
  p<-p + theme_bw()
  print(p)
}
dev.off()

###################################################
## Write to QTL pipeline

wue.total.l<-day_size
wue.total.l<-wue.total.l[,c('genotype', 'treatment', 'dap_i', 'wue_total')]

days<-sort(unique(wue.total.l$dap_i))

ril_wue_ratio_total_qtl<-c()
for(d in 1:length(days)) {
  day<-days[d]
  temp.data<-wue.total.l[wue.total.l$dap_i == day, ]
  colnames(temp.data)[4]<-paste('wue_total', day, sep="_")
  temp.data<-temp.data[,c(1,2,4)]
  if (d == 1) {
    ril_wue_ratio_total_qtl<-temp.data
  }
  
  if (d > 1) {
    ril_wue_ratio_total_qtl<-merge(ril_wue_ratio_total_qtl, temp.data, by=c('genotype', 'treatment'), all=T)
  }
}

# Add misc columns
ril_wue_ratio_total_qtl$Obs<-c(1:nrow(ril_wue_ratio_total_qtl))
ril_wue_ratio_total_qtl$experiment<-rep("BP14", nrow(ril_wue_ratio_total_qtl))
ril_wue_ratio_total_qtl$year<-rep("2014", nrow(ril_wue_ratio_total_qtl))
ril_wue_ratio_total_qtl$plot<-rep("bellweater", nrow(ril_wue_ratio_total_qtl))
ril_wue_ratio_total_qtl$plot_id<-rep("bellweater", nrow(ril_wue_ratio_total_qtl))
ril_wue_ratio_total_qtl$sampling<-rep("bellweater", nrow(ril_wue_ratio_total_qtl))

# Reorder columns
ril_wue_ratio_total_qtl<-ril_wue_ratio_total_qtl[,c(20,21,22,2,23,24,1,25,3:19)]

setwd(wue_results.te.dir)

colnames(ril_wue_ratio_total_qtl)[7]<-c("id")
write.csv(ril_wue_ratio_total_qtl, file="ril_loess_wue_ratio_total_qtl.csv", quote=F, row.names=F)


#################################################################################
# Linear model to calculate TE
#################################################################################

## Use the loess fit values to estimate model fit and model residual (+/- sv_area given the ratio sv_area~water_lost)

setwd(wue_results.te.dir)

## Lets run a loop and fit model sv_area / water used 
## will get model fit and the residual of the model fit

days<-sort(unique(day_size$dap_i))
te_model<-c()
for(d in days){
  temp<-day_size[day_size$dap_i == d,]
  treats<-sort(unique(temp$treatment))
  temp.t<-c()
  for(t in treats){
    temp2<-temp[temp$treatment == t,]
    te_total.mdl<-lm(sv_area~water_lost_total, data=temp2)
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


## Format for te_fit_total QTL pipeline
te_fit.l<-te_model
te_fit.l<-te_fit.l[,c('genotype', 'treatment', 'dap_i', 'te_fit_total')]

days<-sort(unique(te_fit.l$dap_i))

ril_te_fit_total_qtl<-c()
for(d in 1:length(days)) {
  day<-days[d]
  temp.data<-te_fit.l[te_fit.l$dap_i == day, ]
  colnames(temp.data)[4]<-paste('te_fit_total', day, sep="_")
  temp.data<-temp.data[,c(1,2,4)]
  if (d == 1) {
    ril_te_fit_total_qtl<-temp.data
  }
  
  if (d > 1) {
    ril_te_fit_total_qtl<-merge(ril_te_fit_total_qtl, temp.data, by=c('genotype', 'treatment'), all=T)
  }
}

## Add misc columns
ril_te_fit_total_qtl$Obs<-c(1:nrow(ril_te_fit_total_qtl))
ril_te_fit_total_qtl$experiment<-rep("BP14", nrow(ril_te_fit_total_qtl))
ril_te_fit_total_qtl$year<-rep("2014", nrow(ril_te_fit_total_qtl))
ril_te_fit_total_qtl$plot<-rep("bellweater", nrow(ril_te_fit_total_qtl))
ril_te_fit_total_qtl$plot_id<-rep("bellweater", nrow(ril_te_fit_total_qtl))
ril_te_fit_total_qtl$sampling<-rep("bellweater", nrow(ril_te_fit_total_qtl))



## Reorder columns
ril_te_fit_total_qtl<-ril_te_fit_total_qtl[,c(20,21,22,2,23,24,1,25,3:19)]

colnames(ril_te_fit_total_qtl)[7]<-c("id")
write.csv(ril_te_fit_total_qtl, file="ril_loess_te_fit_total_qtl.csv", quote=F, row.names=F)

## Format for te_residual_total QTL pipeline
te_residual.l<-te_model
te_residual.l<-te_residual.l[,c('genotype', 'treatment', 'dap_i', 'te_residual_total')]

days<-sort(unique(te_residual.l$dap_i))

ril_te_residual_total_qtl<-c()
for(d in 1:length(days)) {
  day<-days[d]
  temp.data<-te_residual.l[te_residual.l$dap_i == day, ]
  colnames(temp.data)[4]<-paste('te_residual_total', day, sep="_")
  temp.data<-temp.data[,c(1,2,4)]
  if (d == 1) {
    ril_te_residual_total_qtl<-temp.data
  }
  
  if (d > 1) {
    ril_te_residual_total_qtl<-merge(ril_te_residual_total_qtl, temp.data, by=c('genotype', 'treatment'), all=T)
  }
}

## Add misc columns
ril_te_residual_total_qtl$Obs<-c(1:nrow(ril_te_residual_total_qtl))
ril_te_residual_total_qtl$experiment<-rep("BP14", nrow(ril_te_residual_total_qtl))
ril_te_residual_total_qtl$year<-rep("2014", nrow(ril_te_residual_total_qtl))
ril_te_residual_total_qtl$plot<-rep("bellweater", nrow(ril_te_residual_total_qtl))
ril_te_residual_total_qtl$plot_id<-rep("bellweater", nrow(ril_te_residual_total_qtl))
ril_te_residual_total_qtl$sampling<-rep("bellweater", nrow(ril_te_residual_total_qtl))

## Reorder columns
ril_te_residual_total_qtl<-ril_te_residual_total_qtl[,c(20,21,22,2,23,24,1,25,3:19)]

colnames(ril_te_residual_total_qtl)[7]<-c("id")
write.csv(ril_te_residual_total_qtl, file="ril_loess_te_residual_total_qtl.csv", quote=F, row.names=F)


#################################################
## Now rate statistics
#################################################

setwd(wue_results.te.dir)

## Format te_fit_day for QTL pipeline
te_fit.day.l<-te_model
te_fit.day.l<-te_fit.day.l[,c('genotype', 'treatment', 'dap_i', 'te_fit_day')]

days<-sort(unique(te_fit.day.l$dap_i))

ril_te_fit_day_qtl<-c()
for(d in 1:length(days)) {
  day<-days[d]
  temp.data<-te_fit.day.l[te_fit.day.l$dap_i == day, ]
  colnames(temp.data)[4]<-paste('te_fit_day', day, sep="_")
  temp.data<-temp.data[,c(1,2,4)]
  if (d == 1) {
    ril_te_fit_day_qtl<-temp.data
  }
  
  if (d > 1) {
    ril_te_fit_day_qtl<-merge(ril_te_fit_day_qtl, temp.data, by=c('genotype', 'treatment'), all=T)
  }
}

## Add misc columns
ril_te_fit_day_qtl$Obs<-c(1:nrow(ril_te_fit_day_qtl))
ril_te_fit_day_qtl$experiment<-rep("BP14", nrow(ril_te_fit_day_qtl))
ril_te_fit_day_qtl$year<-rep("2014", nrow(ril_te_fit_day_qtl))
ril_te_fit_day_qtl$plot<-rep("bellweater", nrow(ril_te_fit_day_qtl))
ril_te_fit_day_qtl$plot_id<-rep("bellweater", nrow(ril_te_fit_day_qtl))
ril_te_fit_day_qtl$sampling<-rep("bellweater", nrow(ril_te_fit_day_qtl))

## Reorder columns
ril_te_fit_day_qtl<-ril_te_fit_day_qtl[,c(20,21,22,2,23,24,1,25,3:19)]

colnames(ril_te_fit_day_qtl)[7]<-c("id")
write.csv(ril_te_fit_day_qtl, file="ril_loess_te_fit_day_qtl.csv", quote=F, row.names=F)
 
## Format te_residual_day for QTL pipeline
setwd(wue_results.te.dir)

te_residual.day.l<-te_model
te_residual.day.l<-te_residual.day.l[,c('genotype', 'treatment', 'dap_i', 'te_residual_day')]

days<-sort(unique(te_residual.day.l$dap_i))

ril_te_residual_day_qtl<-c()
for(d in 1:length(days)) {
  day<-days[d]
  temp.data<-te_residual.day.l[te_residual.day.l$dap_i == day, ]
  colnames(temp.data)[4]<-paste('te_residual_day', day, sep="_")
  temp.data<-temp.data[,c(1,2,4)]
  if (d == 1) {
    ril_te_residual_day_qtl<-temp.data
  }
  
  if (d > 1) {
    ril_te_residual_day_qtl<-merge(ril_te_residual_day_qtl, temp.data, by=c('genotype', 'treatment'), all=T)
  }
}

## Add misc columns
ril_te_residual_day_qtl$Obs<-c(1:nrow(ril_te_residual_day_qtl))
ril_te_residual_day_qtl$experiment<-rep("BP14", nrow(ril_te_residual_day_qtl))
ril_te_residual_day_qtl$year<-rep("2014", nrow(ril_te_residual_day_qtl))
ril_te_residual_day_qtl$plot<-rep("bellweater", nrow(ril_te_residual_day_qtl))
ril_te_residual_day_qtl$plot_id<-rep("bellweater", nrow(ril_te_residual_day_qtl))
ril_te_residual_day_qtl$sampling<-rep("bellweater", nrow(ril_te_residual_day_qtl))

## Reorder columns
ril_te_residual_day_qtl<-ril_te_residual_day_qtl[,c(20,21,22,2,23,24,1,25,3:19)]

colnames(ril_te_residual_day_qtl)[7]<-c("id")
write.csv(ril_te_residual_day_qtl, file="ril_loess_te_residual_day_qtl.csv", quote=F, row.names=F)



##########################################################
## Lets examine the model results here
##########################################################
setwd(wue_results.te.dir)

## Residual plot of daily TE
pdf("lm_residual_check_te_fit_day.pdf")

p<-ggplot(te_model, aes(x=te_fit_day, y=te_residual_day, colour=treatment)) + geom_point(size=0.3) + facet_wrap(~dap_i) + scale_color_manual(values=c("orange", "navy")) + theme_bw() + xlab("sv_area ~ water_use_total") + ylab("residual")+ theme(legend.position="none")
q<-p + geom_hline(yintercept=0, linetype=2, color='red')
print(q)

dev.off()

## Plot residual cumulative X sv_area cumulative

pdf("lm_residual_check_te_fit_total.pdf")
p<-ggplot(te_model, aes(x=sv_area, y=te_residual_total, colour=treatment)) + geom_point() + facet_wrap(~dap_i)+ scale_color_manual(values=c("orange", "navy")) + theme_bw()
q<-p + geom_hline(yintercept=0, linetype=2, color='red')
print(q)
dev.off()

## Lets plot these measurements on a line, treatment basis..

## Plot residual total X sv_area totoal
pdf("te_residual_total_vs_sv_area_total.pdf")
p<-ggplot(te_model, aes(x=sv_area, y=te_residual_day, colour=treatment)) + geom_point() + facet_wrap(~dap_i)+ scale_color_manual(values=c("orange", "navy")) + theme_bw()
p
dev.off()

## Plot residual day X sv_area day
pdf("te_residual_day_vs_sv_area_day.pdf")
p<-ggplot(te_model, aes(x=sv_area_day, y=te_residual_day, colour=treatment)) + geom_point() + facet_wrap(~dap_i)+ scale_color_manual(values=c("orange", "navy")) + theme_bw()
p
dev.off()

pdf("te_fit_total_vs_sv_area_total.pdf")
p<-ggplot(te_model, aes(x=sv_area, y=te_fit_total, colour=treatment)) + geom_point() + facet_wrap(~dap_i)+ scale_color_manual(values=c("orange", "navy")) + theme_bw()
p
dev.off()

pdf("te_fit_day_vs_sv_area_day.pdf")
p<-ggplot(te_model, aes(x=sv_area_day, y=te_fit_day, colour=treatment)) + geom_point() + facet_wrap(~dap_i)+ scale_color_manual(values=c("orange", "navy")) + theme_bw()
p
dev.off()



## Lets make some plots of all individuals
setwd(wue_results.te.dir)
genotypes<-sort(unique(te_model$genotype))

pdf("te_model_by_genotype.pdf",onefile = T)

for(g in genotypes){
  temp<-te_model[te_model$genotype == g,]
  ## Plot residual
  p<-ggplot(temp, aes(x=dap_i, y=te_residual_total, colour=treatment)) + geom_line() + scale_color_manual(values=c("orange", "navy")) + theme_bw()
  p<-p + ylab("residual_cumulative") + xlab("Days after planting") + ggtitle(g)
  #print(p)
  q<-ggplot(temp, aes(x=dap_i, y=te_residual_day, colour=treatment)) + geom_line() + scale_color_manual(values=c("orange", "navy")) + theme_bw()
  q<-q+ylab("residual_daily") + xlab("Days after planting")
  #print(q)
  x<-ggplot(temp, aes(x=dap_i, y=te_fit_total, colour=treatment)) + geom_line() + scale_color_manual(values=c("orange", "navy")) + theme_bw()
  x<-x+ylab("fit_cumulative") + xlab("Days after planting") 
  #print(x)
  y<-ggplot(temp, aes(x=dap_i, y=te_fit_day, colour=treatment)) + geom_line() + scale_color_manual(values=c("orange", "navy")) + theme_bw()
  y<-y+ylab("fit_daily") + xlab("Days after planting")
  #print(y)
  grid.arrange(p,q,x,y)
}

dev.off()

save.image('ril_transpiration_efficiency.Rdata')



#############################################################
## Lets look at the ratio of total biomass / total pot water content
## This will give us an idea of how much error there is in watering
#############################################################
setwd(wue_results.te.dir)
load("ril_transpiration_efficiency.Rdata")

## Lets find the average pot weight of each treatment block
## w is name of the data.frame with raw watering data directly from the lemnatec instrument
wet.full.pot.wt<-as.integer(mean(w[w$treatment == 'wet', 'weight.after']))
dry.full.pot.wt<-as.integer(mean(w[w$treatment == 'dry', 'weight.after']))

## Subtract the weight of the carrier (335 grams) and weight of the soil (73 grams) 
carrier.wt<-c(335)
soil.wt<-c(73)
## See methods of Fahlgren et al., 2015

wet.water.weight<-wet.full.pot.wt - carrier.wt - soil.wt
dry.water.weight<-dry.full.pot.wt - carrier.wt - soil.wt

te_model$err.ratio<-rep("NA", nrow(te_model))
te_model[te_model$treatment == 'wet', 'err.ratio']<-(te_model[te_model$treatment == 'wet', 'fw_biomass']/1000)/wet.water.weight
te_model[te_model$treatment == 'dry', 'err.ratio']<-(te_model[te_model$treatment == 'dry', 'fw_biomass']/1000)/dry.water.weight
te_model$err.ratio<-as.numeric(as.character(te_model$err.ratio))

#te_residual$err.ratio<-(te_residual$fw_biomass/1000)/te_residual$water_lost_day


p<-ggplot(te_model, aes(x=dap_i, y=err.ratio, color=treatment)) + geom_smooth()
q<-p + ylab("Biomass (g) / Pot water volume (mL)") + xlab("Days after planting") + scale_color_manual(values=c("orange", "navy")) + theme_bw() + theme(legend.position="none")


pdf('ratio_of_biomass_to_total_pot_water.pdf')
q
dev.off()


#############################################################
## Lets make some plots of the total variance of composite traits
#############################################################

treatments<-sort(unique(te_model$treatment))
days<-sort(unique(te_model$dap_i))
dap<-c()
trmt<-c()
var.te_fit_total<-c()
var.te_residual_total<-c()
var.te_fit_day<-c()
var.te_residual_day<-c()
for(d in days) {
  temp<-te_model[te_model$dap_i == d,]
  for(t in treatments){
  temp2<-temp[temp$treatment == t,]
  temp.var.te_fit_total<-var(temp2$te_fit_total)
  temp.var.te_residual_total<-var(temp2$te_residual_total)
  temp.var.te_fit_day<-var(temp2$te_fit_day)
  temp.var.te_residual_day<-var(temp2$te_residual_day)
  dap<-c(dap,d)
  trmt<-c(trmt, t)
  var.te_fit_total<-c(var.te_fit_total, temp.var.te_fit_total)
  var.te_residual_total<-c(var.te_residual_total, temp.var.te_residual_total)
  var.te_fit_day<-c(var.te_fit_day, temp.var.te_fit_day)
  var.te_residual_day<-c(var.te_residual_day, temp.var.te_residual_day)
  }
}

comp.variance<-cbind(dap, trmt, var.te_fit_total, var.te_residual_total, var.te_fit_day, var.te_residual_day)
comp.variance<-as.data.frame(comp.variance)

comp.variance$dap<-as.numeric(as.character(comp.variance$dap))
comp.variance$var.te_fit_total<-as.numeric(as.character(comp.variance$var.te_fit_total))
comp.variance$var.te_residual_total<-as.numeric(as.character(comp.variance$var.te_residual_total))

p<-ggplot(comp.variance, aes(x=dap, y=var.te_fit_total, color=trmt)) + geom_smooth()
q<-p + geom_smooth(aes(x=dap, y=var.te_residual_total, color=trmt))


rm(list=ls())

