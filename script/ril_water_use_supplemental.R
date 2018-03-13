library(ggplot2)

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

##### END ##### 

## Combine biomass with water data
setwd(wue_results.plant_size.dir)
load('ril_biomass.Rdata')
setwd(wue_results.water_use.dir)
load('ril_water_use.Rdata')

## Lets examine how much water was lost during the first two days
## Add some empties as a line across the bottom
p<-ggplot(data=all_day_water_lost, aes(x=dap_i, y=water_lost_day, colour=treatment)) + geom_smooth() 
p

## Merge water and biomass together
all_day_biomass_water<-merge(ril_loess_sv_area, ril_loess_fw_biomass, by=c("genotype", "treatment","dap_i"))
all_day_biomass_water<-merge(all_day_biomass_water, all_day_water_lost, by=c("genotype", "treatment","dap_i"))

## Calculate fw_biomass added per day
all_day_biomass_water$dap_i<-as.numeric(as.character(all_day_biomass_water$dap_i))
days<-sort(unique(all_day_biomass_water$dap_i))
all_day_size<-c()
for (d in 1:length(days)){
  day<-days[d]
  today<-all_day_biomass_water[all_day_biomass_water$dap_i == day,]
  if (d == 1) {
    genotype<-today$genotype
    treatment<-today$treatment
    dap_i<-today$dap_i
    fw_biomass_day<-today$fw_biomass
    sv_area_day<-today$sv_area
    all_day_size<-cbind(genotype, treatment, dap_i, fw_biomass_day, sv_area_day)
  }
  
  if (d > 1) {
    yest<-day-1
    yesterday<-all_day_biomass_water[all_day_biomass_water$dap_i == yest,]
    yesterday<-yesterday[,-c(3)]
    colnames(yesterday)[c(3,5)]<-c('sv_area.y', 'fw_biomass.y')
    today<-merge(today[,c(1:4,6)], yesterday[,c(1:3,5)], by=c("genotype", "treatment"))
    genotype<-today$genotype 
    treatment<-today$treatment
    dap_i<-today$dap_i
    fw_biomass_day<-today$fw_biomass - today$fw_biomass.y
    sv_area_day<-today$sv_area - today$sv_area.y
    temp<-cbind(genotype, treatment, dap_i, fw_biomass_day, sv_area_day)
    all_day_size<-rbind(all_day_size, temp)
  }
}

all_day_size<-as.data.frame(all_day_size)
all_day_size$dap_i<-as.numeric(as.character(all_day_size$dap_i))
all_day_size$sv_area_day<-as.numeric(as.character(all_day_size$sv_area_day))
all_day_size$fw_biomass_day<-as.numeric(as.character(all_day_size$fw_biomass_day))

all_day_biomass_water<-merge(all_day_biomass_water, all_day_size, by=c("genotype", "treatment","dap_i"))

p<-ggplot(data=all_day_biomass_water, aes(x=dap_i, y=water_lost_total, colour=treatment)) + geom_smooth() 
p

## See how much fw_biomass is added per day
p<-ggplot(data=all_day_size, aes(x=dap_i, y=fw_biomass_day)) + geom_point() + facet_wrap(~treatment)
p

###############################################################
## Calculate the proportional loss of water throughout time
## We will use this plot to diagram the entire experiment

p<-ggplot(all_day_biomass_water, aes(x=dap_i, y=fw_biomass)) + geom_smooth() + facet_wrap(~treatment)
p

p<-ggplot(all_day_biomass_water, aes(x=dap_i, y=fw_biomass_day)) + geom_smooth() + facet_wrap(~treatment)
p

# Here is where we calculate value of the trait relative to the trait maximum at each day
ids<-unique(all_day_biomass_water$genotype)
all_day_biomass_water.r<-c()
for (i in ids){
  temp<-all_day_biomass_water[all_day_biomass_water$genotype == i,]
  temp.w<-temp[temp$treatment == 'wet',]
  if(nrow(temp.w) > 0) {
    fw_max<-max(temp.w$fw_biomass)
    temp.w$fw_biomass.r<-temp.w$fw_biomass/fw_max
    fw_day_max<-max(temp.w$fw_biomass_day)
    temp.w$fw_biomass_day.r<-temp.w$fw_biomass_day/fw_day_max
    fw_rate_max<-max(temp.w$fw_biomass_rate)
    temp.w$fw_rate.r<-temp.w$fw_biomass_rate/fw_rate_max
    sv_max<-max(temp.w$sv_area)
    temp.w$sv_area.r<-temp.w$sv_area/sv_max
    sv_day_max<-max(temp.w$sv_area_day)
    temp.w$sv_area_day.r<-temp.w$sv_area_day/sv_day_max
    sv_rate_max<-max(temp.w$sv_area_rate)
    temp.w$sv_rate.r<-temp.w$sv_area_rate/sv_rate_max
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
    fw_day_max<-max(temp.d$fw_biomass_day)
    temp.d$fw_biomass_day.r<-temp.d$fw_biomass_day/fw_day_max
    fw_rate_max<-max(temp.d$fw_biomass_rate)
    temp.d$fw_rate.r<-temp.d$fw_biomass_rate/fw_rate_max
    sv_max<-max(temp.d$sv_area)
    temp.d$sv_area.r<-temp.d$sv_area/sv_max
    sv_day_max<-max(temp.d$sv_area_day)
    temp.d$sv_area_day.r<-temp.d$sv_area_day/sv_day_max
    sv_rate_max<-max(temp.d$sv_area_rate)
    temp.d$sv_rate.r<-temp.d$sv_area_rate/sv_rate_max
    water_lost_total_max<-max(temp.d$water_lost_total)
    temp.d$water_lost_total.r<-temp.d$water_lost_total/water_lost_total_max
    water_lost_day_max<-max(temp.d$water_lost_day)
    temp.d$water_lost_day.r<-temp.d$water_lost_day/water_lost_day_max
  }
  temp<-rbind(temp.d, temp.w)
  all_day_biomass_water.r<-rbind(all_day_biomass_water.r, temp)
}

## Plot total relative area, total and daily water added
p<-ggplot(data=all_day_biomass_water.r, aes(x=dap_i, y=sv_area.r, colour="dark green")) + geom_smooth() + facet_wrap(~treatment) + scale_color_manual(values=c("dark green")) +ylab("Proportion of maximum value")
g<-p + geom_smooth(aes(x=dap_i, y=water_lost_total.r), colour="blue") 
h <- g+ geom_smooth(aes(x=dap_i, y=water_lost_day.r), colour="light blue") 
h + geom_smooth(aes(x=dap_i, y=sv_area_day.r ), colour="green") 

## This plot re-enforces the idea that there is little plant growth before day 15

########################################################################
# Lets use the relationship between plant size and water use to see 
# how good our predictions are at day 17
########################################################################

## Combine sv_area and total water use estimates at days including 15 and above
day15<-merge(ril_loess_sv_area, ril_loess_water_lost_total, by=c("genotype", "treatment", "dap_i"))

## Calculate the % of maximum sv_area for each plant one each day
ids<-unique(day15$genotype)
day15.r<-c()
for (i in ids){
  temp<-day15[day15$genotype == i,]
  temp.w<-temp[temp$treatment == 'wet',]
  if(nrow(temp.w) > 0) {
    sv_max<-max(temp.w$sv_area)
    temp.w$sv_area.r<-temp.w$sv_area/sv_max
    sv_rate_max<-max(temp.w$sv_area_rate)
    temp.w$sv_rate.r<-temp.w$sv_area_rate/sv_rate_max
    water_lost_total_max<-max(temp.w$water_lost_total)
    temp.w$water_lost_total.r<-temp.w$water_lost_total/water_lost_total_max
  }
  # Now dry 
  temp.d<-temp[temp$treatment == 'dry',]
  if (nrow(temp.d) > 0) {
    sv_max<-max(temp.d$sv_area)
    temp.d$sv_area.r<-temp.d$sv_area/sv_max
    sv_rate_max<-max(temp.d$sv_area_rate)
    temp.d$sv_rate.r<-temp.d$sv_area_rate/sv_rate_max
    water_lost_total_max<-max(temp.d$water_lost_total)
    temp.d$water_lost_total.r<-temp.d$water_lost_total/water_lost_total_max
  }
  temp<-rbind(temp.d, temp.w)
  day15.r<-rbind(day15.r, temp)
}

#########################################################################
## Lets see how good our estimate of total water used before day 17 
## approximates water used during linear growth
#########################################################################

## Lets look at relative sv_area.r vs. days
p<-ggplot(data=day15.r, aes(x=dap_i, y=sv_area.r, colour=treatment)) + geom_smooth() + scale_color_manual(values=c("orange","navy")) + theme_bw()
p
## Looks like slope between 15 - 21 dap is less than is observed 22 - 30 dap

## Lets look at relative sv_area.r vs. water_lost_day
p<-ggplot(data=day15.r, aes(x=sv_area.r, y=water_lost_total, colour=treatment)) + geom_smooth() + scale_color_manual(values=c("orange","navy")) + theme_bw()
p

## Lets investigate slope of water use per relative biomass put on between 17-21 days
e.day15.r<-day15.r[day15.r$dap_i < 22,]
p<-ggplot(data=e.day15.r, aes(x=sv_area.r, y=water_lost_total, colour=treatment)) + geom_smooth() + scale_color_manual(values=c("orange","navy")) + theme_bw()
p
## It looks like the slopes between wet and dry conditions may be different

## Lets see how much sv_area.r is achieved at day 17
p<-ggplot(data=e.day15.r, aes(x=dap_i, y=sv_area.r, colour=treatment)) + geom_point() + scale_color_manual(values=c("orange","navy")) + theme_bw()
p + xlim(15,21)
## Lets see how much water_lost_total is achieved at day 17
p<-ggplot(data=e.day15.r, aes(x=dap_i, y=water_lost_total, colour=treatment)) + geom_point() + scale_color_manual(values=c("orange","navy")) + theme_bw()
p + xlim(15,21)

## Lets try modeling each across total and treatment specific growth rates
e.day17.r<-day15.r[day15.r$dap_i <22 & day15.r$dap_i > 16,]
water_loss_per_biomass.lm<-lm(water_lost_total~sv_area.r, data=e.day17.r)
water_loss_per_biomass.d.lm<-lm(water_lost_total~sv_area.r, data=e.day17.r[e.day17.r$treatment == 'dry',])
water_loss_per_biomass.w.lm<-lm(water_lost_total~sv_area.r, data=e.day17.r[e.day17.r$treatment == 'wet',])

new<-data.frame(sv_area.r = e.day15.r$sv_area.r)
predict.lm(water_loss_per_biomass.lm, newdata=new)

## Lets model the values from day 15 and 16 based upon the values found on days 17 - 21
genos<-sort(unique(e.day15.r$genotype))
treatments<-unique(e.day15.r$treatment)
day15.p<-c()
for (g in genos) {
  temp<-e.day17.r[e.day17.r$genotype == g,]
  temp15<-e.day15.r[e.day15.r$genotype == g,]
  for(t in treatments){
    if(nrow(temp[temp$treatment == t,]) < 1) {next}
    t.temp<-temp[temp$treatment == t,]
    t.temp15<-temp15[temp15$treatment == t,]
    water_area.lm<-lm(water_lost_total~sv_area.r, data=t.temp)
    new<-data.frame(sv_area.r = t.temp15$sv_area.r)
    t.temp15$water_lost_total.p<-predict.lm(water_area.lm, newdata=new)
    day15.p<-rbind(day15.p, t.temp15)
  }
}

## Lets plot the values of total water lost starting on day 15 vs. those predicted at days 15 - 17.
p<-ggplot(data=day15.p, aes(x=water_lost_total, y=water_lost_total.p, colour=treatment)) + geom_point() + scale_color_manual(values=c("orange", "navy")) + ylab("Predicted water lost") + xlab("Calculated water lost")
p + theme_bw()

pdf("projected_te_vs_calculated_te.pdf")
p<-ggplot(data=day15.p, aes(x=dap_i, y=water_lost_total, colour=treatment)) + geom_line(size=0.3) + scale_color_manual(values=c("orange", "navy")) + ylab("Water lost (mL)") + xlab("Days after sewing")
q<-p + geom_line(data=day15.p[day15.p$treatment == 'wet',], aes(x=dap_i, y=water_lost_total.p), colour = 'light blue', size=0.3)
u<-q + geom_line(data=day15.p[day15.p$treatment == 'dry',], aes(x=dap_i, y=water_lost_total.p), colour = 'yellow', size=0.3) + theme_bw()
print(u)
dev.off()

save.image('ril_water_use_supplemental.Rdata')
rm(list=ls())

