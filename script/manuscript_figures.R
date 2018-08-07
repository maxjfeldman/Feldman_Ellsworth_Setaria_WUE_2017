## This is a script to make the figures in the manuscript 
## You can run this after all other scripts have been run


##############################################
# Load in dependencies
##############################################
library(lme4)
library(ggplot2)
library(lattice)
library(lmomco)
library(gplots)
library(matrixStats)
library(grid)
library(gridExtra)
library(pvclust)
library(WGCNA)
library(stringr)
library(VennDiagram)
library(qtl)
library(funqtl)

setwd("/Users/mfeldman/Desktop/Feldman_Ellsworth_Setaria_WUE_2017-master/")


## Tester
#setwd("~/Dropbox/Feldman_Ellsworth_Setaria_WUE_2017/")

## Make the directory of the folder you downloaded the current working directory

home.dir<-getwd()
setwd(home.dir)
load('analysis_fxns.Rdata')


wue_figures.dir<-paste(home.dir, '/figures', sep="")

## If this directory doesn't exist the download is not complete
if (file.exists(wue_figures.dir)){
  print("Directory exists!")
} else {
  dir.create(file.path(wue_figures.dir))
}

wue_figures.main.dir<-paste(wue_figures.dir, '/main', sep="")

## If this directory doesn't exist the download is not complete
if (file.exists(wue_figures.main.dir)){
  print("Directory exists!")
} else {
  dir.create(file.path(wue_figures.main.dir))
}

wue_figures.supplemental.dir<-paste(wue_figures.dir, '/supplemental', sep="")

## If this directory doesn't exist the download is not complete
if (file.exists(wue_figures.supplemental.dir)){
  print("Directory exists!")
} else {
  dir.create(file.path(wue_figures.supplemental.dir))
}


## Make the results directory
wue_results.dir<-paste(home.dir, '/results/', sep="")

## If this directory doesn't exist the download is not complete
if (file.exists(wue_results.dir)){
  print("Directory exists!")
} else {
  dir.create(file.path(wue_results.dir))
}
setwd(wue_results.dir)

## Make a directory for plant size results
wue_results.plant_size.dir<-paste(wue_results.dir, 'plant_size/', sep="")

## If this directory doesn't exist the download is not complete
if (file.exists(wue_results.plant_size.dir)){
  print("Directory exists!")
} else {
  dir.create(file.path(wue_results.plant_size.dir))
}


## Make a directory for water use
wue_results.water_use.dir<-paste(wue_results.dir, 'water_use/', sep="")

## Now one for wue/te results
wue_results.te.dir<-paste(wue_results.dir, 'te/', sep="")

## heritability
wue_results.heritability.dir<-paste(wue_results.dir, 'heritability/', sep="")

## qtl tables
wue_results.qtl.dir<-paste(wue_results.dir, 'qtl_summary/', sep="")

## venn diagrams
wue_results.qtl.venn.dir<-paste(wue_results.qtl.dir, 'venn_diagrams/', sep="")

## clustering diagrams
wue_results.qtl.clustering.dir<-paste(wue_results.qtl.dir, '/qtl_clustering', sep="")

## Lets start process Figures in the order they are in the manuscript

setwd(wue_results.plant_size.dir)
load('ril_biomass.Rdata')

## Do biomass calibration
## Biomass model 
calibrate$treatment2<-paste(calibrate$treatment, "2", sep="_")

## Convert from grams to milligrams
calibrate$fresh_wt<-calibrate$fresh_wt*1000
calibrate$fw.min<-calibrate$fw.min*1000

## Make plot of predicted vs actual biomass based upon minimal model 

setwd(wue_figures.main.dir)
pdf("FIG_1a.pdf")
p<-ggplot(data=calibrate, aes(x=fresh_wt,y=fw.min,colour=treatment)) + geom_point(size=0.75) +  theme_bw()
q<-p +  geom_smooth(method = 'lm', aes(colour = treatment2), size=2) + scale_color_manual(values=c("gold", "orange", "blue","navy")) + ylab("Estimated fresh weight (mg)") + xlab("Fresh weight (mg)") + theme(legend.position="none",axis.text=element_text(size=20), axis.title=element_text(size=24,face="bold")) 
x<-q + annotate("text", x = 27000, y = 41000, label = "R-Sq: 0.86", colour="navy", size=7)
y<-x + annotate("text", x = 20000, y = 5000, label = "R-Sq: 0.74", colour="orange", size=7)
print(y)
dev.off()

## Load in water use data
setwd(wue_results.water_use.dir)
load("ril_water_use_supplemental.Rdata")

setwd(wue_figures.main.dir)

pdf("FIG_1b.pdf")

all_day_biomass_water.r$treatment2<-paste(all_day_biomass_water.r$treatment, "2", sep="_")
p<-ggplot(data=all_day_biomass_water.r, aes(x=dap_i, y=sv_area, colour=treatment)) + geom_smooth() + theme_bw()
p1<-p + theme(panel.background = element_blank(), panel.grid.minor = element_blank(),  axis.ticks.y = element_blank(),legend.position="none",axis.text=element_text(size=20),axis.title=element_text(size=24,face="bold"))
p2<-p1 + geom_smooth(aes(x=dap_i, y=sv_area_day, colour=treatment2), size=2) + scale_color_manual(values=c("light green","Thistle", "dark green", "violet")) + xlim(15,33) + ylab("Plant area (pixel)") + xlab("Days after planting")
print(p2)

dev.off()


# Lets plot some convincing plots about why wating should start at day 15 - 17
empty_pots<-c("Dr27AA001310", "Dr41AA001393", "Dr53AB001468", "Dr65AB001541", "Dr72AA001581")
water.m.all[water.m.all$plantbarcode %in% empty_pots, 'genotype']<-c("EMPTY")
empties<-water.m.all[water.m.all$genotype == "EMPTY",]
empties[empties$treatment == 'dry', 'treatment']<-c('dry_empty')
empties[empties$treatment == 'wet', 'treatment']<-c('wet_empty')

water_2<-water.m.all[water.m.all$genotype != "EMPTY", ]
water_2<-rbind(water.m.all, empties)

setwd(wue_figures.main.dir)

pdf("FIG_1c.pdf")
p<-ggplot(data=water_2, aes(x=dap_i, y=water_lost, colour=treatment)) + geom_smooth(size=2) + scale_colour_manual(values = c("orange","gold","navy","blue"))
q<-p + ylab("Water lost per day (mL)") + xlab("Days after planting") + theme_bw() + xlim(8,33)
x<-q + geom_vline(xintercept=15, linetype=2, color='black', size=1)
y<-x + geom_vline(xintercept=17, linetype=2, color='red', size=1) + theme(legend.position="none", axis.text=element_text(size=20),axis.title=element_text(size=24,face="bold"))
print(y)
dev.off()

vis.wue.all<-merge(vis, water.m.all, by=c("plantbarcode", "cartag", "genotype", "treatment", "dap_i"))
# Lets get value of biomass, sv_area on a per day basis
days<-unique(vis.wue.all$dap_i)

ids<-unique(vis.wue.all$plantbarcode)
vis.wue.all.r<-c()
for (i in ids){
  temp<-vis.wue.all[vis.wue.all$plantbarcode == i,]
  fw_max<-max(temp$fw_biomass)
  temp$fw_biomass.r<-temp$fw_biomass/fw_max
  dw_max<-max(temp$dw_biomass)
  temp$dw_biomass.r<-temp$dw_biomass/dw_max
  sv_max<-max(temp$sv_area)
  temp$sv_area.r<-temp$sv_area/sv_max
  tv_max<-max(temp$tv_area)
  temp$tv_area.r<-temp$tv_area/tv_max
  height_max<-max(temp$height_above_bound)
  temp$height.r<-temp$height_above_bound/height_max
  width_max<-max(temp$width)
  temp$width.r<-temp$width/width_max
  water_lost_max<-max(temp$water_lost_total)
  temp$water_lost_total.r<-temp$water_lost_total/water_lost_max
  water_lost_max<-max(temp$water_lost)
  temp$water_lost_day.r<-temp$water_lost/water_lost_max
  vis.wue.all.r<-rbind(vis.wue.all.r, temp)
}


# Plot relative plant size at day 15 and day 17 to support choice of day
setwd(wue_figures.main.dir)

pdf("FIG_1d.pdf")

p<-ggplot(data=vis.wue.all.r, aes(x=dap_i, y=sv_area.r, colour=treatment)) + geom_smooth(size=2) + scale_colour_manual(values = c("orange","navy")) + theme_bw() + theme(legend.position="none",axis.text=element_text(size=20),axis.title=element_text(size=24,face="bold"))
q<-p + xlim(7,33) + ylab("Proportion of total biomass") + xlab("Days after planting")
x<-q + geom_vline(xintercept=15, linetype=2, color='black',size=1)
y<-x + geom_vline(xintercept=17, linetype=2, color='red',size=1)
print(y)

dev.off()

## Figure S1
## Lets look at the distirbution of these values
setwd(wue_figures.supplemental.dir)


hist(mean_ratio_sv_area$ave_ratio)

max_id<-as.character(mean_ratio_sv_area[mean_ratio_sv_area$ave_ratio == max(mean_ratio_sv_area$ave_ratio),"genotype"])
min_id<-as.character(mean_ratio_sv_area[mean_ratio_sv_area$ave_ratio == min(mean_ratio_sv_area$ave_ratio),"genotype"])

## Lets make a plot of dry/wet ratio over the experiment for A10 and B100
parents<-ratio_sv_area_long[ratio_sv_area_long$id == "A10" | ratio_sv_area_long$id == "B100",]
p<-ggplot(data=parents, aes(x=day, y=ratio, colour=id)) + geom_line() + xlim(17,33)

pdf("FIG_S2.pdf")
## Lets make a plot of dry/wet ratio over the experiment for A10, B100, RIL_099 (max), RIL_010 (min)
parents_max_min<-ratio_sv_area_long[ratio_sv_area_long$id == "A10" | ratio_sv_area_long$id == "B100" | ratio_sv_area_long$id == "RIL_099" | ratio_sv_area_long$id == "RIL_010",]
colnames(parents_max_min)[c(7,9,10)]<-c("Genotype", "Ratio", "Day")
p<-ggplot(data=parents_max_min, aes(x=Day, y=Ratio, colour=Genotype)) + geom_line() + xlim(17,33) + theme_bw() + ylab(c("Ratio (Dry/Wet)")) + xlab("Days after planting")
print(p)

dev.off()

## Figure S2
vis.figs<-vis
vis.figs$fw.all.int = predict.lm(object = fw.all.int, newdata=vis.figs)

setwd(wue_figures.supplemental.dir)

pdf("FIG_S27.pdf")
p<-ggplot(vis.figs, aes(y=fw.all.int, x=dap_i, colour=treatment)) + geom_point(size=0.3)
g<-p+scale_color_manual(values=c("orange", "navy")) + ylab("Fresh weight biomass (g)") + xlab("Days after planting") + theme_bw() + theme(legend.position="none", axis.text=element_text(size=20), axis.title=element_text(size=24,face="bold")) 
q<-g + geom_hline(aes(yintercept=0), color='red', size=2, linetype="dashed") 
print(q)

dev.off()


## Make a plot of WUE ratio on day 15 to show artifact
day15$sv_water_total<-day15$sv_area/day15$water_lost_total
day15[day15$sv_water_total < 1,'sv_water_total']<-c(1)

setwd(wue_figures.supplemental.dir)

pdf("FIG_S29a.pdf")

p<-ggplot(day15, aes(x=dap_i, y=sv_water_total, colour=treatment)) + geom_point(size=0.75) + ylim(-100,15000) +xlim(14.5,33) + scale_colour_manual(values = c("orange","navy")) + theme_bw()
q<-p + xlab("Days after planting") + ylab("Biomass (mg) / total_water_added (mL)") + theme(legend.position="none")

x<-q+annotate("rect",  xmin = 14.5, xmax = 16.5, ymin = -50, ymax = 14000, alpha = .2)
print(x)

dev.off()


e.day17.r<-day15.r[day15.r$dap_i <22 & day15.r$dap_i > 17,]
water_loss_per_biomass.lm<-lm(water_lost_total~sv_area.r, data=e.day17.r)
water_loss_per_biomass.d.lm<-lm(water_lost_total~sv_area.r, data=e.day17.r[e.day17.r$treatment == 'dry',])
water_loss_per_biomass.w.lm<-lm(water_lost_total~sv_area.r, data=e.day17.r[e.day17.r$treatment == 'wet',])

new<-data.frame(sv_area.r = e.day15.r$sv_area.r)
predict.lm(water_loss_per_biomass.lm, newdata=new)

# Lets model the values from day 15 and 16 based upon the values found on days 17 - 21
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


setwd(wue_figures.supplemental.dir)

pdf("FIG_S29b.pdf")

p<-ggplot(data=day15.p, aes(x=dap_i, y=water_lost_total, colour=treatment)) + geom_line(size=2) + scale_color_manual(values=c("orange", "navy")) + ylab("Water lost (mL)") + xlab("Days after planting")
q<-p + geom_line(data=day15.p[day15.p$treatment == 'wet',], aes(x=dap_i, y=water_lost_total.p), colour = 'light blue', size=2)
u<-q + geom_line(data=day15.p[day15.p$treatment == 'dry',], aes(x=dap_i, y=water_lost_total.p), colour = 'yellow', size=2) + theme_bw() + theme(legend.position="none",axis.text=element_text(size=20),axis.title=element_text(size=24,face="bold"))
print(u)

dev.off()


calibrate$fw.min.residual<-residuals(fw.min)

# Make residual plot as supplemental figure
setwd(wue_figures.supplemental.dir)

pdf("FIG_S1.pdf")
p<-ggplot(data=calibrate, aes(x=fresh_wt,y=fw.min.residual,colour=treatment)) + geom_point(size=0.75) +  theme_bw() + scale_color_manual(values=c("orange","navy"))
q<-p + geom_hline(aes(yintercept=0), color='red', size=2, linetype="dashed") + ylab("Residual (mg)") + xlab("Fresh weight (mg)") + theme(legend.position="none",axis.text=element_text(size=20), axis.title=element_text(size=24,face="bold")) 
q
dev.off()



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

setwd(wue_figures.supplemental.dir)

pdf('FIG_S3.pdf')
q
dev.off()

### Lets make some plots of the variance
setwd(wue_results.heritability.dir)
load("ril_heritability.Rdata")

days<-sort(unique(te_model$dap_i))
te_summary_stats<-c()
for(d in days){
  temp<-te_model[te_model$dap_i == d,]
  temp<-temp[complete.cases(temp),]
  treats<-sort(unique(temp$treatment))
  temp.t<-c()
  for(t in treats){
    ## Calculate mean and var within treatment blocks
    temp2<-temp[temp$treatment == t,]
    temp2.total<-temp2[,c('sv_area_total', 'water_lost_total', 'te_fit_total','te_residual_total', 'wue_ratio_total')]
    temp2.day<-temp2[,c('sv_area_day', 'water_lost_day', 'te_fit_day','te_residual_day', 'wue_ratio_day')]
    temp2.total.mean<-colMeans(temp2.total)
    temp2.total.var<-colVars(as.matrix(temp2.total))
    temp2.day.mean<-colMeans(temp2.day)
    temp2.day.var<-colVars(as.matrix(temp2.day))
    temp2.trait.names<-names(temp2.total.mean)
    temp2.dap<-rep(d, length(temp2.total.mean))
    temp2.treat<-rep(t, length(temp2.total.mean))
    temp3<-cbind(temp2.trait.names, temp2.treat, temp2.dap, temp2.total.mean, temp2.total.var, temp2.day.mean, temp2.day.var)
    temp.t<-rbind(temp.t, temp3)
  }
  
  ## Calculate mean and variance across blocks
  temp.total<-temp[,c('sv_area_total', 'water_lost_total', 'te_fit_total','te_residual_total', 'wue_ratio_total')]
  temp.day<-temp[,c('sv_area_day', 'water_lost_day', 'te_fit_day','te_residual_day', 'wue_ratio_day')]
  temp.total.mean<-colMeans(temp.total)
  temp.total.var<-colVars(as.matrix(temp.total))
  temp.day.mean<-colMeans(temp.day)
  temp.day.var<-colVars(as.matrix(temp.day))
  temp.trait.names<-names(temp.total.mean)
  temp.dap<-rep(d, length(temp.total.mean))
  temp.treat<-rep('all', length(temp.total.mean))
  temp.all<-cbind(temp.trait.names, temp.treat, temp.dap, temp.total.mean, temp.total.var, temp.day.mean, temp.day.var)
  
  ## Combine all the summary stats together
  te_summary_stats<-rbind(te_summary_stats, temp.t, temp.all)
}

rownames(te_summary_stats)<-c(1:nrow(te_summary_stats))
colnames(te_summary_stats)<-c('trait', 'treatment', 'dap_i', 'total.mean', 'total.var', 'day.mean', 'day.var')
te_summary_stats<-as.data.frame(te_summary_stats)

for(i in 3:ncol(te_summary_stats)){
  te_summary_stats[,i]<-as.numeric(as.character(te_summary_stats[,i]))
}

## Lets plot these summary statistics
## Mean of large traits total
traits<-unique(te_summary_stats$trait)
for(i in traits){
  temp<-te_summary_stats[te_summary_stats$trait == i,]
  plot.name<-paste(i, '_summary_stats.pdf', sep="")
  pdf(plot.name)
  
  p1<-ggplot(temp, aes(x=dap_i, y=total.mean, color=treatment)) + geom_line() + ggtitle(paste(i, 'total.mean', sep="_"))
  p2<-ggplot(temp, aes(x=dap_i, y=total.var, color=treatment)) + geom_line() + ggtitle(paste(i, 'total.var', sep="_"))
  p3<-ggplot(temp, aes(x=dap_i, y=day.mean, color=treatment)) + geom_line() + ggtitle(paste(i, 'day.mean', sep="_"))
  p4<-ggplot(temp, aes(x=dap_i, y=day.var, color=treatment)) + geom_line() + ggtitle(paste(i, 'day.var', sep="_"))
  
  grid.arrange(p1, p2, p3, p4, ncol = 2, top = i)
  
  dev.off()
}


setwd(wue_figures.supplemental.dir)

pdf("FIG_S28.pdf")
traits<-unique(te_summary_stats$trait)
for(i in traits){
  temp<-te_summary_stats[te_summary_stats$trait == i,]
  
  p1<-ggplot(temp, aes(x=dap_i, y=total.mean, color=treatment)) + geom_line() + ggtitle(paste(i, 'total.mean', sep="_"))
  p2<-ggplot(temp, aes(x=dap_i, y=total.var, color=treatment)) + geom_line() + ggtitle(paste(i, 'total.var', sep="_"))
  p3<-ggplot(temp, aes(x=dap_i, y=day.mean, color=treatment)) + geom_line() + ggtitle(paste(i, 'day.mean', sep="_"))
  p4<-ggplot(temp, aes(x=dap_i, y=day.var, color=treatment)) + geom_line() + ggtitle(paste(i, 'day.var', sep="_"))
  
  grid.arrange(p1, p2, p3, p4, ncol = 2, top = i)
  
}

dev.off()


setwd(wue_results.te.dir)
load("ril_transpiration_analysis.Rdata")

ril_loess_trait.r$treatment2<-paste(ril_loess_trait.r$treatment, 2, sep="_")

days<-sort(unique(ril_loess_trait.r$dap_i))
cor.ac<-c()
cor.dc<-c()
cor.wc<-c()
cor.ad<-c()
cor.dd<-c()
cor.wd<-c()
for(d in days){
  temp<-ril_loess_trait.r[ril_loess_trait.r$dap_i == d,]
  temp.d<-temp[temp$treatment == 'dry',]
  temp.w<-temp[temp$treatment == 'wet',]
  ac<-cor(temp$sv_area, temp$water_lost_total)
  cor.ac<-c(cor.ac, ac)
  dc<-cor(temp.d$sv_area, temp.d$water_lost_total)
  cor.dc<-c(cor.dc, dc)
  wc<-cor(temp.w$sv_area, temp.w$water_lost_total)
  cor.wc<-c(cor.wc, wc)
  
  ad<-cor(temp$sv_area_day, temp$water_lost_day)
  cor.ad<-c(cor.ad, ad)
  dd<-cor(temp.d$sv_area_day, temp.d$water_lost_day)
  cor.dd<-c(cor.dd, dd)
  wd<-cor(temp.w$sv_area_day, temp.w$water_lost_day)
  cor.wd<-c(cor.wd, wd)
}

cor.ac<-cbind(days, cor.ac, rep('all', length(days)), rep('total', length(days)))
colnames(cor.ac)<-c("dap_i", "pcc", "treatment", "time")
cor.dc<-cbind(days, cor.dc, rep('dry', length(days)), rep('total', length(days)))
colnames(cor.dc)<-c("dap_i", "pcc", "treatment", "time")
cor.wc<-cbind(days, cor.wc, rep('wet', length(days)), rep('total', length(days)))
colnames(cor.wc)<-c("dap_i", "pcc", "treatment", "time")

cor.ad<-cbind(days, cor.ad, rep('all', length(days)), rep('day', length(days)))
colnames(cor.ac)<-c("dap_i", "pcc", "treatment", "time")
cor.dd<-cbind(days, cor.dd, rep('dry', length(days)), rep('day', length(days)))
colnames(cor.dd)<-c("dap_i", "pcc", "treatment", "time")
cor.wd<-cbind(days, cor.wd, rep('wet', length(days)), rep('day', length(days)))
colnames(cor.wd)<-c("dap_i", "pcc", "treatment", "time")

biomass_v_water<-rbind(cor.ac, cor.dc, cor.wc, cor.ad, cor.dd, cor.wd)
biomass_v_water<-as.data.frame(biomass_v_water)
biomass_v_water$pcc<-as.numeric(as.character(biomass_v_water$pcc))
biomass_v_water$dap_i<-as.numeric(as.character(biomass_v_water$dap_i))


setwd(wue_figures.main.dir)
pdf("FIG_2a.pdf", width=12, height=3)

p<-ggplot(biomass_v_water[biomass_v_water$time == 'total',], aes(x=dap_i, y=pcc, colour=treatment)) + geom_line(size=2)  + theme_bw() + scale_color_manual(values=c("red", "orange", "navy"))
q<-p + ylab("Pearson's r") + xlab("Days after planting") + ggtitle("Cumulative") + theme(legend.position="none",axis.text=element_text(size=20), axis.title=element_text(size=24,face="bold"),plot.title = element_text(size = 24, face = "bold")) + ylim(0,1)
print(q)

dev.off()

setwd(wue_figures.main.dir)
pdf("FIG_2c.pdf", width=12, height=3)

p<-ggplot(biomass_v_water[biomass_v_water$time == 'day',], aes(x=dap_i, y=pcc, colour=treatment)) + geom_line(size=2)  + theme_bw() + scale_color_manual(values=c("red", "orange", "navy"))
q<-p + ylab("Pearson's r") + xlab("Days after planting") + ggtitle("Rate") + theme(legend.position="none",axis.text=element_text(size=20), axis.title=element_text(size=24,face="bold"),plot.title = element_text(size = 24, face = "bold"))
print(q)

dev.off()


#,axis.text=element_text(size=20), axis.title=element_text(size=24,face="bold"),plot.title = element_text(size = 24, face = "bold")

setwd(wue_figures.main.dir)

pdf("FIG_2b.pdf", width=12, height=3)

p<-ggplot(data=ril_loess_trait.r[ril_loess_trait.r$dap_i %in% c(20, 25, 30),], aes(x=water_lost_total, y=fw_biomass, colour=treatment)) + geom_point(size=2) + facet_wrap(~dap_i) + theme_bw()
q<-p +geom_smooth(method = 'lm', aes(colour = treatment2)) + scale_color_manual(values=c("gold", "orange", "blue","navy")) + ylab("Biomass (mg)") + xlab("Water lost (mL)") + theme(legend.position="none",axis.text=element_text(size=20), axis.title=element_text(size=24,face="bold"),plot.title = element_text(size = 24, face = "bold"))
print(q)
dev.off()


setwd(wue_figures.main.dir)

pdf("FIG_2d.pdf", width=12, height=3)

p<-ggplot(data=ril_loess_trait.r[ril_loess_trait.r$dap_i %in% c(20, 25, 30),], aes(x=water_lost_day, y=fw_biomass_day, colour=treatment)) + geom_point(size=2) + facet_wrap(~dap_i) + theme_bw()
q<-p +geom_smooth(method = 'lm', aes(colour = treatment2)) + scale_color_manual(values=c("gold", "orange", "blue","navy")) + ylab("Biomass (mg)") + xlab("Daily water loss (mL)") + theme(legend.position="none",axis.text=element_text(size=20), axis.title=element_text(size=24,face="bold"),plot.title = element_text(size = 24, face = "bold"))
print(q)

dev.off()


setwd(wue_results.te.dir)
load('ril_transpiration_efficiency.Rdata')


day30<-day_size[day_size$dap_i == 30,]

temp2<-c()
for(t in treats){
  temp<-day30[day30$treatment == t,]
  te_total.mdl<-lm(sv_area~water_lost_total, data=temp)
  te_fit_total<-predict(te_total.mdl)
  te_residual_total<-residuals(te_total.mdl)
  te_day.mdl<-lm(sv_area_day~water_lost_day, data=temp)
  te_fit_day<-predict(te_day.mdl)
  te_residual_day<-residuals(te_day.mdl)
  temp<-cbind(temp, te_fit_total,te_residual_total,te_fit_day,te_residual_day)
  temp2<-rbind(temp2,temp)
}


setwd(wue_figures.main.dir)

#pdf("FIG_3.pdf")
pdf("FIG_3_icon.pdf")
temp2$treatment2<-paste(temp2$treatment, 2, sep="_")
p<-ggplot(temp2, aes(x=water_lost_total, y=sv_area, color=treatment)) + geom_point()
q<-p + theme(legend.position="none",axis.text=element_text(size=20), axis.title=element_text(size=24,face="bold"))
x<-q + geom_segment(aes(xend = water_lost_total, yend = te_fit_total), color="red")
y<-x + geom_point(aes(y = sv_area, color=treatment2))  
z<-x + geom_point(aes(y = te_fit_total, color=treatment2)) + scale_color_manual(values=c("gold","orange", "blue", "navy")) 
#a <- z + theme_bw() + xlab("Water lost (mL)") + ylab("Plant area (pixel)") + ggtitle("30 days after planting") + theme(legend.position="none",axis.text=element_text(size=20), axis.title=element_text(size=24,face="bold"),plot.title = element_text(size = 24, face = "bold"))
a <- z + theme_bw() + xlab("Water lost (mL)") + ylab("Plant area (pixel)") +  theme(legend.position="none",axis.text=element_text(size=20), axis.title=element_text(size=24,face="bold"),plot.title = element_text(size = 24, face = "bold"))

print(a)

dev.off()

#setwd(wue_figures.main.dir)

#pdf("FIG_3.inset.pdf")

#p<-ggplot(temp2, aes(x=water_lost_total, y=sv_area, color=treatment)) + geom_point()
#q<-p + theme(legend.position="none",axis.text=element_text(size=20), axis.title=element_text(size=24,face="bold"))
#x<-q + geom_segment(aes(xend = water_lost_total, yend = te_fit_total), color="red")
#y<-x + geom_point(aes(y = sv_area, color=treatment2)) 
#z<-x + geom_point(aes(y = te_fit_total, color=treatment2)) + scale_color_manual(values=c("gold","orange", "blue", "navy")) 
#a <- z + theme_bw() + xlab("Water lost (mL)") + ylab("Plant area (pixel)") + ggtitle("25 days after planting") + theme(legend.position="none",axis.text=element_text(size=20), axis.title=element_text(size=24,face="bold"),plot.title = element_text(size = 24, face = "bold"))
#b<-a + xlim(600, 900) +ylim(50000,250000)
#print(b)

#dev.off()



setwd(wue_figures.supplemental.dir)

## Plot of the correlation between cumulative WUE ratio and sv_area/biomass
pdf("FIG_S4.pdf")
p<-ggplot(te_model, aes(x=sv_area, y=wue_total, colour=treatment)) + geom_point(size=0.3) + facet_wrap(~dap_i) + scale_color_manual(values=c("orange", "navy")) + theme_bw() + xlab("Plant size (px)") + ylab("WUE Ratio")+ theme(legend.position="none",axis.text.x = element_text(angle = 90))
print(p)
dev.off()

setwd(wue_figures.supplemental.dir)

## Plot of the correlation between daily/rate WUE ratio and sv_area/biomass
pdf("FIG_S5.pdf")
p<-ggplot(te_model, aes(x=sv_area_day, y=wue_day, colour=treatment)) + geom_point(size=0.3) + facet_wrap(~dap_i) + scale_color_manual(values=c("orange", "navy")) + theme_bw() + xlab("Plant size (px)") + ylab("WUE Ratio")+ theme(legend.position="none",axis.text.x = element_text(angle = 90))
print(p)
dev.off()

setwd(wue_figures.supplemental.dir)

## Residual plot of cumulative TE 
pdf("FIG_S6.pdf")
p<-ggplot(te_model, aes(x=te_fit_total, y=te_residual_total, colour=treatment)) + geom_point(size=0.3) + facet_wrap(~dap_i) + scale_color_manual(values=c("orange", "navy")) + theme_bw() + xlab("sv_area ~ water_use_total") + ylab("residual")+ theme(legend.position="none",axis.text.x = element_text(angle = 90))
q<-p + geom_hline(yintercept=0, linetype=2, color='red')
print(q)
dev.off()


setwd(wue_figures.supplemental.dir)
# Residual plot of cumulative TE 
pdf("FIG_S7a.pdf")

p<-ggplot(te_model, aes(x=sv_area, y=te_fit_total, colour=treatment)) + geom_point(size=0.3) + facet_wrap(~dap_i) + scale_color_manual(values=c("orange", "navy")) + theme_bw() + xlab("sv_area") + ylab("model_fit_total")+ theme(legend.position="none",axis.text.x = element_text(angle = 90))
print(p)
dev.off()

setwd(wue_figures.supplemental.dir)

pdf("FIG_S7b.pdf")
p<-ggplot(te_model, aes(x=sv_area_day, y=te_fit_day, colour=treatment)) + geom_point(size=0.3) + facet_wrap(~dap_i) + scale_color_manual(values=c("orange", "navy")) + theme_bw() + xlab("sv_area") + ylab("model_fit_day")+ theme(legend.position="none",axis.text.x = element_text(angle = 90))
print(p)
dev.off()

setwd(wue_figures.supplemental.dir)
# Residual plot of cumulative TE 
pdf("FIG_S8a.pdf")
p<-ggplot(te_model, aes(x=sv_area, y=te_residual_total, colour=treatment)) + geom_point(size=0.3) + facet_wrap(~dap_i) + scale_color_manual(values=c("orange", "navy")) + theme_bw() + xlab("sv_area") + ylab("model_residual_total")+ theme(legend.position="none",axis.text.x = element_text(angle = 90))
print(p)
dev.off()

setwd(wue_figures.supplemental.dir)
pdf("FIG_S8b.pdf")
p<-ggplot(te_model, aes(x=sv_area_day, y=te_residual_day, colour=treatment)) + geom_point(size=0.3) + facet_wrap(~dap_i) + scale_color_manual(values=c("orange", "navy")) + theme_bw() + xlab("sv_area") + ylab("model_residual_day")+ theme(legend.position="none",axis.text.x = element_text(angle = 90))
print(p)
dev.off()


### Lets make some plots of Heritability and % variance explained
setwd(wue_results.heritability.dir)
load("ril_heritability.Rdata")

# Plot heritability for each trait on a cumulative and daily basis
setwd(wue_figures.supplemental.dir)
pdf("FIG_S9.pdf", width=8, height=4)
p<-ggplot(gg.comp.h2, aes(x=dap_i, y=value, color=treatment)) + geom_line(size=2)
g<-p+ scale_color_manual(values=c("red", "orange", "navy")) + ylab("Broad sense H2") + xlab("Days after planting") + facet_wrap(~trait, ncol=5)
d<-g +theme_bw() + ylim(0,1) + theme(legend.position="none", axis.text=element_text(size=20), axis.title=element_text(size=24,face="bold"), plot.title = element_text(size = 24, face = "bold"))
print(d)

p<-ggplot(gg.comp.daily.h2, aes(x=dap_i, y=value, color=treatment)) + geom_line(size=2)
g<-p+ scale_color_manual(values=c("red", "orange", "navy")) + ylab("Broad sense H2") + xlab("Days after planting") + facet_wrap(~trait, ncol=5)
d<-g +theme_bw() + ylim(0,1) + theme(legend.position="none", axis.text=element_text(size=20), axis.title=element_text(size=24,face="bold"), plot.title = element_text(size = 24, face = "bold"))
print(d)

dev.off()

## Now plot % of variance explained
setwd(wue_figures.supplemental.dir)

pdf("FIG_S10.pdf", width = 8, height = 4)
# Plot cumulative
p<-ggplot(gg.prop.var, aes(x=dap_i, y=value, color=type)) + geom_line(size=2)
g<-p+ scale_color_manual(values=c("red", "blue", "purple", "black")) + ylab("% of variance") + xlab("Days after planting") + facet_wrap(~trait)
g<-p + ylab("% of variance") + xlab("Days after planting") + facet_wrap(~trait, ncol=5)
#d<-g +theme_bw() + ylim(0,1) + theme(legend.position="none", axis.text=element_text(size=20), axis.title=element_text(size=24,face="bold"), plot.title = element_text(size = 24, face = "bold"))
d<-g +theme_bw() + ylim(0,1) + theme(axis.text=element_text(size=20), axis.title=element_text(size=24,face="bold"), plot.title = element_text(size = 24, face = "bold"))
print(d)

# Plot daily
p<-ggplot(gg.prop.var.daily, aes(x=dap_i, y=value, color=type)) + geom_line(size=2)
g<-p+ scale_color_manual(values=c("red", "blue", "purple", "black")) + ylab("% of variance") + xlab("Days after planting") + facet_wrap(~trait)
g<-p + ylab("% of variance") + xlab("Days after planting") + facet_wrap(~trait, ncol=5)
#d<-g +theme_bw() + ylim(0,1) + theme(legend.position="none", axis.text=element_text(size=20), axis.title=element_text(size=24,face="bold"), plot.title = element_text(size = 24, face = "bold"))
d<-g +theme_bw() + ylim(0,1) + theme(axis.text=element_text(size=20), axis.title=element_text(size=24,face="bold"), plot.title = element_text(size = 24, face = "bold"))
print(d)

dev.off()


setwd(wue_results.qtl.venn.dir)
load("venn_qtl_analysis.Rdata")

setwd(wue_figures.main.dir)
## Lets plot all QTL found for all cumulative traits

pdf("FIG_4.pdf")
p<-ggplot() + geom_point(data = r.all_total.qtl, aes(x = pos, y = prop.var, shape=plot.char, colour=as.character(plot.col), fill=as.character(plot.col)),size=3, alpha=0.5) + geom_blank(data = blank_data, aes(x = x, y = y)) + facet_wrap(~chr, scales = "free_x") + expand_limits(x = 0) + scale_x_continuous(expand = c(0, 0)) + theme_bw() +   scale_shape_manual(values=c(24,25)) 
print(p + scale_color_manual(values=c("1" = "light green", "2" = "dark green", "3" = "light blue", "4" = "blue", "5" = "yellow", "6" = "orange", "7" = "grey", "8" = "black", "9" = "pink", "10" = "red", "11" = "purple")) + scale_fill_manual(values=c("1" = "green", "2" = "dark green", "3" = "light blue", "4" = "navy", "5" = "yellow", "6" = "orange", "7" = "grey", "8" = "black", "9" = "pink", "10" = "red", "11" = "purple")) + ylab("% Variance")  + xlab("Genome Position") + theme(legend.position = "none"))

dev.off()

setwd(wue_figures.supplemental.dir)

## Plot location of all cumulative QTL condensed
pdf("FIG_S12.pdf")
p<-ggplot() + geom_point(data = t.uni.all_qtl, aes(x = pos, y = prop.var, shape=plot.char, colour=as.character(plot.col), fill=as.character(plot.col)),size=3, alpha=0.5) + geom_blank(data = blank_data, aes(x = x, y = y)) + facet_wrap(~chr, scales = "free_x") + expand_limits(x = 0) + scale_x_continuous(expand = c(0, 0)) + theme_bw() +   scale_shape_manual(values=c(24,25)) 
print(p + scale_color_manual(values=c("1" = "light green", "2" = "dark green", "3" = "light blue", "4" = "blue", "5" = "yellow", "6" = "orange", "7" = "grey", "8" = "black", "9" = "pink", "10" = "red")) + scale_fill_manual(values=c("1" = "light green", "2" = "dark green", "3" = "light blue", "4" = "blue", "5" = "yellow", "6" = "orange", "7" = "grey", "8" = "black", "9" = "pink", "10" = "red")) + ylab("% Variance")  + xlab("Genome Position") + theme(legend.position = "none"))

dev.off()

setwd(wue_figures.supplemental.dir)

## Plot location of all daily QTL
pdf("FIG_S11.pdf")
p<-ggplot() + geom_point(data = r.all_day.qtl, aes(x = pos, y = prop.var, shape=plot.char, colour=as.character(plot.col), fill=as.character(plot.col)),size=3, alpha=0.5) + geom_blank(data = blank_data, aes(x = x, y = y)) + facet_wrap(~chr, scales = "free_x") + expand_limits(x = 0) + scale_x_continuous(expand = c(0, 0)) + theme_bw() +   scale_shape_manual(values=c(24,25)) 
print(p + scale_color_manual(values=c("1" = "light green", "2" = "dark green", "3" = "light blue", "4" = "blue", "5" = "yellow", "6" = "orange", "7" = "grey", "8" = "black", "9" = "pink", "10" = "red")) + scale_fill_manual(values=c("1" = "green", "2" = "dark green", "3" = "light blue", "4" = "navy", "5" = "yellow", "6" = "orange", "7" = "grey", "8" = "black", "9" = "pink", "10" = "red")) + ylab("% Variance")  + xlab("Genome Position") + theme(legend.position = "none"))

dev.off()

setwd(wue_figures.supplemental.dir)

## Plot location of all daily QTL condensed
pdf("FIG_S13.pdf")
p<-ggplot() + geom_point(data = d.uni.all_qtl, aes(x = pos, y = prop.var, shape=plot.char, colour=as.character(plot.col), fill=as.character(plot.col)),size=3, alpha=0.5) + geom_blank(data = blank_data, aes(x = x, y = y)) + facet_wrap(~chr, scales = "free_x") + expand_limits(x = 0) + scale_x_continuous(expand = c(0, 0)) + theme_bw() +   scale_shape_manual(values=c(24,25)) 
print(p + scale_color_manual(values=c("1" = "light green", "2" = "dark green", "3" = "light blue", "4" = "blue", "5" = "yellow", "6" = "orange", "7" = "grey", "8" = "black", "9" = "pink", "10" = "red")) + scale_fill_manual(values=c("1" = "light green", "2" = "dark green", "3" = "light blue", "4" = "blue", "5" = "yellow", "6" = "orange", "7" = "grey", "8" = "black", "9" = "pink", "10" = "red")) + ylab("% Variance")  + xlab("Genome Position") + theme(legend.position = "none"))


dev.off()

setwd(wue_figures.supplemental.dir)

pdf("FIG_S14.pdf")

grid.newpage()
draw.quintuple.venn(area1=length(t.sv_area), 
                    area2=length(t.water_lost), 
                    area3=length(t.wue_ratio),
                    area4=length(t.te_fit), 
                    area5=length(t.te_residual),
                    n12=length(intersect(unique(t.sv_area), unique(t.water_lost))),
                    n13=length(intersect(unique(t.sv_area), unique(t.wue_ratio))),
                    n14=length(intersect(unique(t.sv_area), unique(t.te_fit))),
                    n15=length(intersect(unique(t.sv_area), unique(t.te_residual))),
                    n23=length(intersect(unique(t.water_lost), unique(t.wue_ratio))),
                    n24=length(intersect(unique(t.water_lost), unique(t.te_fit))),
                    n25=length(intersect(unique(t.water_lost), unique(t.te_residual))),
                    n34=length(intersect(unique(t.wue_ratio), unique(t.te_fit))),
                    n35=length(intersect(unique(t.wue_ratio), unique(t.te_residual))),
                    n45=length(intersect(unique(t.te_fit), unique(t.te_residual))),
                    n123=length(intersect(unique(t.sv_area), intersect(unique(t.water_lost), unique(t.wue_ratio)))),
                    n124=length(intersect(unique(t.sv_area), intersect(unique(t.water_lost), unique(t.te_fit)))),
                    n125=length(intersect(unique(t.sv_area), intersect(unique(t.water_lost), unique(t.te_residual)))),
                    n134=length(intersect(unique(t.sv_area), intersect(unique(t.wue_ratio), unique(t.te_fit)))),
                    n135=length(intersect(unique(t.sv_area), intersect(unique(t.wue_ratio), unique(t.te_residual)))),
                    n145=length(intersect(unique(t.sv_area), intersect(unique(t.te_fit), unique(t.te_residual)))),
                    n234=length(intersect(unique(t.water_lost), intersect(unique(t.wue_ratio), unique(t.te_fit)))),
                    n235=length(intersect(unique(t.water_lost), intersect(unique(t.wue_ratio), unique(t.te_residual)))),
                    n245=length(intersect(unique(t.water_lost), intersect(unique(t.te_fit), unique(t.te_residual)))),
                    n345=length(intersect(unique(t.wue_ratio), intersect(unique(t.te_fit), unique(t.te_residual)))),
                    n1234=length(intersect(unique(t.sv_area), intersect(unique(t.water_lost), intersect(unique(t.wue_ratio), unique(t.te_fit))))),
                    n1235=length(intersect(unique(t.sv_area), intersect(unique(t.water_lost), intersect(unique(t.wue_ratio), unique(t.te_residual))))),
                    n1245=length(intersect(unique(t.sv_area), intersect(unique(t.water_lost), intersect(unique(t.te_fit), unique(t.te_residual))))),
                    n1345=length(intersect(unique(t.sv_area), intersect(unique(t.wue_ratio), intersect(unique(t.te_fit), unique(t.te_residual))))),
                    n2345=length(intersect(unique(t.water_lost), intersect(unique(t.wue_ratio), intersect(unique(t.te_fit), unique(t.te_residual))))),
                    n12345=length(intersect(unique(t.sv_area), intersect(unique(t.water_lost), intersect(unique(t.wue_ratio), intersect(unique(t.te_fit), unique(t.te_residual)))))),
                    category=c("Biomass", "Water lost", "WUE ratio", "TE fit", "TE residual"),
                    fill=c("green","light blue", "orange", "dark grey", "red"),
                    cex=2,
                    cat.cex=1.5,
                    cat.just=list(c(0,1), c(-.2,0), c(0,0), c(0,0), c(1,0))
)

dev.off()




########################################################################
# Compare cumulative wet vs dry
########################################################################
setwd(wue_figures.supplemental.dir)

grid.newpage()

pdf("FIG_S15.pdf", width=6, height=6)
draw.pairwise.venn(area1=length(t.d.sv_area),
                   area2=length(t.w.sv_area),
                   cross.area=length(intersect(t.d.sv_area, t.w.sv_area)),
                   c("Plant size (Dry)", "Plant size (Wet)"),
                   fill=c("orange", "navy"),
                   cex=2,
                   cat.just=list(c(0,0),c(1,0)),
                   cat.cex=1.5
)

#dev.off()

sv_area.t.dw_both<-intersect(t.d.sv_area, t.w.sv_area)
sv_area.t.d_only<-t.d.sv_area[!t.d.sv_area %in% t.w.sv_area]
sv_area.t.w_only<-t.w.sv_area[!t.w.sv_area %in% t.d.sv_area]


grid.newpage()
#pdf("water_lost_wet.dry_venn.pdf", width=6, height=6)
draw.pairwise.venn(area1=length(t.d.water_lost),
                   area2=length(t.w.water_lost),
                   cross.area=length(intersect(t.d.water_lost, t.w.water_lost)),
                   c("Water lost (Dry)", "Water lost (Wet)"),
                   fill=c("orange", "navy"),
                   cex=2,
                   cat.just=list(c(0,0),c(1,0)),
                   cat.cex=1.5
)

#dev.off()


water_lost.t.dw_both<-intersect(t.d.water_lost, t.w.water_lost)
water_lost.t.d_only<-t.d.water_lost[!t.d.water_lost %in% t.w.water_lost]
water_lost.t.w_only<-t.w.water_lost[!t.w.water_lost %in% t.d.water_lost]



grid.newpage()
#pdf("wue_ratio_wet.dry_venn.pdf", width=6, height=6)
draw.pairwise.venn(area1=length(t.d.wue_ratio),
                   area2=length(t.w.wue_ratio),
                   cross.area=length(intersect(t.d.wue_ratio, t.w.wue_ratio)),
                   c("WUE ratio (Wet)", "WUE ratio (Dry)"),
                   fill=c("orange", "navy"),
                   cex=2,
                   cat.just=list(c(0,0),c(.75,0)),
                   inverted=TRUE,
                   cat.cex=1.5
)

#dev.off()


wue_ratio.t.dw_both<-intersect(t.d.wue_ratio, t.w.wue_ratio)
wue_ratio.t.d_only<-t.d.wue_ratio[!t.d.wue_ratio %in% t.w.wue_ratio]
wue_ratio.t.w_only<-t.w.wue_ratio[!t.w.wue_ratio %in% t.d.wue_ratio]


grid.newpage()
#pdf("te_fit_wet.dry_venn.pdf", width=6, height=6)
draw.pairwise.venn(area1=length(t.d.te_fit),
                   area2=length(t.w.te_fit),
                   cross.area=length(intersect(t.d.te_fit, t.w.te_fit)),
                   c("WUE fit (Wet)", "WUE fit (Dry)"),
                   fill=c("orange", "navy"),
                   cex=2,
                   cat.just=list(c(0,0),c(1,0)),
                   inverted = TRUE,
                   cat.cex=1.5
)

#dev.off()


te_fit.t.dw_both<-intersect(t.d.te_fit, t.w.te_fit)
te_fit.t.d_only<-t.d.te_fit[!t.d.te_fit %in% t.w.te_fit]
te_fit.t.w_only<-t.w.te_fit[!t.w.te_fit %in% t.d.te_fit]


grid.newpage()
#pdf("te_residual_wet.dry_venn.pdf", width=6, height=6)
draw.pairwise.venn(area1=length(t.d.te_residual),
                   area2=length(t.w.te_residual),
                   cross.area=length(intersect(t.d.te_residual, t.w.te_residual)),
                   c("WUE residual (Dry)", "WUE residual (Wet)"),
                   fill=c("orange", "navy"),
                   cex=2,
                   cat.just=list(c(0,0),c(1,0)),
                   cat.cex=1.5
)

dev.off()


te_residual.t.dw_both<-intersect(t.d.te_residual, t.w.te_residual)
te_residual.t.d_only<-t.d.te_residual[!t.d.te_residual %in% t.w.te_residual]
te_residual.t.w_only<-t.w.te_residual[!t.w.te_residual %in% t.d.te_residual]



########################################################################
## Compare cumulative (total) vs day (rate)
########################################################################

setwd(wue_figures.supplemental.dir)

grid.newpage()
pdf("FIG_S16.pdf")
draw.pairwise.venn(area1=length(t.sv_area),
                   area2=length(d.sv_area),
                   cross.area=length(intersect(t.sv_area, d.sv_area)),
                   c("Biomass total", "Biomass rate"),
                   fill=c("Dark green", "green"),
                   cex=2,
                   cat.cex=1.5
)

#dev.off()

sv_area.td_both<-intersect(t.sv_area, d.sv_area)
sv_area.t_only<-t.sv_area[!t.sv_area %in% d.sv_area]
sv_area.d_only<-d.sv_area[!d.sv_area %in% t.sv_area]


grid.newpage()
#pdf("Figure_X.pdf")
draw.pairwise.venn(area1=length(t.water_lost),
                   area2=length(d.water_lost),
                   cross.area=length(intersect(t.water_lost, d.water_lost)),
                   c("Water lost total", "Water lost rate"),
                   fill=c("Dark blue", "blue"),
                   cex=2,
                   cat.cex=1.5
)

#dev.off()


water_lost.td_both<-intersect(t.water_lost, d.water_lost)
water_lost.t_only<-t.water_lost[!t.water_lost %in% d.water_lost]
water_lost.d_only<-d.water_lost[!d.water_lost %in% t.water_lost]


grid.newpage()
#pdf("Figure_X.pdf")
draw.pairwise.venn(area1=length(t.wue_ratio),
                   area2=length(d.wue_ratio),
                   cross.area=length(intersect(t.wue_ratio, d.wue_ratio)),
                   c("WUE ratio total", "WUE ratio rate"),
                   fill=c("orange", "gold"),
                   cex=2,
                   cat.cex=1.5
)

#dev.off()

wue_ratio.td_both<-intersect(t.wue_ratio, d.wue_ratio)
wue_ratio.t_only<-t.wue_ratio[!t.wue_ratio %in% d.wue_ratio]
wue_ratio.d_only<-d.wue_ratio[!d.wue_ratio %in% t.wue_ratio]


grid.newpage()
#pdf("Figure_X.pdf")
draw.pairwise.venn(area1=length(t.te_fit),
                   area2=length(d.te_fit),
                   cross.area=length(intersect(t.te_fit, d.te_fit)),
                   c("WUE fit total", "WUE fit rate"),
                   fill=c("black", "grey"),
                   cex=2,
                   cat.cex=1.5
)

#dev.off()

te_fit.td_both<-intersect(t.te_fit, d.te_fit)
te_fit.t_only<-t.te_fit[!t.te_fit %in% d.te_fit]
te_fit.d_only<-d.te_fit[!d.te_fit %in% t.te_fit]


grid.newpage()
#pdf("Figure_X.pdf")
draw.pairwise.venn(area1=length(t.te_residual),
                   area2=length(d.te_residual),
                   cross.area=length(intersect(t.te_residual, d.te_residual)),
                   c("WUE residual total", "WUE residual rate"),
                   fill=c("red", "pink"),
                   cex=2,
                   cat.cex=1.5
)

dev.off()


te_residual.td_both<-intersect(t.te_residual, d.te_residual)
te_residual.t_only<-t.te_residual[!t.te_residual %in% d.te_residual]
te_residual.d_only<-d.te_residual[!d.te_residual %in% t.te_residual]




setwd(wue_results.qtl.venn.dir)
load("ril_venn_gxe_qtl_analysis.Rdata")

setwd(wue_figures.supplemental.dir)

## Different types of GxE QTL
pdf("FIG_S17.pdf")
draw.triple.venn(area1=length(comp_diff_snp),
                 area2=length(comp_rel_diff_snp),
                 area3=length(comp_ratio_snp),
                 n12=length(intersect(comp_diff_snp, comp_rel_diff_snp)),
                 n23=length(intersect(comp_rel_diff_snp, comp_ratio_snp)),
                 n13=length(intersect(comp_diff_snp, comp_ratio_snp)),
                 n123=length(intersect(comp_diff_snp, intersect(comp_rel_diff_snp, comp_ratio_snp))),
                 c("Diff", "Relative Diff", "Ratio"),
                 fill=c("red", "blue", "yellow")
)

dev.off()



## venn_cumulative_v_gxe_qtl
setwd(wue_figures.supplemental.dir)

pdf("FIG_S18.pdf")

grid.newpage()
#pdf("Figure_X.pdf")
draw.pairwise.venn(area1=length(c.total.sv_area),
                   area2=length(r.total.sv_area),
                   cross.area=length(intersect(c.total.sv_area, r.total.sv_area)),
                   c("Difference", "Raw"),
                   fill=c("red", "blue")
)

#dev.off()


total.sv_area.raw.comp_both<-intersect(r.total.sv_area, c.total.sv_area)
length(total.sv_area.raw.comp_both)
total.sv_area.raw.comp_both

total.sv_area.raw_only<-r.total.sv_area[!r.total.sv_area %in% c.total.sv_area]
length(total.sv_area.raw_only)
total.sv_area.raw_only


total.sv_area.comp_only<-c.total.sv_area[!c.total.sv_area %in% r.total.sv_area]
length(total.sv_area.comp_only)
total.sv_area.comp_only

length(total.sv_area.raw.comp_both)/length(union(total.sv_area.comp_only, union(total.sv_area.raw_only, total.sv_area.raw.comp_both)))
## 55 %

##############################
## water_lost
##############################

grid.newpage()
#pdf("Figure_X.pdf")
draw.pairwise.venn(area1=length(c.total.water_lost),
                   area2=length(r.total.water_lost),
                   cross.area=length(intersect(c.total.water_lost, r.total.water_lost)),
                   c("Difference", "Raw"),
                   fill=c("red", "blue")
)

#dev.off()


total.water_lost.raw.comp_both<-intersect(r.total.water_lost, c.total.water_lost)
length(total.water_lost.raw.comp_both)
total.water_lost.raw.comp_both

total.water_lost.raw_only<-r.total.water_lost[!r.total.water_lost %in% c.total.water_lost]
length(total.water_lost.raw_only)
total.water_lost.raw_only


total.water_lost.comp_only<-c.total.water_lost[!c.total.water_lost %in% r.total.water_lost]
length(total.water_lost.comp_only)
total.water_lost.comp_only

length(total.water_lost.raw.comp_both)/length(union(total.water_lost.comp_only, union(total.water_lost.raw_only, total.water_lost.raw.comp_both)))
## 80%

##############################
## wue_ratio
##############################

grid.newpage()
#pdf("Figure_X.pdf")
draw.pairwise.venn(area1=length(c.total.wue_ratio),
                   area2=length(r.total.wue_ratio),
                   cross.area=length(intersect(c.total.wue_ratio, r.total.wue_ratio)),
                   c("Difference", "Raw"),
                   fill=c("red", "blue")
)

#dev.off()


total.wue_ratio.raw.comp_both<-intersect(r.total.wue_ratio, c.total.wue_ratio)
length(total.wue_ratio.raw.comp_both)
total.wue_ratio.raw.comp_both

total.wue_ratio.raw_only<-r.total.wue_ratio[!r.total.wue_ratio %in% c.total.wue_ratio]
length(total.wue_ratio.raw_only)
total.wue_ratio.raw_only


total.wue_ratio.comp_only<-c.total.wue_ratio[!c.total.wue_ratio %in% r.total.wue_ratio]
length(total.wue_ratio.comp_only)
total.wue_ratio.comp_only

length(total.wue_ratio.raw.comp_both)/length(union(total.wue_ratio.comp_only, union(total.wue_ratio.raw_only, total.wue_ratio.raw.comp_both)))
## 13 %

##############################
## te_fit
##############################

grid.newpage()
#pdf("Figure_X.pdf")
draw.pairwise.venn(area1=length(c.total.te_fit),
                   area2=length(r.total.te_fit),
                   cross.area=length(intersect(c.total.te_fit, r.total.te_fit)),
                   c("Difference", "Raw"),
                   fill=c("red", "blue")
)

#dev.off()


total.te_fit.raw.comp_both<-intersect(r.total.te_fit, c.total.te_fit)
length(total.te_fit.raw.comp_both)
total.te_fit.raw.comp_both

total.te_fit.raw_only<-r.total.te_fit[!r.total.te_fit %in% c.total.te_fit]
length(total.te_fit.raw_only)
total.te_fit.raw_only


total.te_fit.comp_only<-c.total.te_fit[!c.total.te_fit %in% r.total.te_fit]
length(total.te_fit.comp_only)
total.te_fit.comp_only

length(total.te_fit.raw.comp_both)/length(union(total.te_fit.comp_only, union(total.te_fit.raw_only, total.te_fit.raw.comp_both)))

## 91 %


##############################
## te_residual
##############################

grid.newpage()
#pdf("Figure_X.pdf")
draw.pairwise.venn(area1=length(c.total.te_residual),
                   area2=length(r.total.te_residual),
                   cross.area=length(intersect(c.total.te_residual, r.total.te_residual)),
                   c("Difference", "Raw"),
                   fill=c("red", "blue")
)

#dev.off()

dev.off()

total.te_residual.raw.comp_both<-intersect(r.total.te_residual, c.total.te_residual)
length(total.te_residual.raw.comp_both)
total.te_residual.raw.comp_both

total.te_residual.raw_only<-r.total.te_residual[!r.total.te_residual %in% c.total.te_residual]
length(total.te_residual.raw_only)
total.te_residual.raw_only


total.te_residual.comp_only<-c.total.te_residual[!c.total.te_residual %in% r.total.te_residual]
length(total.te_residual.comp_only)
total.te_residual.comp_only

length(total.te_residual.raw.comp_both)/length(union(total.te_residual.comp_only, union(total.te_residual.raw_only, total.te_residual.raw.comp_both)))
## 22 %


###########################################################################

### QTL shared between diff and cumulative

setwd(wue_figures.supplemental.dir)

## venn_trait_overlap_of_qtl_found_in_both_cumulative_and_gxe

pdf("FIG_S19.pdf")

grid.newpage()
draw.quintuple.venn(area1=length(total.sv_area.raw.comp_both), 
                    area2=length(total.water_lost.raw.comp_both), 
                    area3=length(total.wue_ratio.raw.comp_both),
                    area4=length(total.te_fit.raw.comp_both), 
                    area5=length(total.te_residual.raw.comp_both),
                    n12=length(intersect(unique(total.sv_area.raw.comp_both), unique(total.water_lost.raw.comp_both))),
                    n13=length(intersect(unique(total.sv_area.raw.comp_both), unique(total.wue_ratio.raw.comp_both))),
                    n14=length(intersect(unique(total.sv_area.raw.comp_both), unique(total.te_fit.raw.comp_both))),
                    n15=length(intersect(unique(total.sv_area.raw.comp_both), unique(total.te_residual.raw.comp_both))),
                    n23=length(intersect(unique(total.water_lost.raw.comp_both), unique(total.wue_ratio.raw.comp_both))),
                    n24=length(intersect(unique(total.water_lost.raw.comp_both), unique(total.te_fit.raw.comp_both))),
                    n25=length(intersect(unique(total.water_lost.raw.comp_both), unique(total.te_residual.raw.comp_both))),
                    n34=length(intersect(unique(total.wue_ratio.raw.comp_both), unique(total.te_fit.raw.comp_both))),
                    n35=length(intersect(unique(total.wue_ratio.raw.comp_both), unique(total.te_residual.raw.comp_both))),
                    n45=length(intersect(unique(total.te_fit.raw.comp_both), unique(total.te_residual.raw.comp_both))),
                    n123=length(intersect(unique(total.sv_area.raw.comp_both), intersect(unique(total.water_lost.raw.comp_both), unique(total.wue_ratio.raw.comp_both)))),
                    n124=length(intersect(unique(total.sv_area.raw.comp_both), intersect(unique(total.water_lost.raw.comp_both), unique(total.te_fit.raw.comp_both)))),
                    n125=length(intersect(unique(total.sv_area.raw.comp_both), intersect(unique(total.water_lost.raw.comp_both), unique(total.te_residual.raw.comp_both)))),
                    n134=length(intersect(unique(total.sv_area.raw.comp_both), intersect(unique(total.wue_ratio.raw.comp_both), unique(total.te_fit.raw.comp_both)))),
                    n135=length(intersect(unique(total.sv_area.raw.comp_both), intersect(unique(total.wue_ratio.raw.comp_both), unique(total.te_residual.raw.comp_both)))),
                    n145=length(intersect(unique(total.sv_area.raw.comp_both), intersect(unique(total.te_fit.raw.comp_both), unique(total.te_residual.raw.comp_both)))),
                    n234=length(intersect(unique(total.water_lost.raw.comp_both), intersect(unique(total.wue_ratio.raw.comp_both), unique(total.te_fit.raw.comp_both)))),
                    n235=length(intersect(unique(total.water_lost.raw.comp_both), intersect(unique(total.wue_ratio.raw.comp_both), unique(total.te_residual.raw.comp_both)))),
                    n245=length(intersect(unique(total.water_lost.raw.comp_both), intersect(unique(total.te_fit.raw.comp_both), unique(total.te_residual.raw.comp_both)))),
                    n345=length(intersect(unique(total.wue_ratio.raw.comp_both), intersect(unique(total.te_fit.raw.comp_both), unique(total.te_residual.raw.comp_both)))),
                    n1234=length(intersect(unique(total.sv_area.raw.comp_both), intersect(unique(total.water_lost.raw.comp_both), intersect(unique(total.wue_ratio.raw.comp_both), unique(total.te_fit.raw.comp_both))))),
                    n1235=length(intersect(unique(total.sv_area.raw.comp_both), intersect(unique(total.water_lost.raw.comp_both), intersect(unique(total.wue_ratio.raw.comp_both), unique(total.te_residual.raw.comp_both))))),
                    n1245=length(intersect(unique(total.sv_area.raw.comp_both), intersect(unique(total.water_lost.raw.comp_both), intersect(unique(total.te_fit.raw.comp_both), unique(total.te_residual.raw.comp_both))))),
                    n1345=length(intersect(unique(total.sv_area.raw.comp_both), intersect(unique(total.wue_ratio.raw.comp_both), intersect(unique(total.te_fit.raw.comp_both), unique(total.te_residual.raw.comp_both))))),
                    n2345=length(intersect(unique(total.water_lost.raw.comp_both), intersect(unique(total.wue_ratio.raw.comp_both), intersect(unique(total.te_fit.raw.comp_both), unique(total.te_residual.raw.comp_both))))),
                    n12345=length(intersect(unique(total.sv_area.raw.comp_both), intersect(unique(total.water_lost.raw.comp_both), intersect(unique(total.wue_ratio.raw.comp_both), intersect(unique(total.te_fit.raw.comp_both), unique(total.te_residual.raw.comp_both)))))),
                    category=c("Biomass", "Water lost", "WUE ratio", "TE fit", "TE residual"),
                    fill=c("green","light blue", "orange", "black", "red")
)

dev.off()


## Difference traits only

## venn_trait_overlap_of_qtl_found_only_in_gxe
setwd(wue_figures.supplemental.dir)

pdf("FIG_S20.pdf")
grid.newpage()
draw.quintuple.venn(area1=length(total.sv_area.comp_only), 
                    area2=length(total.water_lost.comp_only), 
                    area3=length(total.wue_ratio.comp_only),
                    area4=length(total.te_fit.comp_only), 
                    area5=length(total.te_residual.comp_only),
                    n12=length(intersect(unique(total.sv_area.comp_only), unique(total.water_lost.comp_only))),
                    n13=length(intersect(unique(total.sv_area.comp_only), unique(total.wue_ratio.comp_only))),
                    n14=length(intersect(unique(total.sv_area.comp_only), unique(total.te_fit.comp_only))),
                    n15=length(intersect(unique(total.sv_area.comp_only), unique(total.te_residual.comp_only))),
                    n23=length(intersect(unique(total.water_lost.comp_only), unique(total.wue_ratio.comp_only))),
                    n24=length(intersect(unique(total.water_lost.comp_only), unique(total.te_fit.comp_only))),
                    n25=length(intersect(unique(total.water_lost.comp_only), unique(total.te_residual.comp_only))),
                    n34=length(intersect(unique(total.wue_ratio.comp_only), unique(total.te_fit.comp_only))),
                    n35=length(intersect(unique(total.wue_ratio.comp_only), unique(total.te_residual.comp_only))),
                    n45=length(intersect(unique(total.te_fit.comp_only), unique(total.te_residual.comp_only))),
                    n123=length(intersect(unique(total.sv_area.comp_only), intersect(unique(total.water_lost.comp_only), unique(total.wue_ratio.comp_only)))),
                    n124=length(intersect(unique(total.sv_area.comp_only), intersect(unique(total.water_lost.comp_only), unique(total.te_fit.comp_only)))),
                    n125=length(intersect(unique(total.sv_area.comp_only), intersect(unique(total.water_lost.comp_only), unique(total.te_residual.comp_only)))),
                    n134=length(intersect(unique(total.sv_area.comp_only), intersect(unique(total.wue_ratio.comp_only), unique(total.te_fit.comp_only)))),
                    n135=length(intersect(unique(total.sv_area.comp_only), intersect(unique(total.wue_ratio.comp_only), unique(total.te_residual.comp_only)))),
                    n145=length(intersect(unique(total.sv_area.comp_only), intersect(unique(total.te_fit.comp_only), unique(total.te_residual.comp_only)))),
                    n234=length(intersect(unique(total.water_lost.comp_only), intersect(unique(total.wue_ratio.comp_only), unique(total.te_fit.comp_only)))),
                    n235=length(intersect(unique(total.water_lost.comp_only), intersect(unique(total.wue_ratio.comp_only), unique(total.te_residual.comp_only)))),
                    n245=length(intersect(unique(total.water_lost.comp_only), intersect(unique(total.te_fit.comp_only), unique(total.te_residual.comp_only)))),
                    n345=length(intersect(unique(total.wue_ratio.comp_only), intersect(unique(total.te_fit.comp_only), unique(total.te_residual.comp_only)))),
                    n1234=length(intersect(unique(total.sv_area.comp_only), intersect(unique(total.water_lost.comp_only), intersect(unique(total.wue_ratio.comp_only), unique(total.te_fit.comp_only))))),
                    n1235=length(intersect(unique(total.sv_area.comp_only), intersect(unique(total.water_lost.comp_only), intersect(unique(total.wue_ratio.comp_only), unique(total.te_residual.comp_only))))),
                    n1245=length(intersect(unique(total.sv_area.comp_only), intersect(unique(total.water_lost.comp_only), intersect(unique(total.te_fit.comp_only), unique(total.te_residual.comp_only))))),
                    n1345=length(intersect(unique(total.sv_area.comp_only), intersect(unique(total.wue_ratio.comp_only), intersect(unique(total.te_fit.comp_only), unique(total.te_residual.comp_only))))),
                    n2345=length(intersect(unique(total.water_lost.comp_only), intersect(unique(total.wue_ratio.comp_only), intersect(unique(total.te_fit.comp_only), unique(total.te_residual.comp_only))))),
                    n12345=length(intersect(unique(total.sv_area.comp_only), intersect(unique(total.water_lost.comp_only), intersect(unique(total.wue_ratio.comp_only), intersect(unique(total.te_fit.comp_only), unique(total.te_residual.comp_only)))))),
                    category=c("Biomass", "Water lost", "WUE ratio", "TE fit", "TE residual"),
                    fill=c("green","light blue", "orange", "black", "red")
)


dev.off()





## Lets make the QTL matrix plots
setwd(wue_results.qtl.dir)
load("qtl_trait_matrix.Rdata")

setwd(wue_figures.main.dir)

cols <- c(
  '0' = "white",
  '1' = "grey80",
  '2' = "black",
  '3' = "pink",
  '4' = "red",
  '5' = "green",
  '6' = "dark green",
  '7' = "blue",
  '8' = "navy",
  '9' = "yellow",
  '10' = "orange"
)

## First cumulative QTL
pdf("FIG_5a.pdf")
heatmap(as.matrix(apply(trait.qtl, 2, as.numeric)), col=cols, scale='none', Rowv = NA, Colv = NA, labRow=row.names(trait.qtl), main="Cumulative trait QTL")
dev.off()

setwd(wue_figures.main.dir)

## Now rate QTL
pdf("FIG_5b.pdf")
heatmap(as.matrix(apply(trait.qtl_rate, 2, as.numeric)), col=cols, scale='none', Rowv = NA, Colv = NA, labRow=row.names(trait.qtl_rate), main="Rate trait QTL")
dev.off()


## Assign the values in the data.frame colors
cols <- c(
  '0' = "white",
  '1' = "black",
  '2' = "red",
  '3' = "dark green",
  '4' = "navy",
  '5' = "orange"
)


setwd(wue_figures.main.dir)

## Plot the locations of the QTL
pdf("FIG_5c.pdf")
heatmap(as.matrix(apply(trait.qtl_diff, 2, as.numeric)), col=cols, scale='none', Rowv = NA, Colv = NA, labRow=row.names(trait.qtl_diff), main="Difference trait QTL")
dev.off()



setwd(wue_results.qtl.dir)
load("qtl_timeseries_plots_rate.Rdata")



## Lets make Figure 6 slod_qtl_cumulative_traits_scanone

setwd(wue_figures.main.dir)

pdf("FIG_6.pdf", width=8, height = 3)
par(mfrow = c(1, 5))  # 1 row and 5 columns

plot(sv_area.dry_out.F, sv_area.wet_out.F, col=c("green", "dark green"), main="Plant size", ylab=c("LOD score"))
abline(h=sv_area.dry_max.perm.F_slod[1], col=c("orange"), lty=2)
abline(h=sv_area.wet_max.perm.F_slod[1], col=c("red"), lty=2)

plot(water_lost.dry_out.F, water_lost.wet_out.F, col=c("blue", "navy"), main="Water lost", ylab=c(""))
abline(h=water_lost.dry_max.perm.F_slod[1], col=c("orange"), lty=2)
abline(h=water_lost.wet_max.perm.F_slod[1], col=c("red"), lty=2)

plot(wue_ratio.dry_out.F, wue_ratio.wet_out.F, col=c("gold", "dark orange"), main="WUE ratio", ylab=c(""))
abline(h=wue_ratio.dry_max.perm.F_slod[1], col=c("orange"), lty=2)
abline(h=wue_ratio.wet_max.perm.F_slod[1], col=c("red"), lty=2)

#plot(te_fit.dry_out.F, te_fit.wet_out.F, col=c("grey80", "black"), main="TE model fit", ylab=c(""))
pdf("WUEfit.pdf")
#plot(te_fit.dry_out.F, te_fit.wet_out.F, col=c("grey80", "black"), ylab=c(""))
plot(te_fit.dry_out.F, te_fit.wet_out.F, col=c("grey80", "black"), ylab=c(""), xlab=c(""))

#abline(h=te_fit.dry_max.perm.F_slod[1], col=c("orange"), lty=2)
abline(h=te_fit.wet_max.perm.F_slod[1], col=c("red"), lty=2)
dev.off()

#plot(te_residual.dry_out.F, te_residual.wet_out.F, col=c("pink", "red"), main="TE model residual", ylab=c(""))
pdf("WUEresidual.pdf")
#plot(te_residual.dry_out.F, te_residual.wet_out.F, col=c("pink", "red"), ylab=c(""))
plot(te_residual.dry_out.F, te_residual.wet_out.F, col=c("pink", "red"), ylab=c(""), xlab=c(""))

#abline(h=te_residual.dry_max.perm.F_slod[1], col=c("orange"), lty=2)
abline(h=te_residual.wet_max.perm.F_slod[1], col=c("red"), lty=2)
dev.off()

dev.off()


## Now plot SLOD for all trait QTL 

setwd(wue_figures.supplemental.dir)

pdf("FIG_S21.pdf", width=8, height = 3)
par(mfrow = c(1, 5))  # 1 row and 5 columns

plotLodProfile(sv_area.dry_slod, col=c("green"), main="Plant size", showallchr=T)
plotLodProfile(sv_area.wet_slod, col=c("dark green"), add=T, showallchr=T)

plotLodProfile(water_lost.dry_slod, col=c("blue"), main="Water lost", showallchr=T, ylab="")
plotLodProfile(water_lost.wet_slod, col=c("navy"), add=T, showallchr=T, ylab="")

plotLodProfile(wue_ratio.wet_slod, col=c("orange"), main="WUE ratio", showallchr=T, ylab="")
plotLodProfile(wue_ratio.dry_slod, col=c("gold"), add=T, showallchr=T, ylab="")

plotLodProfile(te_fit.wet_slod, col=c("black"), main="TE model fit", showallchr=T, ylab="")
plotLodProfile(te_fit.dry_slod, col=c("grey80"), add=T, showallchr=T, ylab="")

plotLodProfile(te_residual.wet_slod, col=c("red"), main="TE model residual", showallchr=T, ylab="")
plotLodProfile(te_residual.dry_slod, col=c("pink"), add=T, showallchr=T, ylab="")

dev.off()


## Lets make Figure S25
## Rate traits using stepwise SLOD
setwd(wue_figures.supplemental.dir)

pdf("FIG_S22.pdf", width=8, height = 3)
par(mfrow = c(1, 5))  # 1 row and 5 columns

plotLodProfile(sv_area.wet_slod, col=c("dark green"), main="Plant size", showallchr=T)
plotLodProfile(sv_area.dry_slod, col=c("green"), add=T, showallchr=T)

plotLodProfile(water_lost.dry_slod, col=c("blue"), main="Water lost", showallchr=T, ylab="")
plotLodProfile(water_lost.wet_slod, col=c("navy"), add=T, showallchr=T, ylab="")

plotLodProfile(wue_ratio.wet_slod, col=c("orange"), main="WUE ratio", showallchr=T, ylab="")
plotLodProfile(wue_ratio.dry_slod, col=c("gold"), add=T, showallchr=T, ylab="")

plotLodProfile(te_fit.wet_slod, col=c("black"), main="TE model fit", showallchr=T, ylab="")
plotLodProfile(te_fit.dry_slod, col=c("grey80"), add=T, showallchr=T, ylab="")

plotLodProfile(te_residual.dry_slod, col=c("pink"), main="TE model residual", showallchr=T, ylab="")
plotLodProfile(te_residual.wet_slod, col=c("red"), add=T, showallchr=T, ylab="")

dev.off()

## Lets make Figure S26
## Rate traits using scanone SLOD
setwd(wue_figures.supplemental.dir)

pdf("FIG_S23.pdf", width=8, height = 3)
par(mfrow = c(1, 5))  # 1 row and 5 columns

plot(sv_area.dry_out.F, sv_area.wet_out.F, col=c("green", "dark green"), main="Plant size", ylab=c("LOD score"))
abline(h=sv_area.dry_max.perm.F_slod[1], col=c("orange"), lty=2)
abline(h=sv_area.wet_max.perm.F_slod[1], col=c("red"), lty=2)

plot(water_lost.dry_out.F, water_lost.wet_out.F, col=c("blue", "navy"), main="Water lost", ylab=c(""))
abline(h=water_lost.dry_max.perm.F_slod[1], col=c("orange"), lty=2)
abline(h=water_lost.wet_max.perm.F_slod[1], col=c("red"), lty=2)

plot(wue_ratio.dry_out.F, wue_ratio.wet_out.F, col=c("gold", "dark orange"), main="WUE ratio", ylab=c(""))
abline(h=wue_ratio.dry_max.perm.F_slod[1], col=c("orange"), lty=2)
abline(h=wue_ratio.wet_max.perm.F_slod[1], col=c("red"), lty=2)

plot(te_fit.dry_out.F, te_fit.wet_out.F, col=c("grey80", "black"), main="TE model fit", ylab=c(""))
abline(h=te_fit.dry_max.perm.F_slod[1], col=c("orange"), lty=2)
abline(h=te_fit.wet_max.perm.F_slod[1], col=c("red"), lty=2)

plot(te_residual.dry_out.F, te_residual.wet_out.F, col=c("pink", "red"), main="TE model residual", ylab=c(""))
abline(h=te_residual.dry_max.perm.F_slod[1], col=c("orange"), lty=2)
abline(h=te_residual.wet_max.perm.F_slod[1], col=c("red"), lty=2)

dev.off()


setwd(wue_results.qtl.clustering.dir)
load("qtl_fixed_fx_size_clustering.Rdata")


traits<-unique(r.all_total.qtl$trait.treatment)
r.all_total.qtl.clust<-c()
for(t in 1:length(traits)){
  t.name<-traits[t]
  if (t == 1) {
    temp<-r.all_total.qtl[r.all_total.qtl$trait.treatment == t.name,]
    temp2<-temp[,c("qtl.name", "signed.prop.var")]
    colnames(temp2)[2]<-t.name
    r.all_total.qtl.clust<-temp2
  }
  
  if (t > 1) {
    temp<-r.all_total.qtl[r.all_total.qtl$trait.treatment == t.name,]
    temp2<-temp[,c("qtl.name", "signed.prop.var")]
    colnames(temp2)[2]<-t.name
    r.all_total.qtl.clust<-merge(r.all_total.qtl.clust, temp2, by = c("qtl.name"))
    
  }
}

r.all_total.qtl.clust.mat<-r.all_total.qtl.clust
row.names(r.all_total.qtl.clust.mat)<-r.all_total.qtl.clust.mat$qtl.name
r.all_total.qtl.clust.mat<-r.all_total.qtl.clust.mat[,-c(1)]


## Lets calculate distance between experiments
dt <- dist(t(r.all_total.qtl.clust.mat), method = "euclidean") # distance matrix
t.fit <- hclust(dt, method="ward.D") 
plot(t.fit, cex=0.8) # display dendogram

## Get the trait names and order of the names from the cluster fit

labs<-t.fit$labels
treatment<-c()
days<-c()
raw_trait<-c()
for(l in 1:length(labs)){
  t<-labs[l]
  fields<-strsplit(as.character(t), "_")
  treat<-fields[[1]][length(fields[[1]])]
  f<-(length(fields[[1]]) - 3)
  day<-fields[[1]][length(fields[[1]])-2]
  r.trait<-fields[[1]][1:f]
  r.trait<-paste(r.trait[1], r.trait[2], r.trait[3], sep="_")
  treatment<-c(treatment, treat)
  raw_trait<-c(raw_trait, r.trait)
  days<-c(days, day)
}

## simplify trait names
raw_trait[raw_trait == "sv_area_total"]<-c("sv_area")
raw_trait[raw_trait == "water_lost_total"]<-c("water_lost")
raw_trait[raw_trait == "wue_total_NA"]<-c("wue")
raw_trait[raw_trait == "te_fit_total"]<-c("fit")
raw_trait[raw_trait == "te_residual_total"]<-c("residual")

trait.cols<-rep('NA', length(raw_trait))
trait.cols[which(raw_trait == 'sv_area')]<-c("green")
trait.cols[which(raw_trait == 'water_lost')]<-c("blue")
trait.cols[which(raw_trait == 'wue')]<-c("orange")
trait.cols[which(raw_trait == 'fit')]<-c("black")
trait.cols[which(raw_trait == 'residual')]<-c("red")


## Add it colors for time (days after planting)
day.cols<-rep('NA', length(days))
days.i<-sort(unique(days))

## For each unique day, add gray scale color to day with white being a small day and black being late day
my_palette <- colorRampPalette(c("white", "black"))(n = length(unique(days.i)))

for(i in 1:length(days.i)){
  pt<-days.i[i]
  col<-my_palette[i]
  day.cols[which(days == pt)]<-col
}

## Add color treatment
treatment.cols<-rep("NA", length(treatment))
treatment.cols[treatment == 'wet']<-c('navy')
treatment.cols[treatment == 'dry']<-c('gold')

## Print the obtained cluster numbers in tree order
all.col<-cbind(day.cols, treatment.cols, trait.cols)

setwd(wue_figures.main.dir)

pdf("FIG_7.pdf")

plotDendroAndColors(t.fit, colors=all.col, groupLabels=c("Days", "Treatment", "Trait"), dendroLabels = FALSE, main="Dendrogram of Traits")

dev.off()

## Lets look at scree plot of clusters

## How many significant clusters?
wss <- (nrow(t(r.all_total.qtl.clust.mat))-1)*sum(apply(t(r.all_total.qtl.clust.mat),2,var))
for (i in 2:20) wss[i] <- sum(kmeans(t(r.all_total.qtl.clust.mat),
                                     centers=i)$withinss)
setwd(wue_figures.supplemental.dir)

pdf("FIG_S24.pdf")
plot(1:20, wss, type="b", xlab="Number of Clusters",
     ylab="Within groups sum of squares")
dev.off()




setwd(wue_results.qtl.dir)
load("fixed_marker_analysis.Rdata")



all_total_rate.qtl<-rbind(sv_area_total_rate.qtl, water_lost_total_rate.qtl, wue_ratio_total_rate.qtl, te_fit_total_rate.qtl, te_residual_total_rate.qtl)
all_total_rate.qtl<-all_total_rate.qtl[grepl("total", all_total_rate.qtl$raw_trait), ]

## Lets change the trait to a factor and place this in the correct order
all_total_rate.qtl$trait.treatment<-as.factor(all_total_rate.qtl$trait.treatment)

## Lets look at the levels of the treatment
print(levels(all_total_rate.qtl$trait.treatment)) 

## Now re-order them. Can you re-order them as a text string rather than column names?
all_total_rate.qtl$trait.treatment = factor(all_total_rate.qtl$trait.treatment,levels(all_total_rate.qtl$trait.treatment)[c(1:2,7:8,9:10,3:4,5:6)])

## Define colors based upon trait order
p.cols<-c("green", "dark green", "blue", "navy","gold", "orange", "grey80", "black", "pink", "red")

## Lets re-order the levels of the base trait name (raw_trait)
all_total_rate.qtl$raw_trait<-as.factor(all_total_rate.qtl$raw_trait)
print(levels(all_total_rate.qtl$raw_trait)) 
all_total_rate.qtl$raw_trait = factor(all_total_rate.qtl$raw_trait,levels(all_total_rate.qtl$raw_trait)[c(1,4,5,2,3)])

## Plot FX of major QTL
p1<-ggplot(all_total_rate.qtl[all_total_rate.qtl$qtl.name == "2@96",], aes(x=day, y=signed.prop.var, colour=trait.treatment)) + geom_line(size=0.3) + facet_wrap(~raw_trait, nrow=1) + scale_colour_manual(values = p.cols) + theme_bw() + labs(title="Major QTL") + ylab("2@96") + theme(legend.position="none") + ylim(c(-15, 30)) + xlab("")
p1<-p1 + geom_hline(yintercept=0, color="brown",linetype="dashed")
p2<-ggplot(all_total_rate.qtl[all_total_rate.qtl$qtl.name == "5@109",], aes(x=day, y=signed.prop.var, colour=trait.treatment)) + geom_line(size=0.3) + facet_wrap(~raw_trait, nrow=1) + scale_colour_manual(values = p.cols) + theme_bw() + ylab("5@109") + theme(legend.position="none") + ylim(c(-15, 30)) + xlab("")
p2<-p2 + geom_hline(yintercept=0, color="brown",linetype="dashed")
p3<-ggplot(all_total_rate.qtl[all_total_rate.qtl$qtl.name == "7@99",], aes(x=day, y=signed.prop.var, colour=trait.treatment)) + geom_line(size=0.3) + facet_wrap(~raw_trait, nrow=1) + scale_colour_manual(values = p.cols) + theme_bw() + ylab("7@99") + theme(legend.position="none") + ylim(c(-15, 30)) + xlab("")
p3<-p3 + geom_hline(yintercept=0, color="brown",linetype="dashed")
p4<-ggplot(all_total_rate.qtl[all_total_rate.qtl$qtl.name == "9@34",], aes(x=day, y=signed.prop.var, colour=trait.treatment)) + geom_line(size=0.3) + facet_wrap(~raw_trait, nrow=1) + scale_colour_manual(values = p.cols) + theme_bw() + ylab("9@34") + xlab("Days after planting") + theme(legend.position="none") + ylim(c(-15, 30)) + xlab("Days after planting")
p4<-p4 + geom_hline(yintercept=0, color="brown",linetype="dashed")

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}



setwd(wue_figures.main.dir)

pdf("FIG_8.pdf")
multiplot(p1, p2, p3, p4, nrow = 4)
dev.off()

## Plot FX of pos.minor.qtl
p1<-ggplot(all_total_rate.qtl[all_total_rate.qtl$qtl.name == "2@113",], aes(x=day, y=signed.prop.var, colour=trait.treatment)) + geom_line(size=0.3) + facet_wrap(~raw_trait, nrow=1) + scale_colour_manual(values = p.cols) + theme_bw() + labs(title="Positive effect QTL") + ylab("2@113") + theme(legend.position="none") + ylim(c(-5, 15)) + xlab("")
p1<-p1 + geom_hline(yintercept=0, color="brown",linetype="dashed")
p2<-ggplot(all_total_rate.qtl[all_total_rate.qtl$qtl.name == "3@48",], aes(x=day, y=signed.prop.var, colour=trait.treatment)) + geom_line(size=0.3) + facet_wrap(~raw_trait, nrow=1) + scale_colour_manual(values = p.cols) + theme_bw() + ylab("3@48") + theme(legend.position="none") + ylim(c(-5, 15)) + xlab("")
p2<-p2 + geom_hline(yintercept=0, color="brown",linetype="dashed")
p3<-ggplot(all_total_rate.qtl[all_total_rate.qtl$qtl.name == "4@52",], aes(x=day, y=signed.prop.var, colour=trait.treatment)) + geom_line(size=0.3) + facet_wrap(~raw_trait, nrow=1) + scale_colour_manual(values = p.cols) + theme_bw() + ylab("4@52") + theme(legend.position="none") + ylim(c(-5, 15)) + xlab("")
p3<-p3 + geom_hline(yintercept=0, color="brown",linetype="dashed")
p4<-ggplot(all_total_rate.qtl[all_total_rate.qtl$qtl.name == "6@65",], aes(x=day, y=signed.prop.var, colour=trait.treatment)) + geom_line(size=0.3) + facet_wrap(~raw_trait, nrow=1) + scale_colour_manual(values = p.cols) + theme_bw() + ylab("6@65") + theme(legend.position="none") + ylim(c(-5, 15)) + xlab("")
p4<-p4 + geom_hline(yintercept=0, color="brown",linetype="dashed")
p5<-ggplot(all_total_rate.qtl[all_total_rate.qtl$qtl.name == "9@127",], aes(x=day, y=signed.prop.var, colour=trait.treatment)) + geom_line(size=0.3) + facet_wrap(~raw_trait, nrow=1) + scale_colour_manual(values = p.cols) + theme_bw() + ylab("9@127") + theme(legend.position="none") + ylim(c(-5, 15)) + xlab("Days after planting")
p5<-p5 + geom_hline(yintercept=0, color="brown",linetype="dashed")

setwd(wue_figures.supplemental.dir)

## Plot QTL where B100 allele gives positive FX
pdf("FIG_S25a.pdf")
multiplot(p1, p2, p3, p4, p5, nrow = 5)
dev.off()

## Plot FX of neg.minor.qtl
p1<-ggplot(all_total_rate.qtl[all_total_rate.qtl$qtl.name == "2@11",], aes(x=day, y=signed.prop.var, colour=trait.treatment)) + geom_line(size=0.3) + facet_wrap(~raw_trait, nrow=1) + scale_colour_manual(values = p.cols) + theme_bw() + labs(title="Negative effect QTL") + ylab("2@11") + theme(legend.position="none") + ylim(c(-15, 5)) + xlab("")
p1<-p1 + geom_hline(yintercept=0, color="brown",linetype="dashed")
p2<-ggplot(all_total_rate.qtl[all_total_rate.qtl$qtl.name == "5@79",], aes(x=day, y=signed.prop.var, colour=trait.treatment)) + geom_line(size=0.3) + facet_wrap(~raw_trait, nrow=1) + scale_colour_manual(values = p.cols) + theme_bw() + ylab("5@79") + theme(legend.position="none") + ylim(c(-15, 5)) + xlab("")
p2<-p2 + geom_hline(yintercept=0, color="brown",linetype="dashed")
p3<-ggplot(all_total_rate.qtl[all_total_rate.qtl$qtl.name == "5@92",], aes(x=day, y=signed.prop.var, colour=trait.treatment)) + geom_line(size=0.3) + facet_wrap(~raw_trait, nrow=1) + scale_colour_manual(values = p.cols) + theme_bw() + ylab("5@95") + theme(legend.position="none") + ylim(c(-15, 5)) + xlab("")
p3<-p3 + geom_hline(yintercept=0, color="brown",linetype="dashed")
p4<-ggplot(all_total_rate.qtl[all_total_rate.qtl$qtl.name == "7@34",], aes(x=day, y=signed.prop.var, colour=trait.treatment)) + geom_line(size=0.3) + facet_wrap(~raw_trait, nrow=1) + scale_colour_manual(values = p.cols) + theme_bw() + ylab("7@34") + xlab("Days after planting") + theme(legend.position="none") + ylim(c(-15, 5)) + xlab("")
p4<-p4 + geom_hline(yintercept=0, color="brown",linetype="dashed")
p5<-ggplot(all_total_rate.qtl[all_total_rate.qtl$qtl.name == "7@51",], aes(x=day, y=signed.prop.var, colour=trait.treatment)) + geom_line(size=0.3) + facet_wrap(~raw_trait, nrow=1) + scale_colour_manual(values = p.cols) + theme_bw() + ylab("7@53") + xlab("Days after planting") + theme(legend.position="none") + ylim(c(-15, 5)) + xlab("Days after planting")
p5<-p5 + geom_hline(yintercept=0, color="brown",linetype="dashed")

setwd(wue_figures.supplemental.dir)
## Plot QTL where B100 allele gives negative FX

pdf("FIG_S25b.pdf")
multiplot(p1, p2, p3, p4, p5, nrow = 5)
dev.off()




## Combine both total and day to get total rate and remove difference QTL
all_total_rate.qtl<-rbind(all_total.qtl, all_day.qtl)
r.all_total_rate.qtl<-all_total_rate.qtl[all_total_rate.qtl$type == 'raw',]

##### CUMULATIVE TRAITS
 
## Make a new category 
r.all_total_rate.qtl$trait.treatment<-paste(r.all_total_rate.qtl$trait, r.all_total_rate.qtl$treatment, sep="_")

## Keep as cumulative 
r.all_total.qtl$trait.treatment<-paste(r.all_total.qtl$trait, r.all_total.qtl$treatment, sep="_")


#### Lets plot % from each parent
days<-unique(r.all_total.qtl$day)
output<-c()
for(d in days){
  temp1<-r.all_total.qtl[r.all_total.qtl$day == d,]
  traits<-unique(temp1$trait.treatment)
  for(t in traits){
    temp2<-temp1[temp1$trait.treatment == t,] 
    a10<-sum(temp2[temp2$additive.fx < 0, 'prop.var'])
    b100<-sum(temp2[temp2$additive.fx > 0, 'prop.var'])
    line1<-c('A10', d, t, a10)
    line2<-c('B100',d, t, b100)
    out<-rbind(line1, line2)
    output<-rbind(output, out)
  }
}

colnames(output)<-c('genotype','dap', 'trait.treatment', 'prop.var')
rownames(output)<-c(1:nrow(output))
output<-as.data.frame(output)
output$prop.var<-as.numeric(as.character(output$prop.var))
output$dap<-as.numeric(as.character(output$dap))

## Lets divide out the trait and the treatment
treatment<-c()
raw_trait<-c()
for (r in 1:nrow(output)) {
  temp<-output[r,]
  fields<-strsplit(as.character(temp$trait.treatment), '_')
  treat<-fields[[1]][length(fields[[1]])]
  f<-(length(fields[[1]]) - 3)
  r.trait<-fields[[1]][1:f]
  r.trait<-paste(r.trait[1], r.trait[2], r.trait[3], sep="_")
  treatment<-c(treatment, treat)
  raw_trait<-c(raw_trait, r.trait)
}

## This is a hack job
raw_trait<-gsub("wue_total_NA", "wue_ratio_total", raw_trait)

## Lets add in additoinal fields
output<-cbind(output, raw_trait)
output<-cbind(output, treatment)

total.output<-output

total.output.dry<-total.output[total.output$treatment == 'dry',]
total.output.wet<-total.output[total.output$treatment == 'wet',]

source("http://peterhaschke.com/Code/multiplot.R")

setwd(wue_figures.main.dir)

pdf("FIG_9.pdf")

# Dry
p1<-ggplot(total.output.dry, aes(dap, prop.var, colour=genotype)) + geom_line() + ggtitle("Dry") + facet_wrap(~raw_trait, nrow=1) + scale_colour_manual(values = c("red", "blue")) + theme_bw() + ylab("% variance explained") + xlab("Days after planting")
# Wet
p2<-ggplot(total.output.wet, aes(dap, prop.var, colour=genotype)) + geom_line() + ggtitle("Wet") + facet_wrap(~raw_trait, nrow=1) + scale_colour_manual(values = c("red", "blue")) + theme_bw() + ylab("% variance explained") + xlab("Days after planting")

multiplot(p1, p2, nrow = 2)

dev.off()




## Lets make the figure
setwd(wue_figures.supplemental.dir)
my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 98)
my_palette[1]<-c("#ffff00")

col_breaks = c(seq(-100,0,length=1),  # for yellow
               seq(0.000000001,.9999999999,length=48), # for blue
               seq(1,1.00000000001, length=1), # for white
               seq(1.000001,5,length=48), # for red
               seq(5.01,1000, length=1)) # extreme red

pdf("FIG_S26.pdf")

heatmap.2(sv_area.ratio,  notecol="black", col=my_palette, Colv="NA", Rowv="NA", trace="none", dendrogram = c("none"), density.info="none", key=F, breaks=col_breaks, symm=F,symkey=F,symbreaks=T, main="Plant size", xlab="Days after planting")

heatmap.2(water_lost.ratio,  notecol="black", col=my_palette, Colv="NA", Rowv="NA", trace="none", dendrogram = c("none"), density.info="none", key=F, breaks=col_breaks, symm=F,symkey=F,symbreaks=T, main="Water lost", xlab="Days after planting")

heatmap.2(wue_ratio.ratio,  notecol="black", col=my_palette, Colv="NA", Rowv="NA", trace="none", dendrogram = c("none"), density.info="none", key=F, breaks=col_breaks, symm=F,symkey=F,symbreaks=T, main="WUE ratio", xlab="Days after planting")

heatmap.2(te_fit.ratio,  notecol="black", col=my_palette, Colv="NA", Rowv="NA", trace="none", dendrogram = c("none"), density.info="none", key=F, breaks=col_breaks, symm=F,symkey=F,symbreaks=T, main="WUE fit", xlab="Days after planting")

heatmap.2(te_residual.ratio,  notecol="black", col=my_palette, Colv="NA", Rowv="NA", trace="none", dendrogram = c("none"), density.info="none", key=F, breaks=col_breaks, symm=F,symkey=F,symbreaks=T, main="WUE residual", xlab="Days after planting")

dev.off()











