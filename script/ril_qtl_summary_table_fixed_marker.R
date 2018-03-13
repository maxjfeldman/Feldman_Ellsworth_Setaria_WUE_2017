
library(ggplot2)
library(gplots)
source("http://peterhaschke.com/Code/multiplot.R")

## The first thing is to set your R session to the base directory you just downloaded from github
setwd()

## Tester
#setwd("~/Dropbox/Feldman_Ellsworth_Setaria_WUE_2017/")

##### CREATE DIRECTORY PATHS ##### 

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
wue_results.dir<-paste(home.dir, '/results', sep="")

## If this directory doesn't exist the download is not complete
if (file.exists(wue_results.dir)){
  print("Directory exists!")
} else {
  dir.create(file.path(wue_results.dir))
}

wue_results.qtl.dir<-paste(wue_results.dir, '/qtl_summary', sep="")

## If this directory doesn't exist the download is not complete
if (file.exists(wue_results.qtl.dir)){
  print("Directory exists!")
} else {
  dir.create(file.path(wue_results.qtl.dir))
}

wue_results.qtl.total_fixed.dir<-paste(wue_results.qtl.dir, '/total_fixed_markers', sep="")

## If this directory doesn't exist the download is not complete
if (file.exists(wue_results.qtl.total_fixed.dir)){
  print("Directory exists!")
} else {
  dir.create(file.path(wue_results.qtl.total_fixed.dir))
}

wue_results.qtl.total_fixed.st.dir<-paste(wue_results.qtl.total_fixed.dir, '/qtl_summary_fixed_marker_dir', sep="")

## If this directory doesn't exist the download is not complete
if (file.exists(wue_results.qtl.total_fixed.st.dir)){
  print("Directory exists!")
} else {
  dir.create(file.path(wue_results.qtl.total_fixed.st.dir))
}

wue_results.qtl.day_fixed.dir<-paste(wue_results.qtl.dir, '/day_fixed_markers', sep="")

## If this directory doesn't exist the download is not complete
if (file.exists(wue_results.qtl.day_fixed.dir)){
  print("Directory exists!")
} else {
  dir.create(file.path(wue_results.qtl.day_fixed.dir))
}

wue_results.qtl.day_fixed.st.dir<-paste(wue_results.qtl.day_fixed.dir, '/qtl_summary_fixed_marker_dir', sep="")

## If this directory doesn't exist the download is not complete
if (file.exists(wue_results.qtl.day_fixed.st.dir)){
  print("Directory exists!")
} else {
  dir.create(file.path(wue_results.qtl.day_fixed.st.dir))
}

##############################################################################################################################
# Load in QTL data from cumulative fixed marker anlaysis
##############################################################################################################################


setwd(wue_results.qtl.total_fixed.st.dir)
sv_area_total.qtl<-read.csv("sv_area_total_3_km_concatenated_summary_table.csv", header=F)
## Remove marker number
sv_area_total.qtl<-sv_area_total.qtl[,-c(1)]
colnames(sv_area_total.qtl)<-c('marker','chr','pos','lod','prop.var','additive.fx','additive.fx_se','trait','treatment','exp','year','type')

water_total.qtl<-read.csv("water_lost_total_3_km_concatenated_summary_table.csv", header=F)
## Remove marker number
water_total.qtl<-water_total.qtl[,-c(1)]
colnames(water_total.qtl)<-c('marker','chr','pos','lod','prop.var','additive.fx','additive.fx_se','trait','treatment','exp','year','type')

wue_ratio_total.qtl<-read.csv("wue_ratio_total_3_km_concatenated_summary_table.csv", header=F)
## Remove marker number
wue_ratio_total.qtl<-wue_ratio_total.qtl[,-c(1)]
colnames(wue_ratio_total.qtl)<-c('marker','chr','pos','lod','prop.var','additive.fx','additive.fx_se','trait','treatment','exp','year','type')

te_fit_total.qtl<-read.csv("te_fit_total_3_km_concatenated_summary_table.csv", header=F)
## Remove marker number
te_fit_total.qtl<-te_fit_total.qtl[,-c(1)]
colnames(te_fit_total.qtl)<-c('marker','chr','pos','lod','prop.var','additive.fx','additive.fx_se','trait','treatment','exp','year','type')

te_residual_total.qtl<-read.csv("te_residual_total_3_km_concatenated_summary_table.csv", header=F)
## Remove marker number
te_residual_total.qtl<-te_residual_total.qtl[,-c(1)]
colnames(te_residual_total.qtl)<-c('marker','chr','pos','lod','prop.var','additive.fx','additive.fx_se','trait','treatment','exp','year','type')

#isotope.qtl<-read.csv("stable_isotope_3_km_concatenated_summary_table.csv", header=F)
## Remove marker number
#isotope.qtl<-isotope.qtl[,-c(1)]
#colnames(isotope.qtl)<-c('marker','chr','pos','lod','prop.var','additive.fx','additive.fx_se','trait','treatment','exp','year','type')


all_total.qtl<-rbind(sv_area_total.qtl, water_total.qtl, wue_ratio_total.qtl, te_fit_total.qtl, te_residual_total.qtl)

## Lets add days after planting to the data.frame from the trait name
all_total.qtl$trait<-as.character(all_total.qtl$trait)
output<-c()
for (r in 1:nrow(all_total.qtl)){
  print(r)
  temp<-all_total.qtl[r,]
  t<-strsplit(temp$trait, "_")
  ## Get the day numer from the trait name 
  e<-length(t[[1]])-1
  day<-t[[1]][e]
  temp<-cbind(temp, day)
  output<-rbind(output, temp)
}

all_total.qtl<-output

all_total.qtl$day<-as.numeric(as.character(all_total.qtl$day))
all_total.qtl<-all_total.qtl[all_total.qtl$day > 16, ]
all_total.qtl$raw_trait<-rep("NA", nrow(all_total.qtl))
## Now lets replace the extended trait name with the actual trait name
all_total.qtl[grep('sv_area_total', all_total.qtl$trait),'raw_trait']<-c('sv_area_total')
all_total.qtl[grep('water_lost_total', all_total.qtl$trait),'raw_trait']<-c('water_lost_total')
all_total.qtl[grep('wue_total', all_total.qtl$trait),'raw_trait']<-c('wue_ratio_total')
all_total.qtl[grep('te_fit_total', all_total.qtl$trait),'raw_trait']<-c('te_fit_total')
all_total.qtl[grep('te_residual_total', all_total.qtl$trait),'raw_trait']<-c('te_residual_total')

## Add name of QTL
all_total.qtl$qtl.name<-paste(all_total.qtl$chr, as.integer(all_total.qtl$pos), sep="@")
## Add signed prop.var
all_total.qtl$signed.prop.var<-all_total.qtl$prop.var
all_total.qtl[all_total.qtl$additive.fx < 0, 'signed.prop.var']<-all_total.qtl[all_total.qtl$additive.fx < 0, 'signed.prop.var'] * -1

## Lets make a field that combines treatment and type
all_total.qtl$treatment.type<-paste(all_total.qtl$treatment, all_total.qtl$type, sep="_")



## Now make the field qtl.name (7@99) a factor and order the factor
keeper_qtl<-c("2@11", "2@96", "2@113", "3@48", "4@52", "5@79", "5@92", "5@109", "6@65", "7@34", "7@51", "7@99","9@34", "9@127")
all_total.qtl$qtl.name<-factor(all_total.qtl$qtl.name, levels=keeper_qtl)

## Lets look at results of biomass
sv_area_total.qtl<-all_total.qtl[all_total.qtl$raw_trait == 'sv_area_total',]

## Get a max and min value to specify y-axis
y.max<-max(sv_area_total.qtl$signed.prop.var)
y.min<-min(sv_area_total.qtl$signed.prop.var)

p<-ggplot(sv_area_total.qtl, aes(x=day, y=signed.prop.var, colour=treatment.type)) + geom_line() + facet_wrap(~qtl.name) + scale_colour_manual(values = c("green", "purple", "red","orange","navy")) + theme_bw() + labs(title="sv_area_total") + ylab("% Additive genetic variance") + xlab("Days after planting") + ylim(y.min,y.max)


## Lets look at results of water_lost_total
water_lost_total.qtl<-all_total.qtl[all_total.qtl$raw_trait == 'water_lost_total',]

## Get a max and min value to specify y-axis
y.max<-max(water_lost_total.qtl$signed.prop.var)
y.min<-min(water_lost_total.qtl$signed.prop.var)

p<-ggplot(water_lost_total.qtl, aes(x=day, y=signed.prop.var, colour=treatment.type)) + geom_line() + facet_wrap(~qtl.name) + scale_colour_manual(values = c("green", "purple", "red","orange","navy")) + theme_bw() + labs(title="water_lost_total") + ylab("% Additive genetic variance") + xlab("Days after planting") + ylim(y.min,y.max)

## Lets look at results of wue ratio
wue_ratio_total.qtl<-all_total.qtl[all_total.qtl$raw_trait == 'wue_ratio_total',]

## Get a max and min value to specify y-axis
y.max<-max(wue_ratio_total.qtl$signed.prop.var)
y.min<-min(wue_ratio_total.qtl$signed.prop.var)

p<-ggplot(wue_ratio_total.qtl, aes(x=day, y=signed.prop.var, colour=treatment.type)) + geom_line() + facet_wrap(~qtl.name) + scale_colour_manual(values = c("green", "purple", "red","orange","navy")) + theme_bw() + labs(title="wue_ratio_total") + ylab("% Additive genetic variance") + xlab("Days after planting") + ylim(y.min,y.max)


## Lets look at results of te_fit
te_fit_total.qtl<-all_total.qtl[all_total.qtl$raw_trait == 'te_fit_total',]

## Get a max and min value to specify y-axis
y.max<-max(te_fit_total.qtl$signed.prop.var)
y.min<-min(te_fit_total.qtl$signed.prop.var)

p<-ggplot(te_fit_total.qtl, aes(x=day, y=signed.prop.var, colour=treatment.type)) + geom_line() + facet_wrap(~qtl.name) + scale_colour_manual(values = c("green", "purple", "red","orange","navy")) + theme_bw() + labs(title="te_fit_total") + ylab("% Additive genetic variance") + xlab("Days after planting") + ylim(y.min,y.max)



## Lets look at results of te_residual
te_residual_total.qtl<-all_total.qtl[all_total.qtl$raw_trait == 'te_residual_total',]

## Get a max and min value to specify y-axis
y.max<-max(te_residual_total.qtl$signed.prop.var)
y.min<-min(te_residual_total.qtl$signed.prop.var)

p<-ggplot(te_residual_total.qtl, aes(x=day, y=signed.prop.var, colour=treatment.type)) + geom_line() + facet_wrap(~qtl.name) + scale_colour_manual(values = c("green", "purple", "red","orange","navy")) + theme_bw() + labs(title="te_residual_total") + ylab("% Additive genetic variance") + xlab("Days after planting") + ylim(y.min,y.max)


############################################################################
## Now compare cumulative traits to one another
############################################################################

r.all_total.qtl<-all_total.qtl[all_total.qtl$type == "raw",]
r.all_total.qtl$trait.treatment<-paste(r.all_total.qtl$raw_trait, r.all_total.qtl$treatment, sep=".")
# Get a max and min value to specify y-axis
y.max<-max(r.all_total.qtl$signed.prop.var)
y.min<-min(r.all_total.qtl$signed.prop.var)

p<-ggplot(r.all_total.qtl, aes(x=day, y=signed.prop.var, colour=trait.treatment)) + geom_line() + facet_wrap(~qtl.name) + scale_colour_manual(values = c("light green", "dark green", "grey","black","pink","red","light blue", "navy", "yellow", "orange"))  + theme_bw() + labs(title="All traits (cumulative)") + ylab("% Additive genetic variance") + xlab("Days after planting") + ylim(y.min,y.max)


## Lets plot difference
d.all_total.qtl<-all_total.qtl[all_total.qtl$type != "raw",]
d.all_total.qtl$trait.treatment<-paste(d.all_total.qtl$raw_trait, d.all_total.qtl$treatment, sep=".")
## Get a max and min value to specify y-axis
y.max<-max(d.all_total.qtl$signed.prop.var)
y.min<-min(d.all_total.qtl$signed.prop.var)

p<-ggplot(d.all_total.qtl, aes(x=day, y=signed.prop.var, colour=trait.treatment)) + geom_line() + facet_wrap(~qtl.name) + scale_colour_manual(values = c("light green", "dark green", "grey","black","pink","red","light blue", "navy", "yellow", "orange"))  + theme_bw() + labs(title="All difference traits (cumulative)") + ylab("% Additive genetic variance") + xlab("totals after planting") + ylim(y.min,y.max)


##############################################################################################################################
## Load in QTL data from daily fixed marker anlaysis
##############################################################################################################################


setwd(wue_results.qtl.day_fixed.st.dir)
sv_area_day.qtl<-read.csv("sv_area_day_3_km_concatenated_summary_table.csv", header=F)
## Remove marker number
sv_area_day.qtl<-sv_area_day.qtl[,-c(1)]
colnames(sv_area_day.qtl)<-c('marker','chr','pos','lod','prop.var','additive.fx','additive.fx_se','trait','treatment','exp','year','type')

water_day.qtl<-read.csv("water_lost_day_3_km_concatenated_summary_table.csv", header=F)
## Remove marker number
water_day.qtl<-water_day.qtl[,-c(1)]
colnames(water_day.qtl)<-c('marker','chr','pos','lod','prop.var','additive.fx','additive.fx_se','trait','treatment','exp','year','type')

wue_ratio_day.qtl<-read.csv("wue_ratio_day_3_km_concatenated_summary_table.csv", header=F)
## Remove marker number
wue_ratio_day.qtl<-wue_ratio_day.qtl[,-c(1)]
colnames(wue_ratio_day.qtl)<-c('marker','chr','pos','lod','prop.var','additive.fx','additive.fx_se','trait','treatment','exp','year','type')

te_fit_day.qtl<-read.csv("te_fit_day_3_km_concatenated_summary_table.csv", header=F)
## Remove marker number
te_fit_day.qtl<-te_fit_day.qtl[,-c(1)]
colnames(te_fit_day.qtl)<-c('marker','chr','pos','lod','prop.var','additive.fx','additive.fx_se','trait','treatment','exp','year','type')

te_residual_day.qtl<-read.csv("te_residual_day_3_km_concatenated_summary_table.csv", header=F)
## Remove marker number
te_residual_day.qtl<-te_residual_day.qtl[,-c(1)]
colnames(te_residual_day.qtl)<-c('marker','chr','pos','lod','prop.var','additive.fx','additive.fx_se','trait','treatment','exp','year','type')



all_day.qtl<-rbind(sv_area_day.qtl, water_day.qtl, wue_ratio_day.qtl, te_fit_day.qtl, te_residual_day.qtl)

## Lets add days after planting to the data.frame from the trait name
all_day.qtl$trait<-as.character(all_day.qtl$trait)
output<-c()
for (r in 1:nrow(all_day.qtl)){
  print(r)
  temp<-all_day.qtl[r,]
  t<-strsplit(temp$trait, "_")
  ## Get the day numer from the trait name 
  e<-length(t[[1]])-1
  day<-t[[1]][e]
  temp<-cbind(temp, day)
  output<-rbind(output, temp)
}

all_day.qtl<-output

all_day.qtl$day<-as.numeric(as.character(all_day.qtl$day))
all_day.qtl<-all_day.qtl[all_day.qtl$day > 16, ]
all_day.qtl$raw_trait<-rep("NA", nrow(all_day.qtl))
## Now lets replace the extended trait name with the actual trait name
all_day.qtl[grep('sv_area_day', all_day.qtl$trait),'raw_trait']<-c('sv_area_day')
all_day.qtl[grep('water_lost', all_day.qtl$trait),'raw_trait']<-c('water_lost_day')
all_day.qtl[grep('wue_day', all_day.qtl$trait),'raw_trait']<-c('wue_ratio_day')
all_day.qtl[grep('te_fit_day', all_day.qtl$trait),'raw_trait']<-c('te_fit_day')
all_day.qtl[grep('te_residual_day', all_day.qtl$trait),'raw_trait']<-c('te_residual_day')

## Add name of QTL
all_day.qtl$qtl.name<-paste(all_day.qtl$chr, as.integer(all_day.qtl$pos), sep="@")

## Now make the field qtl.name (7@99) a factor and order the factor
keeper_qtl<-c("2@11", "2@96", "2@113", "3@48", "4@52", "5@79", "5@92", "5@109", "6@65", "7@34", "7@51", "7@99","9@34", "9@127")
all_day.qtl$qtl.name<-factor(all_day.qtl$qtl.name, levels=keeper_qtl)

## Add signed prop.var
all_day.qtl$signed.prop.var<-all_day.qtl$prop.var
all_day.qtl[all_day.qtl$additive.fx < 0, 'signed.prop.var']<-all_day.qtl[all_day.qtl$additive.fx < 0, 'signed.prop.var'] * -1

## Lets make a field that combines treatment and type
all_day.qtl$treatment.type<-paste(all_day.qtl$treatment, all_day.qtl$type, sep="_")

## Lets look at results of biomass
sv_area_day.qtl<-all_day.qtl[all_day.qtl$raw_trait == 'sv_area_day',]

# Get a max and min value to specify y-axis
y.max<-max(sv_area_day.qtl$signed.prop.var)
y.min<-min(sv_area_day.qtl$signed.prop.var)

p<-ggplot(sv_area_day.qtl, aes(x=day, y=signed.prop.var, colour=treatment.type)) + geom_line() + facet_wrap(~qtl.name) + scale_colour_manual(values = c("green", "purple", "red","orange","navy")) + theme_bw() + labs(title="sv_area_day") + ylab("% Additive genetic variance") + xlab("Days after planting") + ylim(y.min,y.max)

## Lets look at results of water lost
water_lost_day.qtl<-all_day.qtl[all_day.qtl$raw_trait == 'water_lost_day',]

# Get a max and min value to specify y-axis
y.max<-max(water_lost_day.qtl$signed.prop.var)
y.min<-min(water_lost_day.qtl$signed.prop.var)

p<-ggplot(water_lost_day.qtl, aes(x=day, y=signed.prop.var, colour=treatment.type)) + geom_line() + facet_wrap(~qtl.name) + scale_colour_manual(values = c("green", "purple", "red","orange","navy")) + theme_bw() + labs(title="water_lost_day") + ylab("% Additive genetic variance") + xlab("Days after planting") + ylim(y.min,y.max)


## Lets look at results of biomass
wue_ratio_day.qtl<-all_day.qtl[all_day.qtl$raw_trait == 'wue_ratio_day',]

# Get a max and min value to specify y-axis
y.max<-max(wue_ratio_day.qtl$signed.prop.var)
y.min<-min(wue_ratio_day.qtl$signed.prop.var)

p<-ggplot(wue_ratio_day.qtl, aes(x=day, y=signed.prop.var, colour=treatment.type)) + geom_line() + facet_wrap(~qtl.name) + scale_colour_manual(values = c("green", "purple", "red","orange","navy")) + theme_bw() + labs(title="wue_ratio_day") + ylab("% Additive genetic variance") + xlab("Days after planting") + ylim(y.min,y.max)

## Lets look at results of biomass
te_fit_day.qtl<-all_day.qtl[all_day.qtl$raw_trait == 'te_fit_day',]

# Get a max and min value to specify y-axis
y.max<-max(te_fit_day.qtl$signed.prop.var)
y.min<-min(te_fit_day.qtl$signed.prop.var)

p<-ggplot(te_fit_day.qtl, aes(x=day, y=signed.prop.var, colour=treatment.type)) + geom_line() + facet_wrap(~qtl.name) + scale_colour_manual(values = c("green", "purple", "red","orange","navy")) + theme_bw() + labs(title="te_fit_day") + ylab("% Additive genetic variance") + xlab("Days after planting") + ylim(y.min,y.max)


## Lets look at results of biomass
te_residual_day.qtl<-all_day.qtl[all_day.qtl$raw_trait == 'te_residual_day',]

# Get a max and min value to specify y-axis
y.max<-max(te_residual_day.qtl$signed.prop.var)
y.min<-min(te_residual_day.qtl$signed.prop.var)

p<-ggplot(te_residual_day.qtl, aes(x=day, y=signed.prop.var, colour=treatment.type)) + geom_line() + facet_wrap(~qtl.name) + scale_colour_manual(values = c("green", "purple", "red","orange","navy")) + theme_bw() + labs(title="te_residual_day") + ylab("% Additive genetic variance") + xlab("Days after planting") + ylim(y.min,y.max)




############################################################################
# Now compare daily/rate traits to one another
############################################################################

## Lets start with raw
r.all_day.qtl<-all_day.qtl[all_day.qtl$type == "raw",]
r.all_day.qtl$trait.treatment<-paste(r.all_day.qtl$raw_trait, r.all_day.qtl$treatment, sep=".")
# Get a max and min value to specify y-axis
y.max<-max(r.all_day.qtl$signed.prop.var)
y.min<-min(r.all_day.qtl$signed.prop.var)

p<-ggplot(r.all_day.qtl, aes(x=day, y=signed.prop.var, colour=trait.treatment)) + geom_line() + facet_wrap(~qtl.name) + scale_colour_manual(values = c("light green", "dark green", "grey","black","pink","red","light blue", "navy", "yellow", "orange"))  + theme_bw() + labs(title="All traits (cumulative)") + ylab("% Additive genetic variance") + xlab("Days after planting") + ylim(y.min,y.max)


## Lets plot difference
d.all_day.qtl<-all_day.qtl[all_day.qtl$type == "raw",]
d.all_day.qtl$trait.treatment<-paste(d.all_day.qtl$raw_trait, d.all_day.qtl$treatment, sep=".")
# Get a max and min value to specify y-axis
y.max<-max(d.all_day.qtl$signed.prop.var)
y.min<-min(d.all_day.qtl$signed.prop.var)

p<-ggplot(d.all_day.qtl, aes(x=day, y=signed.prop.var, colour=trait.treatment)) + geom_line() + facet_wrap(~qtl.name) + scale_colour_manual(values = c("light green", "dark green", "grey","black","pink","red","light blue", "navy", "yellow", "orange"))  + theme_bw() + labs(title="All difference traits (daily)") + ylab("% Additive genetic variance") + xlab("Days after planting") + ylim(y.min,y.max)


##############################################################################################################################
# Now compare cumulative with daily
##############################################################################################################################

all.qtl<-rbind(all_total.qtl,all_day.qtl)
r.all.qtl<-all.qtl[all.qtl$type == 'raw',]


## Lets look at results of biomass
r.sv_area.qtl<-r.all.qtl[grepl('sv_area', r.all.qtl$raw_trait),]

## Lets add variable treatment trait
r.sv_area.qtl$time<-rep("NA", nrow(r.sv_area.qtl))
r.sv_area.qtl[r.sv_area.qtl$raw_trait == "sv_area_total", "time"]<-c("total")
r.sv_area.qtl[r.sv_area.qtl$raw_trait != "sv_area_total", "time"]<-c("day")
r.sv_area.qtl$treatment.trait<-paste(r.sv_area.qtl$treatment, r.sv_area.qtl$time, sep="_")


## Get a max and min value to specify y-axis
y.max<-max(r.sv_area.qtl$signed.prop.var)
y.min<-min(r.sv_area.qtl$signed.prop.var)

p<-ggplot(r.sv_area.qtl, aes(x=day, y=signed.prop.var, colour=treatment.trait)) + geom_line() + facet_wrap(~qtl.name) + scale_colour_manual(values = c("yellow", "orange", "light blue","navy")) + theme_bw() + labs(title="sv_area") + ylab("% Additive genetic variance") + xlab("Days after planting") + ylim(y.min,y.max)


## Lets look at results of water_lost
r.water_lost.qtl<-r.all.qtl[grepl('water_lost', r.all.qtl$raw_trait),]

## Lets add variable treatment trait
r.water_lost.qtl$time<-rep("NA", nrow(r.water_lost.qtl))
r.water_lost.qtl[r.water_lost.qtl$raw_trait == "water_lost_total", "time"]<-c("total")
r.water_lost.qtl[r.water_lost.qtl$raw_trait != "water_lost_total", "time"]<-c("day")
r.water_lost.qtl$treatment.trait<-paste(r.water_lost.qtl$treatment, r.water_lost.qtl$time, sep="_")


## Get a max and min value to specify y-axis
y.max<-max(r.water_lost.qtl$signed.prop.var)
y.min<-min(r.water_lost.qtl$signed.prop.var)

p<-ggplot(r.water_lost.qtl, aes(x=day, y=signed.prop.var, colour=treatment.trait)) + geom_line() + facet_wrap(~qtl.name) + scale_colour_manual(values = c("yellow", "orange", "light blue","navy")) + theme_bw() + labs(title="water_lost") + ylab("% Additive genetic variance") + xlab("Days after planting") + ylim(y.min,y.max)

## Lets look at results of wue_ratio
r.wue_ratio.qtl<-r.all.qtl[grepl('wue_ratio', r.all.qtl$raw_trait),]

## Lets add variable treatment trait
r.wue_ratio.qtl$time<-rep("NA", nrow(r.wue_ratio.qtl))
r.wue_ratio.qtl[r.wue_ratio.qtl$raw_trait == "wue_ratio_total", "time"]<-c("total")
r.wue_ratio.qtl[r.wue_ratio.qtl$raw_trait != "wue_ratio_total", "time"]<-c("day")
r.wue_ratio.qtl$treatment.trait<-paste(r.wue_ratio.qtl$treatment, r.wue_ratio.qtl$time, sep="_")


## Get a max and min value to specify y-axis
y.max<-max(r.wue_ratio.qtl$signed.prop.var)
y.min<-min(r.wue_ratio.qtl$signed.prop.var)

p<-ggplot(r.wue_ratio.qtl, aes(x=day, y=signed.prop.var, colour=treatment.trait)) + geom_line() + facet_wrap(~qtl.name) + scale_colour_manual(values = c("yellow", "orange", "light blue","navy")) + theme_bw() + labs(title="wue_ratio") + ylab("% Additive genetic variance") + xlab("Days after planting") + ylim(y.min,y.max)


## Lets look at results of te_fit
r.te_fit.qtl<-r.all.qtl[grepl('te_fit', r.all.qtl$raw_trait),]

## Lets add variable treatment trait
r.te_fit.qtl$time<-rep("NA", nrow(r.te_fit.qtl))
r.te_fit.qtl[r.te_fit.qtl$raw_trait == "te_fit_total", "time"]<-c("total")
r.te_fit.qtl[r.te_fit.qtl$raw_trait != "te_fit_total", "time"]<-c("day")
r.te_fit.qtl$treatment.trait<-paste(r.te_fit.qtl$treatment, r.te_fit.qtl$time, sep="_")


## Get a max and min value to specify y-axis
y.max<-max(r.te_fit.qtl$signed.prop.var)
y.min<-min(r.te_fit.qtl$signed.prop.var)

p<-ggplot(r.te_fit.qtl, aes(x=day, y=signed.prop.var, colour=treatment.trait)) + geom_line() + facet_wrap(~qtl.name) + scale_colour_manual(values = c("yellow", "orange", "light blue","navy")) + theme_bw() + labs(title="te_fit") + ylab("% Additive genetic variance") + xlab("Days after planting") + ylim(y.min,y.max)



## Lets look at results of te_residual
r.te_residual.qtl<-r.all.qtl[grepl('te_residual', r.all.qtl$raw_trait),]

## Lets add variable treatment trait
r.te_residual.qtl$time<-rep("NA", nrow(r.te_residual.qtl))
r.te_residual.qtl[r.te_residual.qtl$raw_trait == "te_residual_total", "time"]<-c("total")
r.te_residual.qtl[r.te_residual.qtl$raw_trait != "te_residual_total", "time"]<-c("day")
r.te_residual.qtl$treatment.trait<-paste(r.te_residual.qtl$treatment, r.te_residual.qtl$time, sep="_")



## Get a max and min value to specify y-axis
y.max<-max(r.te_residual.qtl$signed.prop.var)
y.min<-min(r.te_residual.qtl$signed.prop.var)

p<-ggplot(r.te_residual.qtl, aes(x=day, y=signed.prop.var, colour=treatment.trait)) + geom_line() + facet_wrap(~qtl.name) + scale_colour_manual(values = c("yellow", "orange", "light blue","navy")) + theme_bw() + labs(title="te_residual") + ylab("% Additive genetic variance") + xlab("Days after planting") + ylim(y.min,y.max)




setwd(wue_results.qtl.dir)
save.image("fixed_marker_analysis.Rdata")


#######################################################################
## Lets make some other plots using a trait specific color scheme
#######################################################################

## Lets make some FX plots of each trait and its rate of change
## Lets look at results of biomass

## Combine cumulative and daily
sv_area_total_rate.qtl<-rbind(sv_area_total.qtl, sv_area_day.qtl)
sv_area_total_rate.qtl<-sv_area_total_rate.qtl[sv_area_total_rate.qtl$type == 'raw',]

## Get a max and min value to specify y-axis
y.max<-max(sv_area_total_rate.qtl$signed.prop.var)
y.min<-min(sv_area_total_rate.qtl$signed.prop.var)

## Make a new category 
sv_area_total_rate.qtl$trait.treatment<-paste(sv_area_total_rate.qtl$raw_trait, sv_area_total_rate.qtl$treatment, sep="_")

p<-ggplot(sv_area_total_rate.qtl, aes(x=day, y=signed.prop.var, colour=trait.treatment)) + geom_line(size=0.3) + facet_wrap(~qtl.name) + scale_colour_manual(values = c("purple1", "purple4", "green","dark green")) + theme_bw() + labs(title="sv_area") + ylab("% Additive genetic variance") + xlab("Days after planting") + ylim(y.min,y.max)

setwd(wue_results.qtl.dir)

pdf("sv_area_total_rate.qtl_fx.pdf")
print(p)
dev.off()

## Combine cumulative and daily
water_lost_total_rate.qtl<-rbind(water_lost_total.qtl, water_lost_day.qtl)
water_lost_total_rate.qtl<-water_lost_total_rate.qtl[water_lost_total_rate.qtl$type == 'raw',]

## Get a max and min value to specify y-axis
y.max<-max(water_lost_total_rate.qtl$signed.prop.var)
y.min<-min(water_lost_total_rate.qtl$signed.prop.var)

## Make a new category 
water_lost_total_rate.qtl$trait.treatment<-paste(water_lost_total_rate.qtl$raw_trait, water_lost_total_rate.qtl$treatment, sep="_")

p<-ggplot(water_lost_total_rate.qtl, aes(x=day, y=signed.prop.var, colour=trait.treatment)) + geom_line(size=0.3) + facet_wrap(~qtl.name) + scale_colour_manual(values = c("purple1", "purple4", "blue","navy")) + theme_bw() + labs(title="water_lost") + ylab("% Additive genetic variance") + xlab("Days after planting") + ylim(y.min,y.max)

pdf("water_lost_total_rate.qtl_fx.pdf")
print(p)
dev.off()

## Combine cumulative and daily
wue_ratio_total_rate.qtl<-rbind(wue_ratio_total.qtl, wue_ratio_day.qtl)
wue_ratio_total_rate.qtl<-wue_ratio_total_rate.qtl[wue_ratio_total_rate.qtl$type == 'raw',]

## Get a max and min value to specify y-axis
y.max<-max(wue_ratio_total_rate.qtl$signed.prop.var)
y.min<-min(wue_ratio_total_rate.qtl$signed.prop.var)

## Make a new category 
wue_ratio_total_rate.qtl$trait.treatment<-paste(wue_ratio_total_rate.qtl$raw_trait, wue_ratio_total_rate.qtl$treatment, sep="_")

p<-ggplot(wue_ratio_total_rate.qtl, aes(x=day, y=signed.prop.var, colour=trait.treatment)) + geom_line(size=0.3) + facet_wrap(~qtl.name) + scale_colour_manual(values = c("purple1", "purple4", "gold","orange")) + theme_bw() + labs(title="wue_ratio") + ylab("% Additive genetic variance") + xlab("Days after planting") + ylim(y.min,y.max)

pdf("wue_ratio_total_rate.qtl_fx.pdf")
print(p)
dev.off()


## Combine cumulative and daily
te_fit_total_rate.qtl<-rbind(te_fit_total.qtl, te_fit_day.qtl)
te_fit_total_rate.qtl<-te_fit_total_rate.qtl[te_fit_total_rate.qtl$type == 'raw',]

## Get a max and min value to specify y-axis
y.max<-max(te_fit_total_rate.qtl$signed.prop.var)
y.min<-min(te_fit_total_rate.qtl$signed.prop.var)

## Make a new category 
te_fit_total_rate.qtl$trait.treatment<-paste(te_fit_total_rate.qtl$raw_trait, te_fit_total_rate.qtl$treatment, sep="_")

p<-ggplot(te_fit_total_rate.qtl, aes(x=day, y=signed.prop.var, colour=trait.treatment)) + geom_line(size=0.3) + facet_wrap(~qtl.name) + scale_colour_manual(values = c("purple1", "purple4", "grey80","black")) + theme_bw() + labs(title="te_fit") + ylab("% Additive genetic variance") + xlab("Days after planting") + ylim(y.min,y.max)

pdf("te_fit_total_rate.qtl_fx.pdf")
print(p)
dev.off()


## Combine cumulative and daily
te_residual_total_rate.qtl<-rbind(te_residual_total.qtl, te_residual_day.qtl)
te_residual_total_rate.qtl<-te_residual_total_rate.qtl[te_residual_total_rate.qtl$type == 'raw',]

## Get a max and min value to specify y-axis
y.max<-max(te_residual_total_rate.qtl$signed.prop.var)
y.min<-min(te_residual_total_rate.qtl$signed.prop.var)

## Make a new category 
te_residual_total_rate.qtl$trait.treatment<-paste(te_residual_total_rate.qtl$raw_trait, te_residual_total_rate.qtl$treatment, sep="_")

p<-ggplot(te_residual_total_rate.qtl, aes(x=day, y=signed.prop.var, colour=trait.treatment)) + geom_line(size=0.3) + facet_wrap(~qtl.name) + scale_colour_manual(values = c("purple1", "purple4", "pink","red")) + theme_bw() + labs(title="te_residual") + ylab("% Additive genetic variance") + xlab("Days after planting") + ylim(y.min,y.max)

pdf("te_residual_total_rate.qtl_fx.pdf")
print(p)
dev.off()


#######################################################################
# Lets make plot of QTL PVE effects that appears in the manuscript
#######################################################################

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

setwd(wue_results.qtl.dir)
pdf("major_qtl_fx.pdf")
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

setwd(wue_results.qtl.dir)
pdf("positive_minor_qtl_fx.pdf")
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

setwd(wue_results.qtl.dir)
pdf("negative_minor_qtl_fx.pdf")
multiplot(p1, p2, p3, p4, p5, nrow = 5)
dev.off()


#################################################################
## Plot parental contributions
#################################################################

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

setwd(wue_results.qtl.total_fixed.dir)

pdf("contribution_of_parental_alleles_cumulative.pdf")

# Dry
p1<-ggplot(total.output.dry, aes(dap, prop.var, colour=genotype)) + geom_line() + ggtitle("dry") + facet_wrap(~raw_trait, nrow=1) + scale_colour_manual(values = c("red", "blue")) + theme_bw() + ylab("PVE") + xlab("Days after planting")
# Wet
p2<-ggplot(total.output.wet, aes(dap, prop.var, colour=genotype)) + geom_line() + ggtitle("wet") + facet_wrap(~raw_trait, nrow=1) + scale_colour_manual(values = c("red", "blue")) + theme_bw() + ylab("PVE") + xlab("Days after planting")

multiplot(p1, p2, nrow = 2)

dev.off()


##### RATE TRAITS

## Associate treatment and trait
r.all_day.qtl$trait.treatment<-paste(r.all_day.qtl$trait, r.all_day.qtl$treatment, sep="_")

#### Lets plot % from each parent of the rate statistic
days<-unique(r.all_day.qtl$day)
output<-c()
for(d in days){
  temp1<-r.all_day.qtl[r.all_day.qtl$day == d,]
  traits<-unique(temp1$trait.treatment)
  for(t in traits){
    temp2<-temp1[temp1$trait.treatment == t,] 
    a10<-sum(temp2[temp2$additive.fx < 0, 'prop.var'])
    b100<-sum(temp2[temp2$additive.fx > 0, 'prop.var'])
    line1<-c('a10', d, t, a10)
    line2<-c('b100',d, t, b100)
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
  f<-(length(fields[[1]]) - 1)
  r.trait<-fields[[1]][1:f]
  r.trait<-paste(r.trait[1], r.trait[2], r.trait[3], sep="_")
  treatment<-c(treatment, treat)
  raw_trait<-c(raw_trait, r.trait)
}


## Lets add in additoinal fields
output<-cbind(output, raw_trait)
output<-cbind(output, treatment)

output.dry<-output[output$treatment == 'dry',]
output.wet<-output[output$treatment == 'wet',]

setwd(wue_results.qtl.day_fixed.dir)

pdf("contribution_of_parental_alleles_rate.pdf")

# Dry
p1<-ggplot(total.output.dry, aes(dap, prop.var, colour=genotype)) + geom_line() + ggtitle("dry") + facet_wrap(~raw_trait, nrow=1) + scale_colour_manual(values = c("red", "blue")) + theme_bw() + ylab("PVE") + xlab("Days after planting")
# Wet
p2<-ggplot(total.output.wet, aes(dap, prop.var, colour=genotype)) + geom_line() + ggtitle("wet") + facet_wrap(~raw_trait, nrow=1) + scale_colour_manual(values = c("red", "blue")) + theme_bw() + ylab("PVE") + xlab("Days after planting")

multiplot(p1, p2, nrow = 2)

dev.off()


#######################################################################
## Lets exmaine FX across treatment blocks 
## Conditional neutrailty vs antagonistic pleiotropy
#######################################################################

setwd(wue_results.qtl.dir)

## Lets try to make a unique plot across days for each trait

traits<-as.character(unique(r.all_total.qtl$raw_trait))
treatments<-as.character(unique(r.all_total.qtl$treatment))
days<-as.character(unique(r.all_total.qtl$day))

## Start with sv_area (plant size)
sv_area.fx<-r.all_total.qtl[r.all_total.qtl$raw_trait == 'sv_area_total',]
sv_area.ratio.calc<-c()
for(d in days){
  ## Extract wet
  temp1<-sv_area.fx[sv_area.fx$day == d, ]
  temp.w<-temp1[temp1$treatment == 'wet',c("marker", "qtl.name", "day", "signed.prop.var")]
  colnames(temp.w)[4]<-c("signed.prop.var.wet")
  ## Extract dry
  temp.d<-temp1[temp1$treatment == 'dry',c("marker", "qtl.name", "day", "signed.prop.var")]
  colnames(temp.d)[4]<-c("signed.prop.var.dry")
  ## Put in common data.frame
  temp2<-merge(temp.w, temp.d, by=c("marker", "qtl.name", "day"))
  ## Replace very small values with a reasonable value
  temp2[(temp2$signed.prop.var.wet < .2 & temp2$signed.prop.var.wet > 0),"signed.prop.var.wet"]<-c(.2)
  temp2[(temp2$signed.prop.var.wet > -.2 & temp2$signed.prop.var.wet < 0),"signed.prop.var.wet"]<-c(-.2)
  temp2[(temp2$signed.prop.var.dry < .2 & temp2$signed.prop.var.dry > 0),"signed.prop.var.dry"]<-c(.2)
  temp2[(temp2$signed.prop.var.dry > -.2 & temp2$signed.prop.var.dry < 0),"signed.prop.var.dry"]<-c(-.2)
  
  ## Take the ratio between treatment blocks (if ratio > 1 wet fx is greater, if < 1 dry fx is greater)
  temp2$ratio<-temp2$signed.prop.var.wet/temp2$signed.prop.var.dry
  ## If the value is negative and fx are almost zero remove the sign of small fx
  for(row in 1:nrow(temp2)){
    if((temp2[row, "ratio"] < 0) & (abs(temp2[row, "signed.prop.var.wet"]) == 0.2 | (abs(temp2 [row, "signed.prop.var.dry"]) == 0.2))){
      temp2[row, "ratio"]<-abs(temp2[row, "ratio"])
    } 
  }
  ## write to a new data.frame
  sv_area.ratio.calc<-rbind(sv_area.ratio.calc, temp2)
}

## Lets now make it matrix for plotting...
sv_area.ratio<-c()

for (d in days){
  temp1<-sv_area.ratio.calc[sv_area.ratio.calc$day == d, 'ratio']
  sv_area.ratio<-cbind(sv_area.ratio, temp1)
}

## Lets round the values to two decimal places
sv_area.ratio<-round(sv_area.ratio,2)

rownames(sv_area.ratio)<-unique(sv_area.ratio.calc$qtl.name)
colnames(sv_area.ratio)<-days

my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 98)
my_palette[1]<-c("#ffff00")

col_breaks = c(seq(-100,0,length=1),  # for yellow
               seq(0.000000001,.9999999999,length=48), # for blue
               seq(1,1.00000000001, length=1), # for white
               seq(1.000001,5,length=48), # for red
               seq(5.01,1000, length=1)) # extreme red

## Lets re-order the QTL in each data frame to reflect their physical order on chromsomes
sv_area.ratio<-sv_area.ratio[c(1:12,14,13),]

heatmap.2(sv_area.ratio,  notecol="black", col=my_palette, Colv="NA", Rowv="NA", trace="none", dendrogram = c("none"), density.info="none", key=F, breaks=col_breaks, symm=F,symkey=F,symbreaks=T)


## Now process water use 
water_lost.fx<-r.all_total.qtl[r.all_total.qtl$raw_trait == 'water_lost_total',]
water_lost.ratio.calc<-c()
for(d in days){
  ## Extract wet
  temp1<-water_lost.fx[water_lost.fx$day == d, ]
  temp.w<-temp1[temp1$treatment == 'wet',c("marker", "qtl.name", "day", "signed.prop.var")]
  colnames(temp.w)[4]<-c("signed.prop.var.wet")
  ## Extract dry
  temp.d<-temp1[temp1$treatment == 'dry',c("marker", "qtl.name", "day", "signed.prop.var")]
  colnames(temp.d)[4]<-c("signed.prop.var.dry")
  ## Put in common data.frame
  temp2<-merge(temp.w, temp.d, by=c("marker", "qtl.name", "day"))
  ## Replace very small values with a reasonable value
  temp2[(temp2$signed.prop.var.wet < .2 & temp2$signed.prop.var.wet > 0),"signed.prop.var.wet"]<-c(.2)
  temp2[(temp2$signed.prop.var.wet > -.2 & temp2$signed.prop.var.wet < 0),"signed.prop.var.wet"]<-c(-.2)
  temp2[(temp2$signed.prop.var.dry < .2 & temp2$signed.prop.var.dry > 0),"signed.prop.var.dry"]<-c(.2)
  temp2[(temp2$signed.prop.var.dry > -.2 & temp2$signed.prop.var.dry < 0),"signed.prop.var.dry"]<-c(-.2)
  
  ## Take the ratio between treatment blocks (if ratio > 1 wet fx is greater, if < 1 dry fx is greater)
  temp2$ratio<-temp2$signed.prop.var.wet/temp2$signed.prop.var.dry
  ## If the value is negative and fx are almost zero remove the sign of small fx
  for(row in 1:nrow(temp2)){
    if((temp2[row, "ratio"] < 0) & (abs(temp2[row, "signed.prop.var.wet"]) == 0.2 | (abs(temp2 [row, "signed.prop.var.dry"]) == 0.2))){
      temp2[row, "ratio"]<-abs(temp2[row, "ratio"])
    } 
  }
  ## write to a new data.frame
  water_lost.ratio.calc<-rbind(water_lost.ratio.calc, temp2)
}

## Lets now make it matrix for plotting...
water_lost.ratio<-c()

for (d in days){
  temp1<-water_lost.ratio.calc[water_lost.ratio.calc$day == d, 'ratio']
  water_lost.ratio<-cbind(water_lost.ratio, temp1)
}

## Lets round the values to two decimal places
water_lost.ratio<-round(water_lost.ratio,2)

rownames(water_lost.ratio)<-unique(water_lost.ratio.calc$qtl.name)
colnames(water_lost.ratio)<-days

my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 98)
my_palette[1]<-c("#ffff00")

col_breaks = c(seq(-100,0,length=1),  # for yellow
               seq(0.000000001,.9999999999,length=48), # for blue
               seq(1,1.00000000001, length=1), # for white
               seq(1.000001,5,length=48), # for red
               seq(5.01,1000, length=1)) # extreme red

## Lets re-order the QTL in each data frame to reflect their physical order on chromsomes
water_lost.ratio<-water_lost.ratio[c(1:12,14,13),]

heatmap.2(water_lost.ratio,  notecol="black", col=my_palette, Colv="NA", Rowv="NA", trace="none", dendrogram = c("none"), density.info="none", key=F, breaks=col_breaks, symm=F,symkey=F,symbreaks=T)


## Now for WUE ratio
wue_ratio.fx<-r.all_total.qtl[r.all_total.qtl$raw_trait == 'wue_ratio_total',]
wue_ratio.ratio.calc<-c()
for(d in days){
  ## Extract wet
  temp1<-wue_ratio.fx[wue_ratio.fx$day == d, ]
  temp.w<-temp1[temp1$treatment == 'wet',c("marker", "qtl.name", "day", "signed.prop.var")]
  colnames(temp.w)[4]<-c("signed.prop.var.wet")
  ## Extract dry
  temp.d<-temp1[temp1$treatment == 'dry',c("marker", "qtl.name", "day", "signed.prop.var")]
  colnames(temp.d)[4]<-c("signed.prop.var.dry")
  ## Put in common data.frame
  temp2<-merge(temp.w, temp.d, by=c("marker", "qtl.name", "day"))
  ## Replace very small values with a reasonable value
  temp2[(temp2$signed.prop.var.wet < .2 & temp2$signed.prop.var.wet > 0),"signed.prop.var.wet"]<-c(.2)
  temp2[(temp2$signed.prop.var.wet > -.2 & temp2$signed.prop.var.wet < 0),"signed.prop.var.wet"]<-c(-.2)
  temp2[(temp2$signed.prop.var.dry < .2 & temp2$signed.prop.var.dry > 0),"signed.prop.var.dry"]<-c(.2)
  temp2[(temp2$signed.prop.var.dry > -.2 & temp2$signed.prop.var.dry < 0),"signed.prop.var.dry"]<-c(-.2)
  
  ## Take the ratio between treatment blocks (if ratio > 1 wet fx is greater, if < 1 dry fx is greater)
  temp2$ratio<-temp2$signed.prop.var.wet/temp2$signed.prop.var.dry
  ## If the value is negative and fx are almost zero remove the sign of small fx
  for(row in 1:nrow(temp2)){
    if((temp2[row, "ratio"] < 0) & (abs(temp2[row, "signed.prop.var.wet"]) == 0.2 | (abs(temp2 [row, "signed.prop.var.dry"]) == 0.2))){
      temp2[row, "ratio"]<-abs(temp2[row, "ratio"])
    } 
  }
  ## write to a new data.frame
  wue_ratio.ratio.calc<-rbind(wue_ratio.ratio.calc, temp2)
}

## Lets now make it matrix for plotting...
wue_ratio.ratio<-c()

for (d in days){
  temp1<-wue_ratio.ratio.calc[wue_ratio.ratio.calc$day == d, 'ratio']
  wue_ratio.ratio<-cbind(wue_ratio.ratio, temp1)
}

## Lets round the values to two decimal places
wue_ratio.ratio<-round(wue_ratio.ratio,2)

rownames(wue_ratio.ratio)<-unique(wue_ratio.ratio.calc$qtl.name)
colnames(wue_ratio.ratio)<-days

my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 98)
my_palette[1]<-c("#ffff00")

col_breaks = c(seq(-100,0,length=1),  # for yellow
               seq(0.000000001,.9999999999,length=48), # for blue
               seq(1,1.00000000001, length=1), # for white
               seq(1.000001,5,length=48), # for red
               seq(5.01,1000, length=1)) # extreme red

## Lets re-order the QTL in each data frame to reflect their physical order on chromsomes
wue_ratio.ratio<-wue_ratio.ratio[c(1:12,14,13),]

heatmap.2(wue_ratio.ratio,  notecol="black", col=my_palette, Colv="NA", Rowv="NA", trace="none", dendrogram = c("none"), density.info="none", key=F, breaks=col_breaks, symm=F,symkey=F,symbreaks=T)


## Now for TE fit total
te_fit.fx<-r.all_total.qtl[r.all_total.qtl$raw_trait == 'te_fit_total',]
te_fit.ratio.calc<-c()
for(d in days){
  ## Extract wet
  temp1<-te_fit.fx[te_fit.fx$day == d, ]
  temp.w<-temp1[temp1$treatment == 'wet',c("marker", "qtl.name", "day", "signed.prop.var")]
  colnames(temp.w)[4]<-c("signed.prop.var.wet")
  ## Extract dry
  temp.d<-temp1[temp1$treatment == 'dry',c("marker", "qtl.name", "day", "signed.prop.var")]
  colnames(temp.d)[4]<-c("signed.prop.var.dry")
  ## Put in common data.frame
  temp2<-merge(temp.w, temp.d, by=c("marker", "qtl.name", "day"))
  ## Replace very small values with a reasonable value
  temp2[(temp2$signed.prop.var.wet < .2 & temp2$signed.prop.var.wet > 0),"signed.prop.var.wet"]<-c(.2)
  temp2[(temp2$signed.prop.var.wet > -.2 & temp2$signed.prop.var.wet < 0),"signed.prop.var.wet"]<-c(-.2)
  temp2[(temp2$signed.prop.var.dry < .2 & temp2$signed.prop.var.dry > 0),"signed.prop.var.dry"]<-c(.2)
  temp2[(temp2$signed.prop.var.dry > -.2 & temp2$signed.prop.var.dry < 0),"signed.prop.var.dry"]<-c(-.2)
  
  ## Take the ratio between treatment blocks (if ratio > 1 wet fx is greater, if < 1 dry fx is greater)
  temp2$ratio<-temp2$signed.prop.var.wet/temp2$signed.prop.var.dry
  ## If the value is negative and fx are almost zero remove the sign of small fx
  for(row in 1:nrow(temp2)){
    if((temp2[row, "ratio"] < 0) & (abs(temp2[row, "signed.prop.var.wet"]) == 0.2 | (abs(temp2 [row, "signed.prop.var.dry"]) == 0.2))){
      temp2[row, "ratio"]<-abs(temp2[row, "ratio"])
    } 
  }
  ## write to a new data.frame
  te_fit.ratio.calc<-rbind(te_fit.ratio.calc, temp2)
}

## Lets now make it matrix for plotting...
te_fit.ratio<-c()

for (d in days){
  temp1<-te_fit.ratio.calc[te_fit.ratio.calc$day == d, 'ratio']
  te_fit.ratio<-cbind(te_fit.ratio, temp1)
}

## Lets round the values to two decimal places
te_fit.ratio<-round(te_fit.ratio,2)

rownames(te_fit.ratio)<-unique(te_fit.ratio.calc$qtl.name)
colnames(te_fit.ratio)<-days

my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 98)
my_palette[1]<-c("#ffff00")

col_breaks = c(seq(-100,0,length=1),  # for yellow
               seq(0.000000001,.9999999999,length=48), # for blue
               seq(1,1.00000000001, length=1), # for white
               seq(1.000001,5,length=48), # for red
               seq(5.01,1000, length=1)) # extreme red

## Lets re-order the QTL in each data frame to reflect their physical order on chromsomes
te_fit.ratio<-te_fit.ratio[c(1:12,14,13),]
heatmap.2(te_fit.ratio,  notecol="black", col=my_palette, Colv="NA", Rowv="NA", trace="none", dendrogram = c("none"), density.info="none", key=F, breaks=col_breaks, symm=F,symkey=F,symbreaks=T)



## Finally now te_residual
te_residual.fx<-r.all_total.qtl[r.all_total.qtl$raw_trait == 'te_residual_total',]
te_residual.ratio.calc<-c()
for(d in days){
  ## Extract wet
  temp1<-te_residual.fx[te_residual.fx$day == d, ]
  temp.w<-temp1[temp1$treatment == 'wet',c("marker", "qtl.name", "day", "signed.prop.var")]
  colnames(temp.w)[4]<-c("signed.prop.var.wet")
  ## Extract dry
  temp.d<-temp1[temp1$treatment == 'dry',c("marker", "qtl.name", "day", "signed.prop.var")]
  colnames(temp.d)[4]<-c("signed.prop.var.dry")
  ## Put in common data.frame
  temp2<-merge(temp.w, temp.d, by=c("marker", "qtl.name", "day"))
  ## Replace very small values with a reasonable value
  temp2[(temp2$signed.prop.var.wet < .2 & temp2$signed.prop.var.wet > 0),"signed.prop.var.wet"]<-c(.2)
  temp2[(temp2$signed.prop.var.wet > -.2 & temp2$signed.prop.var.wet < 0),"signed.prop.var.wet"]<-c(-.2)
  temp2[(temp2$signed.prop.var.dry < .2 & temp2$signed.prop.var.dry > 0),"signed.prop.var.dry"]<-c(.2)
  temp2[(temp2$signed.prop.var.dry > -.2 & temp2$signed.prop.var.dry < 0),"signed.prop.var.dry"]<-c(-.2)
  
  ## Take the ratio between treatment blocks (if ratio > 1 wet fx is greater, if < 1 dry fx is greater)
  temp2$ratio<-temp2$signed.prop.var.wet/temp2$signed.prop.var.dry
  ## If the value is negative and fx are almost zero remove the sign of small fx
  for(row in 1:nrow(temp2)){
    if((temp2[row, "ratio"] < 0) & (abs(temp2[row, "signed.prop.var.wet"]) == 0.2 | (abs(temp2 [row, "signed.prop.var.dry"]) == 0.2))){
      temp2[row, "ratio"]<-abs(temp2[row, "ratio"])
    } 
  }
  ## write to a new data.frame
  te_residual.ratio.calc<-rbind(te_residual.ratio.calc, temp2)
}

## Lets now make it matrix for plotting...
te_residual.ratio<-c()

for (d in days){
  temp1<-te_residual.ratio.calc[te_residual.ratio.calc$day == d, 'ratio']
  te_residual.ratio<-cbind(te_residual.ratio, temp1)
}

## Lets round the values to two decimal places
te_residual.ratio<-round(te_residual.ratio,2)

rownames(te_residual.ratio)<-unique(te_residual.ratio.calc$qtl.name)
colnames(te_residual.ratio)<-days

my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 98)
my_palette[1]<-c("#ffff00")

col_breaks = c(seq(-100,0,length=1),  # for yellow
               seq(0.000000001,.9999999999,length=48), # for blue
               seq(1,1.00000000001, length=1), # for white
               seq(1.000001,5,length=48), # for red
               seq(5.01,1000, length=1)) # extreme red

## Lets re-order the QTL in each data frame to reflect their physical order on chromsomes
te_residual.ratio<-te_residual.ratio[c(1:12,14,13),]

heatmap.2(te_residual.ratio,  notecol="black", col=my_palette, Colv="NA", Rowv="NA", trace="none", dendrogram = c("none"), density.info="none", key=F, breaks=col_breaks, symm=F,symkey=F,symbreaks=T)

## Lets make the figure
setwd(wue_results.qtl.dir)

pdf("relative_fx_of_qtl_across_treatment_blocks.pdf")

heatmap.2(sv_area.ratio,  notecol="black", col=my_palette, Colv="NA", Rowv="NA", trace="none", dendrogram = c("none"), density.info="none", key=F, breaks=col_breaks, symm=F,symkey=F,symbreaks=T, main="Plant size", xlab="Days after planting")

heatmap.2(water_lost.ratio,  notecol="black", col=my_palette, Colv="NA", Rowv="NA", trace="none", dendrogram = c("none"), density.info="none", key=F, breaks=col_breaks, symm=F,symkey=F,symbreaks=T, main="Water lost", xlab="Days after planting")

heatmap.2(wue_ratio.ratio,  notecol="black", col=my_palette, Colv="NA", Rowv="NA", trace="none", dendrogram = c("none"), density.info="none", key=F, breaks=col_breaks, symm=F,symkey=F,symbreaks=T, main="WUE ratio", xlab="Days after planting")

heatmap.2(te_fit.ratio,  notecol="black", col=my_palette, Colv="NA", Rowv="NA", trace="none", dendrogram = c("none"), density.info="none", key=F, breaks=col_breaks, symm=F,symkey=F,symbreaks=T, main="TE fit", xlab="Days after planting")

heatmap.2(te_residual.ratio,  notecol="black", col=my_palette, Colv="NA", Rowv="NA", trace="none", dendrogram = c("none"), density.info="none", key=F, breaks=col_breaks, symm=F,symkey=F,symbreaks=T, main="TE residual", xlab="Days after planting")

dev.off()


setwd(wue_results.qtl.dir)
save.image("fixed_marker_analysis.Rdata")


rm(list=ls())
