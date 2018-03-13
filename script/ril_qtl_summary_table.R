## This script loads in the results of the QTL analysis
## and condenses/unifies the results by collapsing all SNPs
## into a common SNP/QTL location

## This is done by searching for the most significant QTL on each chromosome (LOD)
## And calling all other SNPs within a 10 cM radius as the same QTL


library(ggplot2)


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

## Make a directory for isotope results
wue_results.isotope.dir<-paste(wue_results.dir, 'isotope/', sep="")

## If this directory doesn't exist the download is not complete
if (file.exists(wue_results.isotope.dir)){
  print("Directory exists!")
} else {
  dir.create(file.path(wue_results.isotope.dir))
}

## Make a directory path for the QTL summary tables
wue_results.qtl.dir<-paste(wue_results.dir, 'qtl_summary/', sep="")

## Plant size
wue_results.qtl.plant_size.dir<-paste(wue_results.qtl.dir, 'plant_size/', sep="")

## Cumulative
wue_results.qtl.plant_size.total.dir<-paste(wue_results.qtl.plant_size.dir, 'sv_area_total/', sep="")
wue_results.qtl.plant_size.total.dry.dir<-paste(wue_results.qtl.plant_size.total.dir, 'dry/', sep="")
wue_results.qtl.plant_size.total.wet.dir<-paste(wue_results.qtl.plant_size.total.dir, 'wet/', sep="")

## Rate
wue_results.qtl.plant_size.day.dir<-paste(wue_results.qtl.plant_size.dir, 'sv_area_day/', sep="")
wue_results.qtl.plant_size.day.dry.dir<-paste(wue_results.qtl.plant_size.day.dir, 'dry/', sep="")
wue_results.qtl.plant_size.day.wet.dir<-paste(wue_results.qtl.plant_size.day.dir, 'wet/', sep="")

## Water Use
wue_results.qtl.water_lost.dir<-paste(wue_results.qtl.dir, 'water_lost/', sep="")

## Cumulative
wue_results.qtl.water_lost.total.dir<-paste(wue_results.qtl.water_lost.dir, 'water_lost_total/', sep="")
wue_results.qtl.water_lost.total.dry.dir<-paste(wue_results.qtl.water_lost.total.dir, 'dry/', sep="")
wue_results.qtl.water_lost.total.wet.dir<-paste(wue_results.qtl.water_lost.total.dir, 'wet/', sep="")

## Rate
wue_results.qtl.water_lost.day.dir<-paste(wue_results.qtl.water_lost.dir, 'water_lost_day/', sep="")
wue_results.qtl.water_lost.day.dry.dir<-paste(wue_results.qtl.water_lost.day.dir, 'dry/', sep="")
wue_results.qtl.water_lost.day.wet.dir<-paste(wue_results.qtl.water_lost.day.dir, 'wet/', sep="")

## WUE ratio
wue_results.qtl.wue_ratio.dir<-paste(wue_results.qtl.dir, 'wue_ratio/', sep="")

## Cumulative
wue_results.qtl.wue_ratio.total.dir<-paste(wue_results.qtl.wue_ratio.dir, 'wue_ratio_total/', sep="")
wue_results.qtl.wue_ratio.total.dry.dir<-paste(wue_results.qtl.wue_ratio.total.dir, 'dry/', sep="")
wue_results.qtl.wue_ratio.total.wet.dir<-paste(wue_results.qtl.wue_ratio.total.dir, 'wet/', sep="")

## Rate
wue_results.qtl.wue_ratio.day.dir<-paste(wue_results.qtl.wue_ratio.dir, 'wue_ratio_day/', sep="")
wue_results.qtl.wue_ratio.day.dry.dir<-paste(wue_results.qtl.wue_ratio.day.dir, 'dry/', sep="")
wue_results.qtl.wue_ratio.day.wet.dir<-paste(wue_results.qtl.wue_ratio.day.dir, 'wet/', sep="")

## TE fit
wue_results.qtl.te_fit.dir<-paste(wue_results.qtl.dir, 'te_fit/', sep="")

## Cumulative
wue_results.qtl.te_fit.total.dir<-paste(wue_results.qtl.te_fit.dir, 'te_fit_total/', sep="")
wue_results.qtl.te_fit.total.dry.dir<-paste(wue_results.qtl.te_fit.total.dir, 'dry/', sep="")
wue_results.qtl.te_fit.total.wet.dir<-paste(wue_results.qtl.te_fit.total.dir, 'wet/', sep="")

## Rate
wue_results.qtl.te_fit.day.dir<-paste(wue_results.qtl.te_fit.dir, 'te_fit_day/', sep="")
wue_results.qtl.te_fit.day.dry.dir<-paste(wue_results.qtl.te_fit.day.dir, 'dry/', sep="")
wue_results.qtl.te_fit.day.wet.dir<-paste(wue_results.qtl.te_fit.day.dir, 'wet/', sep="")

## TE fit MA

## Cumulative
wue_results.qtl.te_ma_fit.total.dir<-paste(wue_results.qtl.te_fit.dir, 'te_ma_fit_total/', sep="")
wue_results.qtl.te_ma_fit.total.dry.dir<-paste(wue_results.qtl.te_ma_fit.total.dir, 'dry/', sep="")
wue_results.qtl.te_ma_fit.total.wet.dir<-paste(wue_results.qtl.te_ma_fit.total.dir, 'wet/', sep="")

## Rate
wue_results.qtl.te_ma_fit.day.dir<-paste(wue_results.qtl.te_fit.dir, 'te_ma_fit_day/', sep="")
wue_results.qtl.te_ma_fit.day.dry.dir<-paste(wue_results.qtl.te_ma_fit.day.dir, 'dry/', sep="")
wue_results.qtl.te_ma_fit.day.wet.dir<-paste(wue_results.qtl.te_ma_fit.day.dir, 'wet/', sep="")

## TE fit RMA

## Cumulative
wue_results.qtl.te_rma_fit.total.dir<-paste(wue_results.qtl.te_fit.dir, 'te_rma_fit_total/', sep="")
wue_results.qtl.te_rma_fit.total.dry.dir<-paste(wue_results.qtl.te_rma_fit.total.dir, 'dry/', sep="")
wue_results.qtl.te_rma_fit.total.wet.dir<-paste(wue_results.qtl.te_rma_fit.total.dir, 'wet/', sep="")

## Rate
## DID NOT ANALYZE RATE FOR RMA

## TE fit SMA

## Cumulative
wue_results.qtl.te_sma_fit.total.dir<-paste(wue_results.qtl.te_fit.dir, 'te_sma_fit_total/', sep="")
wue_results.qtl.te_sma_fit.total.dry.dir<-paste(wue_results.qtl.te_sma_fit.total.dir, 'dry/', sep="")
wue_results.qtl.te_sma_fit.total.wet.dir<-paste(wue_results.qtl.te_sma_fit.total.dir, 'wet/', sep="")

## Rate
## DID NOT ANALYZE RATE FOR SMA

## TE residual
wue_results.qtl.te_residual.dir<-paste(wue_results.qtl.dir, 'te_residual/', sep="")

## Cumulative
wue_results.qtl.te_residual.total.dir<-paste(wue_results.qtl.te_residual.dir, 'te_residual_total/', sep="")
wue_results.qtl.te_residual.total.dry.dir<-paste(wue_results.qtl.te_residual.total.dir, 'dry/', sep="")
wue_results.qtl.te_residual.total.wet.dir<-paste(wue_results.qtl.te_residual.total.dir, 'wet/', sep="")

## Rate
wue_results.qtl.te_residual.day.dir<-paste(wue_results.qtl.te_residual.dir, 'te_residual_day/', sep="")
wue_results.qtl.te_residual.day.dry.dir<-paste(wue_results.qtl.te_residual.day.dir, 'dry/', sep="")
wue_results.qtl.te_residual.day.wet.dir<-paste(wue_results.qtl.te_residual.day.dir, 'wet/', sep="")

## TE residual MA

## Cumulative
wue_results.qtl.te_ma_residual.total.dir<-paste(wue_results.qtl.te_residual.dir, 'te_ma_residual_total/', sep="")
wue_results.qtl.te_ma_residual.total.dry.dir<-paste(wue_results.qtl.te_ma_residual.total.dir, 'dry/', sep="")
wue_results.qtl.te_ma_residual.total.wet.dir<-paste(wue_results.qtl.te_ma_residual.total.dir, 'wet/', sep="")

## Rate
wue_results.qtl.te_ma_residual.day.dir<-paste(wue_results.qtl.te_residual.dir, 'te_ma_residual_day/', sep="")
wue_results.qtl.te_ma_residual.day.dry.dir<-paste(wue_results.qtl.te_ma_residual.day.dir, 'dry/', sep="")
wue_results.qtl.te_ma_residual.day.wet.dir<-paste(wue_results.qtl.te_ma_residual.day.dir, 'wet/', sep="")

## TE residual RMA

## Cumulative
wue_results.qtl.te_rma_residual.total.dir<-paste(wue_results.qtl.te_residual.dir, 'te_rma_residual_total/', sep="")
wue_results.qtl.te_rma_residual.total.dry.dir<-paste(wue_results.qtl.te_rma_residual.total.dir, 'dry/', sep="")
wue_results.qtl.te_rma_residual.total.wet.dir<-paste(wue_results.qtl.te_rma_residual.total.dir, 'wet/', sep="")

## Rate
## DID NOT ANALYZE RATE FOR RMA

## TE residual SMA

## Cumulative
wue_results.qtl.te_sma_residual.total.dir<-paste(wue_results.qtl.te_residual.dir, 'te_sma_residual_total/', sep="")
wue_results.qtl.te_sma_residual.total.dry.dir<-paste(wue_results.qtl.te_sma_residual.total.dir, 'dry/', sep="")
wue_results.qtl.te_sma_residual.total.wet.dir<-paste(wue_results.qtl.te_sma_residual.total.dir, 'wet/', sep="")

## Rate
## DID NOT ANALYZE RATE FOR SMA

## Isotope
wue_results.qtl.isotope.dir<-paste(wue_results.qtl.dir, 'isotope/', sep="")


##############################################################################################################################
# Load in QTL data from cumulative datasets
##############################################################################################################################

setwd(wue_results.qtl.te_fit.total.dir)
te_fit_total.qtl<-read.csv("te_fit_total_2_concatenated_summary_table.csv", header=F)
colnames(te_fit_total.qtl)<-c('marker','chr','pos','lod','prop.var','additive.fx','additive.fx_se','L.CI_marker','L.CI_chr','L.CI_pos','L.CI_lod','R.CI_marker','R.CI_chr','R.CI_pos','R.CI_lod','trait','treatment','exp','year','type')

setwd(wue_results.qtl.te_residual.total.dir)
te_residual_total.qtl<-read.csv("te_residual_total_2_concatenated_summary_table.csv", header=F)
colnames(te_residual_total.qtl)<-c('marker','chr','pos','lod','prop.var','additive.fx','additive.fx_se','L.CI_marker','L.CI_chr','L.CI_pos','L.CI_lod','R.CI_marker','R.CI_chr','R.CI_pos','R.CI_lod','trait','treatment','exp','year','type')


setwd(wue_results.qtl.te_ma_fit.total.dir)
te_ma_fit_total.qtl<-read.csv("te_ma_fit_total_2_concatenated_summary_table.csv", header=F)
colnames(te_ma_fit_total.qtl)<-c('marker','chr','pos','lod','prop.var','additive.fx','additive.fx_se','L.CI_marker','L.CI_chr','L.CI_pos','L.CI_lod','R.CI_marker','R.CI_chr','R.CI_pos','R.CI_lod','trait','treatment','exp','year','type')

setwd(wue_results.qtl.te_ma_residual.total.dir)
te_ma_residual_total.qtl<-read.csv("te_ma_residual_total_2_concatenated_summary_table.csv", header=F)
colnames(te_ma_residual_total.qtl)<-c('marker','chr','pos','lod','prop.var','additive.fx','additive.fx_se','L.CI_marker','L.CI_chr','L.CI_pos','L.CI_lod','R.CI_marker','R.CI_chr','R.CI_pos','R.CI_lod','trait','treatment','exp','year','type')


setwd(wue_results.qtl.wue_ratio.total.dir)
wue_ratio_total.qtl<-read.csv("wue_ratio_total_2_concatenated_summary_table.csv", header=F)
colnames(wue_ratio_total.qtl)<-c('marker','chr','pos','lod','prop.var','additive.fx','additive.fx_se','L.CI_marker','L.CI_chr','L.CI_pos','L.CI_lod','R.CI_marker','R.CI_chr','R.CI_pos','R.CI_lod','trait','treatment','exp','year','type')

setwd(wue_results.qtl.plant_size.total.dir)
sv_area_total.qtl<-read.csv("sv_area_total_2_concatenated_summary_table.csv", header=F)
colnames(sv_area_total.qtl)<-c('marker','chr','pos','lod','prop.var','additive.fx','additive.fx_se','L.CI_marker','L.CI_chr','L.CI_pos','L.CI_lod','R.CI_marker','R.CI_chr','R.CI_pos','R.CI_lod','trait','treatment','exp','year','type')

## Lets remove days before day 17
dap_i<-c()
for(i in 1:nrow(sv_area_total.qtl)){
  x<-strsplit(as.character(sv_area_total.qtl[i,'trait']), '_')[[1]]
  d<-x[length(x)-1]
  dap_i<-c(dap_i, d)
}

## Convert dap_i vector to numeric
dap_i<-as.numeric(dap_i)

sv_area_total.qtl<-cbind(sv_area_total.qtl, dap_i)
colnames(sv_area_total.qtl)
sv_area_total.qtl<-sv_area_total.qtl[sv_area_total.qtl$dap_i > 16,]
sv_area_total.qtl<-sv_area_total.qtl[,-c(21)]

## Water use size
setwd(wue_results.qtl.water_lost.total.dir)
water_total.qtl<-read.csv("water_lost_total_2_concatenated_summary_table.csv", header=F)
colnames(water_total.qtl)<-c('marker','chr','pos','lod','prop.var','additive.fx','additive.fx_se','L.CI_marker','L.CI_chr','L.CI_pos','L.CI_lod','R.CI_marker','R.CI_chr','R.CI_pos','R.CI_lod','trait','treatment','exp','year','type')

## Lets not include the isotopes
#setwd(wue_results.qtl.isotope.dir)
#isotope.qtl<-read.csv("stable_isotope4_concatenated_summary_table.csv", header=F)
#colnames(isotope.qtl)<-c('marker','chr','pos','lod','prop.var','additive.fx','additive.fx_se','L.CI_marker','L.CI_chr','L.CI_pos','L.CI_lod','R.CI_marker','R.CI_chr','R.CI_pos','R.CI_lod','trait','treatment','exp','year','type')

##############################################################################################################################
## Load in genetic map for plotting
##############################################################################################################################

setwd(wue_data.dir)
map<-read.table("GBS_map_A10xB100_v0.96.csv", header=T, sep=",")
chrs<-t(map[1,2:ncol(map)])
pos<-t(map[2,2:ncol(map)])

genome<-cbind(chrs, pos)
colnames(genome)<-c('chrs', 'pos')
genome<-as.data.frame(genome)
genome$chrs<-as.numeric(as.character(genome$chrs))
genome$pos<-as.numeric(as.character(genome$pos))

c1.max<-max(genome[genome$chrs == 1, 'pos'])
c2.max<-max(genome[genome$chrs == 2, 'pos'])
c3.max<-max(genome[genome$chrs == 3, 'pos'])
c4.max<-max(genome[genome$chrs == 4, 'pos'])
c5.max<-max(genome[genome$chrs == 5, 'pos'])
c6.max<-max(genome[genome$chrs == 6, 'pos'])
c7.max<-max(genome[genome$chrs == 7, 'pos'])
c8.max<-max(genome[genome$chrs == 8, 'pos'])
c9.max<-max(genome[genome$chrs == 9, 'pos'])

blank_data<-data.frame(chr=c('1','1','2','2','3','3','4','4','5','5','6','6','7','7','8','8','9','9'), x=c(0,c1.max,0,c2.max,0,c3.max,0,c4.max,0,c5.max,0,c6.max,0,c7.max,0,c8.max,0,c9.max), y=0)


##############################################################################################################################
## Lets condense all QTL calculated on a cumulative basis
##############################################################################################################################

## Remove all QTL detected from the difference (GXE) traits
## lm fit
r.te_fit_total.qtl <- te_fit_total.qtl[te_fit_total.qtl$type == 'raw',]
r.uni.te_fit_total.qtl <- unify_marker(r.te_fit_total.qtl)
r.cond.te_fit_total.qtl <- condense_qtl(r.te_fit_total.qtl)

uni.te_fit_total.qtl <- unify_marker(te_fit_total.qtl)
cond.te_fit_total.qtl <- condense_qtl(te_fit_total.qtl)

setwd(wue_results.qtl.te_fit.dir)

make_qtl_common_plot(r.te_fit_total.qtl, "r.te_fit_total.qtl.pdf")
make_qtl_common_plot(r.uni.te_fit_total.qtl, "r.uni.te_fit_total.qtl.pdf")

make_qtl_common_plot_diff(te_fit_total.qtl, "te_fit_total.qtl.pdf")
make_qtl_common_plot_diff(uni.te_fit_total.qtl, "uni.te_fit_total.qtl.pdf")

## lm residual

r.te_residual_total.qtl <- te_residual_total.qtl[te_residual_total.qtl$type == 'raw',]
r.uni.te_residual_total.qtl <- unify_marker(r.te_residual_total.qtl)
r.cond.te_residual_total.qtl <- condense_qtl(r.te_residual_total.qtl)

uni.te_residual_total.qtl <- unify_marker(te_residual_total.qtl)
cond.te_residual_total.qtl <- condense_qtl(te_residual_total.qtl)

setwd(wue_results.qtl.te_residual.dir)

make_qtl_common_plot(r.te_residual_total.qtl, "r.te_residual_total.qtl.pdf")
make_qtl_common_plot(r.uni.te_residual_total.qtl, "r.uni.te_residual_total.qtl.pdf")

make_qtl_common_plot_diff(te_residual_total.qtl, "te_residual_total.qtl.pdf")
make_qtl_common_plot_diff(uni.te_residual_total.qtl, "uni.te_residual_total.qtl.pdf")

## MAJOR AXIS
## MA fit total
r.te_ma_fit_total.qtl <- te_ma_fit_total.qtl[te_ma_fit_total.qtl$type == 'raw',]
r.uni.te_ma_fit_total.qtl <- unify_marker(r.te_ma_fit_total.qtl)
r.cond.te_ma_fit_total.qtl <- condense_qtl(r.te_ma_fit_total.qtl)

uni.te_ma_fit_total.qtl <- unify_marker(te_ma_fit_total.qtl)
cond.te_ma_fit_total.qtl <- condense_qtl(te_ma_fit_total.qtl)

setwd(wue_results.qtl.te_fit.dir)
make_qtl_common_plot(r.te_ma_fit_total.qtl, "r.te_ma_fit_total.qtl.pdf")
make_qtl_common_plot(r.uni.te_ma_fit_total.qtl, "r.uni.te_ma_fit_total.qtl.pdf")

make_qtl_common_plot_diff(te_ma_fit_total.qtl, "te_ma_fit_total.qtl.pdf")
make_qtl_common_plot_diff(uni.te_ma_fit_total.qtl, "uni.te_ma_fit_total.qtl.pdf")

# MA residual total
r.te_ma_residual_total.qtl <- te_ma_residual_total.qtl[te_ma_residual_total.qtl$type == 'raw',]
r.uni.te_ma_residual_total.qtl <- unify_marker(r.te_ma_residual_total.qtl)
r.cond.te_ma_residual_total.qtl <- condense_qtl(r.te_ma_residual_total.qtl)

uni.te_ma_residual_total.qtl <- unify_marker(te_ma_residual_total.qtl)
cond.te_ma_residual_total.qtl <- condense_qtl(te_ma_residual_total.qtl)

setwd(wue_results.qtl.te_fit.dir)

make_qtl_common_plot(r.te_ma_residual_total.qtl, "r.te_ma_residual_total.qtl.pdf")
make_qtl_common_plot(r.uni.te_ma_residual_total.qtl, "r.uni.te_ma_residual_total.qtl.pdf")

make_qtl_common_plot_diff(te_ma_residual_total.qtl, "te_ma_residual_total.qtl.pdf")
make_qtl_common_plot_diff(uni.te_ma_residual_total.qtl, "uni.te_ma_residual_total.qtl.pdf")


## WUE
r.wue_ratio_total.qtl <- wue_ratio_total.qtl[wue_ratio_total.qtl$type == 'raw',]
r.uni.wue_ratio_total.qtl <- unify_marker(r.wue_ratio_total.qtl)
r.cond.wue_ratio_total.qtl <- condense_qtl(r.wue_ratio_total.qtl)

uni.wue_ratio_total.qtl <- unify_marker(wue_ratio_total.qtl)
cond.wue_ratio_total.qtl <- condense_qtl(wue_ratio_total.qtl)

setwd(wue_results.qtl.wue_ratio.dir)

make_qtl_common_plot(r.wue_ratio_total.qtl, "r.wue_ratio_total.qtl.pdf")
make_qtl_common_plot(r.uni.wue_ratio_total.qtl, "r.uni.wue_ratio_total.qtl.pdf")

make_qtl_common_plot_diff(wue_ratio_total.qtl, "wue_ratio_total.qtl.pdf")
make_qtl_common_plot_diff(uni.wue_ratio_total.qtl, "uni.wue_ratio_total.qtl.pdf")

## sv area total
r.sv_area_total.qtl<-sv_area_total.qtl[sv_area_total.qtl$type == 'raw',]
r.uni.sv_area_total.qtl<-unify_marker(r.sv_area_total.qtl)
r.cond.sv_area_total.qtl<-condense_qtl(r.sv_area_total.qtl)

uni.sv_area_total.qtl<-unify_marker(sv_area_total.qtl)
cond.sv_area_total.qtl<-condense_qtl(sv_area_total.qtl)

setwd(wue_results.qtl.plant_size.dir)

make_qtl_common_plot(r.sv_area_total.qtl, "r.sv_area_total.qtl.pdf")
make_qtl_common_plot(r.uni.sv_area_total.qtl, "r.uni.sv_area_total.qtl.pdf")

make_qtl_common_plot_diff(sv_area_total.qtl, "sv_area_total.qtl.pdf")
make_qtl_common_plot_diff(uni.sv_area_total.qtl, "uni.sv_area_total.qtl.pdf")

## water total
r.water_total.qtl<-water_total.qtl[water_total.qtl$type == 'raw',]
r.uni.water_total.qtl<-unify_marker(r.water_total.qtl)
r.cond.water_total.qtl<-condense_qtl(r.water_total.qtl)

uni.water_total.qtl<-unify_marker(water_total.qtl)
cond.water_total.qtl<-condense_qtl(water_total.qtl)

setwd(wue_results.qtl.water_lost.dir)

make_qtl_common_plot(r.water_total.qtl, "r.water_total.qtl.pdf")
make_qtl_common_plot(r.uni.water_total.qtl, "r.uni.water_total.qtl.pdf")

make_qtl_common_plot_diff(water_total.qtl, "water_total.qtl.pdf")
make_qtl_common_plot_diff(uni.water_total.qtl, "uni.water_total.qtl.pdf")

## isotope
#r.isotope.qtl <- isotope.qtl[isotope.qtl$type == 'raw',]
#r.uni.isotope.qtl <- unify_marker(r.isotope.qtl)
#r.cond.isotope.qtl <- condense_qtl(r.isotope.qtl)

#setwd(wue_results.qtl.isotope.dir)

#make_qtl_common_plot(r.isotope.qtl, "r.isotope.qtl.pdf")
#make_qtl_common_plot(r.uni.isotope.qtl, "r.uni.isotope.qtl.pdf")

#make_qtl_common_plot_diff(isotope.qtl, "isotope.qtl.pdf")
#make_qtl_common_plot_diff(r.uni.isotope.qtl, "uni.isotope.qtl.pdf")



## Lets combine all traits
## Note that this does not include daily rates

## Lets not include the isotope data 
#all_total.qtl <- rbind(te_fit_total.qtl, te_residual_total.qtl, wue_ratio_total.qtl, sv_area_total.qtl, water_total.qtl, isotope.qtl)
all_total.qtl <- rbind(te_fit_total.qtl, te_residual_total.qtl, wue_ratio_total.qtl, sv_area_total.qtl, water_total.qtl)

r.all_total.qtl <- all_total.qtl[all_total.qtl$type == 'raw',]

## Used two different ways to summarize data
r.uni.all_total.qtl<-unify_marker(r.all_total.qtl)
r.cond.all_total.qtl<-condense_qtl(r.all_total.qtl)

uni.all_total.qtl<-unify_marker(all_total.qtl)
cond.all_total.qtl<-condense_qtl(all_total.qtl)

r.uni.all_total.qtl<-uni.all_total.qtl[uni.all_total.qtl$type == 'raw',]

# We get 23 unique QTL for all traits
cond.all_total.qtl$qtl_count<-as.numeric(as.character(cond.all_total.qtl$qtl_count))

# Lets see how many times each QTL is observed
sort(table(uni.all_total.qtl$marker))
sort(table(r.uni.all_total.qtl$marker))


##############################################################################################################################
## Load in QTL data from daily datasets
##############################################################################################################################

## TE fit
setwd(wue_results.qtl.te_fit.day.dir)
te_fit_day.qtl<-read.csv("te_fit_day_2_concatenated_summary_table.csv", header=F)
colnames(te_fit_day.qtl)<-c('marker','chr','pos','lod','prop.var','additive.fx','additive.fx_se','L.CI_marker','L.CI_chr','L.CI_pos','L.CI_lod','R.CI_marker','R.CI_chr','R.CI_pos','R.CI_lod','trait','treatment','exp','year','type')

## TE MA fit
setwd(wue_results.qtl.te_ma_fit.day.dir)
te_ma_fit_day.qtl<-read.csv("te_ma_fit_day_2_concatenated_summary_table.csv", header=F)
colnames(te_ma_fit_day.qtl)<-c('marker','chr','pos','lod','prop.var','additive.fx','additive.fx_se','L.CI_marker','L.CI_chr','L.CI_pos','L.CI_lod','R.CI_marker','R.CI_chr','R.CI_pos','R.CI_lod','trait','treatment','exp','year','type')

## TE MA residual
setwd(wue_results.qtl.te_ma_residual.day.dir)
te_ma_residual_day.qtl<-read.csv("te_ma_residual_day_2_concatenated_summary_table.csv", header=F)
colnames(te_ma_residual_day.qtl)<-c('marker','chr','pos','lod','prop.var','additive.fx','additive.fx_se','L.CI_marker','L.CI_chr','L.CI_pos','L.CI_lod','R.CI_marker','R.CI_chr','R.CI_pos','R.CI_lod','trait','treatment','exp','year','type')

## TE residual
setwd(wue_results.qtl.te_residual.day.dir)
te_residual_day.qtl<-read.csv("te_residual_day_2_concatenated_summary_table.csv", header=F)
colnames(te_residual_day.qtl)<-c('marker','chr','pos','lod','prop.var','additive.fx','additive.fx_se','L.CI_marker','L.CI_chr','L.CI_pos','L.CI_lod','R.CI_marker','R.CI_chr','R.CI_pos','R.CI_lod','trait','treatment','exp','year','type')

## WUE ratio
setwd(wue_results.qtl.wue_ratio.day.dir)
wue_ratio_day.qtl<-read.csv("wue_ratio_day_2_concatenated_summary_table.csv", header=F)
colnames(wue_ratio_day.qtl)<-c('marker','chr','pos','lod','prop.var','additive.fx','additive.fx_se','L.CI_marker','L.CI_chr','L.CI_pos','L.CI_lod','R.CI_marker','R.CI_chr','R.CI_pos','R.CI_lod','trait','treatment','exp','year','type')

## Plant size
setwd(wue_results.qtl.plant_size.day.dir)
sv_area_day.qtl<-read.csv("sv_area_day_2_concatenated_summary_table.csv", header=F)
colnames(sv_area_day.qtl)<-c('marker','chr','pos','lod','prop.var','additive.fx','additive.fx_se','L.CI_marker','L.CI_chr','L.CI_pos','L.CI_lod','R.CI_marker','R.CI_chr','R.CI_pos','R.CI_lod','trait','treatment','exp','year','type')

## Water lost
setwd(wue_results.qtl.water_lost.day.dir)
water_day.qtl<-read.csv("water_lost_day_2_concatenated_summary_table.csv", header=F)
colnames(water_day.qtl)<-c('marker','chr','pos','lod','prop.var','additive.fx','additive.fx_se','L.CI_marker','L.CI_chr','L.CI_pos','L.CI_lod','R.CI_marker','R.CI_chr','R.CI_pos','R.CI_lod','trait','treatment','exp','year','type')


# Lets try condensing each trait set independently
# Remove all QTL detected from the difference (GXE) traits
r.te_fit_day.qtl<-te_fit_day.qtl[te_fit_day.qtl$type == 'raw',]
r.uni.te_fit_day.qtl<-unify_marker(r.te_fit_day.qtl)
r.cond.te_fit_day.qtl<-condense_qtl(r.te_fit_day.qtl)

uni.te_fit_day.qtl<-unify_marker(te_fit_day.qtl)
cond.te_fit_day.qtl<-condense_qtl(te_fit_day.qtl)

setwd(wue_results.qtl.te_fit.dir)

make_qtl_common_plot(r.te_fit_day.qtl, "r.te_fit_day.qtl.pdf")
make_qtl_common_plot(r.uni.te_fit_day.qtl, "r.uni.te_fit_day.qtl.pdf")

make_qtl_common_plot_diff(te_fit_day.qtl, "te_fit_day.qtl.pdf")
make_qtl_common_plot_diff(uni.te_fit_day.qtl, "uni.te_fit_day.qtl.pdf")

r.te_residual_day.qtl<-te_residual_day.qtl[te_residual_day.qtl$type == 'raw',]
r.uni.te_residual_day.qtl<-unify_marker(r.te_residual_day.qtl)
r.cond.te_residual_day.qtl<-condense_qtl(r.te_residual_day.qtl)

uni.te_residual_day.qtl<-unify_marker(te_residual_day.qtl)
cond.te_residual_day.qtl<-condense_qtl(te_residual_day.qtl)

setwd(wue_results.qtl.te_residual.dir)

make_qtl_common_plot(r.te_residual_day.qtl, "r.te_residual_day.qtl.pdf")
make_qtl_common_plot(r.uni.te_residual_day.qtl, "r.uni.te_residual_day.qtl.pdf")

make_qtl_common_plot_diff(te_residual_day.qtl, "te_residual_day.qtl.pdf")
make_qtl_common_plot_diff(uni.te_residual_day.qtl, "uni.te_residual_day.qtl.pdf")

### MAJOR AXIS
r.te_ma_fit_day.qtl<-te_fit_day.qtl[te_ma_fit_day.qtl$type == 'raw',]
r.uni.te_ma_fit_day.qtl<-unify_marker(r.te_ma_fit_day.qtl)
r.cond.te_ma_fit_day.qtl<-condense_qtl(r.te_ma_fit_day.qtl)

uni.te_ma_fit_day.qtl<-unify_marker(te_ma_fit_day.qtl)
cond.te_ma_fit_day.qtl<-condense_qtl(te_ma_fit_day.qtl)

setwd(wue_results.qtl.te_fit.dir)

make_qtl_common_plot(r.te_ma_fit_day.qtl, "r.te_ma_fit_day.qtl.pdf")
make_qtl_common_plot(r.uni.te_ma_fit_day.qtl, "r.uni.te_ma_fit_day.qtl.pdf")

make_qtl_common_plot_diff(te_ma_fit_day.qtl, "te_ma_fit_day.qtl.pdf")
make_qtl_common_plot_diff(uni.te_ma_fit_day.qtl, "uni.te_ma_fit_day.qtl.pdf")


r.te_ma_residual_day.qtl<-te_ma_residual_day.qtl[te_ma_residual_day.qtl$type == 'raw',]
r.uni.te_ma_residual_day.qtl<-unify_marker(r.te_ma_residual_day.qtl)
r.cond.te_ma_residual_day.qtl<-condense_qtl(r.te_ma_residual_day.qtl)

uni.te_ma_residual_day.qtl<-unify_marker(te_ma_residual_day.qtl)
cond.te_ma_residual_day.qtl<-condense_qtl(te_ma_residual_day.qtl)

setwd(wue_results.qtl.te_residual.dir)

make_qtl_common_plot(r.te_ma_residual_day.qtl, "r.te_ma_residual_day.qtl.pdf")
make_qtl_common_plot(r.uni.te_ma_residual_day.qtl, "r.uni.te_ma_residual_day.qtl.pdf")

make_qtl_common_plot_diff(te_ma_residual_day.qtl, "te_ma_residual_day.qtl.pdf")
make_qtl_common_plot_diff(uni.te_ma_residual_day.qtl, "uni.te_ma_residual_day.qtl.pdf")


r.sv_area_day.qtl<-sv_area_day.qtl[sv_area_day.qtl$type == 'raw',]
r.uni.sv_area_day.qtl<-unify_marker(r.sv_area_day.qtl)
r.cond.sv_area_day.qtl<-condense_qtl(r.sv_area_day.qtl)

uni.sv_area_day.qtl<-unify_marker(sv_area_day.qtl)
cond.sv_area_day.qtl<-condense_qtl(sv_area_day.qtl)

setwd(wue_results.qtl.plant_size.dir)

make_qtl_common_plot(r.sv_area_day.qtl, "r.sv_area_day.qtl.pdf")
make_qtl_common_plot(r.uni.sv_area_day.qtl, "r.uni.sv_area_day.qtl.pdf")

make_qtl_common_plot_diff(sv_area_day.qtl, "sv_area_day.qtl.pdf")
make_qtl_common_plot_diff(uni.sv_area_day.qtl, "uni.sv_area_day.qtl.pdf")


r.water_day.qtl<-water_day.qtl[water_day.qtl$type == 'raw',]
r.uni.water_day.qtl<-unify_marker(r.water_day.qtl)
r.cond.water_day.qtl<-condense_qtl(r.water_day.qtl)

uni.water_day.qtl<-unify_marker(water_day.qtl)
cond.water_day.qtl<-condense_qtl(water_day.qtl)

setwd(wue_results.qtl.water_lost.dir)

make_qtl_common_plot(r.water_day.qtl, "r.water_day.qtl.pdf")
make_qtl_common_plot(r.uni.water_day.qtl, "r.uni.water_day.qtl.pdf")

make_qtl_common_plot_diff(water_day.qtl, "water_day.qtl.pdf")
make_qtl_common_plot_diff(uni.water_day.qtl, "uni.water_day.qtl.pdf")


##############################################################################################################################
## Combine all traits
##############################################################################################################################


## Lets combine all traits and condense them together
#all_day.qtl<-rbind(te_fit_day.qtl, te_residual_day.qtl, wue_ratio_day.qtl, sv_area_day.qtl, water_day.qtl, isotope.qtl)

## Lets not include the isotope data
all_day.qtl<-rbind(te_fit_day.qtl, te_residual_day.qtl, wue_ratio_day.qtl, sv_area_day.qtl, water_day.qtl)

r.all_day.qtl<-all_day.qtl[all_day.qtl$type == 'raw',]

## Get a list of traits which we do not want to include
#traits_not_te<-c("gN_m2_BP14", "d13C_z_BP14", "gC_m2_BP14", "gN_m2_z_BP14", "gC_m2_z_BP14")

#r.all_day.qtl<-r.all_day.qtl[!r.all_day.qtl$trait %in% traits_not_te, ]
r.all_day.qtl$marker<-as.character(r.all_day.qtl$marker)

#all_day.qtl<-all_day.qtl[!all_day.qtl$trait %in% traits_not_te, ]
all_day.qtl$marker<-as.character(all_day.qtl$marker)

##############################################################################################################################
# Lets try combining cumulative and daily QTL results together and then unify/condensing
##############################################################################################################################

## Raw traits only (no GxE)
r.all_day.qtl<-r.all_day.qtl[r.all_day.qtl$trait != 'd13C_BP14',]
r.all.qtl<-rbind(r.all_total.qtl, r.all_day.qtl)

## "Unify" the raw traits (no GxE)
r.uni.all.qtl<-unify_marker(r.all.qtl)
r.cond.all.qtl<-condense_qtl(r.all.qtl)

all_day.qtl<-all_day.qtl[all_day.qtl$trait != 'd13C_BP14',]
all.qtl<-rbind(all_total.qtl, all_day.qtl)


all.qtl.with_diff_ratio<-all.qtl
## Lets remove the GxE ratio QTL because they can be very noisy.
all.qtl<-all.qtl[all.qtl$type != 'comp_ratio',]

uni.all.qtl<-unify_marker(all.qtl)
cond.all.qtl<-condense_qtl(all.qtl)

uni.all.qtl.with_diff_ratio<-unify_marker(all.qtl.with_diff_ratio)
cond.all.qtl.with_diff_ratio<-unify_marker(all.qtl.with_diff_ratio)

setwd(wue_results.qtl.dir)

make_qtl_common_plot(r.all.qtl, "r.all.qtl.pdf")
make_qtl_common_plot_diff(uni.all.qtl, "uni.all.qtl.pdf")

all_day.qtl<-all_day.qtl[all_day.qtl$trait != 'd13C_BP14',]
all.qtl<-rbind(all_total.qtl, all_day.qtl)

#traits_not_te<-c("gN_m2_BP14", "d13C_z_BP14", "gC_m2_BP14", "gN_m2_z_BP14", "gC_m2_z_BP14","CN_BP14", "CN_z_BP14","d13C_BP14")
#all.qtl<-all.qtl[!all.qtl$trait %in% traits_not_te, ]

setwd(wue_results.qtl.dir)

make_qtl_common_plot(r.all.qtl, "r.all.qtl.pdf")
make_qtl_common_plot(r.uni.all.qtl, "r.uni.all.qtl.pdf")

setwd(wue_results.qtl.dir)

make_qtl_common_plot_diff(all.qtl, "all.qtl.pdf")
make_qtl_common_plot_diff(uni.all.qtl.with_diff_ratio, "uni.all.qtl.with_diff_ratio.pdf")

r.uni.all.qtl<-uni.all.qtl[uni.all.qtl$type == 'raw',]

## Lets see how many times each QTL is observed
length(sort(table(r.uni.all.qtl$marker)))
sort(table(r.uni.all.qtl$marker))



## Lets ignore those QTL detected 10 times or greater
found_X10<-names(sort(table(r.uni.all.qtl$marker))[sort(table(r.uni.all.qtl$marker)) > 9])
unique(r.uni.all.qtl[r.uni.all.qtl$marker %in% found_X10,'marker'])


## We get 14 QTL locations
known_qtl_markers<-r.uni.all.qtl[r.uni.all.qtl$marker %in% found_X10,c(1,2,3,4)]
known_qtl_markers<-aggregate(known_qtl_markers, by=list(known_qtl_markers$marker), FUN=mean)
known_qtl_markers<-known_qtl_markers[,-c(2)]
colnames(known_qtl_markers)<-c('marker', 'chr', 'pos', 'mean_lod')

## Write out these markers for the fixed model
write.csv(known_qtl_markers, file="14_markers_for_fixed_model.csv", quote=F, row.names=F)

## Save image for next step
save.image("qtl_summary_4_qtl_trait_matrix.Rdata")

rm(list=ls())


