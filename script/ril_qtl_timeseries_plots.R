library(qtl)
library(funqtl)



###########################################################################
# QTL Thru time... SLOD, MLOD, etc.
###########################################################################
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



##### PLANT SIZE

## Plant size dry and wet

setwd(wue_results.qtl.plant_size.total.dry.dir)
load('timeseries_dry_cross.object.raw.Rdata')

# Re-name scanone and stepwise QTL results
sv_area.dry_out.F<-out.F
sv_area.dry_slod<-refqtlslod_dry_slod
sv_area.dry_max.perm.F_slod<-max.perm.F

setwd(wue_results.qtl.plant_size.total.wet.dir)
load('timeseries_wet_cross.object.raw.Rdata')

# Re-name scanone and stepwise QTL results
sv_area.wet_out.F<-out.F
sv_area.wet_slod<-refqtlslod_wet_slod
sv_area.wet_max.perm.F_slod<-max.perm.F

##### PLANT WATER USE

## Plant water use
setwd(wue_results.qtl.water_lost.total.dry.dir)
load('timeseries_dry_cross.object.raw.Rdata')

# Re-name scanone and stepwise QTL results
water_lost.dry_out.F<-out.F
water_lost.dry_slod<-refqtlslod_dry_slod
water_lost.dry_max.perm.F_slod<-max.perm.F

setwd(wue_results.qtl.water_lost.total.wet.dir)
load('timeseries_wet_cross.object.raw.Rdata')

# Re-name scanone and stepwise QTL results
water_lost.wet_out.F<-out.F
water_lost.wet_slod<-refqtlslod_wet_slod
water_lost.wet_max.perm.F_slod<-max.perm.F

##### WUE RATIO

## Plant WUE ratio
setwd(wue_results.qtl.wue_ratio.total.dry.dir)
load('timeseries_dry_cross.object.raw.Rdata')


# Re-name scanone and stepwise QTL results
wue_ratio.dry_out.F<-out.F
wue_ratio.dry_slod<-refqtlslod_dry_slod
wue_ratio.dry_max.perm.F_slod<-max.perm.F

setwd(wue_results.qtl.wue_ratio.total.wet.dir)
load('timeseries_wet_cross.object.raw.Rdata')

# Re-name scanone and stepwise QTL results
wue_ratio.wet_out.F<-out.F
wue_ratio.wet_slod<-refqtlslod_wet_slod
wue_ratio.wet_max.perm.F_slod<-max.perm.F


##### TE MODEL FIT

## TE model fit
setwd(wue_results.qtl.te_fit.total.dry.dir)
load('timeseries_dry_cross.object.raw.Rdata')


# Re-name scanone and stepwise QTL results
te_fit.dry_out.F<-out.F
te_fit.dry_slod<-refqtlslod_dry_slod
te_fit.dry_max.perm.F_slod<-max.perm.F

setwd(wue_results.qtl.te_fit.total.wet.dir)
load('timeseries_wet_cross.object.raw.Rdata')

# Re-name scanone and stepwise QTL results
te_fit.wet_out.F<-out.F
te_fit.wet_slod<-refqtlslod_wet_slod
te_fit.wet_max.perm.F_slod<-max.perm.F


##### TE MODEL RESIDUAL

## TE model residual
setwd(wue_results.qtl.te_residual.total.dry.dir)
load('timeseries_dry_cross.object.raw.Rdata')


# Re-name scanone and stepwise QTL results
te_residual.dry_out.F<-out.F
te_residual.dry_slod<-refqtlslod_dry_slod
te_residual.dry_max.perm.F_slod<-max.perm.F

setwd(wue_results.qtl.te_residual.total.wet.dir)
load('timeseries_wet_cross.object.raw.Rdata')

# Re-name scanone and stepwise QTL results
te_residual.wet_out.F<-out.F
te_residual.wet_slod<-refqtlslod_wet_slod
te_residual.wet_max.perm.F_slod<-max.perm.F

#######################################################################
# Lets make the figures
#######################################################################

## Lets make Figure 6

setwd(wue_results.qtl.dir)

pdf("slod_qtl_cumulative_traits_scanone.pdf", width=8, height = 3)
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


## Lets make Figure S20

setwd(wue_results.qtl.dir)
pdf("slod_qtl_cumulative_traits_stepwise.pdf", width=8, height = 3)
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


setwd(wue_results.qtl.dir)
save.image("qtl_timeseries_plots_cumulative.Rdata")


############### NOTE WE ARE OVER WRITING THE NAMES OF THESE VARIABLES!!!


############ NOW LETS WORK WITH DAILY/RATE
setwd(wue_results.qtl.plant_size.day.dry.dir)
load('timeseries_dry_cross.object.raw.Rdata')

# Re-name scanone and stepwise QTL results
sv_area.dry_out.F<-out.F
sv_area.dry_slod<-refqtlslod_dry_slod
sv_area.dry_max.perm.F_slod<-max.perm.F

setwd(wue_results.qtl.plant_size.day.wet.dir)
load('timeseries_wet_cross.object.raw.Rdata')

# Re-name scanone and stepwise QTL results
sv_area.wet_out.F<-out.F
sv_area.wet_slod<-refqtlslod_wet_slod
sv_area.wet_max.perm.F_slod<-max.perm.F


##### PLANT WATER USE

## Plant water use
setwd(wue_results.qtl.water_lost.day.dry.dir)
load('timeseries_dry_cross.object.raw.Rdata')

# Re-name scanone and stepwise QTL results
water_lost.dry_out.F<-out.F
water_lost.dry_slod<-refqtlslod_dry_slod
water_lost.dry_max.perm.F_slod<-max.perm.F

setwd(wue_results.qtl.water_lost.day.wet.dir)
load('timeseries_wet_cross.object.raw.Rdata')

# Re-name scanone and stepwise QTL results
water_lost.wet_out.F<-out.F
water_lost.wet_slod<-refqtlslod_wet_slod
water_lost.wet_max.perm.F_slod<-max.perm.F

##### WUE RATIO

## Plant WUE ratio
setwd(wue_results.qtl.wue_ratio.day.dry.dir)
load('timeseries_dry_cross.object.raw.Rdata')

# Re-name scanone and stepwise QTL results
wue_ratio.dry_out.F<-out.F
wue_ratio.dry_slod<-refqtlslod_dry_slod
wue_ratio.dry_max.perm.F_slod<-max.perm.F

setwd(wue_results.qtl.wue_ratio.day.wet.dir)
load('timeseries_wet_cross.object.raw.Rdata')

# Re-name scanone and stepwise QTL results
wue_ratio.wet_out.F<-out.F
wue_ratio.wet_slod<-refqtlslod_wet_slod
wue_ratio.wet_max.perm.F_slod<-max.perm.F


##### TE MODEL FIT

## TE model fit
setwd(wue_results.qtl.te_fit.day.dry.dir)
load('timeseries_dry_cross.object.raw.Rdata')


# Re-name scanone and stepwise QTL results
te_fit.dry_out.F<-out.F
te_fit.dry_slod<-refqtlslod_dry_slod
te_fit.dry_max.perm.F_slod<-max.perm.F

setwd(wue_results.qtl.te_fit.day.wet.dir)
load('timeseries_wet_cross.object.raw.Rdata')

# Re-name scanone and stepwise QTL results
te_fit.wet_out.F<-out.F
te_fit.wet_slod<-refqtlslod_wet_slod
te_fit.wet_max.perm.F_slod<-max.perm.F


##### TE MODEL RESIDUAL

## TE model residual
setwd(wue_results.qtl.te_residual.day.dry.dir)
load('timeseries_dry_cross.object.raw.Rdata')

# Re-name scanone and stepwise QTL results
te_residual.dry_out.F<-out.F
te_residual.dry_slod<-refqtlslod_dry_slod
te_residual.dry_max.perm.F_slod<-max.perm.F

setwd(wue_results.qtl.te_residual.day.wet.dir)
load('timeseries_wet_cross.object.raw.Rdata')

# Re-name scanone and stepwise QTL results
te_residual.wet_out.F<-out.F
te_residual.wet_slod<-refqtlslod_wet_slod
te_residual.wet_max.perm.F_slod<-max.perm.F


## Lets make Figure S21

setwd(wue_results.qtl.dir)
pdf("slod_qtl_rate_traits_stepwise.pdf", width=8, height = 3)
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

## Lets make Figure S22

setwd(wue_results.qtl.dir)
pdf("slod_qtl_rate_traits_scanone.pdf", width=8, height = 3)
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


setwd(wue_results.qtl.dir)
save.image("qtl_timeseries_plots_rate.Rdata")

rm(list=ls())


