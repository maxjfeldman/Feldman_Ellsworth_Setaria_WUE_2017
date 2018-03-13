## This is a script to download the files needed to perform the analysis

## The first thing is to set your R session to the base directory you just downloaded from github
## insert path below...
setwd()

## Tester
#setwd("~/Dropbox/Feldman_Ellsworth_Setaria_WUE_2017/")

##### CREATE DIRECTORY PATHS ##### 

## Make the directory of the folder you downloaded the current working directory
home.dir<-getwd()
setwd(home.dir)
load('analysis_fxns.Rdata')
## Lets make a path to the download

download.dir<-paste(home.dir, '/feldman_ellsworth_setaria_wue_2017_data', sep="")

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




## Download bigger R workspace files from zenodo repository
download.file('https://zenodo.org/record/1117062/files/feldman_ellsworth_setaria_wue_2017_data.zip',
              'setaria_ril_wue_data_zenodo.zip', mode='wb', method='curl')

## Unzip this directory
unzip('setaria_ril_wue_data_zenodo.zip')
#unzip('setaria_ril_wue_figshare_data.zip')

## Now lets transfer the files in the zenodo download to the correct locations in the git download

## Move vis.data file to the data directory
visdata<-paste(download.dir, "/vis_setaria_ril.csv", sep="")
file.copy(from=visdata, to=wue_data.dir)

##### PLANT SIZE

## Move plant_size .Rdata files to the qtl_summary/plant_size directory
download.qtl_summary.dir<-paste(download.dir, "/qtl_summary", sep="")
download.qtl_summary.plant_size.dir<-paste(download.qtl_summary.dir, "/plant_size", sep="")
download.qtl_summary.plant_size.total.dir<-paste(download.qtl_summary.plant_size.dir, "/sv_area_total", sep="")
download.qtl_summary.plant_size.day.dir<-paste(download.qtl_summary.plant_size.dir, "/sv_area_day", sep="")

##### CUMULATIVE

## Cumulative area dry
sv_area_total.dry.qtl.file<-paste(download.qtl_summary.plant_size.total.dir, "/dry/timeseries_dry_cross.object.raw.Rdata", sep="")
## Copy file
file.copy(from=sv_area_total.dry.qtl.file, to=wue_results.qtl.plant_size.total.dry.dir, overwrite = TRUE)


## Cumulative area wet
sv_area_total.wet.qtl.file<-paste(download.qtl_summary.plant_size.total.dir, "/wet/timeseries_wet_cross.object.raw.Rdata", sep="")
## Copy file
file.copy(from=sv_area_total.wet.qtl.file, to=wue_results.qtl.plant_size.total.wet.dir, overwrite = TRUE)

##### RATE

## Rate area dry
sv_area_day.dry.qtl.file<-paste(download.qtl_summary.plant_size.day.dir, "/dry/timeseries_dry_cross.object.raw.Rdata", sep="")
## Copy file
file.copy(from=sv_area_day.dry.qtl.file, to=wue_results.qtl.plant_size.day.dry.dir, overwrite = TRUE)


## Rate area wet
sv_area_day.wet.qtl.file<-paste(download.qtl_summary.plant_size.day.dir, "/wet/timeseries_wet_cross.object.raw.Rdata", sep="")
## Copy file
file.copy(from=sv_area_day.wet.qtl.file, to=wue_results.qtl.plant_size.day.wet.dir, overwrite = TRUE)


##### WATER LOST

## Move water_lost .Rdata files to the qtl_summary/water_lost directory
download.qtl_summary.water_lost.dir<-paste(download.qtl_summary.dir, "/water_lost", sep="")
download.qtl_summary.water_lost.total.dir<-paste(download.qtl_summary.water_lost.dir, "/water_lost_total", sep="")
download.qtl_summary.water_lost.day.dir<-paste(download.qtl_summary.water_lost.dir, "/water_lost_day", sep="")

##### CUMULATIVE

## Cumulative water lost dry
water_lost_total.dry.qtl.file<-paste(download.qtl_summary.water_lost.total.dir, "/dry/timeseries_dry_cross.object.raw.Rdata", sep="")
## Copy file
file.copy(from=water_lost_total.dry.qtl.file, to=wue_results.qtl.water_lost.total.dry.dir, overwrite = TRUE)


## Cumulative water lost wet
water_lost_total.wet.qtl.file<-paste(download.qtl_summary.water_lost.total.dir, "/wet/timeseries_wet_cross.object.raw.Rdata", sep="")
## Copy file
file.copy(from=water_lost_total.wet.qtl.file, to=wue_results.qtl.water_lost.total.wet.dir, overwrite = TRUE)

##### RATE

## Rate water lost dry
water_lost_day.dry.qtl.file<-paste(download.qtl_summary.water_lost.day.dir, "/dry/timeseries_dry_cross.object.raw.Rdata", sep="")
## Copy file
file.copy(from=water_lost_day.dry.qtl.file, to=wue_results.qtl.water_lost.day.dry.dir, overwrite = TRUE)


## Rate water lost wet
water_lost_day.wet.qtl.file<-paste(download.qtl_summary.water_lost.day.dir, "/wet/timeseries_wet_cross.object.raw.Rdata", sep="")
## Copy file
file.copy(from=water_lost_day.wet.qtl.file, to=wue_results.qtl.water_lost.day.wet.dir, overwrite = TRUE)


##### WUE RATIO

## Move wue ratio .Rdata files to the qtl_summary/wue_ratio directory
download.qtl_summary.wue_ratio.dir<-paste(download.qtl_summary.dir, "/wue_ratio", sep="")
download.qtl_summary.wue_ratio.total.dir<-paste(download.qtl_summary.wue_ratio.dir, "/wue_ratio_total", sep="")
download.qtl_summary.wue_ratio.day.dir<-paste(download.qtl_summary.wue_ratio.dir, "/wue_ratio_day", sep="")

##### CUMULATIVE

## Cumulative wue ratio dry
wue_ratio_total.dry.qtl.file<-paste(download.qtl_summary.wue_ratio.total.dir, "/dry/timeseries_dry_cross.object.raw.Rdata", sep="")
## Copy file
file.copy(from=wue_ratio_total.dry.qtl.file, to=wue_results.qtl.wue_ratio.total.dry.dir, overwrite = TRUE)


## Cumulative wue ratio wet
wue_ratio_total.wet.qtl.file<-paste(download.qtl_summary.wue_ratio.total.dir, "/wet/timeseries_wet_cross.object.raw.Rdata", sep="")
## Copy file
file.copy(from=wue_ratio_total.wet.qtl.file, to=wue_results.qtl.wue_ratio.total.wet.dir, overwrite = TRUE)

##### RATE

## Rate wue ratio dry
wue_ratio_day.dry.qtl.file<-paste(download.qtl_summary.wue_ratio.day.dir, "/dry/timeseries_dry_cross.object.raw.Rdata", sep="")
## Copy file
file.copy(from=wue_ratio_day.dry.qtl.file, to=wue_results.qtl.wue_ratio.day.dry.dir, overwrite = TRUE)


## Rate wue ratio wet
wue_ratio_day.wet.qtl.file<-paste(download.qtl_summary.wue_ratio.day.dir, "/wet/timeseries_wet_cross.object.raw.Rdata", sep="")
## Copy file
file.copy(from=wue_ratio_day.wet.qtl.file, to=wue_results.qtl.wue_ratio.day.wet.dir, overwrite = TRUE)




##### TE_FIT

## Move te_fit .Rdata files to the qtl_summary/te_fit directory
download.qtl_summary.te_fit.dir<-paste(download.qtl_summary.dir, "/te_fit", sep="")
download.qtl_summary.te_fit.total.dir<-paste(download.qtl_summary.te_fit.dir, "/te_fit_total", sep="")
download.qtl_summary.te_fit.day.dir<-paste(download.qtl_summary.te_fit.dir, "/te_fit_day", sep="")

##### CUMULATIVE

## Cumulative te_fit dry
te_fit_total.dry.qtl.file<-paste(download.qtl_summary.te_fit.total.dir, "/dry/timeseries_dry_cross.object.raw.Rdata", sep="")
## Copy file
file.copy(from=te_fit_total.dry.qtl.file, to=wue_results.qtl.te_fit.total.dry.dir, overwrite = TRUE)


## Cumulative te_fit wet
te_fit_total.wet.qtl.file<-paste(download.qtl_summary.te_fit.total.dir, "/wet/timeseries_wet_cross.object.raw.Rdata", sep="")
## Copy file
file.copy(from=te_fit_total.wet.qtl.file, to=wue_results.qtl.te_fit.total.wet.dir, overwrite = TRUE)

##### RATE

## Rate te_fit dry
te_fit_day.dry.qtl.file<-paste(download.qtl_summary.te_fit.day.dir, "/dry/timeseries_dry_cross.object.raw.Rdata", sep="")
## Copy file
file.copy(from=te_fit_day.dry.qtl.file, to=wue_results.qtl.te_fit.day.dry.dir, overwrite = TRUE)


## Rate te_fit wet
te_fit_day.wet.qtl.file<-paste(download.qtl_summary.te_fit.day.dir, "/wet/timeseries_wet_cross.object.raw.Rdata", sep="")
## Copy file
file.copy(from=te_fit_day.wet.qtl.file, to=wue_results.qtl.te_fit.day.wet.dir, overwrite = TRUE)




##### TE_RESIDUAL

## Move te_residual .Rdata files to the qtl_summary/te_residual directory
download.qtl_summary.te_residual.dir<-paste(download.qtl_summary.dir, "/te_residual", sep="")
download.qtl_summary.te_residual.total.dir<-paste(download.qtl_summary.te_residual.dir, "/te_residual_total", sep="")
download.qtl_summary.te_residual.day.dir<-paste(download.qtl_summary.te_residual.dir, "/te_residual_day", sep="")

##### CUMULATIVE

## Cumulative te_residual dry
te_residual_total.dry.qtl.file<-paste(download.qtl_summary.te_residual.total.dir, "/dry/timeseries_dry_cross.object.raw.Rdata", sep="")
## Copy file
file.copy(from=te_residual_total.dry.qtl.file, to=wue_results.qtl.te_residual.total.dry.dir, overwrite = TRUE)


## Cumulative te_residual wet
te_residual_total.wet.qtl.file<-paste(download.qtl_summary.te_residual.total.dir, "/wet/timeseries_wet_cross.object.raw.Rdata", sep="")
## Copy file
file.copy(from=te_residual_total.wet.qtl.file, to=wue_results.qtl.te_residual.total.wet.dir, overwrite = TRUE)

##### RATE

## Rate te_residual dry
te_residual_day.dry.qtl.file<-paste(download.qtl_summary.te_residual.day.dir, "/dry/timeseries_dry_cross.object.raw.Rdata", sep="")
## Copy file
file.copy(from=te_residual_day.dry.qtl.file, to=wue_results.qtl.te_residual.day.dry.dir, overwrite = TRUE)


## Rate te_residual wet
te_residual_day.wet.qtl.file<-paste(download.qtl_summary.te_residual.day.dir, "/wet/timeseries_wet_cross.object.raw.Rdata", sep="")
## Copy file
file.copy(from=te_residual_day.wet.qtl.file, to=wue_results.qtl.te_residual.day.wet.dir, overwrite = TRUE)


