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

setwd(wue_data.dir)

## Add planting data
planting_date = as.POSIXct("2014-1-13")
data<-read.table("vis_setaria_ril.csv", sep=";", header=T)
data<-data[data$plantbarcode != '0370',]

## Convert to date, make that value an integer
data$date<-as.POSIXct(data$timestamp, origin = "1970-01-01")
data$dap <- as.numeric(data$date - planting_date)
data$dap_i<-as.integer(data$dap)

## Read a key that associates phenotyper ID (aka 'plant_id') to RIL id
key<-read.csv('phe2ril2.csv', header=F)
colnames(key)<-c('phe', 'ril')

## Make the RIL name same as in QTL map file
key$ril<-sprintf("%03d", key$ril)
key$ril<-paste('RIL_', key$ril, sep="")

## Add a new (empty column called genotype)
data$genotype<-NA
## Add genotypes 
for(i in 1:nrow(key)) {
  phe_name<-key[i,1]
  ril_name<-key[i,2]
  data$genotype[grep(paste('Dr', phe_name, 'A', sep=""), data$plantbarcode)]<-ril_name
}

## Add parental genotypes
data$genotype[grep('Dp1', data$plantbarcode)]<-c("A10")
data$genotype[grep('Dp2', data$plantbarcode)]<-c("B100")

## Add treatment
data$treatment<-c('none')
data$treatment[grep("AA", data$plantbarcode)]<-c('wet')
data$treatment[grep("AB", data$plantbarcode)]<-c('dry')

## Lets remove all the coordinate based landmark features and other non-informative fields
data<-data[,c(1:4,8:12,13:14,18:41,47:55,70:71)]

## Decided to not remove images where the plant object is out of bound
## This likely will contribute useful data whereas for plant height it would not

#in_bound_flag<-data[,c(9,11,22,45)]
#in_bound_flag<-unique(in_bound_flag)

## Do calibration for area and height 
zoom.lm = lm(zoom.camera ~ zoom, data=data.frame(zoom=c(1,6000), zoom.camera=c(1,6)))
summary(zoom.lm)

wue_data_zoom_cal_file<-paste(wue_data.dir, "/zoom_calibration_data.txt", sep="")

## Get calibration data
if (!file.exists(wue_data_zoom_cal_file)) {
  download.file('http://files.figshare.com/2084101/zoom_calibration_data.txt',
                'zoom_calibration_data.txt')
}


## Note there are not calibration points for all of the top-view zooms present in this experiment
## We cannot use these values in our predictive model

## Read in thsi table and translate pixel to cm
z.data = read.table(file=wue_data_zoom_cal_file, sep="\t", header=TRUE)
z.data$px_cm = z.data$length_px / z.data$length_cm
z.data$zoom.camera = predict(object = zoom.lm, newdata=z.data)

## Define zooms
data$zoom = sub('z', '', data$zoom)
data$zoom<-as.numeric(data$zoom)
data$sv.zoom.camera = predict(object = zoom.lm, newdata=data)
data$tv.zoom.camera = predict(object = zoom.lm, newdata=data)

## Look for model that works 
## Exponential first
area.coef = coef(nls(log(rel_area) ~ log(a * exp(b * zoom.camera)),
                     z.data, start = c(a = 1, b = 0.01)))
area.coef = data.frame(a = area.coef[1], b = area.coef[2])
area.nls = nls(rel_area ~ a * exp(b * zoom.camera),
               data = z.data, start = c(a = area.coef$a, b = area.coef$b))
summary(area.nls)

## Polynomial
area.pol = lm(rel_area ~ zoom.camera + I(zoom.camera^2), z.data)
summary(area.pol)

## Compare models
AIC(area.nls, area.pol)
## Exponential is best (AIC score)

## Now do the same thing for non-area features
## Exponential
len.coef = coef(nls(log(px_cm) ~ log(a * exp(b * zoom.camera)),
                    z.data[z.data$camera == 'VIS SV',], start = c(a = 1, b = 0.01)))
len.coef = data.frame(a = len.coef[1], b = len.coef[2])
len.nls = nls(px_cm ~ a * exp(b * zoom.camera),
              data = z.data[z.data$camera == 'VIS SV',],
              start = c(a = len.coef$a, b = len.coef$b))
summary(len.nls)

## Polynomial
len.poly = lm(px_cm ~ zoom.camera + I(zoom.camera^2),
              data=z.data[z.data$camera == 'VIS SV',])
summary(len.poly)

AIC(len.nls, len.poly)


## Create zoom adjusted values for side view and top view
vis.data <- data
sv.zoom.camera <- vis.data[vis.data$camera == 'SV', c(4, 5, 8:10, 45, 47)]
colnames(sv.zoom.camera)[7] <- c('zoom.camera')
tv.zoom.camera <- vis.data[vis.data$camera == 'TV', c(4, 5, 8:10, 45, 48)]
colnames(tv.zoom.camera)[7] <- c('zoom.camera')
zoom.camera <- rbind(sv.zoom.camera, tv.zoom.camera)
vis.data <- merge(vis.data, zoom.camera, by = c("frame", "cartag", "zoom","plantbarcode","camera","dap_i")) ##THIS DOESN'T WORK BECAUSE THESE COLUMNS ARE FOUND IN ZOOM.CAMERA

vis.data$sv_rel_area = predict(object = area.nls, newdata = vis.data)
vis.data$tv_rel_area = predict(object = area.nls, newdata = vis.data)

################### Calculate Object Lengths

vis.data$px_cm = predict(object = len.poly, newdata=vis.data)
vis.data$perimeter = vis.data$perimeter / vis.data$px_cm
vis.data$width = vis.data$width / vis.data$px_cm
vis.data$height = vis.data$height / vis.data$px_cm
vis.data$longest_axis = vis.data$longest_axis / vis.data$px_cm
vis.data$ellipse_major_axis = vis.data$ellipse_major_axis / vis.data$px_cm
vis.data$ellipse_minor_axis = vis.data$ellipse_minor_axis / vis.data$px_cm
vis.data$height_above_bound = vis.data$height_above_bound / vis.data$px_cm
vis.data$height_below_bound = vis.data$height_below_bound / vis.data$px_cm



measured_biomass_data<-paste(wue_data.dir, "/measured_biomass_setaria_ril.csv", sep="")


# Load in biomass data
biomass<-read.csv(measured_biomass_data)

## For each car find the sampling day for model calibration
barcode<-biomass$plantbarcode
# Get the a data.frame with plantbarcode and days after planting
barcode_day<-vis.data[vis.data$plantbarcode %in% barcode,c(4,6)]

# Get the max day for each barcode
uni.barcode<-as.character(unique(barcode_day$plantbarcode))

max.days<-c()
for (u in uni.barcode){
  temp<-barcode_day[barcode_day$plantbarcode == u,]
  max.d<-temp[temp$dap_i == max(temp$dap_i),]
  max.days<-rbind(max.days, max.d)
}

## Looks like plants were sampled late
max.days<-unique(max.days)

## Now lets get the biomass with final date
biomass<-merge(biomass, max.days, by=c('plantbarcode'))
biomass$plantbarcode<-as.character(biomass$plantbarcode)

## The cartags don't line up between datasets. Likey due to error while recording cartags manually.
biomass<-biomass[,c(1:4,6)]
cal.table<-merge(vis.data, biomass, by=c('plantbarcode', 'dap_i'))

## Lets get average zoom adjusted area
cal.sv<-cal.table[cal.table$camera == 'SV',]
cal.sv$sv_area<-cal.sv$area/cal.sv$sv_rel_area
## Lets remove all plants that have area of zero
cal.sv<-cal.sv[cal.sv$area != 0,]
cal.ag<-cal.sv[,c(1,2,4,5,6,12:56)]
## Get the average of all 4 rotations
cal.sv.ag<-aggregate(cal.ag[,c(7:16,18:39,41:50)], by=list(cal.ag$plantbarcode, cal.ag$dap_i, cal.ag$treatment), mean)
colnames(cal.sv.ag)[1:3]<-c("plantbarcode", "dap_i", "treatment")

cal.tv<-cal.table[cal.table$camera == 'TV',]
cal.tv$tv_area<-cal.tv$area/cal.tv$tv_rel_area
## Lets remove all plants that have area of zero
cal.tv<-cal.tv[cal.tv$area != 0,]
cal.tv<-cal.tv[,c(1,2,4,5,12,14:29,37,46,48,51,56)]
colnames(cal.tv)[c(6:22)]<-paste(colnames(cal.tv)[c(6:22)], 'tv', sep="_")

## Now lets combine the sideview and topview to model area
calibrate<-merge(cal.sv, cal.tv, by=c('plantbarcode', 'dap_i', 'cartag', 'zoom', 'treatment','genotype'))

############# FRESH WEIGHT
## See what happens when you predict biomass with all attributes no interactions
fw.all = lm(fresh_wt ~ sv_area + perimeter + height_above_bound + width + solidity, calibrate)
summary(fw.all)

## Try stepwise addition and subtraction
fw.all.step = stepAIC(fw.all, direction="both")
summary(fw.all.step)
anova(fw.all, fw.all.step)
AIC(fw.all,fw.all.step)

## Results are identical
## All components are predictive (width is boarderline)

## Plot predicted vs. actual
calibrate$fw.all = predict.lm(object = fw.all, newdata=calibrate)
p<-ggplot(data=calibrate, aes(x=fresh_wt,y=fw.all,colour=treatment)) + geom_point()
p

## Try the same thing as above but with all possible interactions
fw.all.int = lm(fresh_wt ~ sv_area * perimeter * height_above_bound * width * solidity, calibrate)
summary(fw.all.int)

## Use the stepwise method to get best fit model
fw.all.int.step = stepAIC(fw.all.int, direction="both")
summary(fw.all.int.step)
anova(fw.all.int,fw.all.int.step)
AIC(fw.all.int, fw.all.int.step)

## Fit looks good but it is a over fit of the data
calibrate$fw.all.int.step = predict.lm(object = fw.all.int.step, newdata=calibrate)
p<-ggplot(data=calibrate, aes(x=fresh_wt,y=fw.all.int.step,colour=treatment)) + geom_point()
p


# Try a simpler version of the model chose by AIC
fw.aic = lm(fresh_wt ~ sv_area + perimeter + height_above_bound + width + solidity + sv_area:height_above_bound + sv_area:width + sv_area:perimeter + height_above_bound:width, calibrate)
summary(fw.aic)
fw.aic.step = stepAIC(fw.aic, direction="both")
summary(fw.aic.step)
anova(fw.aic, fw.aic.step)
AIC(fw.aic, fw.aic.step)

# Looks better but still over fits the data
calibrate$fw.aic = predict.lm(object = fw.aic.step, newdata=calibrate)
p<-ggplot(data=calibrate, aes(x=fresh_wt,y=fw.aic,colour=treatment)) + geom_point()
p

# Examine how the minimal model projects the data.
fw.min = lm(fresh_wt ~ sv_area, calibrate)
summary(fw.min)
anova(fw.aic.step, fw.min)
AIC(fw.aic.step, fw.min)
calibrate$fw.min = predict.lm(object = fw.aic, newdata=calibrate)
p<-ggplot(data=calibrate, aes(x=fresh_wt,y=fw.all.int.step,colour=treatment)) + geom_point()
p
# Relatively high support for the minimal model

# Patrick suggested we try a slightly more complex model. It looks to fit the data much better 
# and doesn't exhibit artifacts of over-fitting.
fw.basic = lm(fresh_wt ~ sv_area + perimeter + height_above_bound, calibrate)
AIC(fw.min, fw.basic)
calibrate$fw.basic = predict.lm(object = fw.basic, newdata=calibrate)
p<-ggplot(data=calibrate, aes(x=fresh_wt,y=fw.basic,colour=treatment)) + geom_point()
p

# A colleague recommended adding treatment as a factor to the model
# It improves the fit this is unacceptable because treatment does not have an influence on biomass early in experiment
fw.basic.t = lm(fresh_wt ~ sv_area + perimeter + height_above_bound + treatment, calibrate)
calibrate$fw.basic.t = predict.lm(object = fw.basic.t, newdata=calibrate)
p<-ggplot(data=calibrate, aes(x=fresh_wt,y=fw.basic.t,colour=treatment)) + geom_point()
p
# Will use fw.basic for fresh weight (so as not to have treatment influence values before treatment)
AIC(fw.basic, fw.basic.t)


############# DRY WEIGHT
# Now lets do the same as above for the dry weight measurements (dry_wt)

# See what happens when you predict biomass with all attributes no interactions
dw.all = lm(dry_wt ~ sv_area + solidity + perimeter + width + height_above_bound, calibrate)
dw.all.step = stepAIC(dw.all, direction="both")
summary(dw.all)
summary(dw.all.step)
AIC(dw.all,dw.all.step)

calibrate$dw.all.step = predict.lm(object = dw.all.step, newdata=calibrate)
p<-ggplot(data=calibrate, aes(x=dry_wt,y=dw.all.step, colour=treatment)) + geom_point()

dw.all.int = lm(dry_wt ~ sv_area * height_above_bound * width, calibrate)
dw.all.int.step = stepAIC(dw.all.int, direction="both")
summary(dw.all.int.step)
AIC(dw.all.int, dw.all.int.step)

calibrate$dw.all.int.step = predict.lm(object = dw.all.int.step, newdata=calibrate)
p<-ggplot(data=calibrate, aes(x=dry_wt,y=dw.all.int.step,colour=treatment)) + geom_point()

AIC(dw.all, dw.all.int.step)

dw.aic = lm(dry_wt ~ sv_area * height_above_bound + width, calibrate)
summary(dw.aic)
anova(dw.all.int.step, dw.aic)
AIC(dw.all.int.step, dw.aic)
calibrate$dw.aic = predict.lm(object = dw.aic, newdata=calibrate)
p<-ggplot(data=calibrate, aes(x=dry_wt,y=dw.aic,colour=treatment)) + geom_point()

dw.aic.step = stepAIC(dw.aic, direction="both")
summary(dw.aic.step)
AIC(dw.aic, dw.aic.step)

dw.min = lm(dry_wt ~ sv_area, calibrate)
summary(dw.min)

anova(dw.aic, dw.min)
AIC(dw.aic, dw.min)
calibrate$dw.min = predict.lm(object = dw.aic, newdata=calibrate)
p<-ggplot(data=calibrate, aes(x=fresh_wt,y=dw.min,colour=treatment)) + geom_point()

# Will use dw.aic for dry weight
dw.basic = lm(dry_wt ~ sv_area + width + height_above_bound, calibrate)
calibrate$dw.basic = predict.lm(object = dw.basic, newdata=calibrate)
p<-ggplot(data=calibrate, aes(x=fresh_wt,y=dw.basic,colour=treatment)) + geom_point()
p
AIC(dw.basic, dw.min)

# Will use dw.basic for dry weight
dw.basic.t = lm(dry_wt ~ sv_area + width + height_above_bound + treatment, calibrate)
calibrate$dw.basic.t = predict.lm(object = dw.basic.t, newdata=calibrate)
p<-ggplot(data=calibrate, aes(x=fresh_wt,y=dw.basic.t,colour=treatment)) + geom_point()
p
AIC(dw.basic, dw.basic.t)

######### Use model to calculate biomass

# Lets get average zoom adjusted area
vis.sv<-vis.data[vis.data$camera == 'SV',]
vis.sv$sv_area<-vis.sv$area/vis.sv$sv_rel_area
vis.sv$hull.area<-vis.sv$hull.area/vis.sv$sv_rel_area

# Lets remove all plants that have area of zero
vis.sv<-vis.sv[vis.sv$area != 0,]
# Remove all images where plants are out of bounds
vis.sv$in_bounds<-as.character(vis.sv$in_bounds) 
# Chose not to remove values which were out of bounds
#vis.sv<-vis.sv[vis.sv$in_bounds == 'True', ]
vis.ag<-vis.sv[,c(1,2,4,5,6,12:53)]

# Get the average of all 4 rotations
vis.sv.ag<-aggregate(vis.ag[,c(7:16,18:39,41:47)], by=list(vis.ag$plantbarcode, vis.ag$genotype, vis.ag$dap_i, vis.ag$treatment), mean)
colnames(vis.sv.ag)[1:4]<-c("plantbarcode", "genotype", "dap_i", "treatment")
vis.sv.ag<-vis.sv.ag[,-c(38,39,41,42)]

vis.tv<-vis.data[vis.data$camera == 'TV',]
vis.tv$tv_area<-vis.tv$area/vis.tv$tv_rel_area
# Lets remove all plants that have area of zero
vis.tv<-vis.tv[vis.tv$area != 0,]
vis.tv<-vis.tv[,c(1,2,4,5,6,12,14:29,37,46,48,51,53)]
colnames(vis.tv)[c(7:23)]<-paste(colnames(vis.tv)[c(7:23)], 'tv', sep="_")

# Merge top view and side view
vis<-merge(vis.sv.ag, vis.tv, by=c('plantbarcode', 'dap_i', 'treatment','genotype'))

# Use models to predict fresh weight and dry weight biomass

# Fresh weight biomass (using model fw.min) other methods lead to large artifacts
vis$fw_biomass = predict.lm(object = fw.basic, newdata=vis)

# Dry weight biomass (using model dw.min) other methods lead to large artifacts
vis$dw_biomass = predict.lm(object = dw.basic, newdata=vis)

# Lets do some diagnosic plotting
p<-ggplot(data=vis, aes(x=fw_biomass, y = dw_biomass, colour=treatment)) + geom_point(size=0.3)
p<-ggplot(data=vis, aes(x=fw_biomass, y = sv_area, colour=treatment)) + geom_point(size=0.3)
p<-ggplot(data=vis, aes(x=dw_biomass, y = sv_area, colour=treatment)) + geom_point(size=0.3)

p<-ggplot(data=vis, aes(x=dap_i, y = fw_biomass, colour=treatment)) + geom_point(size=0.3)
p<-ggplot(data=vis, aes(x=dap_i, y = dw_biomass, colour=treatment)) + geom_point(size=0.3)

# Too many negative values for both. Take the minimum value + 0.1 and add it to all
fw_min<-min(vis$fw_biomass) - 0.1
vis$fw_biomass<-vis$fw_biomass - fw_min
dw_min<-min(vis$dw_biomass) - 0.1
vis$dw_biomass<-vis$dw_biomass - dw_min

p<-ggplot(data=vis, aes(x=dap_i, y = fw_biomass, colour=treatment)) + geom_point(size=0.3)
p<-ggplot(data=vis, aes(x=dap_i, y = dw_biomass, colour=treatment)) + geom_point(size=0.3)


p<-ggplot(data=vis, aes(x=dap_i, y = height_above_bound, colour=treatment)) + geom_point(size=0.3)

# Need to zoom in 17 - 32 dap_i
close_up<-vis[vis$dap_i > 16,]
plot(close_up$height_above_bound~close_up$dap_i, ylim=c(0,25), pch=20, cex=0.3)
close_up<-close_up[close_up$plantbarcode != 'Dr27AA001310',]
close_up<-close_up[close_up$plantbarcode != 'Dr41AA001393',]
close_up<-close_up[close_up$plantbarcode != 'Dr53AB001468',]
close_up<-close_up[close_up$plantbarcode != 'Dr72AA001581',]
close_up<-close_up[close_up$plantbarcode != 'Dr65AB001541',]

filter_1<-close_up[close_up$height_above_bound > 5,]
# Get empty pots
empty_pots<-c("Dr27AA001310", "Dr41AA001393", "Dr53AB001468", "Dr65AB001541", "Dr72AA001581")

# Now remove all of those images from RIL mean that are no good
vis<-vis[!vis$plantbarcode %in% empty_pots, ]
plot(vis$height_above_bound~vis$dap_i, pch=20, cex=0.3, ylab="height", xlab="day")

# Lets convert biomass to milligrams (mg) from grams (g)
vis$fw_biomass<-vis$fw_biomass*1000
vis$dw_biomass<-vis$dw_biomass*1000

plot(vis$fw_biomass~vis$dap_i, pch=20, cex=0.3, ylab="fw_biomass (mg)", xlab="day")
plot(vis$dw_biomass~vis$dap_i, pch=20, cex=0.3, ylab="dw_biomass (mg)", xlab="day")

plot(vis$fw_biomass~vis$dap_i, pch=20, cex=0.3, ylab="fw_biomass (mg)", xlab="day", ylim=c(0,3000))
vis[!vis$fw_biomass < 3000 & vis$dap_i >29,]

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

setwd(wue_results.plant_size.dir)
save.image('ril_biomass.Rdata')


############################################
# Now go for a loess fit for biomass
############################################

genos<-as.character(sort(unique(vis$genotype)))
treatments<-unique(vis$treatment)
dap_i<-unique(vis$dap_i)
dap_i<-sort(as.numeric(as.character(dap_i)))
times = seq(from = 8, to = 33, by=0.1)

report.loess.values(vis, 'fw_biomass', genos, treatments, dap_i, min(dap_i), max(dap_i), "ril_loess_estimates_fw_biomass_total.pdf")

# Give the output an appropriate name
ril_loess_fw_biomass<-ril_loess_model_fit
# Rename the growth report
fw_biomass_report<-growth_rate_report


# Format for QTL pipeline
biomass.l<-ril_loess_fw_biomass

days<-sort(unique(biomass.l$dap_i))
colnames(biomass.l)
colnames(biomass.l)[1]<-c("genotype")
ril_fw_biomass_qtl<-c()
for(d in 1:length(days)) {
  day<-days[d]
  temp.data<-biomass.l[biomass.l$dap_i == day, ]
  colnames(temp.data)[4]<-paste('fw_biomass_total', day, sep="_")
  temp.data<-temp.data[,c(1,2,4)]
  if (d == 1) {
    ril_fw_biomass_qtl<-temp.data
  }
  
  if (d > 1) {
    ril_fw_biomass_qtl<-merge(ril_fw_biomass_qtl, temp.data, by=c('genotype', 'treatment'), all=T)
  }
}

# Add misc columns
ril_fw_biomass_qtl$Obs<-c(1:nrow(ril_fw_biomass_qtl))
ril_fw_biomass_qtl$experiment<-rep("BP14", nrow(ril_fw_biomass_qtl))
ril_fw_biomass_qtl$year<-rep("2014", nrow(ril_fw_biomass_qtl))
ril_fw_biomass_qtl$plot<-rep("bellweater", nrow(ril_fw_biomass_qtl))
ril_fw_biomass_qtl$plot_id<-rep("bellweater", nrow(ril_fw_biomass_qtl))
ril_fw_biomass_qtl$sampling<-rep("bellweater", nrow(ril_fw_biomass_qtl))

ril_fw_biomass_qtl<-ril_fw_biomass_qtl[,c(29,30,31,2,32,33,1,34,3:28)]
colnames(ril_fw_biomass_qtl)[7]<-c("id")
## Write to file
write.csv(ril_fw_biomass_qtl, file="ril_loess_fw_biomass_qtl.csv", quote=F, row.names=F)


############################################
# Now go for a loess fit for sv_area
############################################

genos<-as.character(sort(unique(vis$genotype)))
treatments<-unique(vis$treatment)
dap_i<-unique(vis$dap_i)
dap_i<-sort(as.numeric(as.character(dap_i)))
times = seq(from = 8, to = 33, by=0.1)
report.loess.values(vis, 'sv_area', genos, treatments, dap_i, min(dap_i), max(dap_i), "ril_loess_estimates_sv_area_total.pdf")

# Remove fine grain variables
ril_loess_sv_area<-ril_loess_model_fit
# Rename the growth report
sv_area_report<-growth_rate_report

# Format for QTL pipeline
biomass.l<-ril_loess_sv_area

days<-sort(unique(biomass.l$dap_i))
colnames(biomass.l)
colnames(biomass.l)[1]<-c("genotype")
ril_sv_area_qtl<-c()
for(d in 1:length(days)) {
  day<-days[d]
  temp.data<-biomass.l[biomass.l$dap_i == day, ]
  colnames(temp.data)[4]<-paste('sv_area_total', day, sep="_")
  temp.data<-temp.data[,c(1,2,4)]
  if (d == 1) {
    ril_sv_area_qtl<-temp.data
  }
  
  if (d > 1) {
    ril_sv_area_qtl<-merge(ril_sv_area_qtl, temp.data, by=c('genotype', 'treatment'), all=T)
  }
}

# Add misc columns
ril_sv_area_qtl$Obs<-c(1:nrow(ril_sv_area_qtl))
ril_sv_area_qtl$experiment<-rep("BP14", nrow(ril_sv_area_qtl))
ril_sv_area_qtl$year<-rep("2014", nrow(ril_sv_area_qtl))
ril_sv_area_qtl$plot<-rep("bellweater", nrow(ril_sv_area_qtl))
ril_sv_area_qtl$plot_id<-rep("bellweater", nrow(ril_sv_area_qtl))
ril_sv_area_qtl$sampling<-rep("bellweater", nrow(ril_sv_area_qtl))

ril_sv_area_qtl<-ril_sv_area_qtl[,c(29,30,31,2,32,33,1,34,3:28)]
colnames(ril_sv_area_qtl)[7]<-c("id")
write.csv(ril_sv_area_qtl, file="ril_loess_sv_area_total_qtl.csv", quote=F, row.names=F)

save.image('ril_biomass.Rdata')
rm(list=ls())

