## Lets make some Venn diagrams of what we saw in the ril_qtl_trait_matrix_plot
## ril_qtl_trait_matrix_plot script

## This serves as a quick way to identify trends and acts as a sanity check


################### 
library(ggplot2)
library(stringr)
library(VennDiagram)

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

## Make a directory path for the QTL summary tables
wue_results.qtl.dir<-paste(wue_results.dir, 'qtl_summary/', sep="")

## Make the directory to store venn diagrams in
wue_results.qtl.venn.dir<-paste(wue_results.qtl.dir, 'venn_diagrams/', sep="")

## If this directory doesn't exist the download is not complete
if (file.exists(wue_results.qtl.venn.dir)){
  print("Directory exists!")
} else {
  dir.create(file.path(wue_results.qtl.venn.dir))
}


##################################################################


setwd(wue_results.qtl.dir)
load("qtl_trait_matrix.Rdata")

## Lets get a list of files to delete later
files<-ls()

## Get cumulative/total QTL
t.uni.all_qtl<-r.uni.all.qtl[grepl('_total', r.uni.all.qtl$trait),]

## Get daily/rate QTL
d.uni.all_qtl<-r.uni.all.qtl[!grepl('_total', r.uni.all.qtl$trait),]

## Get QTL unique to each trait
t.sv_area<-unique(t.uni.all_qtl[grepl("sv_area", t.uni.all_qtl$trait),'qtl.name'])
## 18
t.water_lost<-unique(t.uni.all_qtl[grepl("water_lost", t.uni.all_qtl$trait),'qtl.name'])
## 8
t.wue_ratio<-unique(t.uni.all_qtl[grepl("wue_total", t.uni.all_qtl$trait),'qtl.name'])
## 12
t.te_fit<-unique(t.uni.all_qtl[grepl("te_fit", t.uni.all_qtl$trait),'qtl.name'])
## 11
t.te_residual<-unique(t.uni.all_qtl[grepl("te_residual", t.uni.all_qtl$trait),'qtl.name'])
## 10

all_qtl.t.nr<-unique(c(t.sv_area, t.water_lost, t.wue_ratio, t.te_fit, t.te_residual))

###############################################
## Lets look at QTL over lap for all traits
###############################################

########################################################################
## Cumulative (total)
########################################################################

setwd(wue_results.qtl.venn.dir)

pdf("venn_all_cumulative_traits.pdf")

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


intersect(unique(t.sv_area), intersect(unique(t.water_lost), intersect(unique(t.wue_ratio), intersect(unique(t.te_fit), unique(t.te_residual)))))
## 1 QTL found in all traits: 2@96

intersect(unique(t.sv_area), intersect(unique(t.water_lost), intersect(unique(t.wue_ratio), unique(t.te_fit))))
## 4 QTL found in everything but TE residual: 2@96 3@48 5@109 7@99 7@34 9@34 

intersect(unique(t.sv_area), intersect(unique(t.wue_ratio), unique(t.te_residual)))
##  QTL found for TE residual, WUE ratio and Biomass: 2@96, 2@113, 2@11, 5@92, 5@79, 9@127

intersect(unique(t.water_lost), intersect(unique(t.sv_area), unique(t.te_fit)))
## QTL found for water lost, sv area, te fit: 6@65, 7@51

intersect(unique(t.sv_area), intersect(unique(t.wue_ratio), unique(t.te_residual)))

no_biomass_u.qtl<-unique(c(t.water_lost, t.wue_ratio, t.te_fit, t.te_residual))
biomass_u.qtl<-unique(t.sv_area)
biomass_u.qtl[!biomass_u.qtl %in% no_biomass_u.qtl]
## QTL unique to biomass 3@21, 5@119, 6@80, 9@138   

no_te_residual_u.qtl<-unique(c(t.water_lost, t.wue_ratio, t.te_fit, t.sv_area))
te_residual_u.qtl<-unique(t.te_residual)
te_residual_u.qtl[!te_residual_u.qtl %in% no_te_residual_u.qtl]
## QTL unique to te_residual 2@82 3@77, 6@47 

no_wue_ratio_u.qtl<-unique(c(t.water_lost, t.te_fit, t.te_residual, t.sv_area))
wue_ratio_u.qtl<-unique(t.wue_ratio)
wue_ratio_u.qtl[!wue_ratio_u.qtl %in% no_wue_ratio_u.qtl]
## QTL unique to wue_ratio 0

no_te_fit_u.qtl<-unique(c(t.water_lost, t.wue_ratio, t.te_residual, t.sv_area))
te_fit_u.qtl<-unique(t.te_fit)
te_fit_u.qtl[!te_fit_u.qtl %in% no_te_fit_u.qtl]
## QTL unique to te_fit is 5@39

## How many QTL for each trait?
length(t.sv_area)
length(t.water_lost)
length(t.wue_ratio)
length(t.te_fit)
length(t.te_residual)


########################################################################
# Daily (rate)
########################################################################

## Get QTL unique to each trait
d.sv_area<-unique(d.uni.all_qtl[grepl("sv_area", d.uni.all_qtl$trait),'qtl.name'])
d.water_lost<-unique(d.uni.all_qtl[grepl("water_lost", d.uni.all_qtl$trait),'qtl.name'])
d.wue_ratio<-unique(d.uni.all_qtl[grepl("wue_day", d.uni.all_qtl$trait),'qtl.name'])
d.te_fit<-unique(d.uni.all_qtl[grepl("te_fit", d.uni.all_qtl$trait),'qtl.name'])
d.te_residual<-unique(d.uni.all_qtl[grepl("te_residual", d.uni.all_qtl$trait),'qtl.name'])

all_qtl.d.nr<-unique(c(d.sv_area, d.water_lost, d.wue_ratio, d.te_fit, d.te_residual))


grid.newpage()
draw.quintuple.venn(area1=length(d.sv_area), 
                    area2=length(d.water_lost), 
                    area3=length(d.wue_ratio),
                    area4=length(d.te_fit), 
                    area5=length(d.te_residual),
                    n12=length(intersect(unique(d.sv_area), unique(d.water_lost))),
                    n13=length(intersect(unique(d.sv_area), unique(d.wue_ratio))),
                    n14=length(intersect(unique(d.sv_area), unique(d.te_fit))),
                    n15=length(intersect(unique(d.sv_area), unique(d.te_residual))),
                    n23=length(intersect(unique(d.water_lost), unique(d.wue_ratio))),
                    n24=length(intersect(unique(d.water_lost), unique(d.te_fit))),
                    n25=length(intersect(unique(d.water_lost), unique(d.te_residual))),
                    n34=length(intersect(unique(d.wue_ratio), unique(d.te_fit))),
                    n35=length(intersect(unique(d.wue_ratio), unique(d.te_residual))),
                    n45=length(intersect(unique(d.te_fit), unique(d.te_residual))),
                    n123=length(intersect(unique(d.sv_area), intersect(unique(d.water_lost), unique(d.wue_ratio)))),
                    n124=length(intersect(unique(d.sv_area), intersect(unique(d.water_lost), unique(d.te_fit)))),
                    n125=length(intersect(unique(d.sv_area), intersect(unique(d.water_lost), unique(d.te_residual)))),
                    n134=length(intersect(unique(d.sv_area), intersect(unique(d.wue_ratio), unique(d.te_fit)))),
                    n135=length(intersect(unique(d.sv_area), intersect(unique(d.wue_ratio), unique(d.te_residual)))),
                    n145=length(intersect(unique(d.sv_area), intersect(unique(d.te_fit), unique(d.te_residual)))),
                    n234=length(intersect(unique(d.water_lost), intersect(unique(d.wue_ratio), unique(d.te_fit)))),
                    n235=length(intersect(unique(d.water_lost), intersect(unique(d.wue_ratio), unique(d.te_residual)))),
                    n245=length(intersect(unique(d.water_lost), intersect(unique(d.te_fit), unique(d.te_residual)))),
                    n345=length(intersect(unique(d.wue_ratio), intersect(unique(d.te_fit), unique(d.te_residual)))),
                    n1234=length(intersect(unique(d.sv_area), intersect(unique(d.water_lost), intersect(unique(d.wue_ratio), unique(d.te_fit))))),
                    n1235=length(intersect(unique(d.sv_area), intersect(unique(d.water_lost), intersect(unique(d.wue_ratio), unique(d.te_residual))))),
                    n1245=length(intersect(unique(d.sv_area), intersect(unique(d.water_lost), intersect(unique(d.te_fit), unique(d.te_residual))))),
                    n1345=length(intersect(unique(d.sv_area), intersect(unique(d.wue_ratio), intersect(unique(d.te_fit), unique(d.te_residual))))),
                    n2345=length(intersect(unique(d.water_lost), intersect(unique(d.wue_ratio), intersect(unique(d.te_fit), unique(d.te_residual))))),
                    n12345=length(intersect(unique(d.sv_area), intersect(unique(d.water_lost), intersect(unique(d.wue_ratio), intersect(unique(d.te_fit), unique(d.te_residual)))))),
                    category=c("Biomass", "Water lost", "WUE ratio", "TE fit", "TE residual"),
                    fill=c("green","light blue", "dark blue", "orange", "red"),
                    cex=2,
                    cat.cex=1.5
)




intersect(unique(d.sv_area), intersect(unique(d.water_lost), intersect(unique(d.wue_ratio), intersect(unique(d.te_fit), unique(d.te_residual)))))
## Only QTL found in all traits: 5@92, 5@79, 9@34

intersect(unique(d.sv_area), intersect(unique(d.water_lost), intersect(unique(d.wue_ratio), unique(d.te_fit))))
## 6 QTL found in everything but TE residual: 2@96, 3@92, 5@109, 3@79, 6@65, 7@99, 9@34

intersect(unique(d.sv_area), intersect(unique(d.wue_ratio), unique(d.te_residual)))
## 3 QTL found for TE residual, WUE ratio and Biomass: 3@48, 2@11, 9@127

no_biomass_u.qtl<-unique(c(d.water_lost, d.wue_ratio, d.te_fit, d.te_residual))
biomass_u.qtl<-unique(d.sv_area)
biomass_u.qtl[!biomass_u.qtl %in% no_biomass_u.qtl]
## QTL unique to biomass 5@119, 9@138   

no_te_residual_u.qtl<-unique(c(d.water_lost, d.wue_ratio, d.te_fit, d.sv_area))
te_residual_u.qtl<-unique(d.te_residual)
te_residual_u.qtl[!te_residual_u.qtl %in% no_te_residual_u.qtl]
## QTL unique to te_residual 1@54, 3@61, 3@77

no_wue_ratio_u.qtl<-unique(c(d.water_lost, d.te_fit, d.te_residual, d.sv_area))
wue_ratio_u.qtl<-unique(d.wue_ratio)
wue_ratio_u.qtl[!wue_ratio_u.qtl %in% no_wue_ratio_u.qtl]
## QTL unique to te_fit is 9@87 2@82


## How many QTL for each trait?
length(d.sv_area)
length(d.water_lost)
length(d.wue_ratio)
length(d.te_fit)
length(d.te_residual)



## Lets examine the intersection between rate and cumulative QTL
union(all_qtl.t.nr, all_qtl.d.nr)
## 28

intersect(all_qtl.t.nr, all_qtl.d.nr)
## 22

## Lets get the QTL specific to cumulative and rate

## QTL unique to cumulative traits
all_qtl.t.nr[!all_qtl.t.nr %in% all_qtl.d.nr]

## QTL unique to rate traits
all_qtl.d.nr[!all_qtl.d.nr %in% all_qtl.t.nr]


########################################################################
# Compare cumulative (total) vs day (rate)
########################################################################


grid.newpage()
pdf("venn_cumulative_v_rate_qtl.pdf")
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
                   c("TE fit total", "TE fit rate"),
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
                   c("TE residual total", "TE residual rate"),
                   fill=c("red", "pink"),
                   cex=2,
                   cat.cex=1.5
)

dev.off()


te_residual.td_both<-intersect(t.te_residual, d.te_residual)
te_residual.t_only<-t.te_residual[!t.te_residual %in% d.te_residual]
te_residual.d_only<-d.te_residual[!d.te_residual %in% t.te_residual]

########################################################################
# Compare treatments dry (cumulative/total) vs wet (cumulative/total)
########################################################################

## Lets make a venn diagram for each trait cumulative wet vs dry
## Get QTL unique to each trait in each condition
t.w.sv_area<-unique(t.uni.all_qtl[grepl("sv_area", t.uni.all_qtl$trait) & t.uni.all_qtl$treatment == 'wet','qtl.name'])
t.w.water_lost<-unique(t.uni.all_qtl[grepl("water_lost", t.uni.all_qtl$trait) & t.uni.all_qtl$treatment == 'wet','qtl.name'])
t.w.wue_ratio<-unique(t.uni.all_qtl[grepl("wue_total", t.uni.all_qtl$trait) & t.uni.all_qtl$treatment == 'wet','qtl.name'])
t.w.te_fit<-unique(t.uni.all_qtl[grepl("te_fit", t.uni.all_qtl$trait) & t.uni.all_qtl$treatment == 'wet','qtl.name'])
t.w.te_residual<-unique(t.uni.all_qtl[grepl("te_residual", t.uni.all_qtl$trait) & t.uni.all_qtl$treatment == 'wet','qtl.name'])

t.d.sv_area<-unique(t.uni.all_qtl[grepl("sv_area", t.uni.all_qtl$trait) & t.uni.all_qtl$treatment == 'dry','qtl.name'])
t.d.water_lost<-unique(t.uni.all_qtl[grepl("water_lost", t.uni.all_qtl$trait) & t.uni.all_qtl$treatment == 'dry','qtl.name'])
t.d.wue_ratio<-unique(t.uni.all_qtl[grepl("wue_total", t.uni.all_qtl$trait) & t.uni.all_qtl$treatment == 'dry','qtl.name'])
t.d.te_fit<-unique(t.uni.all_qtl[grepl("te_fit", t.uni.all_qtl$trait) & t.uni.all_qtl$treatment == 'dry','qtl.name'])
t.d.te_residual<-unique(t.uni.all_qtl[grepl("te_residual", t.uni.all_qtl$trait) & t.uni.all_qtl$treatment == 'dry','qtl.name'])

## Make venns of dry vs wet

grid.newpage()

pdf("venn_dry_v_wet_cumulative_qtl.pdf")
draw.pairwise.venn(area1=length(t.d.sv_area),
                   area2=length(t.w.sv_area),
                   cross.area=length(intersect(t.d.sv_area, t.w.sv_area)),
                   c("Biomass dry", "Biomass wet"),
                   fill=c("orange", "navy"),
                   cex=2,
                   cat.cex=1.5
)

#dev.off()

sv_area.t.dw_both<-intersect(t.d.sv_area, t.w.sv_area)
sv_area.t.d_only<-t.d.sv_area[!t.d.sv_area %in% t.w.sv_area]
sv_area.t.w_only<-t.w.sv_area[!t.w.sv_area %in% t.d.sv_area]


grid.newpage()
#pdf("Figure_X.pdf")
draw.pairwise.venn(area1=length(t.d.water_lost),
                   area2=length(t.w.water_lost),
                   cross.area=length(intersect(t.d.water_lost, t.w.water_lost)),
                   c("Water lost dry", "Water lost wet"),
                   fill=c("orange", "navy"),
                   cex=2,
                   cat.cex=1.5
)

#dev.off()


water_lost.t.dw_both<-intersect(t.d.water_lost, t.w.water_lost)
water_lost.t.d_only<-t.d.water_lost[!t.d.water_lost %in% t.w.water_lost]
water_lost.t.w_only<-t.w.water_lost[!t.w.water_lost %in% t.d.water_lost]



grid.newpage()
#pdf("Figure_X.pdf")
draw.pairwise.venn(area1=length(t.d.wue_ratio),
                   area2=length(t.w.wue_ratio),
                   cross.area=length(intersect(t.d.wue_ratio, t.w.wue_ratio)),
                   c("WUE ratio dry", "WUE ratio wet"),
                   fill=c("orange", "navy"),
                   cex=2,
                   cat.cex=1.5
)

#dev.off()


wue_ratio.t.dw_both<-intersect(t.d.wue_ratio, t.w.wue_ratio)
wue_ratio.t.d_only<-t.d.wue_ratio[!t.d.wue_ratio %in% t.w.wue_ratio]
wue_ratio.t.w_only<-t.w.wue_ratio[!t.w.wue_ratio %in% t.d.wue_ratio]


grid.newpage()
#pdf("Figure_X.pdf")
draw.pairwise.venn(area1=length(t.d.te_fit),
                   area2=length(t.w.te_fit),
                   cross.area=length(intersect(t.d.te_fit, t.w.te_fit)),
                   c("TE fit dry", "TE fit wet"),
                   fill=c("orange", "navy"),
                   cex=2,
                   cat.cex=1.5
)

dev.off()


te_fit.t.dw_both<-intersect(t.d.te_fit, t.w.te_fit)
te_fit.t.d_only<-t.d.te_fit[!t.d.te_fit %in% t.w.te_fit]
te_fit.t.w_only<-t.w.te_fit[!t.w.te_fit %in% t.d.te_fit]


grid.newpage()
#pdf("Figure_X.pdf")
draw.pairwise.venn(area1=length(t.d.te_residual),
                   area2=length(t.w.te_residual),
                   cross.area=length(intersect(t.d.te_residual, t.w.te_residual)),
                   c("TE residual dry", "TE residual wet"),
                   fill=c("orange", "navy"),
                   cex=2,
                   cat.cex=1.5
)

dev.off()


te_residual.t.dw_both<-intersect(t.d.te_residual, t.w.te_residual)
te_residual.t.d_only<-t.d.te_residual[!t.d.te_residual %in% t.w.te_residual]
te_residual.t.w_only<-t.w.te_residual[!t.w.te_residual %in% t.d.te_residual]


########################################################################
# Compare treatments dry (rate/day) vs wet (rate/day)
########################################################################

## Lets make a venn diagram for each trait rate wet vs dry
## Get QTL unique to each trait in each condition
d.w.sv_area<-unique(d.uni.all_qtl[grepl("sv_area", d.uni.all_qtl$trait) & d.uni.all_qtl$treatment == 'wet','qtl.name'])
d.w.water_lost<-unique(d.uni.all_qtl[grepl("water_lost", d.uni.all_qtl$trait) & d.uni.all_qtl$treatment == 'wet','qtl.name'])
d.w.wue_ratio<-unique(d.uni.all_qtl[grepl("wue_total", d.uni.all_qtl$trait) & d.uni.all_qtl$treatment == 'wet','qtl.name'])
d.w.te_fit<-unique(d.uni.all_qtl[grepl("te_fit", d.uni.all_qtl$trait) & d.uni.all_qtl$treatment == 'wet','qtl.name'])
d.w.te_residual<-unique(d.uni.all_qtl[grepl("te_residual", d.uni.all_qtl$trait) & d.uni.all_qtl$treatment == 'wet','qtl.name'])

d.d.sv_area<-unique(d.uni.all_qtl[grepl("sv_area", d.uni.all_qtl$trait) & d.uni.all_qtl$treatment == 'dry','qtl.name'])
d.d.water_lost<-unique(d.uni.all_qtl[grepl("water_lost", d.uni.all_qtl$trait) & d.uni.all_qtl$treatment == 'dry','qtl.name'])
d.d.wue_ratio<-unique(d.uni.all_qtl[grepl("wue_total", d.uni.all_qtl$trait) & d.uni.all_qtl$treatment == 'dry','qtl.name'])
d.d.te_fit<-unique(d.uni.all_qtl[grepl("te_fit", d.uni.all_qtl$trait) & d.uni.all_qtl$treatment == 'dry','qtl.name'])
d.d.te_residual<-unique(d.uni.all_qtl[grepl("te_residual", d.uni.all_qtl$trait) & d.uni.all_qtl$treatment == 'dry','qtl.name'])



grid.newpage()
#pdf("Figure_X.pdf")
draw.pairwise.venn(area1=length(d.d.sv_area),
                   area2=length(d.w.sv_area),
                   cross.area=length(intersect(d.d.sv_area, d.w.sv_area)),
                   c("Biomass dry", "Biomass wet"),
                   fill=c("orange", "navy"),
                   cex=2,
                   cat.cex=1.5
)

#dev.off()

sv_area.d.dw_both<-intersect(d.d.sv_area, d.w.sv_area)
sv_area.d.d_only<-d.d.sv_area[!d.d.sv_area %in% d.w.sv_area]
sv_area.d.w_only<-d.w.sv_area[!d.w.sv_area %in% d.d.sv_area]


grid.newpage()
#pdf("Figure_X.pdf")
draw.pairwise.venn(area1=length(d.d.water_lost),
                   area2=length(d.w.water_lost),
                   cross.area=length(intersect(d.d.water_lost, d.w.water_lost)),
                   c("Water lost dry", "Water lost wet"),
                   fill=c("orange", "navy"),
                   cex=2,
                   cat.cex=1.5
)

#dev.off()

water_lost.d.dw_both<-intersect(d.d.water_lost, d.w.water_lost)
water_lost.d.d_only<-d.d.water_lost[!d.d.water_lost %in% d.w.water_lost]
water_lost.d.w_only<-d.w.water_lost[!d.w.water_lost %in% d.d.water_lost]

grid.newpage()
#pdf("Figure_X.pdf")
draw.pairwise.venn(area1=length(d.d.wue_ratio),
                   area2=length(d.w.wue_ratio),
                   cross.area=length(intersect(d.d.wue_ratio, d.w.wue_ratio)),
                   c("WUE ratio dry", "WUE ratio wet"),
                   fill=c("orange", "navy"),
                   cex=2,
                   cat.cex=1.5
)

#dev.off()

wue_ratio.d.dw_both<-intersect(d.d.wue_ratio, d.w.wue_ratio)
wue_ratio.d.d_only<-d.d.wue_ratio[!d.d.wue_ratio %in% d.w.wue_ratio]
wue_ratio.d.w_only<-d.w.wue_ratio[!d.w.wue_ratio %in% d.d.wue_ratio]

grid.newpage()
#pdf("Figure_X.pdf")
draw.pairwise.venn(area1=length(d.d.te_fit),
                   area2=length(d.w.te_fit),
                   cross.area=length(intersect(d.d.te_fit, d.w.te_fit)),
                   c("TE fit dry", "TE fit wet"),
                   fill=c("orange", "navy"),
                   cex=2,
                   cat.cex=1.5
)

#dev.off()

te_fit.d.dw_both<-intersect(d.d.te_fit, d.w.te_fit)
te_fit.d.d_only<-d.d.te_fit[!d.d.te_fit %in% d.w.te_fit]
te_fit.d.w_only<-d.w.te_fit[!d.w.te_fit %in% d.d.te_fit]

grid.newpage()
#pdf("Figure_X.pdf")
draw.pairwise.venn(area1=length(d.d.te_residual),
                   area2=length(d.w.te_residual),
                   cross.area=length(intersect(d.d.te_residual, d.w.te_residual)),
                   c("TE residual dry", "TE residual wet"),
                   fill=c("orange", "navy"),
                   cex=2,
                   cat.cex=1.5
)

dev.off()

te_residual.d.dw_both<-intersect(d.d.te_residual, d.w.te_residual)
te_residual.d.d_only<-d.d.te_residual[!d.d.te_residual %in% d.w.te_residual]
te_residual.d.w_only<-d.w.te_residual[!d.w.te_residual %in% d.d.te_residual]

########################################################################
# Lets combine these QTL subsets into a list (by trait)
########################################################################

sv_area.qtl<-list(sv_area.td_both, sv_area.t.d_only, sv_area.t.w_only, sv_area.t.dw_both, sv_area.t.d_only, sv_area.t.w_only, sv_area.d.dw_both, sv_area.d.d_only, sv_area.d.w_only)
names(sv_area.qtl)<-c('sv_area.td_both', 'sv_area.t.d_only', 'sv_area.t.w_only', 'sv_area.t.dw_both', 'sv_area.t.d_only', 'sv_area.t.w_only', 'sv_area.d.dw_both', 'sv_area.d.d_only', 'sv_area.d.w_only')

water_lost.qtl<-list(water_lost.td_both, water_lost.t.d_only, water_lost.t.w_only, water_lost.t.dw_both, water_lost.t.d_only, water_lost.t.w_only, water_lost.d.dw_both, water_lost.d.d_only, water_lost.d.w_only)
names(water_lost.qtl)<-c('water_lost.td_both', 'water_lost.t.d_only', 'water_lost.t.w_only', 'water_lost.t.dw_both', 'water_lost.t.d_only', 'water_lost.t.w_only', 'water_lost.d.dw_both', 'water_lost.d.d_only', 'water_lost.d.w_only')

wue_ratio.qtl<-list(wue_ratio.td_both, wue_ratio.t.d_only, wue_ratio.t.w_only, wue_ratio.t.dw_both, wue_ratio.t.d_only, wue_ratio.t.w_only, wue_ratio.d.dw_both, wue_ratio.d.d_only, wue_ratio.d.w_only)
names(wue_ratio.qtl)<-c('wue_ratio.td_both', 'wue_ratio.t.d_only', 'wue_ratio.t.w_only', 'wue_ratio.t.dw_both', 'wue_ratio.t.d_only', 'wue_ratio.t.w_only', 'wue_ratio.d.dw_both', 'wue_ratio.d.d_only', 'wue_ratio.d.w_only')

te_fit.qtl<-list(te_fit.td_both, te_fit.t.d_only, te_fit.t.w_only, te_fit.t.dw_both, te_fit.t.d_only, te_fit.t.w_only, te_fit.d.dw_both, te_fit.d.d_only, te_fit.d.w_only)
names(te_fit.qtl)<-c('te_fit.td_both', 'te_fit.t.d_only', 'te_fit.t.w_only', 'te_fit.t.dw_both', 'te_fit.t.d_only', 'te_fit.t.w_only', 'te_fit.d.dw_both', 'te_fit.d.d_only', 'te_fit.d.w_only')

te_residual.qtl<-list(te_residual.td_both, te_residual.t.d_only, te_residual.t.w_only, te_residual.t.dw_both, te_residual.t.d_only, te_residual.t.w_only, te_residual.d.dw_both, te_residual.d.d_only, te_residual.d.w_only)
names(te_residual.qtl)<-c('te_residual.td_both', 'te_residual.t.d_only', 'te_residual.t.w_only', 'te_residual.t.dw_both', 'te_residual.t.d_only', 'te_residual.t.w_only', 'te_residual.d.dw_both', 'te_residual.d.d_only', 'te_residual.d.w_only')



##################################################################
# Lets get count of all markers and QTL for both cumulative and rate traits
##################################################################

## remove isotope traits

## isotope traits already removed
#traits_not_te_2<-c("gN_m2_BP14", "d13C_z_BP14", "gC_m2_BP14", "gN_m2_z_BP14", "gC_m2_z_BP14", "CN_BP14", "CN_z_BP14")
#all_total.qtl<-all_total.qtl[!all_total.qtl$trait %in% traits_not_te_2, ]
## how many unique QTL positions for raw traits?
r.all_total.qtl<-all_total.qtl[all_total.qtl$type == "raw",]
length(unique(r.all_total.qtl$marker))
## 86 SNP

## How many collapsed QTL positions?
length(unique(t.uni.all_qtl$qtl.name))
## 23 QTL

## now daily/rate

## isotope traits already removed
#all_day.qtl<-all_day.qtl[!all_day.qtl$trait %in% traits_not_te, ]
## how many unique QTL positions?
r.all_day.qtl<-all_day.qtl[all_day.qtl$type == "raw",]
length(unique(r.all_day.qtl$marker))
## 106

## How many collapsed QTL positions?
length(unique(d.uni.all_qtl$qtl.name))
## 27

## Unique in total
unique(t.uni.all_qtl$qtl.name)[!unique(t.uni.all_qtl$qtl.name) %in% unique(d.uni.all_qtl$qtl.name)]
## 6@47

## Unique in rate
unique(d.uni.all_qtl$qtl.name)[!unique(d.uni.all_qtl$qtl.name) %in% unique(t.uni.all_qtl$qtl.name)]
## 1@54 2@58 3@4 3@61 8@35

grid.newpage()
#pdf("X")
draw.pairwise.venn(area1=length(unique(t.uni.all_qtl$qtl.name)),
                   area2=length(unique(d.uni.all_qtl$qtl.name)),
                   cross.area=length(intersect(unique(t.uni.all_qtl$qtl.name), unique(d.uni.all_qtl$qtl.name))),
                   c("Cumulative", "Rate"),
                   fill=c("green", "red"),
                   cex=2,
                   cat.cex=1.5
)

#dev.off()


trait.rate_both<-intersect(unique(t.uni.all_qtl$qtl.name), unique(d.uni.all_qtl$qtl.name))
trait_only<-unique(t.uni.all_qtl$qtl.name)[!unique(t.uni.all_qtl$qtl.name) %in% unique(d.uni.all_qtl$qtl.name)]
rate_only<-unique(d.uni.all_qtl$qtl.name)[!unique(d.uni.all_qtl$qtl.name) %in% unique(t.uni.all_qtl$qtl.name)]



##################################################
## Lets compare the QTL detected by rate statistics across treatment blocks
##################################################

## Plant size
sv_area.d.dw_both
sv_area.t.dw_both

sv_area.d.d_only
sv_area.t.d_only

sv_area.d.w_only
sv_area.t.w_only

## Water loss
water_lost.d.dw_both
water_lost.t.dw_both

water_lost.d.d_only
water_lost.t.d_only

water_lost.d.w_only
water_lost.t.w_only

## WUE Ratio
wue_ratio.d.dw_both
wue_ratio.t.dw_both

wue_ratio.d.d_only
wue_ratio.t.d_only

wue_ratio.d.w_only
wue_ratio.t.w_only

## TE fit
te_fit.d.dw_both
te_fit.t.dw_both

te_fit.d.d_only
te_fit.t.d_only

te_fit.d.w_only
te_fit.t.w_only


## TE residual
te_residual.d.dw_both
te_residual.t.dw_both

te_residual.d.d_only
te_residual.t.d_only

te_residual.d.w_only
te_residual.t.w_only

setwd(wue_results.qtl.dir)

## Write out all individual SNP locations
## This is Table S2

write.csv(all.qtl, file="all_snp_locations.csv", quote=F, row.names=F)

## Write out all condensed QTL
## This is Table S3
#write.csv(uni.all_qtl, file="TABLE_S3.csv", quote=F, row.names=F)

write.csv(uni.all.qtl, file="all_condensed_qtl_locations.csv", quote=F, row.names=F)


##############################################################################################################################
# Make the plots of all QTL identified in different types of contrasts
##############################################################################################################################

######################################
# All QTL in found positions (no gxe)
######################################

## isotope traits have been removed
#r.isotope.qtl<-all.qtl[all.qtl$trait == 'd13C_BP14' & all.qtl$type == 'raw',]
#r.all_total.qtl<-r.all_total.qtl[!r.all_total.qtl$trait %in% c("d13C_BP14", "CN_BP14", "gC_m2_BP14", "gN_m2_BP14", "gN_m2_z_BP14","d13C_z_BP14 ", "gC_m2_z_BP14","CN_z_BP14"),]
#r.all_total.qtl[r.all_total.qtl$trait == 'd13C_BP14',]


## Here we extract trait names
raw_trait<-c()
for(i in 1:nrow(r.all_total.qtl)){
  fulltrait<-as.character(r.all_total.qtl[i,'trait'])
  temp<-strsplit(fulltrait, '_')
  trait.name<-temp[[1]][1:(length(temp[[1]])-2)]
  trait.name<-paste(trait.name, collapse="_")
  raw_trait<-c(raw_trait, trait.name)
}
#raw_trait<-raw_trait[!raw_trait %in% c("d13C", "CN", "gC_m2", "gN_m2")]
r.all_total.qtl$raw_trait<-raw_trait
unique(raw_trait)


r.all_total.qtl$chr<-factor(r.all_total.qtl$chr, levels=c(1,2,3,4,5,6,7,8,9))

fx.size<-r.all_total.qtl$additive.fx
fx.size<-as.numeric(as.character(fx.size))

plot.char<-c()
for(i in 1:length(fx.size)){
  if (fx.size[i] > 0) {plot.char<-c(plot.char, '24')}
  if (fx.size[i] < 0) {plot.char<-c(plot.char, '25')}
}

r.all_total.qtl$plot.char<-plot.char
r.all_total.qtl$plot.char<-as.factor(r.all_total.qtl$plot.char)

r.all_total.qtl$raw_trait.treat<-paste(r.all_total.qtl$raw_trait, r.all_total.qtl$treatment, sep=".")

## Set factor levels
r.all_total.qtl$raw_trait.treat<-factor(r.all_total.qtl$raw_trait.treat, levels=c('sv_area_total.dry','sv_area_total.wet','water_lost_total.dry','water_lost_total.wet','wue_total.dry','wue_total.wet','te_fit_total.dry','te_fit_total.wet','te_residual_total.dry','te_residual_total.wet'))

raw_trait.treat<-as.character(r.all_total.qtl$raw_trait.treat)
raw_trait.treat_name<-unique(raw_trait.treat)

## Assign factor levels
raw_trait.treat_name<-c('sv_area_total.dry','sv_area_total.wet','water_lost_total.dry','water_lost_total.wet','wue_total.dry','wue_total.wet','te_fit_total.dry','te_fit_total.wet','te_residual_total.dry','te_residual_total.wet')


plot.col<-c()
for(i in 1:length(raw_trait.treat)){
  logical<-raw_trait.treat[i] == raw_trait.treat_name
  col<-which(logical, arr.ind=TRUE)
  plot.col<-c(plot.col, col)
}

r.all_total.qtl$plot.col<-plot.col


setwd(wue_results.qtl.dir)

## This is Figure_4

pdf("cumulative_snp_locations_no_gxe.pdf")
p<-ggplot() + geom_point(data = r.all_total.qtl, aes(x = pos, y = prop.var, shape=plot.char, colour=as.character(plot.col), fill=as.character(plot.col)),size=3, alpha=0.5) + geom_blank(data = blank_data, aes(x = x, y = y)) + facet_wrap(~chr, scales = "free_x") + expand_limits(x = 0) + scale_x_continuous(expand = c(0, 0)) + theme_bw() +   scale_shape_manual(values=c(24,25)) 
print(p + scale_color_manual(values=c("1" = "light green", "2" = "dark green", "3" = "light blue", "4" = "blue", "5" = "yellow", "6" = "orange", "7" = "grey", "8" = "black", "9" = "pink", "10" = "red", "11" = "purple")) + scale_fill_manual(values=c("1" = "green", "2" = "dark green", "3" = "light blue", "4" = "navy", "5" = "yellow", "6" = "orange", "7" = "grey", "8" = "black", "9" = "pink", "10" = "red", "11" = "purple")) + ylab("% Variance")  + xlab("Genome Position") + theme(legend.position = "none"))

dev.off()

## Make FIG_4 again with stable isotope data plotted on it

#fx.size.i<-r.isotope.qtl$additive.fx
#fx.size.i<-as.numeric(as.character(fx.size.i))
#plot.char.i<-c()
#for(i in 1:length(fx.size.i)){
#  if (fx.size.i[i] > 0) {plot.char.i<-c(plot.char.i, '24')}
#  if (fx.size.i[i] < 0) {plot.char.i<-c(plot.char.i, '25')}
#}

#r.isotope.qtl$plot.char<-plot.char.i
#r.isotope.qtl$plot.char<-as.factor(r.isotope.qtl$plot.char)

#pdf("FIG_4_with_d13C.pdf")
#p<-ggplot() + geom_point(data = r.all_total.qtl, aes(x = pos, y = prop.var, shape=plot.char, colour=as.character(plot.col), fill=as.character(plot.col)),size=3, alpha=0.5) + geom_blank(data = blank_data, aes(x = x, y = y)) + facet_wrap(~chr, scales = "free_x") + expand_limits(x = 0) + scale_x_continuous(expand = c(0, 0)) + theme_bw() +   scale_shape_manual(values=c(24,25)) 
#q<-p + scale_color_manual(values=c("1" = "light green", "2" = "dark green", "3" = "light blue", "4" = "blue", "5" = "yellow", "6" = "orange", "7" = "grey", "8" = "black", "9" = "pink", "10" = "red", "11" = "purple")) + scale_fill_manual(values=c("1" = "green", "2" = "dark green", "3" = "light blue", "4" = "navy", "5" = "yellow", "6" = "orange", "7" = "grey", "8" = "black", "9" = "pink", "10" = "red", "11" = "purple")) + ylab("% Variance")  + xlab("Genome Position") + theme(legend.position = "none")
#y<-q+geom_point(data=r.isotope.qtl, aes(x=pos, y=prop.var, shape=plot.char), colour="purple", fill="purple",size=5, alpha=0.5)
#y
#dev.off()


######################################
# All QTL in found positions
######################################
t.uni.all_qtl<-uni.all_qtl[grepl('_total', uni.all_qtl$trait),]
#t.uni.isotope.qtl<-uni.all_qtl[uni.all_qtl$trait == 'd13C_BP14',]

## Here we extract trait names
raw_trait<-c()
for(i in 1:nrow(t.uni.all_qtl)){
  fulltrait<-as.character(t.uni.all_qtl[i,'trait'])
  temp<-strsplit(fulltrait, '_')
  trait.name<-temp[[1]][1:(length(temp[[1]])-2)]
  trait.name<-paste(trait.name, collapse="_")
  raw_trait<-c(raw_trait, trait.name)
}

t.uni.all_qtl$raw_trait<-raw_trait
unique(raw_trait)


t.uni.all_qtl$chr<-factor(t.uni.all_qtl$chr, levels=c(1,2,3,4,5,6,7,8,9))

fx.size<-t.uni.all_qtl$additive.fx
fx.size<-as.numeric(as.character(fx.size))

plot.char<-c()
for(i in 1:length(fx.size)){
  if (fx.size[i] > 0) {plot.char<-c(plot.char, '24')}
  if (fx.size[i] < 0) {plot.char<-c(plot.char, '25')}
}

t.uni.all_qtl$plot.char<-plot.char
t.uni.all_qtl$plot.char<-as.factor(t.uni.all_qtl$plot.char)

t.uni.all_qtl$raw_trait.treat<-paste(t.uni.all_qtl$raw_trait, t.uni.all_qtl$treatment, sep=".")

## Set factor levels
t.uni.all_qtl$raw_trait.treat<-factor(t.uni.all_qtl$raw_trait.treat, levels=c('sv_area_total.dry','sv_area_total.wet','water_lost_total.dry','water_lost_total.wet','wue_total.dry','wue_total.wet','te_fit_total.dry','te_fit_total.wet','te_residual_total.dry','te_residual_total.wet'))

raw_trait.treat<-as.character(t.uni.all_qtl$raw_trait.treat)

raw_trait.treat_name<-unique(raw_trait.treat)
raw_trait.treat_name<-c('sv_area_total.dry','sv_area_total.wet','water_lost_total.dry','water_lost_total.wet','wue_total.dry','wue_total.wet','te_fit_total.dry','te_fit_total.wet','te_residual_total.dry','te_residual_total.wet')


plot.col<-c()
for(i in 1:length(raw_trait.treat)){
  logical<-raw_trait.treat[i] == raw_trait.treat_name
  col<-which(logical, arr.ind=TRUE)
  plot.col<-c(plot.col, col)
}

t.uni.all_qtl$plot.col<-plot.col


setwd(wue_results.qtl.dir)


pdf("cumulative_qtl_locations_no_gxe_condensed.pdf")
p<-ggplot() + geom_point(data = t.uni.all_qtl, aes(x = pos, y = prop.var, shape=plot.char, colour=as.character(plot.col), fill=as.character(plot.col)),size=3, alpha=0.5) + geom_blank(data = blank_data, aes(x = x, y = y)) + facet_wrap(~chr, scales = "free_x") + expand_limits(x = 0) + scale_x_continuous(expand = c(0, 0)) + theme_bw() +   scale_shape_manual(values=c(24,25)) 
print(p + scale_color_manual(values=c("1" = "light green", "2" = "dark green", "3" = "light blue", "4" = "blue", "5" = "yellow", "6" = "orange", "7" = "grey", "8" = "black", "9" = "pink", "10" = "red")) + scale_fill_manual(values=c("1" = "light green", "2" = "dark green", "3" = "light blue", "4" = "blue", "5" = "yellow", "6" = "orange", "7" = "grey", "8" = "black", "9" = "pink", "10" = "red")) + ylab("% Variance")  + xlab("Genome Position") + theme(legend.position = "none"))

dev.off()

## Lets make another plot with d13C

#fx.size.i<-t.uni.isotope.qtl$additive.fx
#fx.size.i<-as.numeric(as.character(fx.size.i))

#plot.char.i<-c()
#for(i in 1:length(fx.size.i)){
#  if (fx.size.i[i] > 0) {plot.char.i<-c(plot.char.i, '24')}
#  if (fx.size.i[i] < 0) {plot.char.i<-c(plot.char.i, '25')}
#}

#t.uni.isotope.qtl$plot.char<-plot.char.i
#t.uni.isotope.qtl$plot.char<-as.factor(t.uni.isotope.qtl$plot.char)


## This is Fig_S12 with d13C included
#pdf("FIG_S12_with_d13C.pdf")
#p<-ggplot() + geom_point(data = t.uni.all_qtl, aes(x = pos, y = prop.var, shape=plot.char, colour=as.character(plot.col), fill=as.character(plot.col)),size=3, alpha=0.5) + geom_blank(data = blank_data, aes(x = x, y = y)) + facet_wrap(~chr, scales = "free_x") + expand_limits(x = 0) + scale_x_continuous(expand = c(0, 0)) + theme_bw() +   scale_shape_manual(values=c(24,25)) 
#q<-p + scale_color_manual(values=c("1" = "light green", "2" = "dark green", "3" = "light blue", "4" = "blue", "5" = "yellow", "6" = "orange", "7" = "grey", "8" = "black", "9" = "pink", "10" = "red")) + scale_fill_manual(values=c("1" = "light green", "2" = "dark green", "3" = "light blue", "4" = "blue", "5" = "yellow", "6" = "orange", "7" = "grey", "8" = "black", "9" = "pink", "10" = "red")) + ylab("% Variance")  + xlab("Genome Position") + theme(legend.position = "none")
#y<-q+geom_point(data=t.uni.isotope.qtl, aes(x=pos, y=prop.var, shape=plot.char), colour="purple", fill="purple",size=5, alpha=0.5)
#y

#dev.off()


r.all_day.qtl<-all_day.qtl[all_day.qtl$type == "raw",]

#####################################################################
# Now for the daily/rate statistics
#####################################################################

## Here we extract trait names
raw_trait<-c()
for(i in 1:nrow(r.all_day.qtl)){
  fulltrait<-as.character(r.all_day.qtl[i,'trait'])
  temp<-strsplit(fulltrait, '_')
  trait.name<-temp[[1]][1:(length(temp[[1]])-2)]
  trait.name<-paste(trait.name, collapse="_")
  raw_trait<-c(raw_trait, trait.name)
}

r.all_day.qtl$raw_trait<-raw_trait
unique(raw_trait)


r.all_day.qtl$chr<-factor(r.all_day.qtl$chr, levels=c(1,2,3,4,5,6,7,8,9))

fx.size<-r.all_day.qtl$additive.fx
fx.size<-as.numeric(as.character(fx.size))

plot.char<-c()
for(i in 1:length(fx.size)){
  if (fx.size[i] > 0) {plot.char<-c(plot.char, '24')}
  if (fx.size[i] < 0) {plot.char<-c(plot.char, '25')}
}

r.all_day.qtl$plot.char<-plot.char
r.all_day.qtl$plot.char<-as.factor(r.all_day.qtl$plot.char)

r.all_day.qtl$raw_trait.treat<-paste(r.all_day.qtl$raw_trait, r.all_day.qtl$treatment, sep=".")

## Set factor levels
r.all_day.qtl$raw_trait.treat<-factor(r.all_day.qtl$raw_trait.treat, levels=c('sv_area_day.dry','sv_area_day.wet','water_lost.dry','water_lost.wet','wue_day.dry','wue_day.wet','te_fit_day.dry','te_fit_day.wet','te_residual_day.dry','te_residual_day.wet'))

raw_trait.treat<-as.character(r.all_day.qtl$raw_trait.treat)
raw_trait.treat_name<-unique(raw_trait.treat)

## Assign factor levels
raw_trait.treat_name<-c('sv_area_day.dry','sv_area_day.wet','water_lost.dry','water_lost.wet','wue_day.dry','wue_day.wet','te_fit_day.dry','te_fit_day.wet','te_residual_day.dry','te_residual_day.wet')

plot.col<-c()
for(i in 1:length(raw_trait.treat)){
  logical<-raw_trait.treat[i] == raw_trait.treat_name
  col<-which(logical, arr.ind=TRUE)
  plot.col<-c(plot.col, col)
}

r.all_day.qtl$plot.col<-plot.col



## This is Figure_S13
pdf("rate_snp_locations_no_gxe.pdf")
p<-ggplot() + geom_point(data = r.all_day.qtl, aes(x = pos, y = prop.var, shape=plot.char, colour=as.character(plot.col), fill=as.character(plot.col)),size=3, alpha=0.5) + geom_blank(data = blank_data, aes(x = x, y = y)) + facet_wrap(~chr, scales = "free_x") + expand_limits(x = 0) + scale_x_continuous(expand = c(0, 0)) + theme_bw() +   scale_shape_manual(values=c(24,25)) 
print(p + scale_color_manual(values=c("1" = "light green", "2" = "dark green", "3" = "light blue", "4" = "blue", "5" = "yellow", "6" = "orange", "7" = "grey", "8" = "black", "9" = "pink", "10" = "red")) + scale_fill_manual(values=c("1" = "green", "2" = "dark green", "3" = "light blue", "4" = "navy", "5" = "yellow", "6" = "orange", "7" = "grey", "8" = "black", "9" = "pink", "10" = "red")) + ylab("% Variance")  + xlab("Genome Position") + theme(legend.position = "none"))

dev.off()


######################################
# All QTL in found positions
######################################


## Here we extract trait names
raw_trait<-c()
for(i in 1:nrow(d.uni.all_qtl)){
  fulltrait<-as.character(d.uni.all_qtl[i,'trait'])
  temp<-strsplit(fulltrait, '_')
  trait.name<-temp[[1]][1:(length(temp[[1]])-2)]
  trait.name<-paste(trait.name, collapse="_")
  raw_trait<-c(raw_trait, trait.name)
}

d.uni.all_qtl$raw_trait<-raw_trait
unique(raw_trait)


d.uni.all_qtl$chr<-factor(d.uni.all_qtl$chr, levels=c(1,2,3,4,5,6,7,8,9))

fx.size<-d.uni.all_qtl$additive.fx
fx.size<-as.numeric(as.character(fx.size))

plot.char<-c()
for(i in 1:length(fx.size)){
  if (fx.size[i] > 0) {plot.char<-c(plot.char, '24')}
  if (fx.size[i] < 0) {plot.char<-c(plot.char, '25')}
}

d.uni.all_qtl$plot.char<-plot.char
d.uni.all_qtl$plot.char<-as.factor(d.uni.all_qtl$plot.char)

d.uni.all_qtl$raw_trait.treat<-paste(d.uni.all_qtl$raw_trait, d.uni.all_qtl$treatment, sep=".")

## Set factor levels
d.uni.all_qtl$raw_trait.treat<-factor(d.uni.all_qtl$raw_trait.treat, levels=c('sv_area_day.dry','sv_area_day.wet','water_lost.dry','water_lost.wet','wue_day.dry','wue_day.wet','te_fit_day.dry','te_fit_day.wet','te_residual_day.dry','te_residual_day.wet'))

raw_trait.treat<-as.character(d.uni.all_qtl$raw_trait.treat)

raw_trait.treat_name<-unique(raw_trait.treat)
raw_trait.treat_name<-c('sv_area_day.dry','sv_area_day.wet','water_lost.dry','water_lost.wet','wue_day.dry','wue_day.wet','te_fit_day.dry','te_fit_day.wet','te_residual_day.dry','te_residual_day.wet')


plot.col<-c()
for(i in 1:length(raw_trait.treat)){
  logical<-raw_trait.treat[i] == raw_trait.treat_name
  col<-which(logical, arr.ind=TRUE)
  plot.col<-c(plot.col, col)
}

d.uni.all_qtl$plot.col<-plot.col

## This is Fig_S14
pdf("rate_qtl_locations_no_gxe_condensed.pdf")
p<-ggplot() + geom_point(data = d.uni.all_qtl, aes(x = pos, y = prop.var, shape=plot.char, colour=as.character(plot.col), fill=as.character(plot.col)),size=3, alpha=0.5) + geom_blank(data = blank_data, aes(x = x, y = y)) + facet_wrap(~chr, scales = "free_x") + expand_limits(x = 0) + scale_x_continuous(expand = c(0, 0)) + theme_bw() +   scale_shape_manual(values=c(24,25)) 
print(p + scale_color_manual(values=c("1" = "light green", "2" = "dark green", "3" = "light blue", "4" = "blue", "5" = "yellow", "6" = "orange", "7" = "grey", "8" = "black", "9" = "pink", "10" = "red")) + scale_fill_manual(values=c("1" = "light green", "2" = "dark green", "3" = "light blue", "4" = "blue", "5" = "yellow", "6" = "orange", "7" = "grey", "8" = "black", "9" = "pink", "10" = "red")) + ylab("% Variance")  + xlab("Genome Position") + theme(legend.position = "none"))


dev.off()

setwd(wue_results.qtl.venn.dir)

save.image("venn_qtl_analysis.Rdata")

rm(list=ls())


