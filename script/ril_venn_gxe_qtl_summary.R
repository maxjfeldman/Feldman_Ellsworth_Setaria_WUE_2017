
##################################################################
library(ggplot2)
library(stringr)
library(VennDiagram)


##################################################################
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


setwd(wue_results.qtl.venn.dir)
load("venn_qtl_analysis.Rdata")

##################################################################
## We will start by analyzing all the GxE QTL 
## trait difference, trait relative difference and trait ratio
##################################################################

## Lets examine the GxE SNPs and QTL 

## How many SNPs for each type of GxE trait
t.all.qtl.with_diff_ratio<-all.qtl.with_diff_ratio[grepl("_total", all.qtl.with_diff_ratio$trait),]

## Mathematical difference
length(unique(t.all.qtl.with_diff_ratio[t.all.qtl.with_diff_ratio$type == 'comp_diff','marker']))
## 43

## Relative difference
length(unique(t.all.qtl.with_diff_ratio[t.all.qtl.with_diff_ratio$type == 'comp_rel_diff','marker']))
## 40

## Trait ratio
length(unique(t.all.qtl.with_diff_ratio[t.all.qtl.with_diff_ratio$type == 'comp_ratio','marker']))
## 104

## All GxE formulations
length(unique(t.all.qtl.with_diff_ratio[all.qtl.with_diff_ratio$type == 'comp_diff' | all.qtl.with_diff_ratio$type == 'comp_rel_diff' | all.qtl.with_diff_ratio$type == 'comp_ratio','marker']))
## 148


comp_diff_snp<-unique(t.all.qtl.with_diff_ratio[t.all.qtl.with_diff_ratio$type == 'comp_diff','marker'])
comp_rel_diff_snp<-unique(t.all.qtl.with_diff_ratio[t.all.qtl.with_diff_ratio$type == 'comp_rel_diff','marker'])
comp_ratio_snp<-unique(t.all.qtl.with_diff_ratio[t.all.qtl.with_diff_ratio$type == 'comp_ratio','marker'])


setwd(wue_results.qtl.venn.dir)

grid.newpage()

## FIG_S19

pdf("venn_comparing_gxe_trait_type_qtl_numbers.pdf")
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

## First lets unify them into QTL
uni.all.qtl.with_diff_ratio<-unify_marker(all.qtl.with_diff_ratio)
cond.all.qtl.with_diff_ratio<-condense_qtl(all.qtl.with_diff_ratio)


## Extract day after planting from trait field
dap_i<-c()
for(i in 1:nrow(uni.all.qtl.with_diff_ratio)){
  x<-strsplit(as.character(uni.all.qtl.with_diff_ratio[i,'trait']), '_')[[1]]
  d<-x[length(x)-1]
  dap_i<-c(dap_i, d)
}

## Convert dap_i vector to numeric
dap_i<-as.numeric(dap_i)
## Stable isotopes were only measured at a single time point so day value is 'NA'
## Lets replace that value with 30 DAP (should be no isotope QTL left in the set by this point)
dap_i[is.na(dap_i)]<-c(30)

## Add DAP
uni.all.qtl.with_diff_ratio<-cbind(uni.all.qtl.with_diff_ratio, dap_i)
uni.all.qtl.with_diff_ratio$dap_i<-as.numeric(as.character(uni.all.qtl.with_diff_ratio$dap_i))
uni.all.qtl.with_diff_ratio$qtl.name<-paste(uni.all.qtl.with_diff_ratio$chr, as.integer(uni.all.qtl.with_diff_ratio$pos), sep="@")

## Remove early days (already removed in prior script)
uni.all.qtl.with_diff_ratio<-uni.all.qtl.with_diff_ratio[uni.all.qtl.with_diff_ratio$dap_i > 16,]

## Lets take the comparison/GxE QTL out into a new dataframe

c.uni.all.qtl.with_diff_ratio<-uni.all.qtl.with_diff_ratio[uni.all.qtl.with_diff_ratio$type != 'raw',]


## Extract the simplest trait names
full_trait_names<-c.uni.all.qtl.with_diff_ratio$trait
raw_trait<-c()
for(i in 1:length(full_trait_names)){
  name<-full_trait_names[i]
  fields<-strsplit(as.character(name), "_")
  r.trait<-fields[[1]][1:3]
  r.trait<-paste(r.trait[1], r.trait[2], r.trait[3], sep="_")
  raw_trait<-c(raw_trait, r.trait)
}

## Simplify trait names
raw_trait[grepl("sv_area", raw_trait)]<-c("sv_area")
raw_trait[grepl("water_lost", raw_trait)]<-c("water_lost")
raw_trait[grepl("wue", raw_trait)]<-c("wue")
raw_trait[grepl("te_fit", raw_trait)]<-c("fit")
raw_trait[grepl("te_residual", raw_trait)]<-c("residual")

## Add this base trait name into the table
c.uni.all.qtl.with_diff_ratio<-cbind(c.uni.all.qtl.with_diff_ratio, raw_trait)

## Lets examine the cumulative/total traits
c.uni.all.qtl.with_diff_ratio.total<-c.uni.all.qtl.with_diff_ratio[grepl("total", c.uni.all.qtl.with_diff_ratio$trait),]



## How many unique QTL for comparison traits
length(unique(c.uni.all.qtl.with_diff_ratio.total$marker))
## 22

## Lets list them...
unique(c.uni.all.qtl.with_diff_ratio.total$marker)


## How many QTL for each type of comparison trait
length(unique(c.uni.all.qtl.with_diff_ratio.total[c.uni.all.qtl.with_diff_ratio.total$type == 'comp_diff','marker']))
## 22
length(unique(c.uni.all.qtl.with_diff_ratio.total[c.uni.all.qtl.with_diff_ratio.total$type == 'comp_rel_diff','marker']))
## 20
length(unique(c.uni.all.qtl.with_diff_ratio.total[c.uni.all.qtl.with_diff_ratio.total$type == 'comp_ratio','marker']))
## 38

comp_diff_qtl<-unique(c.uni.all.qtl.with_diff_ratio.total[c.uni.all.qtl.with_diff_ratio.total$type == 'comp_diff','marker'])
comp_rel_diff_qtl<-unique(c.uni.all.qtl.with_diff_ratio.total[c.uni.all.qtl.with_diff_ratio.total$type == 'comp_rel_diff','marker'])
comp_ratio_qtl<-unique(c.uni.all.qtl.with_diff_ratio.total[c.uni.all.qtl.with_diff_ratio.total$type == 'comp_ratio','marker'])



grid.newpage()
#pdf("Figure_X.pdf")
draw.triple.venn(area1=length(comp_diff_qtl),
                 area2=length(comp_rel_diff_qtl),
                 area3=length(comp_ratio_qtl),
                 n12=length(intersect(comp_diff_qtl, comp_rel_diff_qtl)),
                 n23=length(intersect(comp_rel_diff_qtl, comp_ratio_qtl)),
                 n13=length(intersect(comp_diff_qtl, comp_ratio_qtl)),
                 n123=length(intersect(comp_diff_qtl, intersect(comp_rel_diff_qtl, comp_ratio_qtl))),
                 c("Diff", "Relative Diff", "Ratio"),
                 fill=c("red", "blue", "yellow")
)

#dev.off()



c.uni.sv_area.qtl.with_diff_ratio.total<-c.uni.all.qtl.with_diff_ratio.total[grepl("sv_area", c.uni.all.qtl.with_diff_ratio.total$trait),]
p<-ggplot(c.uni.sv_area.qtl.with_diff_ratio.total, aes(x=dap_i, y=lod, colour=treatment)) + geom_line() + facet_wrap(~qtl.name)

## Water
c.uni.water.qtl.with_diff_ratio.total<-c.uni.all.qtl.with_diff_ratio.total[grepl("water_lost", c.uni.all.qtl.with_diff_ratio.total$trait),]


## Lets compare the uniqueness of QTL between biomass and water (Cumulative/total)
intersect(unique(c.uni.sv_area.qtl.with_diff_ratio.total$qtl.name), unique(c.uni.water.qtl.with_diff_ratio.total$qtl.name))
unique(c.uni.sv_area.qtl.with_diff_ratio.total$qtl.name)[!unique(c.uni.sv_area.qtl.with_diff_ratio.total$qtl.name) %in% unique(c.uni.water.qtl.with_diff_ratio.total$qtl.name)]
unique(c.uni.water.qtl.with_diff_ratio.total$qtl.name)[!unique(c.uni.water.qtl.with_diff_ratio.total$qtl.name) %in% unique(c.uni.sv_area.qtl.with_diff_ratio.total$qtl.name)]



################################################################
# Break out daily (rate) meaurements
################################################################

## Get rate/day out of unified set
c.uni.all.qtl.with_diff_ratio.day<-c.uni.all.qtl.with_diff_ratio[!grepl('_total', c.uni.all.qtl.with_diff_ratio$trait),]
unique(c.uni.all.qtl.with_diff_ratio.day$qtl.name)
length(unique(c.uni.all.qtl.with_diff_ratio.day$qtl.name))
# 46

## How many QTL for each type of comparison trait
length(unique(c.uni.all.qtl.with_diff_ratio.day[c.uni.all.qtl.with_diff_ratio.day$type == 'comp_diff','marker']))
## 25
length(unique(c.uni.all.qtl.with_diff_ratio.day[c.uni.all.qtl.with_diff_ratio.day$type == 'comp_rel_diff','marker']))
## 25
length(unique(c.uni.all.qtl.with_diff_ratio.day[c.uni.all.qtl.with_diff_ratio.day$type == 'comp_ratio','marker']))
## 36



## Lets look at the over lap of QTL between types

########################################################################
# Cumulative/total
########################################################################

## Biomass
temp<-c.uni.all.qtl.with_diff_ratio.total[c.uni.all.qtl.with_diff_ratio.total$raw_trait == 'sv_area',]
sv.area_diff.qtl<-unique(temp[temp$type == 'comp_diff','qtl.name'])
sv.area_rel_diff.qtl<-unique(temp[temp$type == 'comp_rel_diff','qtl.name'])
sv.area_ratio.qtl<-unique(temp[temp$type == 'comp_ratio','qtl.name'])

setwd(wue_results.qtl.venn.dir)

#pdf("X.pdf")
grid.newpage()
g<-draw.triple.venn(area1=length(sv.area_diff.qtl),
                    area2=length(sv.area_rel_diff.qtl),
                    area3=length(sv.area_ratio.qtl),
                    n12=length(intersect(sv.area_diff.qtl, sv.area_rel_diff.qtl)),
                    n23=length(intersect(sv.area_rel_diff.qtl, sv.area_ratio.qtl)),
                    n13=length(intersect(sv.area_diff.qtl, sv.area_ratio.qtl)),
                    n123=length(intersect(sv.area_diff.qtl, intersect(sv.area_rel_diff.qtl, sv.area_ratio.qtl))),
                    c("Diff", "Relative Diff", "Ratio"),
                    fill=c("red", "blue", "yellow")
)
require(gridExtra)
grid.arrange(gTree(children=g), top="Plant area")

#dev.off()

## Water lost
temp<-c.uni.all.qtl.with_diff_ratio.total[c.uni.all.qtl.with_diff_ratio.total$raw_trait == 'water_lost',]
water_lost_diff.qtl<-unique(temp[temp$type == 'comp_diff','qtl.name'])
water_lost_rel_diff.qtl<-unique(temp[temp$type == 'comp_rel_diff','qtl.name'])
water_lost_ratio.qtl<-unique(temp[temp$type == 'comp_ratio','qtl.name'])


grid.newpage()
#pdf("Figure_X.pdf")
g<-draw.triple.venn(area1=length(water_lost_diff.qtl),
                    area2=length(water_lost_rel_diff.qtl),
                    area3=length(water_lost_ratio.qtl),
                    n12=length(intersect(water_lost_diff.qtl, water_lost_rel_diff.qtl)),
                    n23=length(intersect(water_lost_rel_diff.qtl, water_lost_ratio.qtl)),
                    n13=length(intersect(water_lost_diff.qtl, water_lost_ratio.qtl)),
                    n123=length(intersect(water_lost_diff.qtl, intersect(water_lost_rel_diff.qtl, water_lost_ratio.qtl))),
                    c("Diff", "Relative Diff", "Ratio"),
                    fill=c("red", "blue", "yellow")
)


require(gridExtra)
grid.arrange(gTree(children=g), top="Water Lost")

#dev.off()


## WUE ratio
temp<-c.uni.all.qtl.with_diff_ratio.total[c.uni.all.qtl.with_diff_ratio.total$raw_trait == 'wue',]
wue_ratio_diff.qtl<-unique(temp[temp$type == 'comp_diff','qtl.name'])
wue_ratio_rel_diff.qtl<-unique(temp[temp$type == 'comp_rel_diff','qtl.name'])
wue_ratio_ratio.qtl<-unique(temp[temp$type == 'comp_ratio','qtl.name'])


grid.newpage()
#pdf("Figure_X.pdf")
g<-draw.triple.venn(area1=length(wue_ratio_diff.qtl),
                    area2=length(wue_ratio_rel_diff.qtl),
                    area3=length(wue_ratio_ratio.qtl),
                    n12=length(intersect(wue_ratio_diff.qtl, wue_ratio_rel_diff.qtl)),
                    n23=length(intersect(wue_ratio_rel_diff.qtl, wue_ratio_ratio.qtl)),
                    n13=length(intersect(wue_ratio_diff.qtl, wue_ratio_ratio.qtl)),
                    n123=length(intersect(wue_ratio_diff.qtl, intersect(wue_ratio_rel_diff.qtl, wue_ratio_ratio.qtl))),
                    c("Diff", "Relative Diff", "Ratio"),
                    fill=c("red", "blue", "yellow")
)

require(gridExtra)
grid.arrange(gTree(children=g), top="WUE Ratio")
#dev.off()

## TE fit
temp<-c.uni.all.qtl.with_diff_ratio.total[c.uni.all.qtl.with_diff_ratio.total$raw_trait == 'fit',]
te_fit_diff.qtl<-unique(temp[temp$type == 'comp_diff','qtl.name'])
te_fit_rel_diff.qtl<-unique(temp[temp$type == 'comp_rel_diff','qtl.name'])
te_fit_ratio.qtl<-unique(temp[temp$type == 'comp_ratio','qtl.name'])


grid.newpage()
#pdf("Figure_X.pdf")
g<-draw.triple.venn(area1=length(te_fit_diff.qtl),
                    area2=length(te_fit_rel_diff.qtl),
                    area3=length(te_fit_ratio.qtl),
                    n12=length(intersect(te_fit_diff.qtl, te_fit_rel_diff.qtl)),
                    n23=length(intersect(te_fit_rel_diff.qtl, te_fit_ratio.qtl)),
                    n13=length(intersect(te_fit_diff.qtl, te_fit_ratio.qtl)),
                    n123=length(intersect(te_fit_diff.qtl, intersect(te_fit_rel_diff.qtl, te_fit_ratio.qtl))),
                    c("Diff", "Relative Diff", "Ratio"),
                    fill=c("red", "blue", "yellow")
)
require(gridExtra)
grid.arrange(gTree(children=g), top="TE_model_fit")

#dev.off()

## TE residual
temp<-c.uni.all.qtl.with_diff_ratio.total[c.uni.all.qtl.with_diff_ratio.total$raw_trait == 'residual',]
te_residual_diff.qtl<-unique(temp[temp$type == 'comp_diff','qtl.name'])
te_residual_rel_diff.qtl<-unique(temp[temp$type == 'comp_rel_diff','qtl.name'])
te_residual_ratio.qtl<-unique(temp[temp$type == 'comp_ratio','qtl.name'])


grid.newpage()
#pdf("Figure_X.pdf")
g<-draw.triple.venn(area1=length(te_residual_diff.qtl),
                    area2=length(te_residual_rel_diff.qtl),
                    area3=length(te_residual_ratio.qtl),
                    n12=length(intersect(te_residual_diff.qtl, te_residual_rel_diff.qtl)),
                    n23=length(intersect(te_residual_rel_diff.qtl, te_residual_ratio.qtl)),
                    n13=length(intersect(te_residual_diff.qtl, te_residual_ratio.qtl)),
                    n123=length(intersect(te_residual_diff.qtl, intersect(te_residual_rel_diff.qtl, te_residual_ratio.qtl))),
                    c("Diff", "Relative Diff", "Ratio"),
                    fill=c("red", "blue", "yellow")
)
require(gridExtra)
grid.arrange(gTree(children=g), top="TE_model_residual")


#dev.off()


########################################################################
# Daily/rate
########################################################################


## Biomass
temp<-c.uni.all.qtl.with_diff_ratio.day[c.uni.all.qtl.with_diff_ratio.day$raw_trait == 'sv_area',]
sv.area_diff.qtl<-unique(temp[temp$type == 'comp_diff','qtl.name'])
sv.area_rel_diff.qtl<-unique(temp[temp$type == 'comp_rel_diff','qtl.name'])
sv.area_ratio.qtl<-unique(temp[temp$type == 'comp_ratio','qtl.name'])


#pdf("FIG_X.pdf")
grid.newpage()
draw.triple.venn(area1=length(sv.area_diff.qtl),
                 area2=length(sv.area_rel_diff.qtl),
                 area3=length(sv.area_ratio.qtl),
                 n12=length(intersect(sv.area_diff.qtl, sv.area_rel_diff.qtl)),
                 n23=length(intersect(sv.area_rel_diff.qtl, sv.area_ratio.qtl)),
                 n13=length(intersect(sv.area_diff.qtl, sv.area_ratio.qtl)),
                 n123=length(intersect(sv.area_diff.qtl, intersect(sv.area_rel_diff.qtl, sv.area_ratio.qtl))),
                 c("Diff", "Relative Diff", "Ratio"),
                 fill=c("red", "blue", "yellow")
)

#dev.off()

## Water lost
temp<-c.uni.all.qtl.with_diff_ratio.day[c.uni.all.qtl.with_diff_ratio.day$raw_trait == 'water_lost',]
temp<-temp[!grepl('_total', temp$trait),]
water_lost_diff.qtl<-unique(temp[temp$type == 'comp_diff','qtl.name'])
water_lost_rel_diff.qtl<-unique(temp[temp$type == 'comp_rel_diff','qtl.name'])
water_lost_ratio.qtl<-unique(temp[temp$type == 'comp_ratio','qtl.name'])


grid.newpage()
#pdf("Figure_X.pdf")
draw.triple.venn(area1=length(water_lost_diff.qtl),
                 area2=length(water_lost_rel_diff.qtl),
                 area3=length(water_lost_ratio.qtl),
                 n12=length(intersect(water_lost_diff.qtl, water_lost_rel_diff.qtl)),
                 n23=length(intersect(water_lost_rel_diff.qtl, water_lost_ratio.qtl)),
                 n13=length(intersect(water_lost_diff.qtl, water_lost_ratio.qtl)),
                 n123=length(intersect(water_lost_diff.qtl, intersect(water_lost_rel_diff.qtl, water_lost_ratio.qtl))),
                 c("Diff", "Relative Diff", "Ratio"),
                 fill=c("red", "blue", "yellow")
)

#dev.off()


## WUE ratio
temp<-c.uni.all.qtl.with_diff_ratio.day[c.uni.all.qtl.with_diff_ratio.day$raw_trait == 'wue',]
wue_ratio_diff.qtl<-unique(temp[temp$type == 'comp_diff','qtl.name'])
wue_ratio_rel_diff.qtl<-unique(temp[temp$type == 'comp_rel_diff','qtl.name'])
wue_ratio_ratio.qtl<-unique(temp[temp$type == 'comp_ratio','qtl.name'])


grid.newpage()
#pdf("Figure_X.pdf")
draw.triple.venn(area1=length(wue_ratio_diff.qtl),
                 area2=length(wue_ratio_rel_diff.qtl),
                 area3=length(wue_ratio_ratio.qtl),
                 n12=length(intersect(wue_ratio_diff.qtl, wue_ratio_rel_diff.qtl)),
                 n23=length(intersect(wue_ratio_rel_diff.qtl, wue_ratio_ratio.qtl)),
                 n13=length(intersect(wue_ratio_diff.qtl, wue_ratio_ratio.qtl)),
                 n123=length(intersect(wue_ratio_diff.qtl, intersect(wue_ratio_rel_diff.qtl, wue_ratio_ratio.qtl))),
                 c("Diff", "Relative Diff", "Ratio"),
                 fill=c("red", "blue", "yellow")
)

#dev.off()

## TE fit
temp<-c.uni.all.qtl.with_diff_ratio.day[c.uni.all.qtl.with_diff_ratio.day$raw_trait == 'fit',]
te_fit_diff.qtl<-unique(temp[temp$type == 'comp_diff','qtl.name'])
te_fit_rel_diff.qtl<-unique(temp[temp$type == 'comp_rel_diff','qtl.name'])
te_fit_ratio.qtl<-unique(temp[temp$type == 'comp_ratio','qtl.name'])


grid.newpage()
#pdf("Figure_X.pdf")
draw.triple.venn(area1=length(te_fit_diff.qtl),
                 area2=length(te_fit_rel_diff.qtl),
                 area3=length(te_fit_ratio.qtl),
                 n12=length(intersect(te_fit_diff.qtl, te_fit_rel_diff.qtl)),
                 n23=length(intersect(te_fit_rel_diff.qtl, te_fit_ratio.qtl)),
                 n13=length(intersect(te_fit_diff.qtl, te_fit_ratio.qtl)),
                 n123=length(intersect(te_fit_diff.qtl, intersect(te_fit_rel_diff.qtl, te_fit_ratio.qtl))),
                 c("Diff", "Relative Diff", "Ratio"),
                 fill=c("red", "blue", "yellow")
)

#dev.off()

## TE residual
temp<-c.uni.all.qtl.with_diff_ratio.day[c.uni.all.qtl.with_diff_ratio.day$raw_trait == 'residual',]
te_residual_diff.qtl<-unique(temp[temp$type == 'comp_diff','qtl.name'])
te_residual_rel_diff.qtl<-unique(temp[temp$type == 'comp_rel_diff','qtl.name'])
te_residual_ratio.qtl<-unique(temp[temp$type == 'comp_ratio','qtl.name'])


grid.newpage()
#pdf("Figure_X.pdf")
draw.triple.venn(area1=length(te_residual_diff.qtl),
                 area2=length(te_residual_rel_diff.qtl),
                 area3=length(te_residual_ratio.qtl),
                 n12=length(intersect(te_residual_diff.qtl, te_residual_rel_diff.qtl)),
                 n23=length(intersect(te_residual_rel_diff.qtl, te_residual_ratio.qtl)),
                 n13=length(intersect(te_residual_diff.qtl, te_residual_ratio.qtl)),
                 n123=length(intersect(te_residual_diff.qtl, intersect(te_residual_rel_diff.qtl, te_residual_ratio.qtl))),
                 c("Diff", "Relative Diff", "Ratio"),
                 fill=c("red", "blue", "yellow")
)

#dev.off()


################################################################
## The GxE trait ratio QTL look too variable/noisy
## Here we will focus on comparing the GxE QTL that correspond to 
## Only the difference and relative difference traits
################################################################

## How many SNPs for GxE traits?

## Mathematical difference
length(unique(t.all.qtl.with_diff_ratio[t.all.qtl.with_diff_ratio$type == 'comp_diff','marker']))
## 43

## Relative difference
length(unique(t.all.qtl.with_diff_ratio[t.all.qtl.with_diff_ratio$type == 'comp_rel_diff','marker']))
## 40


## How many unique QTL for comparison traits
length(unique(c.uni.all.total.qtl$marker))
## 22
unique(c.uni.all.total.qtl$marker)

## How many QTL for each type of comparison trait
length(unique(c.uni.all.total.qtl[c.uni.all.total.qtl$type == 'comp_diff','marker']))
## 20
length(unique(c.uni.all.total.qtl[c.uni.all.total.qtl$type == 'comp_rel_diff','marker']))
## 18


comp_diff_qtl<-unique(c.uni.all.total.qtl[c.uni.all.total.qtl$type == 'comp_diff','marker'])
comp_rel_diff_qtl<-unique(c.uni.all.total.qtl[c.uni.all.total.qtl$type == 'comp_rel_diff','marker'])

## Lets get some counts of the type of difference QTL detected
length(unique(comp_diff_qtl))
## comp_diff 20
length(unique(comp_rel_diff_qtl))
## comp_rel_diff 18

## Plot Venn Diagram of the QTL
setwd(wue_results.qtl.venn.dir)


grid.newpage()
#pdf("venn_gxe_qtl_types_cumulative.pdf")
draw.pairwise.venn(area1=length(comp_diff_qtl),
                   area2=length(comp_rel_diff_qtl),
                   cross.area=length(intersect(comp_diff_qtl, comp_rel_diff_qtl)),
                   c("Diff", "Relative Diff"),
                   fill=c("red", "blue")
)

## Pretty similar results, share 16 QTL in common


c.uni.sv_area.qtl<-c.uni.all.total.qtl[grepl("sv_area", c.uni.all.total.qtl$trait),]
p<-ggplot(c.uni.sv_area.qtl, aes(x=dap_i, y=lod, colour=treatment)) + geom_line() + facet_wrap(~qtl.name)

## Water
c.uni.water.qtl<-c.uni.all.total.qtl[grepl("water_lost", c.uni.all.total.qtl$trait),]


## Lets compare the uniqueness of QTL between biomass and water (Cumulative/total)
intersect(unique(c.uni.sv_area.qtl$qtl.name), unique(c.uni.water.qtl$qtl.name))
unique(c.uni.sv_area.qtl$qtl.name)[!unique(c.uni.sv_area.qtl$qtl.name) %in% unique(c.uni.water.qtl$qtl.name)]
unique(c.uni.water.qtl$qtl.name)[!unique(c.uni.water.qtl$qtl.name) %in% unique(c.uni.sv_area.qtl$qtl.name)]


################################################################
# Break out daily (rate) meaurements
################################################################

## Get rate/day out of unified set
c.uni.all.day.qtl<-c.uni.all.qtl[!grepl('_total', c.uni.all.qtl$trait),]
unique(c.uni.all.day.qtl$qtl.name)
length(unique(c.uni.all.day.qtl$qtl.name))
# 26

## How many QTL for each type of comparison trait
length(unique(c.uni.all.day.qtl[c.uni.all.day.qtl$type == 'comp_diff','marker']))
## 22
length(unique(c.uni.all.day.qtl[c.uni.all.day.qtl$type == 'comp_rel_diff','marker']))
## 21


## Lets look at the over lap of QTL between types of difference formulation

########################################################################
# Cumulative/total
########################################################################

## Biomass
temp<-c.uni.all.total.qtl[c.uni.all.total.qtl$raw_trait == 'sv_area',]
sv.area_diff.qtl<-unique(temp[temp$type == 'comp_diff','qtl.name'])
sv.area_rel_diff.qtl<-unique(temp[temp$type == 'comp_rel_diff','qtl.name'])


grid.newpage()
#pdf("venn_gxe_qtl_types_day.pdf")
draw.pairwise.venn(area1=length(comp_diff_qtl),
                   area2=length(comp_rel_diff_qtl),
                   cross.area=length(intersect(comp_diff_qtl, comp_rel_diff_qtl)),
                   c("Diff", "Relative Diff"),
                   fill=c("red", "blue")
)

#dev.off()

## Plotting similarity between GxE trait formulations
#pdf("X.pdf")


grid.newpage()

g<-draw.pairwise.venn(area1=length(sv.area_diff.qtl),
                   area2=length(sv.area_rel_diff.qtl),
                   cross.area=length(intersect(sv.area_diff.qtl, sv.area_rel_diff.qtl)),
                   c("Diff", "Relative Diff"),
                   fill=c("red", "blue")
)

require(gridExtra)
grid.arrange(gTree(children=g), top="Plant area")

#dev.off()

## Water lost
temp<-c.uni.all.total.qtl[c.uni.all.total.qtl$raw_trait == 'water_lost',]
water_lost_diff.qtl<-unique(temp[temp$type == 'comp_diff','qtl.name'])
water_lost_rel_diff.qtl<-unique(temp[temp$type == 'comp_rel_diff','qtl.name'])
water_lost_ratio.qtl<-unique(temp[temp$type == 'comp_ratio','qtl.name'])


grid.newpage()
#pdf("Figure_X.pdf")

g<-draw.pairwise.venn(area1=length(water_lost_diff.qtl),
                      area2=length(water_lost_rel_diff.qtl),
                      cross.area=length(intersect(water_lost_rel_diff.qtl, water_lost_rel_diff.qtl)),
                      c("Diff", "Relative Diff"),
                      fill=c("red", "blue")
)

require(gridExtra)
grid.arrange(gTree(children=g), top="Water Lost")

#dev.off()


## WUE ratio
temp<-c.uni.all.total.qtl[c.uni.all.total.qtl$raw_trait == 'wue',]
wue_ratio_diff.qtl<-unique(temp[temp$type == 'comp_diff','qtl.name'])
wue_ratio_rel_diff.qtl<-unique(temp[temp$type == 'comp_rel_diff','qtl.name'])
wue_ratio_ratio.qtl<-unique(temp[temp$type == 'comp_ratio','qtl.name'])


grid.newpage()
#pdf("Figure_X.pdf")

g<-draw.pairwise.venn(area1=length(wue_ratio_diff.qtl),
                      area2=length(wue_ratio_rel_diff.qtl),
                      cross.area=length(intersect(wue_ratio_diff.qtl, wue_ratio_rel_diff.qtl)),
                      c("Diff", "Relative Diff"),
                      fill=c("red", "blue")
)

require(gridExtra)
grid.arrange(gTree(children=g), top="WUE Ratio")
#dev.off()

## TE fit
temp<-c.uni.all.total.qtl[c.uni.all.total.qtl$raw_trait == 'fit',]
te_fit_diff.qtl<-unique(temp[temp$type == 'comp_diff','qtl.name'])
te_fit_rel_diff.qtl<-unique(temp[temp$type == 'comp_rel_diff','qtl.name'])
te_fit_ratio.qtl<-unique(temp[temp$type == 'comp_ratio','qtl.name'])


grid.newpage()
#pdf("Figure_X.pdf")
g<-draw.pairwise.venn(area1=length(te_fit_diff.qtl),
                      area2=length(te_fit_rel_diff.qtl),
                      cross.area=length(intersect(te_fit_diff.qtl, te_fit_rel_diff.qtl)),
                      c("Diff", "Relative Diff"),
                      fill=c("red", "blue")
)

require(gridExtra)
grid.arrange(gTree(children=g), top="TE_model_fit")

#dev.off()

## TE residual
temp<-c.uni.all.total.qtl[c.uni.all.total.qtl$raw_trait == 'residual',]
te_residual_diff.qtl<-unique(temp[temp$type == 'comp_diff','qtl.name'])
te_residual_rel_diff.qtl<-unique(temp[temp$type == 'comp_rel_diff','qtl.name'])
te_residual_ratio.qtl<-unique(temp[temp$type == 'comp_ratio','qtl.name'])


grid.newpage()
#pdf("Figure_X.pdf")

g<-draw.pairwise.venn(area1=length(te_residual_diff.qtl),
                      area2=length(te_residual_rel_diff.qtl),
                      cross.area=length(intersect(te_residual_diff.qtl, te_residual_rel_diff.qtl)),
                      c("Diff", "Relative Diff"),
                      fill=c("red", "blue")
)

require(gridExtra)
grid.arrange(gTree(children=g), top="TE_model_residual")


dev.off()


########################################################################
# Daily/rate
########################################################################
## Here we extract trait names
raw_trait<-c()
for(i in 1:nrow(c.uni.all.day.qtl)){
  fulltrait<-as.character(c.uni.all.day.qtl[i,'trait'])
  temp<-strsplit(fulltrait, '_')
  trait.name<-temp[[1]][1:(length(temp[[1]])-2)]
  trait.name<-paste(trait.name, collapse="_")
  raw_trait<-c(raw_trait, trait.name)
}

## simplify trait names
raw_trait[grepl("sv_area", raw_trait)]<-c("sv_area")
raw_trait[grepl("water_lost", raw_trait)]<-c("water_lost")
raw_trait[grepl("wue", raw_trait)]<-c("wue")
raw_trait[grepl("te_fit", raw_trait)]<-c("fit")
raw_trait[grepl("te_residual", raw_trait)]<-c("residual")


c.uni.all.day.qtl$raw_trait<-raw_trait
unique(raw_trait)

## Biomass
temp<-c.uni.all.day.qtl[c.uni.all.day.qtl$raw_trait == 'sv_area',]
sv.area_diff.qtl<-unique(temp[temp$type == 'comp_diff','qtl.name'])
sv.area_rel_diff.qtl<-unique(temp[temp$type == 'comp_rel_diff','qtl.name'])


#pdf("FIG_X.pdf")
grid.newpage()

draw.pairwise.venn(area1=length(sv.area_diff.qtl),
                   area2=length(sv.area_rel_diff.qtl),
                   cross.area=length(intersect(sv.area_diff.qtl, sv.area_rel_diff.qtl)),
                   c("Diff", "Relative Diff"),
                   fill=c("red", "blue")
)
#dev.off()

## Water lost
temp<-c.uni.all.day.qtl[c.uni.all.day.qtl$raw_trait == 'water_lost',]
temp<-temp[!grepl('_total', temp$trait),]
water_lost_diff.qtl<-unique(temp[temp$type == 'comp_diff','qtl.name'])
water_lost_rel_diff.qtl<-unique(temp[temp$type == 'comp_rel_diff','qtl.name'])

grid.newpage()
#pdf("Figure_X.pdf")

draw.pairwise.venn(area1=length(water_lost_diff.qtl),
                   area2=length(water_lost_rel_diff.qtl),
                   cross.area=length(intersect(water_lost_diff.qtl, water_lost_rel_diff.qtl)),
                   c("Diff", "Relative Diff"),
                   fill=c("red", "blue")
)
#dev.off()


## WUE ratio
temp<-c.uni.all.day.qtl[c.uni.all.day.qtl$raw_trait == 'wue',]
wue_ratio_diff.qtl<-unique(temp[temp$type == 'comp_diff','qtl.name'])
wue_ratio_rel_diff.qtl<-unique(temp[temp$type == 'comp_rel_diff','qtl.name'])
wue_ratio_ratio.qtl<-unique(temp[temp$type == 'comp_ratio','qtl.name'])


grid.newpage()
#pdf("Figure_X.pdf")

draw.pairwise.venn(area1=length(wue_ratio_diff.qtl),
                   area2=length(wue_ratio_rel_diff.qtl),
                   cross.area=length(intersect(wue_ratio_diff.qtl, wue_ratio_rel_diff.qtl)),
                   c("Diff", "Relative Diff"),
                   fill=c("red", "blue")
)
#dev.off()

## TE fit
temp<-c.uni.all.day.qtl[c.uni.all.day.qtl$raw_trait == 'fit',]
te_fit_diff.qtl<-unique(temp[temp$type == 'comp_diff','qtl.name'])
te_fit_rel_diff.qtl<-unique(temp[temp$type == 'comp_rel_diff','qtl.name'])
te_fit_ratio.qtl<-unique(temp[temp$type == 'comp_ratio','qtl.name'])


grid.newpage()
#pdf("Figure_X.pdf")

draw.pairwise.venn(area1=length(te_fit_diff.qtl),
                   area2=length(te_fit_rel_diff.qtl),
                   cross.area=length(intersect(te_fit_diff.qtl, te_fit_rel_diff.qtl)),
                   c("Diff", "Relative Diff"),
                   fill=c("red", "blue")
)
#dev.off()

## TE residual
temp<-c.uni.all.day.qtl[c.uni.all.day.qtl$raw_trait == 'residual',]
te_residual_diff.qtl<-unique(temp[temp$type == 'comp_diff','qtl.name'])
te_residual_rel_diff.qtl<-unique(temp[temp$type == 'comp_rel_diff','qtl.name'])
te_residual_ratio.qtl<-unique(temp[temp$type == 'comp_ratio','qtl.name'])


grid.newpage()
#pdf("Figure_X.pdf")

draw.pairwise.venn(area1=length(te_residual_diff.qtl),
                   area2=length(te_residual_rel_diff.qtl),
                   cross.area=length(intersect(te_residual_diff.qtl, te_residual_rel_diff.qtl)),
                   c("Diff", "Relative Diff"),
                   fill=c("red", "blue")
)
#dev.off()



########################################################################
## Lets compare QTL presence between raw and comparison
########################################################################

## Both Total/Cumulative and Daily/Rate
## Comparison
c.uni.all.qtl_names<-unique(c.uni.all.qtl$qtl.name)

## Raw
r.uni.all.qtl_names<-unique(r.uni.all.qtl$qtl.name)

grid.newpage()
#pdf("Figure_X.pdf")
draw.pairwise.venn(area1=length(c.uni.all.qtl_names),
                   area2=length(r.uni.all.qtl_names),
                   cross.area=length(intersect(c.uni.all.qtl_names, r.uni.all.qtl_names)),
                   c("Difference", "Raw"),
                   fill=c("red", "blue")
)

#dev.off()

uni.all.raw.comp_both<-intersect(c.uni.all.qtl_names, r.uni.all.qtl_names)
length(uni.all.raw.comp_both)
uni.all.raw.comp_both

uni.all.raw_only<-r.uni.all.qtl_names[!r.uni.all.qtl_names %in% c.uni.all.qtl_names]
length(uni.all.raw_only)
uni.all.raw_only


uni.all.comp_only<-c.uni.all.qtl_names[!c.uni.all.qtl_names %in% r.uni.all.qtl_names]
length(uni.all.comp_only)
uni.all.comp_only

length(uni.all.raw_only)/length(uni.all.raw.comp_both)

########################################################################
## Total/Cumulative
########################################################################

c.uni.all.total.qtl_names<-unique(c.uni.all.total.qtl$qtl.name)
r.uni.all.total.qtl_names<-unique(r.uni.all.total.qtl$qtl.name)

grid.newpage()
#pdf("Figure_X.pdf")
draw.pairwise.venn(area1=length(c.uni.all.total.qtl_names),
                   area2=length(r.uni.all.total.qtl_names),
                   cross.area=length(intersect(c.uni.all.total.qtl_names, r.uni.all.total.qtl_names)),
                   c("Difference", "Raw"),
                   fill=c("red", "blue")
)

#dev.off()


uni.all.total.raw.comp_both<-intersect(c.uni.all.total.qtl_names, r.uni.all.total.qtl_names)
length(uni.all.total.raw.comp_both)
uni.all.total.raw.comp_both

uni.all.total.raw_only<-r.uni.all.total.qtl_names[!r.uni.all.total.qtl_names %in% c.uni.all.total.qtl_names]
length(uni.all.total.raw_only)
uni.all.total.raw_only

uni.all.total.comp_only<-c.uni.all.total.qtl_names[!c.uni.all.total.qtl_names %in% r.uni.all.total.qtl_names]
length(uni.all.total.comp_only)
uni.all.total.comp_only


########################################################################
## Daily/rate
########################################################################

c.uni.all.day.qtl_names<-unique(c.uni.all.day.qtl$qtl.name)
r.uni.all.day.qtl_names<-unique(r.uni.all.day.qtl$qtl.name)

grid.newpage()
#pdf("Figure_X.pdf")
draw.pairwise.venn(area1=length(c.uni.all.day.qtl_names),
                   area2=length(r.uni.all.day.qtl_names),
                   cross.area=length(intersect(c.uni.all.day.qtl_names, r.uni.all.day.qtl_names)),
                   c("Difference", "Raw"),
                   fill=c("red", "blue")
)

#dev.off()


uni.all.day.raw.comp_both<-intersect(c.uni.all.day.qtl_names, r.uni.all.day.qtl_names)
length(uni.all.day.raw.comp_both)
uni.all.day.raw.comp_both


uni.all.day.raw_only<-r.uni.all.day.qtl_names[!r.uni.all.day.qtl_names %in% c.uni.all.day.qtl_names]
length(uni.all.day.raw_only)
uni.all.day.raw_only

uni.all.day.comp_only<-c.uni.all.day.qtl_names[!c.uni.all.day.qtl_names %in% r.uni.all.day.qtl_names]
length(uni.all.day.comp_only)
uni.all.day.comp_only


########################################################################
## Lets compare total/cumulative raw vs comparison for each trait
########################################################################

## Get QTL unique to each trait for cumulative raw traits
r.total.sv_area<-unique(r.uni.all.total.qtl[grepl("sv_area", r.uni.all.total.qtl$trait),'qtl.name'])
r.total.water_lost<-unique(r.uni.all.total.qtl[grepl("water_lost", r.uni.all.total.qtl$trait),'qtl.name'])
r.total.wue_ratio<-unique(r.uni.all.total.qtl[grepl("wue_total", r.uni.all.total.qtl$trait),'qtl.name'])
r.total.te_fit<-unique(r.uni.all.total.qtl[grepl("te_fit", r.uni.all.total.qtl$trait),'qtl.name'])
r.total.te_residual<-unique(r.uni.all.total.qtl[grepl("te_residual", r.uni.all.total.qtl$trait),'qtl.name'])

r.total.all.qtl<-unique(c(r.total.sv_area, r.total.water_lost, r.total.wue_ratio, r.total.te_fit, r.total.te_residual))

## Get QTL unique to each trait for cumulative raw traits
c.total.sv_area<-unique(c.uni.all.total.qtl[grepl("sv_area", c.uni.all.total.qtl$trait),'qtl.name'])
c.total.water_lost<-unique(c.uni.all.total.qtl[grepl("water_lost", c.uni.all.total.qtl$trait),'qtl.name'])
c.total.wue_ratio<-unique(c.uni.all.total.qtl[grepl("wue_total", c.uni.all.total.qtl$trait),'qtl.name'])
c.total.te_fit<-unique(c.uni.all.total.qtl[grepl("te_fit", c.uni.all.total.qtl$trait),'qtl.name'])
c.total.te_residual<-unique(c.uni.all.total.qtl[grepl("te_residual", c.uni.all.total.qtl$trait),'qtl.name'])

c.total.all.qtl<-unique(c(c.total.sv_area, c.total.water_lost, c.total.wue_ratio, c.total.te_fit, c.total.te_residual))



##############################
## sv_area
##############################

## FIG S20
setwd(wue_results.qtl.venn.dir)

pdf("venn_cumulative_v_gxe_qtl.pdf")

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

setwd(wue_results.qtl.venn.dir)

## FIG S21

pdf("venn_trait_overlap_of_qtl_found_in_both_cumulative_and_gxe.pdf")

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

## QTL shared by all traits but TE residual
intersect(unique(total.sv_area.raw.comp_both), intersect(unique(total.water_lost.raw.comp_both), intersect(unique(total.wue_ratio.raw.comp_both), unique(total.te_fit.raw.comp_both))))
## 7@99

## QTL shared by all traits but WUE ratio
intersect(unique(total.sv_area.raw.comp_both), intersect(unique(total.water_lost.raw.comp_both), intersect(unique(total.te_fit.raw.comp_both), unique(total.te_residual.raw.comp_both))))## 7@99
## 2@96

## QTL shared by all traits but water loss
intersect(unique(total.sv_area.raw.comp_both), intersect(unique(total.wue_ratio.raw.comp_both), intersect(unique(total.te_fit.raw.comp_both), unique(total.te_residual.raw.comp_both))))
## 0

## QTL shared by sv_area, water_lost, and TE fit
intersect(unique(total.sv_area.raw.comp_both), intersect(unique(total.water_lost.raw.comp_both), unique(total.te_fit.raw.comp_both)))
## "2@96"  "3@48"  "5@109" "7@99"  "7@51"  "7@34"  "9@34"

## QTL shared by sv_area, water_lost, TE fit and WUE ratio
intersect(unique(total.sv_area.raw.comp_both), intersect(unique(total.water_lost.raw.comp_both), intersect(unique(total.wue_ratio.raw.comp_both), unique(total.te_fit.raw.comp_both))))

## QTL shared by sv_area, WUE_ratio, TE_fit, TE residual
intersect(unique(total.sv_area.raw.comp_both), intersect(unique(total.wue_ratio.raw.comp_both), intersect(unique(total.te_fit.raw.comp_both), unique(total.te_residual.raw.comp_both))))

## QTL shared by sv_area and WUE ratio
intersect(unique(total.sv_area.comp_only), unique(total.wue_ratio.comp_only))

###########################################################################

## Difference traits only

## FIG S22

pdf("venn_trait_overlap_of_qtl_found_only_in_gxe.pdf")
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
## For Fig S19

## QTL shared by all traits
intersect(unique(total.sv_area.comp_only), intersect(unique(total.water_lost.comp_only), intersect(unique(total.wue_ratio.comp_only), intersect(unique(total.te_fit.comp_only), unique(total.te_residual.comp_only)))))
## 0

## QTL shared by all traits except TE residual
intersect(unique(total.sv_area.raw.comp_both), intersect(unique(total.water_lost.raw.comp_both), intersect(unique(total.wue_ratio.raw.comp_both), unique(total.te_fit.raw.comp_both))))
## 7@99

## QTL shared by sv_area, water_lost, and TE fit
intersect(unique(total.sv_area.comp_only), intersect(unique(total.water_lost.comp_only), unique(total.te_fit.comp_only)))

## QTL shared by water loss and TE fit
intersect(unique(total.water_lost.comp_only), unique(total.te_fit.comp_only))

## QTL shared by plant size and wue ratio
intersect(unique(total.sv_area.comp_only), unique(total.wue_ratio.comp_only))
##  5@15

intersect(unique(total.te_fit.comp_only), unique(total.te_residual.comp_only))



########################################################################
## Lets compare rate/daily raw vs comparison for each trait
########################################################################

## Get QTL unique to each trait for cumulative raw traits
r.day.sv_area<-unique(r.uni.all.day.qtl[grepl("sv_area", r.uni.all.day.qtl$trait),'qtl.name'])
r.day.water_lost<-unique(r.uni.all.day.qtl[grepl("water_lost", r.uni.all.day.qtl$trait),'qtl.name'])
r.day.wue_ratio<-unique(r.uni.all.day.qtl[grepl("wue_day", r.uni.all.day.qtl$trait),'qtl.name'])
r.day.te_fit<-unique(r.uni.all.day.qtl[grepl("te_fit", r.uni.all.day.qtl$trait),'qtl.name'])
r.day.te_residual<-unique(r.uni.all.day.qtl[grepl("te_residual", r.uni.all.day.qtl$trait),'qtl.name'])

r.day.all.qtl<-unique(c(r.day.sv_area, r.day.water_lost, r.day.wue_ratio, r.day.te_fit, r.day.te_residual))

## Get QTL unique to each trait for cumulative raw traits
c.day.sv_area<-unique(c.uni.all.day.qtl[grepl("sv_area", c.uni.all.day.qtl$trait),'qtl.name'])
c.day.water_lost<-unique(c.uni.all.day.qtl[grepl("water_lost", c.uni.all.day.qtl$trait),'qtl.name'])
c.day.wue_ratio<-unique(c.uni.all.day.qtl[grepl("wue_day", c.uni.all.day.qtl$trait),'qtl.name'])
c.day.te_fit<-unique(c.uni.all.day.qtl[grepl("te_fit", c.uni.all.day.qtl$trait),'qtl.name'])
c.day.te_residual<-unique(c.uni.all.day.qtl[grepl("te_residual", c.uni.all.day.qtl$trait),'qtl.name'])

c.day.all.qtl<-unique(c(c.day.sv_area, c.day.water_lost, c.day.wue_ratio, c.day.te_fit, c.day.te_residual))



##############################
## sv_area
##############################

grid.newpage()
#pdf("Figure_X.pdf")
draw.pairwise.venn(area1=length(c.day.sv_area),
                   area2=length(r.day.sv_area),
                   cross.area=length(intersect(c.day.sv_area, r.day.sv_area)),
                   c("Difference", "Raw"),
                   fill=c("red", "blue")
)

#dev.off()


day.sv_area.raw.comp_both<-intersect(r.day.sv_area, c.day.sv_area)
length(day.sv_area.raw.comp_both)
day.sv_area.raw.comp_both

day.sv_area.raw_only<-r.day.sv_area[!r.day.sv_area %in% c.day.sv_area]
length(day.sv_area.raw_only)
day.sv_area.raw_only


day.sv_area.comp_only<-c.day.sv_area[!c.day.sv_area %in% r.day.sv_area]
length(day.sv_area.comp_only)
day.sv_area.comp_only


##############################
## water_lost
##############################

grid.newpage()
#pdf("Figure_X.pdf")
draw.pairwise.venn(area1=length(c.day.water_lost),
                   area2=length(r.day.water_lost),
                   cross.area=length(intersect(c.day.water_lost, r.day.water_lost)),
                   c("Difference", "Raw"),
                   fill=c("red", "blue")
)

#dev.off()


day.water_lost.raw.comp_both<-intersect(r.day.water_lost, c.day.water_lost)
length(day.water_lost.raw.comp_both)
day.water_lost.raw.comp_both

day.water_lost.raw_only<-r.day.water_lost[!r.day.water_lost %in% c.day.water_lost]
length(day.water_lost.raw_only)
day.water_lost.raw_only


day.water_lost.comp_only<-c.day.water_lost[!c.day.water_lost %in% r.day.water_lost]
length(day.water_lost.comp_only)
day.water_lost.comp_only


##############################
## wue_ratio
##############################

grid.newpage()
#pdf("Figure_X.pdf")
draw.pairwise.venn(area1=length(c.day.wue_ratio),
                   area2=length(r.day.wue_ratio),
                   cross.area=length(intersect(c.day.wue_ratio, r.day.wue_ratio)),
                   c("Difference", "Raw"),
                   fill=c("red", "blue")
)

#dev.off()


day.wue_ratio.raw.comp_both<-intersect(r.day.wue_ratio, c.day.wue_ratio)
length(day.wue_ratio.raw.comp_both)
day.wue_ratio.raw.comp_both

day.wue_ratio.raw_only<-r.day.wue_ratio[!r.day.wue_ratio %in% c.day.wue_ratio]
length(day.wue_ratio.raw_only)
day.wue_ratio.raw_only


day.wue_ratio.comp_only<-c.day.wue_ratio[!c.day.wue_ratio %in% r.day.wue_ratio]
length(day.wue_ratio.comp_only)
day.wue_ratio.comp_only


##############################
## te_fit
##############################

grid.newpage()
#pdf("Figure_X.pdf")
draw.pairwise.venn(area1=length(c.day.te_fit),
                   area2=length(r.day.te_fit),
                   cross.area=length(intersect(c.day.te_fit, r.day.te_fit)),
                   c("Difference", "Raw"),
                   fill=c("red", "blue")
)

#dev.off()


day.te_fit.raw.comp_both<-intersect(r.day.te_fit, c.day.te_fit)
length(day.te_fit.raw.comp_both)
day.te_fit.raw.comp_both

day.te_fit.raw_only<-r.day.te_fit[!r.day.te_fit %in% c.day.te_fit]
length(day.te_fit.raw_only)
day.te_fit.raw_only


day.te_fit.comp_only<-c.day.te_fit[!c.day.te_fit %in% r.day.te_fit]
length(day.te_fit.comp_only)
day.te_fit.comp_only


##############################
## te_residual
##############################

grid.newpage()
#pdf("Figure_X.pdf")
draw.pairwise.venn(area1=length(c.day.te_residual),
                   area2=length(r.day.te_residual),
                   cross.area=length(intersect(c.day.te_residual, r.day.te_residual)),
                   c("Difference", "Raw"),
                   fill=c("red", "blue")
)

#dev.off()


day.te_residual.raw.comp_both<-intersect(r.day.te_residual, c.day.te_residual)
length(day.te_residual.raw.comp_both)
day.te_residual.raw.comp_both

day.te_residual.raw_only<-r.day.te_residual[!r.day.te_residual %in% c.day.te_residual]
length(day.te_residual.raw_only)
day.te_residual.raw_only


day.te_residual.comp_only<-c.day.te_residual[!c.day.te_residual %in% r.day.te_residual]
length(day.te_residual.comp_only)
day.te_residual.comp_only



































###########################################################################

### QTL shared between diff and daily/rate trait

## FIG S18?

#pdf("X.pdf")

grid.newpage()
draw.quintuple.venn(area1=length(day.sv_area.raw.comp_both), 
                    area2=length(day.water_lost.raw.comp_both), 
                    area3=length(day.wue_ratio.raw.comp_both),
                    area4=length(day.te_fit.raw.comp_both), 
                    area5=length(day.te_residual.raw.comp_both),
                    n12=length(intersect(unique(day.sv_area.raw.comp_both), unique(day.water_lost.raw.comp_both))),
                    n13=length(intersect(unique(day.sv_area.raw.comp_both), unique(day.wue_ratio.raw.comp_both))),
                    n14=length(intersect(unique(day.sv_area.raw.comp_both), unique(day.te_fit.raw.comp_both))),
                    n15=length(intersect(unique(day.sv_area.raw.comp_both), unique(day.te_residual.raw.comp_both))),
                    n23=length(intersect(unique(day.water_lost.raw.comp_both), unique(day.wue_ratio.raw.comp_both))),
                    n24=length(intersect(unique(day.water_lost.raw.comp_both), unique(day.te_fit.raw.comp_both))),
                    n25=length(intersect(unique(day.water_lost.raw.comp_both), unique(day.te_residual.raw.comp_both))),
                    n34=length(intersect(unique(day.wue_ratio.raw.comp_both), unique(day.te_fit.raw.comp_both))),
                    n35=length(intersect(unique(day.wue_ratio.raw.comp_both), unique(day.te_residual.raw.comp_both))),
                    n45=length(intersect(unique(day.te_fit.raw.comp_both), unique(day.te_residual.raw.comp_both))),
                    n123=length(intersect(unique(day.sv_area.raw.comp_both), intersect(unique(day.water_lost.raw.comp_both), unique(day.wue_ratio.raw.comp_both)))),
                    n124=length(intersect(unique(day.sv_area.raw.comp_both), intersect(unique(day.water_lost.raw.comp_both), unique(day.te_fit.raw.comp_both)))),
                    n125=length(intersect(unique(day.sv_area.raw.comp_both), intersect(unique(day.water_lost.raw.comp_both), unique(day.te_residual.raw.comp_both)))),
                    n134=length(intersect(unique(day.sv_area.raw.comp_both), intersect(unique(day.wue_ratio.raw.comp_both), unique(day.te_fit.raw.comp_both)))),
                    n135=length(intersect(unique(day.sv_area.raw.comp_both), intersect(unique(day.wue_ratio.raw.comp_both), unique(day.te_residual.raw.comp_both)))),
                    n145=length(intersect(unique(day.sv_area.raw.comp_both), intersect(unique(day.te_fit.raw.comp_both), unique(day.te_residual.raw.comp_both)))),
                    n234=length(intersect(unique(day.water_lost.raw.comp_both), intersect(unique(day.wue_ratio.raw.comp_both), unique(day.te_fit.raw.comp_both)))),
                    n235=length(intersect(unique(day.water_lost.raw.comp_both), intersect(unique(day.wue_ratio.raw.comp_both), unique(day.te_residual.raw.comp_both)))),
                    n245=length(intersect(unique(day.water_lost.raw.comp_both), intersect(unique(day.te_fit.raw.comp_both), unique(day.te_residual.raw.comp_both)))),
                    n345=length(intersect(unique(day.wue_ratio.raw.comp_both), intersect(unique(day.te_fit.raw.comp_both), unique(day.te_residual.raw.comp_both)))),
                    n1234=length(intersect(unique(day.sv_area.raw.comp_both), intersect(unique(day.water_lost.raw.comp_both), intersect(unique(day.wue_ratio.raw.comp_both), unique(day.te_fit.raw.comp_both))))),
                    n1235=length(intersect(unique(day.sv_area.raw.comp_both), intersect(unique(day.water_lost.raw.comp_both), intersect(unique(day.wue_ratio.raw.comp_both), unique(day.te_residual.raw.comp_both))))),
                    n1245=length(intersect(unique(day.sv_area.raw.comp_both), intersect(unique(day.water_lost.raw.comp_both), intersect(unique(day.te_fit.raw.comp_both), unique(day.te_residual.raw.comp_both))))),
                    n1345=length(intersect(unique(day.sv_area.raw.comp_both), intersect(unique(day.wue_ratio.raw.comp_both), intersect(unique(day.te_fit.raw.comp_both), unique(day.te_residual.raw.comp_both))))),
                    n2345=length(intersect(unique(day.water_lost.raw.comp_both), intersect(unique(day.wue_ratio.raw.comp_both), intersect(unique(day.te_fit.raw.comp_both), unique(day.te_residual.raw.comp_both))))),
                    n12345=length(intersect(unique(day.sv_area.raw.comp_both), intersect(unique(day.water_lost.raw.comp_both), intersect(unique(day.wue_ratio.raw.comp_both), intersect(unique(day.te_fit.raw.comp_both), unique(day.te_residual.raw.comp_both)))))),
                    category=c("Biomass", "Water lost", "WUE ratio", "TE fit", "TE residual"),
                    fill=c("green","light blue", "orange", "black", "red")
)


#dev.off()

## QTL shared by all traits but TE residual
intersect(unique(day.sv_area.raw.comp_both), intersect(unique(day.water_lost.raw.comp_both), intersect(unique(day.wue_ratio.raw.comp_both), unique(day.te_fit.raw.comp_both))))
## 2@96 7@99

## QTL shared by all traits but WUE ratio
intersect(unique(day.sv_area.raw.comp_both), intersect(unique(day.water_lost.raw.comp_both), intersect(unique(day.te_fit.raw.comp_both), unique(day.te_residual.raw.comp_both))))
## 5@79 9@34

## QTL shared by all traits but water loss
intersect(unique(day.sv_area.raw.comp_both), intersect(unique(day.wue_ratio.raw.comp_both), intersect(unique(day.te_fit.raw.comp_both), unique(day.te_residual.raw.comp_both))))
## 0

## QTL shared by sv_area, water_lost, and TE fit
intersect(unique(day.sv_area.raw.comp_both), intersect(unique(day.water_lost.raw.comp_both), unique(day.te_fit.raw.comp_both)))
## "2@96"  "5@79"  "5@109" "7@99"  "7@51"  "7@34"  "9@34"

## QTL shared by sv_area, water_lost, TE fit and WUE ratio
intersect(unique(day.sv_area.raw.comp_both), intersect(unique(day.water_lost.raw.comp_both), intersect(unique(day.wue_ratio.raw.comp_both), unique(day.te_fit.raw.comp_both))))
## 2@96 7@99

## QTL shared by sv_area, WUE_ratio, TE_fit, TE residual
intersect(unique(day.sv_area.raw.comp_both), intersect(unique(day.wue_ratio.raw.comp_both), intersect(unique(day.te_fit.raw.comp_both), unique(day.te_residual.raw.comp_both))))

## QTL shared by sv_area and WUE ratio
intersect(unique(day.sv_area.comp_only), unique(day.wue_ratio.comp_only))

###########################################################################

## Daily/Rate GxE traits only


#pdf("X.pdf")
grid.newpage()
draw.quintuple.venn(area1=length(day.sv_area.comp_only), 
                    area2=length(day.water_lost.comp_only), 
                    area3=length(day.wue_ratio.comp_only),
                    area4=length(day.te_fit.comp_only), 
                    area5=length(day.te_residual.comp_only),
                    n12=length(intersect(unique(day.sv_area.comp_only), unique(day.water_lost.comp_only))),
                    n13=length(intersect(unique(day.sv_area.comp_only), unique(day.wue_ratio.comp_only))),
                    n14=length(intersect(unique(day.sv_area.comp_only), unique(day.te_fit.comp_only))),
                    n15=length(intersect(unique(day.sv_area.comp_only), unique(day.te_residual.comp_only))),
                    n23=length(intersect(unique(day.water_lost.comp_only), unique(day.wue_ratio.comp_only))),
                    n24=length(intersect(unique(day.water_lost.comp_only), unique(day.te_fit.comp_only))),
                    n25=length(intersect(unique(day.water_lost.comp_only), unique(day.te_residual.comp_only))),
                    n34=length(intersect(unique(day.wue_ratio.comp_only), unique(day.te_fit.comp_only))),
                    n35=length(intersect(unique(day.wue_ratio.comp_only), unique(day.te_residual.comp_only))),
                    n45=length(intersect(unique(day.te_fit.comp_only), unique(day.te_residual.comp_only))),
                    n123=length(intersect(unique(day.sv_area.comp_only), intersect(unique(day.water_lost.comp_only), unique(day.wue_ratio.comp_only)))),
                    n124=length(intersect(unique(day.sv_area.comp_only), intersect(unique(day.water_lost.comp_only), unique(day.te_fit.comp_only)))),
                    n125=length(intersect(unique(day.sv_area.comp_only), intersect(unique(day.water_lost.comp_only), unique(day.te_residual.comp_only)))),
                    n134=length(intersect(unique(day.sv_area.comp_only), intersect(unique(day.wue_ratio.comp_only), unique(day.te_fit.comp_only)))),
                    n135=length(intersect(unique(day.sv_area.comp_only), intersect(unique(day.wue_ratio.comp_only), unique(day.te_residual.comp_only)))),
                    n145=length(intersect(unique(day.sv_area.comp_only), intersect(unique(day.te_fit.comp_only), unique(day.te_residual.comp_only)))),
                    n234=length(intersect(unique(day.water_lost.comp_only), intersect(unique(day.wue_ratio.comp_only), unique(day.te_fit.comp_only)))),
                    n235=length(intersect(unique(day.water_lost.comp_only), intersect(unique(day.wue_ratio.comp_only), unique(day.te_residual.comp_only)))),
                    n245=length(intersect(unique(day.water_lost.comp_only), intersect(unique(day.te_fit.comp_only), unique(day.te_residual.comp_only)))),
                    n345=length(intersect(unique(day.wue_ratio.comp_only), intersect(unique(day.te_fit.comp_only), unique(day.te_residual.comp_only)))),
                    n1234=length(intersect(unique(day.sv_area.comp_only), intersect(unique(day.water_lost.comp_only), intersect(unique(day.wue_ratio.comp_only), unique(day.te_fit.comp_only))))),
                    n1235=length(intersect(unique(day.sv_area.comp_only), intersect(unique(day.water_lost.comp_only), intersect(unique(day.wue_ratio.comp_only), unique(day.te_residual.comp_only))))),
                    n1245=length(intersect(unique(day.sv_area.comp_only), intersect(unique(day.water_lost.comp_only), intersect(unique(day.te_fit.comp_only), unique(day.te_residual.comp_only))))),
                    n1345=length(intersect(unique(day.sv_area.comp_only), intersect(unique(day.wue_ratio.comp_only), intersect(unique(day.te_fit.comp_only), unique(day.te_residual.comp_only))))),
                    n2345=length(intersect(unique(day.water_lost.comp_only), intersect(unique(day.wue_ratio.comp_only), intersect(unique(day.te_fit.comp_only), unique(day.te_residual.comp_only))))),
                    n12345=length(intersect(unique(day.sv_area.comp_only), intersect(unique(day.water_lost.comp_only), intersect(unique(day.wue_ratio.comp_only), intersect(unique(day.te_fit.comp_only), unique(day.te_residual.comp_only)))))),
                    category=c("Biomass", "Water lost", "WUE ratio", "TE fit", "TE residual"),
                    fill=c("green","light blue", "orange", "black", "red")
)


dev.off()

## QTL shared by sv area and TE fit
intersect(unique(day.sv_area.comp_only), unique(day.te_fit.comp_only))
## 7@13

## QTL shared by water loss and TE fit
intersect(unique(day.water_lost.comp_only), unique(day.te_fit.comp_only))
## 3@48

## QTL shared by plant size and wue ratio
intersect(unique(day.sv_area.comp_only), unique(day.wue_ratio.comp_only))
##  0


## Get only the daily gxe QTL


## Make a unique data frame for difference and days
d.uni.all.qtl<-c.uni.all.qtl[grepl('_day', c.uni.all.qtl$trait),]

## Lets get out water_lost (day), field doesn't have day in name
temp<-c.uni.all.qtl[grepl('water_lost_', c.uni.all.qtl$trait),]
temp<-temp[!(grepl('_total', temp$trait)),]
d.uni.all.qtl<-rbind(d.uni.all.qtl, temp)



## Get QTL unique to each trait
d.sv_area<-unique(d.uni.all.qtl[grepl("sv_area", d.uni.all.qtl$trait),'qtl.name'])
d.water_lost<-unique(d.uni.all.qtl[grepl("water_lost", d.uni.all.qtl$trait),'qtl.name'])
d.wue_ratio<-unique(d.uni.all.qtl[grepl("wue_day", d.uni.all.qtl$trait),'qtl.name'])
d.te_fit<-unique(d.uni.all.qtl[grepl("te_fit", d.uni.all.qtl$trait),'qtl.name'])
d.te_residual<-unique(d.uni.all.qtl[grepl("te_residual", d.uni.all.qtl$trait),'qtl.name'])

d.all_qtl.nr<-unique(c(d.sv_area, d.water_lost, d.wue_ratio, d.te_fit, d.te_residual))


##############################
# sv_area
##############################

grid.newpage()
#pdf("Figure_X.pdf")
draw.pairwise.venn(area1=length(d.sv_area),
                   area2=length(day.sv_area.raw.comp_both),
                   cross.area=length(intersect(d.sv_area, day.sv_area.raw.comp_both)),
                   c("Difference", "Raw"),
                   fill=c("red", "blue")
)

#dev.off()


d.sv_area.raw.comp_both<-intersect(day.sv_area.raw.comp_both, d.sv_area)
length(d.sv_area.raw.comp_both)
d.sv_area.raw.comp_both

d.sv_area.raw_only<-day.sv_area.raw.comp_both[!day.sv_area.raw.comp_both %in% d.sv_area]
length(d.sv_area.raw_only)
d.sv_area.raw_only

d.sv_area.comp_only<-d.sv_area[!d.sv_area %in% day.sv_area.raw.comp_both]
length(d.sv_area.comp_only)
d.sv_area.comp_only



##############################
# water_lost
##############################

grid.newpage()
#pdf("Figure_X.pdf")
draw.pairwise.venn(area1=length(d.water_lost),
                   area2=length(day.water_lost.raw.comp_both),
                   cross.area=length(intersect(d.water_lost, day.water_lost.raw.comp_both)),
                   c("Difference", "Raw"),
                   fill=c("red", "blue")
)

#dev.off()


d.water_lost.raw.comp_both<-intersect(day.water_lost.raw.comp_both, d.water_lost)
length(d.water_lost.raw.comp_both)
d.water_lost.raw.comp_both

d.water_lost.raw_only<-day.water_lost.raw.comp_both[!day.water_lost.raw.comp_both %in% d.water_lost]
length(d.water_lost.raw_only)
d.water_lost.raw_only

d.water_lost.comp_only<-d.water_lost[!d.water_lost %in% day.water_lost.raw.comp_both]
length(d.water_lost.comp_only)
d.water_lost.comp_only

##############################
# wue_ratio
##############################

grid.newpage()
#pdf("Figure_X.pdf")
draw.pairwise.venn(area1=length(d.wue_ratio),
                   area2=length(day.wue_ratio.raw.comp_both),
                   cross.area=length(intersect(d.wue_ratio, day.wue_ratio.raw.comp_both)),
                   c("Difference", "Raw"),
                   fill=c("red", "blue")
)

#dev.off()


d.wue_ratio.raw.comp_both<-intersect(day.wue_ratio.raw.comp_both, d.wue_ratio)
length(d.wue_ratio.raw.comp_both)
d.wue_ratio.raw.comp_both

d.wue_ratio.raw_only<-day.wue_ratio.raw.comp_both[!day.wue_ratio.raw.comp_both %in% d.wue_ratio]
length(d.wue_ratio.raw_only)
d.wue_ratio.raw_only

d.wue_ratio.comp_only<-d.wue_ratio[!d.wue_ratio %in% day.wue_ratio.raw.comp_both]
length(d.wue_ratio.comp_only)
d.wue_ratio.comp_only


##############################
# te_fit
##############################

grid.newpage()
#pdf("Figure_X.pdf")
draw.pairwise.venn(area1=length(d.te_fit),
                   area2=length(day.te_fit.raw.comp_both),
                   cross.area=length(intersect(d.te_fit, day.te_fit.raw.comp_both)),
                   c("Difference", "Raw"),
                   fill=c("red", "blue")
)

#dev.off()


d.te_fit.raw.comp_both<-intersect(day.te_fit.raw.comp_both, d.te_fit)
length(d.te_fit.raw.comp_both)
d.te_fit.raw.comp_both

d.te_fit.raw_only<-day.te_fit.raw.comp_both[!day.te_fit.raw.comp_both %in% d.te_fit]
length(d.te_fit.raw_only)
d.te_fit.raw_only

d.te_fit.comp_only<-d.te_fit[!d.te_fit %in% day.te_fit.raw.comp_both]
length(d.te_fit.comp_only)
d.te_fit.comp_only


##############################
# te_residual
##############################

grid.newpage()
#pdf("Figure_X.pdf")
draw.pairwise.venn(area1=length(d.te_residual),
                   area2=length(day.te_residual.raw.comp_both),
                   cross.area=length(intersect(d.te_residual, day.te_residual.raw.comp_both)),
                   c("Difference", "Raw"),
                   fill=c("red", "blue")
)

#dev.off()


d.te_residual.raw.comp_both<-intersect(day.te_residual.raw.comp_both, d.te_residual)
length(d.te_residual.raw.comp_both)
d.te_residual.raw.comp_both

d.te_residual.raw_only<-day.te_residual.raw.comp_both[!day.te_residual.raw.comp_both %in% d.te_residual]
length(d.te_residual.raw_only)
d.te_residual.raw_only

d.te_residual.comp_only<-d.te_residual[!d.te_residual %in% day.te_residual.raw.comp_both]
length(d.te_residual.comp_only)
d.te_residual.comp_only





### QTL shared between diff and cumulative 
grid.newpage()
draw.quintuple.venn(area1=length(day.sv_area.raw.comp_both), 
                    area2=length(day.water_lost.raw.comp_both), 
                    area3=length(day.wue_ratio.raw.comp_both),
                    area4=length(day.te_fit.raw.comp_both), 
                    area5=length(day.te_residual.raw.comp_both),
                    n12=length(intersect(unique(day.sv_area.raw.comp_both), unique(day.water_lost.raw.comp_both))),
                    n13=length(intersect(unique(day.sv_area.raw.comp_both), unique(day.wue_ratio.raw.comp_both))),
                    n14=length(intersect(unique(day.sv_area.raw.comp_both), unique(day.te_fit.raw.comp_both))),
                    n15=length(intersect(unique(day.sv_area.raw.comp_both), unique(day.te_residual.raw.comp_both))),
                    n23=length(intersect(unique(day.water_lost.raw.comp_both), unique(day.wue_ratio.raw.comp_both))),
                    n24=length(intersect(unique(day.water_lost.raw.comp_both), unique(day.te_fit.raw.comp_both))),
                    n25=length(intersect(unique(day.water_lost.raw.comp_both), unique(day.te_residual.raw.comp_both))),
                    n34=length(intersect(unique(day.wue_ratio.raw.comp_both), unique(day.te_fit.raw.comp_both))),
                    n35=length(intersect(unique(day.wue_ratio.raw.comp_both), unique(day.te_residual.raw.comp_both))),
                    n45=length(intersect(unique(day.te_fit.raw.comp_both), unique(day.te_residual.raw.comp_both))),
                    n123=length(intersect(unique(day.sv_area.raw.comp_both), intersect(unique(day.water_lost.raw.comp_both), unique(day.wue_ratio.raw.comp_both)))),
                    n124=length(intersect(unique(day.sv_area.raw.comp_both), intersect(unique(day.water_lost.raw.comp_both), unique(day.te_fit.raw.comp_both)))),
                    n125=length(intersect(unique(day.sv_area.raw.comp_both), intersect(unique(day.water_lost.raw.comp_both), unique(day.te_residual.raw.comp_both)))),
                    n134=length(intersect(unique(day.sv_area.raw.comp_both), intersect(unique(day.wue_ratio.raw.comp_both), unique(day.te_fit.raw.comp_both)))),
                    n135=length(intersect(unique(day.sv_area.raw.comp_both), intersect(unique(day.wue_ratio.raw.comp_both), unique(day.te_residual.raw.comp_both)))),
                    n145=length(intersect(unique(day.sv_area.raw.comp_both), intersect(unique(day.te_fit.raw.comp_both), unique(day.te_residual.raw.comp_both)))),
                    n234=length(intersect(unique(day.water_lost.raw.comp_both), intersect(unique(day.wue_ratio.raw.comp_both), unique(day.te_fit.raw.comp_both)))),
                    n235=length(intersect(unique(day.water_lost.raw.comp_both), intersect(unique(day.wue_ratio.raw.comp_both), unique(day.te_residual.raw.comp_both)))),
                    n245=length(intersect(unique(day.water_lost.raw.comp_both), intersect(unique(day.te_fit.raw.comp_both), unique(day.te_residual.raw.comp_both)))),
                    n345=length(intersect(unique(day.wue_ratio.raw.comp_both), intersect(unique(day.te_fit.raw.comp_both), unique(day.te_residual.raw.comp_both)))),
                    n1234=length(intersect(unique(day.sv_area.raw.comp_both), intersect(unique(day.water_lost.raw.comp_both), intersect(unique(day.wue_ratio.raw.comp_both), unique(day.te_fit.raw.comp_both))))),
                    n1235=length(intersect(unique(day.sv_area.raw.comp_both), intersect(unique(day.water_lost.raw.comp_both), intersect(unique(day.wue_ratio.raw.comp_both), unique(day.te_residual.raw.comp_both))))),
                    n1245=length(intersect(unique(day.sv_area.raw.comp_both), intersect(unique(day.water_lost.raw.comp_both), intersect(unique(day.te_fit.raw.comp_both), unique(day.te_residual.raw.comp_both))))),
                    n1345=length(intersect(unique(day.sv_area.raw.comp_both), intersect(unique(day.wue_ratio.raw.comp_both), intersect(unique(day.te_fit.raw.comp_both), unique(day.te_residual.raw.comp_both))))),
                    n2345=length(intersect(unique(day.water_lost.raw.comp_both), intersect(unique(day.wue_ratio.raw.comp_both), intersect(unique(day.te_fit.raw.comp_both), unique(day.te_residual.raw.comp_both))))),
                    n12345=length(intersect(unique(day.sv_area.raw.comp_both), intersect(unique(day.water_lost.raw.comp_both), intersect(unique(day.wue_ratio.raw.comp_both), intersect(unique(day.te_fit.raw.comp_both), unique(day.te_residual.raw.comp_both)))))),
                    category=c("Biomass", "Water lost", "WUE ratio", "TE fit", "TE residual"),
                    fill=c("green","light blue", "orange", "red", "black")
)


## QTL shared by all traits
intersect(unique(d.sv_area.raw.comp_both), intersect(unique(d.water_lost.raw.comp_both), intersect(unique(d.wue_ratio.raw.comp_both), intersect(unique(d.te_fit.raw.comp_both), unique(d.te_residual.raw.comp_both)))))
## 5@95


## QTL shared by all traits except TE residual
intersect(unique(d.sv_area.raw.comp_both), intersect(unique(d.water_lost.raw.comp_both), intersect(unique(d.wue_ratio.raw.comp_both), unique(d.te_fit.raw.comp_both))))
## 2@86, 2@96, 5@95, 7@99

## QTL shared by all traits but WUE ratio
intersect(unique(d.sv_area.raw.comp_both), intersect(unique(d.water_lost.raw.comp_both), intersect(unique(d.te_fit.raw.comp_both), unique(d.te_residual.raw.comp_both))))
## 5@95, 5@79, 9@36

## QTL shared by biomass, water lost and te fit
intersect(unique(d.sv_area.raw.comp_both), intersect(unique(d.water_lost.raw.comp_both), unique(d.te_fit.raw.comp_both)))
## 2@86, 2@96, 5@95, 5@79, 7@99, 7@53, 7@34, 9@36

## QTL shared by water lost and TE fit
intersect(unique(d.water_lost.raw.comp_both), unique(d.te_fit.raw.comp_both))
## 2@86,  2@96,  5@95,  5@106, 5@79,  5@39,  6@73,  7@99,  7@53,  7@34 

## QTL shared by TE residual and WUE ratio
intersect(unique(d.wue_ratio.raw.comp_both), unique(d.te_residual.raw.comp_both))
## 5@95,  9@127



## Difference traits only
grid.newpage()
draw.quintuple.venn(area1=length(day.sv_area.comp_only), 
                    area2=length(day.water_lost.comp_only), 
                    area3=length(day.wue_ratio.comp_only),
                    area4=length(day.te_fit.comp_only), 
                    area5=length(day.te_residual.comp_only),
                    n12=length(intersect(unique(day.sv_area.comp_only), unique(day.water_lost.comp_only))),
                    n13=length(intersect(unique(day.sv_area.comp_only), unique(day.wue_ratio.comp_only))),
                    n14=length(intersect(unique(day.sv_area.comp_only), unique(day.te_fit.comp_only))),
                    n15=length(intersect(unique(day.sv_area.comp_only), unique(day.te_residual.comp_only))),
                    n23=length(intersect(unique(day.water_lost.comp_only), unique(day.wue_ratio.comp_only))),
                    n24=length(intersect(unique(day.water_lost.comp_only), unique(day.te_fit.comp_only))),
                    n25=length(intersect(unique(day.water_lost.comp_only), unique(day.te_residual.comp_only))),
                    n34=length(intersect(unique(day.wue_ratio.comp_only), unique(day.te_fit.comp_only))),
                    n35=length(intersect(unique(day.wue_ratio.comp_only), unique(day.te_residual.comp_only))),
                    n45=length(intersect(unique(day.te_fit.comp_only), unique(day.te_residual.comp_only))),
                    n123=length(intersect(unique(day.sv_area.comp_only), intersect(unique(day.water_lost.comp_only), unique(day.wue_ratio.comp_only)))),
                    n124=length(intersect(unique(day.sv_area.comp_only), intersect(unique(day.water_lost.comp_only), unique(day.te_fit.comp_only)))),
                    n125=length(intersect(unique(day.sv_area.comp_only), intersect(unique(day.water_lost.comp_only), unique(day.te_residual.comp_only)))),
                    n134=length(intersect(unique(day.sv_area.comp_only), intersect(unique(day.wue_ratio.comp_only), unique(day.te_fit.comp_only)))),
                    n135=length(intersect(unique(day.sv_area.comp_only), intersect(unique(day.wue_ratio.comp_only), unique(day.te_residual.comp_only)))),
                    n145=length(intersect(unique(day.sv_area.comp_only), intersect(unique(day.te_fit.comp_only), unique(day.te_residual.comp_only)))),
                    n234=length(intersect(unique(day.water_lost.comp_only), intersect(unique(day.wue_ratio.comp_only), unique(day.te_fit.comp_only)))),
                    n235=length(intersect(unique(day.water_lost.comp_only), intersect(unique(day.wue_ratio.comp_only), unique(day.te_residual.comp_only)))),
                    n245=length(intersect(unique(day.water_lost.comp_only), intersect(unique(day.te_fit.comp_only), unique(day.te_residual.comp_only)))),
                    n345=length(intersect(unique(day.wue_ratio.comp_only), intersect(unique(day.te_fit.comp_only), unique(day.te_residual.comp_only)))),
                    n1234=length(intersect(unique(day.sv_area.comp_only), intersect(unique(day.water_lost.comp_only), intersect(unique(day.wue_ratio.comp_only), unique(day.te_fit.comp_only))))),
                    n1235=length(intersect(unique(day.sv_area.comp_only), intersect(unique(day.water_lost.comp_only), intersect(unique(day.wue_ratio.comp_only), unique(day.te_residual.comp_only))))),
                    n1245=length(intersect(unique(day.sv_area.comp_only), intersect(unique(day.water_lost.comp_only), intersect(unique(day.te_fit.comp_only), unique(day.te_residual.comp_only))))),
                    n1345=length(intersect(unique(day.sv_area.comp_only), intersect(unique(day.wue_ratio.comp_only), intersect(unique(day.te_fit.comp_only), unique(day.te_residual.comp_only))))),
                    n2345=length(intersect(unique(day.water_lost.comp_only), intersect(unique(day.wue_ratio.comp_only), intersect(unique(day.te_fit.comp_only), unique(day.te_residual.comp_only))))),
                    n12345=length(intersect(unique(day.sv_area.comp_only), intersect(unique(day.water_lost.comp_only), intersect(unique(day.wue_ratio.comp_only), intersect(unique(day.te_fit.comp_only), unique(day.te_residual.comp_only)))))),
                    category=c("Biomass", "Water lost", "WUE ratio", "TE fit", "TE residual"),
                    fill=c("green","light blue", "orange", "black", "red")
)


setwd(wue_results.qtl.venn.dir)
save.image("ril_venn_gxe_qtl_analysis.Rdata")

rm(list=ls())

