################### 
library(ggplot2)
library(stringr)
library(VennDiagram)
library(lme4)
library(lattice)

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

setwd(wue_results.qtl.dir)


load("qtl_summary_4_qtl_trait_matrix.Rdata")

## Lets get some summary statistics of how many SNPs observed for all traits

## Number of SNP for cumulative traits with no comparision/GxE traits
length(unique(r.all_total.qtl$marker))
# 86 unique SNPs

## Number of SNP for cumulative traits with no comparision/GxE traits
length(unique(r.all_day.qtl$marker))
# 106 unique SNPs


## Extract day after planting from trait field
dap_i<-c()
for(i in 1:nrow(uni.all.qtl)){
  x<-strsplit(as.character(uni.all.qtl[i,'trait']), '_')[[1]]
  d<-x[length(x)-1]
  dap_i<-c(dap_i, d)
}

## Convert dap_i vector to numeric
dap_i<-as.numeric(dap_i)
## Stable isotopes were only measured at a single time point so day value is 'NA'
## Lets replace that value with 30 DAP (should be no isotope QTL left in the set by this point)
dap_i[is.na(dap_i)]<-c(30)

## Add DAP
uni.all.qtl<-cbind(uni.all.qtl, dap_i)
uni.all.qtl$dap_i<-as.numeric(as.character(uni.all.qtl$dap_i))
uni.all.qtl$qtl.name<-paste(uni.all.qtl$chr, as.integer(uni.all.qtl$pos), sep="@")

## Remove early days (already removed in prior script)
uni.all.qtl<-uni.all.qtl[uni.all.qtl$dap_i > 16,]

###########################################################
## Lets start by plotting the location of all cumulative QTL

## Get only raw (no GxE traits)
r.uni.all.qtl<-uni.all.qtl[uni.all.qtl$type == 'raw',]

## Take only the total/cumulative traits and not the daily/rate ones
r.uni.all.total.qtl<-r.uni.all.qtl[grepl("total", r.uni.all.qtl$trait),]

## Extract the simplest trait names
full_trait_names<-r.uni.all.total.qtl$trait
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
r.uni.all.total.qtl<-cbind(r.uni.all.total.qtl, raw_trait)

## Lets combine the base trait and treatment
r.uni.all.total.qtl$trait.treatment<-paste(r.uni.all.total.qtl$raw_trait, r.uni.all.total.qtl$treatment, sep=".")

unique(r.uni.all.total.qtl[r.uni.all.total.qtl$trait.treatment == 'residual.dry', 'qtl.name'])

## Re-order the QTLs based upon order not character or factor
r.t.qtls<-unique(r.uni.all.total.qtl$qtl.name)
r.t.qtls<-r.t.qtls[c(4,3,1,2,7,5,6,8,12,11,9,10,13,16,14,15,19,18,17,20,23,21,22)]

## Use the line below if the GxE trait ratio QTL are included
#r.t.qtls<-r.t.qtls[c(1:4,7,6,5,8,9,13,12,10,11,14:17,21,20,18,19,22,25,23,24)]

trait.treatments<-sort(unique(r.uni.all.total.qtl$trait.treatment))

## Identify which QTL are found for which trait
trait.qtl<-c()
for(t in trait.treatments){
  loci<-unique(r.uni.all.total.qtl[r.uni.all.total.qtl$trait.treatment == t, 'qtl.name'])
  trait.loci<-c(t,r.t.qtls %in% loci)
  trait.qtl<-rbind(trait.qtl, trait.loci)
}


colnames(trait.qtl)<-c('qtl.name',r.t.qtls)
rownames(trait.qtl)<-trait.qtl[,1]
trait.qtl<-trait.qtl[,-c(1)]

trait.qtl<-t(trait.qtl)
trait.qtl<-as.data.frame(trait.qtl)

## Code the presenece (nubmer) or absence (0) for each QTL
## We will use this coding to color the matrix plot

trait.qtl$fit.dry<-as.character(trait.qtl$fit.dry)
trait.qtl[trait.qtl$fit.dry == "FALSE",'fit.dry']<-c(0)
trait.qtl[trait.qtl$fit.dry == "TRUE",'fit.dry']<-c(1)

trait.qtl$fit.wet<-as.character(trait.qtl$fit.wet)
trait.qtl[trait.qtl$fit.wet == "FALSE",'fit.wet']<-c(0)
trait.qtl[trait.qtl$fit.wet == "TRUE",'fit.wet']<-c(2)


trait.qtl$residual.dry<-as.character(trait.qtl$residual.dry)
trait.qtl[trait.qtl$residual.dry == "FALSE",'residual.dry']<-c(0)
trait.qtl[trait.qtl$residual.dry == "TRUE",'residual.dry']<-c(3)

trait.qtl$residual.wet<-as.character(trait.qtl$residual.wet)
trait.qtl[trait.qtl$residual.wet == "FALSE",'residual.wet']<-c(0)
trait.qtl[trait.qtl$residual.wet == "TRUE",'residual.wet']<-c(4)


trait.qtl$sv_area.dry<-as.character(trait.qtl$sv_area.dry)
trait.qtl[trait.qtl$sv_area.dry == "FALSE",'sv_area.dry']<-c(0)
trait.qtl[trait.qtl$sv_area.dry == "TRUE",'sv_area.dry']<-c(5)

trait.qtl$sv_area.wet<-as.character(trait.qtl$sv_area.wet)
trait.qtl[trait.qtl$sv_area.wet == "FALSE",'sv_area.wet']<-c(0)
trait.qtl[trait.qtl$sv_area.wet == "TRUE",'sv_area.wet']<-c(6)


trait.qtl$water_lost.dry<-as.character(trait.qtl$water_lost.dry)
trait.qtl[trait.qtl$water_lost.dry == "FALSE",'water_lost.dry']<-c(0)
trait.qtl[trait.qtl$water_lost.dry == "TRUE",'water_lost.dry']<-c(7)

trait.qtl$water_lost.wet<-as.character(trait.qtl$water_lost.wet)
trait.qtl[trait.qtl$water_lost.wet == "FALSE",'water_lost.wet']<-c(0)
trait.qtl[trait.qtl$water_lost.wet == "TRUE",'water_lost.wet']<-c(8)


trait.qtl$wue.dry<-as.character(trait.qtl$wue.dry)
trait.qtl[trait.qtl$wue.dry == "FALSE",'wue.dry']<-c(0)
trait.qtl[trait.qtl$wue.dry == "TRUE",'wue.dry']<-c(9)

trait.qtl$wue.wet<-as.character(trait.qtl$wue.wet)
trait.qtl[trait.qtl$wue.wet == "FALSE",'wue.wet']<-c(0)
trait.qtl[trait.qtl$wue.wet == "TRUE",'wue.wet']<-c(10)

## Lets order the traits in a logical way
trait.qtl<-trait.qtl[,c(5:10,1:4)]


## Add colors to the trait cateogry the reflect the assignments in the manuscript
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

trait.qtl<-t(trait.qtl)
for(c in 1:ncol(trait.qtl)){
  trait.qtl[,c]<-as.numeric(as.character(trait.qtl[,c]))
}

## Plot the presence and absence of QTL for the total/cumulative traits

setwd(wue_results.qtl.dir)

pdf("cumulative_trait_qtl_grid.pdf")
heatmap(as.matrix(apply(trait.qtl, 2, as.numeric)), col=cols, scale='none', Rowv = NA, Colv = NA, labRow=row.names(trait.qtl), main="Cumulative trait QTL")
dev.off()


###########################################################
## Lets now do the same thing for the rate traits
## Lets start by plotting the location of all cumulative QTL
r.uni.all.qtl<-uni.all.qtl[uni.all.qtl$type == 'raw',]
r.uni.all.day.qtl<-r.uni.all.qtl[!grepl("total", r.uni.all.qtl$trait),]

full_trait_names<-r.uni.all.day.qtl$trait
raw_trait<-c()
for(i in 1:length(full_trait_names)){
  name<-full_trait_names[i]
  fields<-strsplit(as.character(name), "_")
  r.trait<-fields[[1]][1:3]
  r.trait<-paste(r.trait[1], r.trait[2], r.trait[3], sep="_")
  raw_trait<-c(raw_trait, r.trait)
}

## simplify trait names
raw_trait[grepl("sv_area", raw_trait)]<-c("sv_area")
raw_trait[grepl("water_lost", raw_trait)]<-c("water_lost")
raw_trait[grepl("wue", raw_trait)]<-c("wue")
raw_trait[grepl("te_fit", raw_trait)]<-c("fit")
raw_trait[grepl("te_residual", raw_trait)]<-c("residual")

## Add this base trait name into the table
r.uni.all.day.qtl<-cbind(r.uni.all.day.qtl, raw_trait)

## Lets combine the base trait and treatment
r.uni.all.day.qtl$trait.treatment<-paste(r.uni.all.day.qtl$raw_trait, r.uni.all.day.qtl$treatment, sep=".")

## Order the QTL to reflect their physical order not factor or character string order
r.d.qtls<-unique(r.uni.all.day.qtl$qtl.name)
r.d.qtls<-r.d.qtls[c(1,6,4,5,2,3,8,11,7,10,9,12,16,15,13,14,17:19,22,21,20,23:24,27,25,26)]

## Use the line below if the GxE trait ratio QTL are included
#r.d.qtls<-r.d.qtls[c(1,2,4,3,5,6,9,11,8,7,10,13,12,17,16,14,15,18,19:21,25,24,22,23,26:28,31,29,30)]

trait.treatments<-sort(unique(as.character((r.uni.all.day.qtl$trait.treatment))))

trait.qtl_rate<-c()
for(t in trait.treatments){
  loci<-unique(r.uni.all.day.qtl[r.uni.all.day.qtl$trait.treatment == t, 'qtl.name'])
  trait.loci<-c(t,r.d.qtls %in% loci)
  trait.qtl_rate<-rbind(trait.qtl_rate, trait.loci)
}


colnames(trait.qtl_rate)<-c('qtl.name',r.d.qtls)
rownames(trait.qtl_rate)<-trait.qtl_rate[,1]
trait.qtl_rate<-trait.qtl_rate[,-c(1)]

trait.qtl_rate<-t(trait.qtl_rate)
trait.qtl_rate<-as.data.frame(trait.qtl_rate)

trait.qtl_rate$fit.dry<-as.character(trait.qtl_rate$fit.dry)
trait.qtl_rate[trait.qtl_rate$fit.dry == "FALSE",'fit.dry']<-c(0)
trait.qtl_rate[trait.qtl_rate$fit.dry == "TRUE",'fit.dry']<-c(1)

trait.qtl_rate$fit.wet<-as.character(trait.qtl_rate$fit.wet)
trait.qtl_rate[trait.qtl_rate$fit.wet == "FALSE",'fit.wet']<-c(0)
trait.qtl_rate[trait.qtl_rate$fit.wet == "TRUE",'fit.wet']<-c(2)


trait.qtl_rate$residual.dry<-as.character(trait.qtl_rate$residual.dry)
trait.qtl_rate[trait.qtl_rate$residual.dry == "FALSE",'residual.dry']<-c(0)
trait.qtl_rate[trait.qtl_rate$residual.dry == "TRUE",'residual.dry']<-c(3)

trait.qtl_rate$residual.wet<-as.character(trait.qtl_rate$residual.wet)
trait.qtl_rate[trait.qtl_rate$residual.wet == "FALSE",'residual.wet']<-c(0)
trait.qtl_rate[trait.qtl_rate$residual.wet == "TRUE",'residual.wet']<-c(4)


trait.qtl_rate$sv_area.dry<-as.character(trait.qtl_rate$sv_area.dry)
trait.qtl_rate[trait.qtl_rate$sv_area.dry == "FALSE",'sv_area.dry']<-c(0)
trait.qtl_rate[trait.qtl_rate$sv_area.dry == "TRUE",'sv_area.dry']<-c(5)

trait.qtl_rate$sv_area.wet<-as.character(trait.qtl_rate$sv_area.wet)
trait.qtl_rate[trait.qtl_rate$sv_area.wet == "FALSE",'sv_area.wet']<-c(0)
trait.qtl_rate[trait.qtl_rate$sv_area.wet == "TRUE",'sv_area.wet']<-c(6)


trait.qtl_rate$water_lost.dry<-as.character(trait.qtl_rate$water_lost.dry)
trait.qtl_rate[trait.qtl_rate$water_lost.dry == "FALSE",'water_lost.dry']<-c(0)
trait.qtl_rate[trait.qtl_rate$water_lost.dry == "TRUE",'water_lost.dry']<-c(7)

trait.qtl_rate$water_lost.wet<-as.character(trait.qtl_rate$water_lost.wet)
trait.qtl_rate[trait.qtl_rate$water_lost.wet == "FALSE",'water_lost.wet']<-c(0)
trait.qtl_rate[trait.qtl_rate$water_lost.wet == "TRUE",'water_lost.wet']<-c(8)


trait.qtl_rate$wue.dry<-as.character(trait.qtl_rate$wue.dry)
trait.qtl_rate[trait.qtl_rate$wue.dry == "FALSE",'wue.dry']<-c(0)
trait.qtl_rate[trait.qtl_rate$wue.dry == "TRUE",'wue.dry']<-c(9)

trait.qtl_rate$wue.wet<-as.character(trait.qtl_rate$wue.wet)
trait.qtl_rate[trait.qtl_rate$wue.wet == "FALSE",'wue.wet']<-c(0)
trait.qtl_rate[trait.qtl_rate$wue.wet == "TRUE",'wue.wet']<-c(10)

trait.qtl_rate<-trait.qtl_rate[,c(5:10,1:4)]

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

trait.qtl_rate<-t(trait.qtl_rate)
for(c in 1:ncol(trait.qtl_rate)){
  trait.qtl_rate[,c]<-as.numeric(as.character(trait.qtl_rate[,c]))
}

pdf("rate_trait_qtl_grid.pdf")
heatmap(as.matrix(apply(trait.qtl_rate, 2, as.numeric)), col=cols, scale='none', Rowv = NA, Colv = NA, labRow=row.names(trait.qtl_rate), main="Rate trait QTL")
dev.off()


##################################################################
## Lets try plotting the GxE/comparison/difference QTL in the same fashion
##################################################################

## c. tands for "comparison"
## Take only comparison traits from the all.qtl data.frame
c.all.qtl<-all.qtl[all.qtl$type != 'raw', ]

## Lets do the same thing for the 
c.uni.all.qtl<-uni.all.qtl[uni.all.qtl$type != 'raw', ]

## Get number of QTL
unique(c.uni.all.qtl$marker)
length(unique(c.uni.all.qtl$marker))
## 30 unique QTL

## Lets get a count of total number of SNPs
length(unique(c.all.qtl[c.all.qtl$type == 'comp_diff', 'marker']))
length(unique(c.all.qtl[c.all.qtl$type == 'comp_rel_diff', 'marker']))

## comp_ratio has been removed
#length(unique(c.all.qtl[c.all.qtl$type == 'comp_ratio', 'marker']))

## Lets see how many SNPs for cumulative (total) and daily (rate)

## Start with comp/GxE total SNPs
t.all.qtl<-c.all.qtl[grepl('_total',c.all.qtl$trait), ]
unique(t.all.qtl$marker)
length(unique(t.all.qtl$marker))
## 52
length(unique(t.all.qtl[t.all.qtl$type == 'comp_diff', 'marker']))
## 43
length(unique(t.all.qtl[t.all.qtl$type == 'comp_rel_diff', 'marker']))
## 40

## comp_ratio has been removed
#length(unique(t.all.qtl[t.all.qtl$type == 'comp_ratio', 'marker']))
## 143

## Get daily/rate SNPs for the cumulative traits
d.all.qtl<-c.all.qtl[!grepl('_total',c.all.qtl$trait), ]
unique(t.all.qtl$marker)
length(unique(d.all.qtl$marker))
## 65
length(unique(d.all.qtl[d.all.qtl$type == 'comp_diff', 'marker']))
## 53
length(unique(d.all.qtl[d.all.qtl$type == 'comp_rel_diff', 'marker']))
## 54

## comp_ratio has been removed
#length(unique(d.all.qtl[d.all.qtl$type == 'comp_ratio', 'marker']))
## 131

## Get out a data.frame that contains only "unified" cumluative QTL
## For the comparison/GxE traits
c.uni.all.total.qtl<-c.uni.all.qtl[grepl('_total', c.uni.all.qtl$trait),]

dap_i<-c()
for(i in 1:nrow(c.uni.all.total.qtl)){
  x<-strsplit(as.character(c.uni.all.total.qtl[i,'trait']), '_')[[1]]
  d<-x[length(x)-1]
  dap_i<-c(dap_i, d)
}

## Convert dap_i vector to numeric
dap_i<-as.numeric(dap_i)

## Stable isotopes were only measured at a single time point so day value is 'NA'
## Lets replace that value with 30 DAP
## There should be no stable isotope data left in this data.frame
dap_i[is.na(dap_i)]<-c(30)

c.uni.all.total.qtl<-cbind(c.uni.all.total.qtl, dap_i)
c.uni.all.total.qtl$dap_i<-as.numeric(as.character(c.uni.all.total.qtl$dap_i))
c.uni.all.total.qtl$qtl.name<-paste(c.uni.all.total.qtl$chr, as.integer(c.uni.all.total.qtl$pos), sep="@")

## Data from before 16 DAP has been removed previously
c.uni.all.total.qtl<-c.uni.all.total.qtl[c.uni.all.total.qtl$dap_i > 16,]

## How many unique QTL for cumulative/total comparison/GxE traits
length(unique(c.uni.all.total.qtl$marker))
## 22
unique(c.uni.all.total.qtl$marker)

## How many QTL for each type of comparison trait
length(unique(c.uni.all.total.qtl[c.uni.all.total.qtl$type == 'comp_diff','marker']))
## 20
length(unique(c.uni.all.total.qtl[c.uni.all.total.qtl$type == 'comp_rel_diff','marker']))
## 18

## comp_ratio traits have been removed
#length(unique(c.uni.all.total.qtl[c.uni.all.total.qtl$type == 'comp_ratio','marker']))
## 14

comp_diff_qtl<-unique(c.uni.all.qtl[c.uni.all.qtl$type == 'comp_diff','marker'])
comp_rel_diff_qtl<-unique(c.uni.all.qtl[c.uni.all.qtl$type == 'comp_rel_diff','marker'])

## comp_ratio traits have been removed
#comp_ratio_qtl<-unique(c.uni.all.qtl[c.uni.all.qtl$type == 'comp_ratio','marker'])

## Lets get some counts of the type of difference QTL detected
table(c.uni.all.qtl$type)
## comp_diff 564
## comp_rel_diff 546
## comp_ratio 0



###########################################################
## Lets extract the trait names from the cumulative/total comparison/GxE traits

full_trait_names<-c.uni.all.total.qtl$trait
raw_trait<-c()
for(i in 1:length(full_trait_names)){
  name<-full_trait_names[i]
  fields<-strsplit(as.character(name), "_")
  r.trait<-fields[[1]][1:3]
  r.trait<-paste(r.trait[1], r.trait[2], r.trait[3], sep="_")
  raw_trait<-c(raw_trait, r.trait)
}

## simplify trait names
raw_trait[grepl("sv_area", raw_trait)]<-c("sv_area")
raw_trait[grepl("water_lost", raw_trait)]<-c("water_lost")
raw_trait[grepl("wue", raw_trait)]<-c("wue")
raw_trait[grepl("te_fit", raw_trait)]<-c("fit")
raw_trait[grepl("te_residual", raw_trait)]<-c("residual")

## Add this base trait name into the table
c.uni.all.total.qtl<-cbind(c.uni.all.total.qtl, raw_trait)


diff_qtls<-unique(c.uni.all.total.qtl$qtl.name)

## order the QTLs by their location in the genome not factor level or character order
diff_qtls<-diff_qtls[c(1,5,4,2,3,8,6,7,14,13,12,9,10,11,15,19,18,17,20,16,21,22)]

traits<-unique(raw_trait)

## Get the QTL locations associated with each trait and put them into a data.frame
trait.qtl_diff<-c()
for(t in traits){
  loci<-unique(c.uni.all.total.qtl[c.uni.all.total.qtl$raw_trait == t, 'qtl.name'])
  trait.loci<-c(t,diff_qtls %in% loci)
  trait.qtl_diff<-rbind(trait.qtl_diff, trait.loci)
}


colnames(trait.qtl_diff)<-c('qtl.name',diff_qtls)
rownames(trait.qtl_diff)<-trait.qtl_diff[,1]
trait.qtl_diff<-trait.qtl_diff[,-c(1)]

trait.qtl_diff<-t(trait.qtl_diff)
trait.qtl_diff<-as.data.frame(trait.qtl_diff)


## Recode the QTL locations as indicator variables
trait.qtl_diff$fit<-as.character(trait.qtl_diff$fit)
trait.qtl_diff[trait.qtl_diff$fit == "FALSE",'fit']<-c(0)
trait.qtl_diff[trait.qtl_diff$fit == "TRUE",'fit']<-c(1)

trait.qtl_diff$residual<-as.character(trait.qtl_diff$residual)
trait.qtl_diff[trait.qtl_diff$residual == "FALSE",'residual']<-c(0)
trait.qtl_diff[trait.qtl_diff$residual == "TRUE",'residual']<-c(2)


trait.qtl_diff$sv_area<-as.character(trait.qtl_diff$sv_area)
trait.qtl_diff[trait.qtl_diff$sv_area == "FALSE",'sv_area']<-c(0)
trait.qtl_diff[trait.qtl_diff$sv_area == "TRUE",'sv_area']<-c(3)


trait.qtl_diff$water_lost<-as.character(trait.qtl_diff$water_lost)
trait.qtl_diff[trait.qtl_diff$water_lost == "FALSE",'water_lost']<-c(0)
trait.qtl_diff[trait.qtl_diff$water_lost == "TRUE",'water_lost']<-c(4)


trait.qtl_diff$wue<-as.character(trait.qtl_diff$wue)
trait.qtl_diff[trait.qtl_diff$wue == "FALSE",'wue']<-c(0)
trait.qtl_diff[trait.qtl_diff$wue == "TRUE",'wue']<-c(5)


## Re-order the traits to look nice
trait.qtl_diff<-trait.qtl_diff[,c(3,4,5,2,1)]

## Assign the values in the data.frame colors
cols <- c(
  '0' = "white",
  '1' = "black",
  '2' = "red",
  '3' = "dark green",
  '4' = "navy",
  '5' = "orange"
)

trait.qtl_diff<-t(trait.qtl_diff)
for(c in 1:ncol(trait.qtl_diff)){
  trait.qtl_diff[,c]<-as.numeric(as.character(trait.qtl_diff[,c]))
}

## Plot the locations of the QTL
pdf("diff_trait_qtl_grid.pdf")
heatmap(as.matrix(apply(trait.qtl_diff, 2, as.numeric)), col=cols, scale='none', Rowv = NA, Colv = NA, labRow=row.names(trait.qtl_diff), main="Difference trait QTL")
dev.off()


## Save.image
setwd(wue_results.qtl.dir)
save.image("qtl_trait_matrix.Rdata")





