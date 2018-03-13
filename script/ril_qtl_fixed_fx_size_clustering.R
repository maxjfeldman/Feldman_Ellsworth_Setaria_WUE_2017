library(ggplot2)
library(gplots)
library(pvclust)
library(WGCNA)

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

wue_results.qtl.clustering.dir<-paste(wue_results.qtl.dir, '/qtl_clustering', sep="")

## If this directory doesn't exist the download is not complete
if (file.exists(wue_results.qtl.clustering.dir)){
  print("Directory exists!")
} else {
  dir.create(file.path(wue_results.qtl.clustering.dir))
}


setwd(wue_results.qtl.dir)
load("fixed_marker_analysis.Rdata")

setwd(wue_results.qtl.clustering.dir)

## Combine both total and day to get total rate and remove difference QTL
all_total_rate.qtl<-rbind(all_total.qtl, all_day.qtl)
r.all_total_rate.qtl<-all_total_rate.qtl[all_total_rate.qtl$type == 'raw',]

######################################################################################
## Clustering across both dry and wet together (cumulative traits)
######################################################################################

## Make a new category 
r.all_total_rate.qtl$trait.treatment<-paste(r.all_total_rate.qtl$trait, r.all_total_rate.qtl$treatment, sep="_")

## Keep as cumulative 
r.all_total.qtl$trait.treatment<-paste(r.all_total.qtl$trait, r.all_total.qtl$treatment, sep="_")

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

setwd(wue_results.qtl.clustering.dir)

pdf("clustering_of_cumulative_traits_wardsD.pdf")

plotDendroAndColors(t.fit, colors=all.col, groupLabels=c("Days", "Treatment", "Trait"), dendroLabels = FALSE, main="Dendrogram of Traits")

dev.off()

## How many significant clusters?
wss <- (nrow(t(r.all_total.qtl.clust.mat))-1)*sum(apply(t(r.all_total.qtl.clust.mat),2,var))
for (i in 2:13) wss[i] <- sum(kmeans(t(r.all_total.qtl.clust.mat),
                                     centers=i)$withinss)


setwd(wue_results.qtl.clustering.dir)

pdf("scree_of_cumulative_trait_clusters.pdf")
plot(1:13, wss, type="b", xlab="Number of Clusters",
     ylab="Within groups sum of squares")
dev.off()

plot(t.fit) # display dendogram
rect.hclust(t.fit, k=3, border="red")


## Lets try using the pvclust algorithm
##PV_r.all_total.qtl.clust.mat<-pvclust(r.all_total.qtl.clust.mat, method.dist=c("euclidean"), method.hclust=c("ward.D"))
##plot(PV_r.all_total.qtl.clust.mat)


r.all_day.qtl<-all_day.qtl[all_day.qtl$type == 'raw',]


## Make a new category 
r.all_day.qtl$trait.treatment<-paste(r.all_day.qtl$trait, r.all_day.qtl$treatment, sep="_")


traits<-unique(r.all_day.qtl$trait.treatment)
r.all_day.qtl.clust<-c()
for(t in 1:length(traits)){
  t.name<-traits[t]
  if (t == 1) {
    temp<-r.all_day.qtl[r.all_day.qtl$trait.treatment == t.name,]
    temp2<-temp[,c("qtl.name", "signed.prop.var")]
    colnames(temp2)[2]<-t.name
    r.all_day.qtl.clust<-temp2
  }
  
  if (t > 1) {
    temp<-r.all_day.qtl[r.all_day.qtl$trait.treatment == t.name,]
    temp2<-temp[,c("qtl.name", "signed.prop.var")]
    colnames(temp2)[2]<-t.name
    r.all_day.qtl.clust<-merge(r.all_day.qtl.clust, temp2, by = c("qtl.name"))
    
  }
}

r.all_day.qtl.clust.mat<-r.all_day.qtl.clust
row.names(r.all_day.qtl.clust.mat)<-r.all_day.qtl.clust.mat$qtl.name
r.all_day.qtl.clust.mat<-r.all_day.qtl.clust.mat[,-c(1)]


## Lets calculate distance between experiments
dd <- dist(t(r.all_day.qtl.clust.mat), method = "euclidean") # distance matrix
d.fit <- hclust(dd, method="ward.D") 
plot(d.fit, cex=0.8) # display dendogram

## Get the trait names and order of the names from the cluster fit

labs<-d.fit$labels
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
raw_trait[raw_trait == "sv_area_day"]<-c("sv_area")
raw_trait[raw_trait == "water_lost_NA"]<-c("water_lost")
raw_trait[raw_trait == "wue_day_NA"]<-c("wue")
raw_trait[raw_trait == "te_fit_day"]<-c("fit")
raw_trait[raw_trait == "te_residual_day"]<-c("residual")

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

setwd(wue_results.qtl.clustering.dir)

pdf("clustering_of_all_traits_wardsD.pdf")

plotDendroAndColors(d.fit, colors=all.col, groupLabels=c("Days", "Treatment", "Trait"), dendroLabels = FALSE, main="Dendrogram of Traits")

dev.off()

## How many significant clusters?
wss <- (nrow(t(r.all_day.qtl.clust.mat))-1)*sum(apply(t(r.all_day.qtl.clust.mat),2,var))
for (i in 2:13) wss[i] <- sum(kmeans(t(r.all_day.qtl.clust.mat),
                                     centers=i)$withinss)


setwd(wue_results.qtl.clustering.dir)

pdf("scree_of_rate_trait_clusters.pdf")
plot(1:13, wss, type="b", xlab="Number of Clusters",
     ylab="Within groups sum of squares")
dev.off()

plot(d.fit) # display dendogram
rect.hclust(d.fit, k=3, border="red")

setwd(wue_results.qtl.clustering.dir)
save.image("qtl_fixed_fx_size_clustering.Rdata")


#######################################################################
## Script below contains no Figures in manuscript
#######################################################################

## Start with euclidean distance between QTL
d <- dist(r.all_total.qtl.clust.mat, method = "euclidean") # distance matrix
fit <- hclust(d, method="ward.D") 
plot(fit) # display dendogram

wss <- (nrow(r.all_total.qtl.clust.mat)-1)*sum(apply(r.all_total.qtl.clust.mat,2,var))
for (i in 2:10) wss[i] <- sum(kmeans(r.all_total.qtl.clust.mat,
                                     centers=i)$withinss)
plot(1:10, wss, type="b", xlab="Number of Clusters",
     ylab="Within groups sum of squares")

## Lets try using the pvclust algorithm
#PV_wr.all_total.qtl.clust.mat<-pvclust(t(wr.all_total.qtl.clust.mat), method.dist=c("euclidean"), method.hclust=c("ward.D"))
#plot(PV_wr.all_total.qtl.clust.mat)

plot(fit) # display dendogram
rect.hclust(fit, k=4, border="red")

## Here we are making a plot of the euclidean distance between QTL
nHalf <- 100

Min <- min(r.all_total.qtl.clust.mat)
Max <- max(r.all_total.qtl.clust.mat)
Thresh <- 0

## Make vector of colors for values below threshold
rc1 <- colorRampPalette(colors = c("blue", "white"), space="Lab")(nHalf)    
## Make vector of colors for values above threshold
rc2 <- colorRampPalette(colors = c("white", "red"), space="Lab")(nHalf)
rampcols <- c(rc1, rc2)

rb1 <- seq(Min, Thresh, length.out=nHalf+1)
rb2 <- seq(Thresh, Max, length.out=nHalf+1)[-1]
rampbreaks <- c(rb1, rb2)

#pdf("all_total_both.qtl.clust.pdf")
heatmap.2(as.matrix(r.all_total.qtl.clust.mat), distfun = function(x) dist(r.all_total.qtl.clust.mat,method = 'euclidean'), hclustfun = function(x) hclust(x,method = 'ward.D'), col=rampcols, breaks=rampbreaks, dendrogram=c("row"), Colv="NA", trace="none", density.info="none", key=F,margins=c(10,5))
#dev.off()


########## Lets try signed correlation distance now

## First calculate correlation between genes, then calculate a distance matrix based upon these

## By keeping values signed, this will identify general genetic trends that influence trait

cor.r.all_total.qtl.clust.mat<- cor(t(r.all_total.qtl.clust.mat))
cd <- dist(cor.r.all_total.qtl.clust.mat, method = "euclidean") # distance matrix

## Use wards method to cluster
c.fit <- hclust(cd, method="ward.D") 
plot(c.fit) # display dendogram

## Lets see how many clusters to expect using scree plot 
## wss is within sum of squares
c.wss <- (nrow(cor.r.all_total.qtl.clust.mat)-1)*sum(apply(cor.r.all_total.qtl.clust.mat,2,var))
for (i in 2:10) c.wss[i] <- sum(kmeans(cor.r.all_total.qtl.clust.mat,
                                       centers=i)$withinss)
plot(1:10, c.wss, type="b", xlab="Number of Clusters",
     ylab="Within groups sum of squares")

## Lets try using the pvclust algorithm
#PV_cor.wr.all_total.qtl.clust.mat<-pvclust(cor.wr.all_total.qtl.clust.mat, method.hclust=c("ward.D"))
#plot(PV_cor.wr.all_total.qtl.clust.mat)

groups <- cutree(fit, k=6) # cut tree into 6 clusters
# draw dendogram with red borders around the 6 clusters 
c.fit <- hclust(cd, method="ward.D") 
plot(c.fit) # display dendogram
rect.hclust(c.fit, k=6, border="red")


## Lets make heatmap of correlation distance between QTL 
## By keeping values signed, this will identify general genetic trends that influence trait
#pdf("all_total_both_signed.cor_qtl.clust.pdf")
heatmap.2(as.matrix(r.all_total.qtl.clust.mat),  col=rampcols, breaks=rampbreaks, distfun = function(x) dist(cor.r.all_total.qtl.clust.mat, method = 'euclidean'), hclustfun = function(x) hclust(x,method = 'ward.D'), dendrogram=c("row"), Colv="NA", trace="none", density.info="none", key=F,margins=c(10,5))
#dev.off()

########## Lets try unsigned correlation distance now

## Here we are looking for similar temporal trends in fx size regardless of directional effect
## Make negative values positive
p.r.all_total.qtl.clust.mat<-abs(r.all_total.qtl.clust.mat)
p.cor.r.all_total.qtl.clust.mat<- cor(t(p.r.all_total.qtl.clust.mat))
p.cd <- dist(p.cor.r.all_total.qtl.clust.mat, method = "euclidean") # distance matrix

## Use wards method to cluster
p.c.fit <- hclust(p.cd, method="ward.D") 
plot(p.c.fit) # display dendogram

## Lets see how many clusters to expect using scree plot 
## wss is within sum of squares
p.c.wss <- (nrow(p.cor.r.all_total.qtl.clust.mat)-1)*sum(apply(p.cor.r.all_total.qtl.clust.mat,2,var))
for (i in 2:10) p.c.wss[i] <- sum(kmeans(p.cor.r.all_total.qtl.clust.mat,
                                         centers=i)$withinss)
plot(1:10, p.c.wss, type="b", xlab="Number of Clusters",
     ylab="Within groups sum of squares")

## Lets try using the pvclust algorithm
#PV_cor.wr.all_total.qtl.clust.mat<-pvclust(cor.wr.all_total.qtl.clust.mat, method.hclust=c("ward.D"))
#plot(PV_cor.wr.all_total.qtl.clust.mat)

groups <- cutree(p.c.fit, k=6) # cut tree into 6 clusters
# draw dendogram with red borders around the 6 clusters 
p.c.fit <- hclust(cd, method="ward.D") 
plot(p.c.fit) # display dendogram
rect.hclust(p.c.fit, k=6, border="red")


## Lets make heatmap of correlation distance between QTL 
## Here we are looking for similar temporal trends in fx size regardless of directional effect
#pdf("all_total_both_unsigned.cor_qtl.clust.pdf")
heatmap.2(as.matrix(p.r.all_total.qtl.clust.mat),  col=rampcols, breaks=rampbreaks, distfun = function(x) dist(p.cor.r.all_total.qtl.clust.mat, method = 'euclidean'), hclustfun = function(x) hclust(x,method = 'ward.D'), dendrogram=c("row"), Colv="NA", trace="none", density.info="none", key=F,margins=c(10,5))
#dev.off()

####################################
## Lets split up into QTL of positive and negative effects and do the same thing
####################################
qtl.ave<-apply(r.all_total.qtl.clust.mat, 1, mean)
pos.fx<-names(qtl.ave[qtl.ave > 0])
neg.fx<-names(qtl.ave[qtl.ave < 0])

pos.r.all_total.qtl<-r.all_total.qtl.clust.mat[rownames(r.all_total.qtl.clust.mat) %in% pos.fx,]

cor.pos.r.all_total.qtl.clust.mat<- cor(t(pos.r.all_total.qtl))

cd <- dist(cor.pos.r.all_total.qtl.clust.mat, method = "euclidean") # distance matrix

## Use wards method to cluster
c.fit <- hclust(cd, method="ward.D") 
plot(c.fit) # display dendogram

## Lets see how many clusters to expect using scree plot 
## wss is within sum of squares
c.wss <- (nrow(cor.pos.r.all_total.qtl.clust.mat)-1)*sum(apply(cor.pos.r.all_total.qtl.clust.mat,2,var))
for (i in 2:6) c.wss[i] <- sum(kmeans(cor.pos.r.all_total.qtl.clust.mat,
                                      centers=i)$withinss)
plot(1:6, c.wss, type="b", xlab="Number of Clusters",
     ylab="Within groups sum of squares")

groups <- cutree(fit, k=4) # cut tree into 4 clusters
# draw dendogram with red borders around the 4 clusters 
c.fit <- hclust(cd, method="ward.D") 
plot(c.fit) # display dendogram
rect.hclust(c.fit, k=4, border="red")

my_palette <- colorRampPalette(c("white", "pink"))(n = 50)

## Lets make heatmap of correlation distance between QTL (Not as informative)
#pdf("all_total_both_pos.fx.cor_qtl.clust.pdf")
heatmap.2(as.matrix(pos.r.all_total.qtl), distfun = function(x) dist(cor.pos.r.all_total.qtl.clust.mat, method = 'euclidean'), hclustfun = function(x) hclust(x,method = 'ward.D'), col=my_palette, dendrogram=c("row"), Colv="NA", trace="none", density.info="none", key=F,margins=c(10,5))
#dev.off()

## Negative 
neg.r.all_total.qtl<-r.all_total.qtl.clust.mat[rownames(r.all_total.qtl.clust.mat) %in% neg.fx,]

cor.neg.r.all_total.qtl.clust.mat<- cor(t(neg.r.all_total.qtl))

cd <- dist(cor.neg.r.all_total.qtl.clust.mat, method = "euclidean") # distance matrix

## Use wards method to cluster
c.fit <- hclust(cd, method="ward.D") 
plot(c.fit) # display dendogram

## Lets see how many clusters to expect using scree plot 
## wss is within sum of squares
c.wss <- (nrow(cor.neg.r.all_total.qtl.clust.mat)-1)*sum(apply(cor.neg.r.all_total.qtl.clust.mat,2,var))
for (i in 2:6) c.wss[i] <- sum(kmeans(cor.neg.r.all_total.qtl.clust.mat,
                                      centers=i)$withinss)
plot(1:6, c.wss, type="b", xlab="Number of Clusters",
     ylab="Within groups sum of squares")

groups <- cutree(fit, k=3) # cut tree into 4 clusters
# draw dendogram with red borders around the 4 clusters 
c.fit <- hclust(cd, method="ward.D") 
plot(c.fit) # display dendogram
rect.hclust(c.fit, k=3, border="red")

my_palette <- colorRampPalette(c("blue", "white"))(n = 50)
## Lets make heatmap of correlation distance between QTL (Not as informative)
#pdf("all_total_both_neg.fx.cor_qtl.clust.pdf")
heatmap.2(as.matrix(neg.r.all_total.qtl), distfun = function(x) dist(cor.neg.r.all_total.qtl.clust.mat, method = 'euclidean'), hclustfun = function(x) hclust(x,method = 'ward.D'), col=my_palette, dendrogram=c("row"), Colv="NA", trace="none", density.info="none", key=F,margins=c(10,5))
#dev.off()


######################################################################################
## Clustering across both dry and wet together (rate traits)
######################################################################################

## Keep as cumulative 
r.all_day.qtl$trait.treatment<-paste(r.all_day.qtl$trait, r.all_day.qtl$treatment, sep="_")


traits<-unique(r.all_day.qtl$trait.treatment)
r.all_day.qtl.clust<-c()
for(t in 1:length(traits)){
  t.name<-traits[t]
  if (t == 1) {
    temp<-r.all_day.qtl[r.all_day.qtl$trait.treatment == t.name,]
    temp2<-temp[,c("qtl.name", "signed.prop.var")]
    colnames(temp2)[2]<-t.name
    r.all_day.qtl.clust<-temp2
  }
  
  if (t > 1) {
    temp<-r.all_day.qtl[r.all_day.qtl$trait.treatment == t.name,]
    temp2<-temp[,c("qtl.name", "signed.prop.var")]
    colnames(temp2)[2]<-t.name
    r.all_day.qtl.clust<-merge(r.all_day.qtl.clust, temp2, by = c("qtl.name"))
    
  }
}

r.all_day.qtl.clust.mat<-r.all_day.qtl.clust
row.names(r.all_day.qtl.clust.mat)<-r.all_day.qtl.clust.mat$qtl.name
r.all_day.qtl.clust.mat<-r.all_day.qtl.clust.mat[,-c(1)]


## Start with euclidean distance between QTL
d <- dist(r.all_day.qtl.clust.mat, method = "euclidean") # distance matrix
fit <- hclust(d, method="ward.D") 
plot(fit) # display dendogram

wss <- (nrow(r.all_day.qtl.clust.mat)-1)*sum(apply(r.all_day.qtl.clust.mat,2,var))
for (i in 2:10) wss[i] <- sum(kmeans(r.all_day.qtl.clust.mat,
                                     centers=i)$withinss)
plot(1:10, wss, type="b", xlab="Number of Clusters",
     ylab="Within groups sum of squares")

## Lets try using the pvclust algorithm
#PV_wr.all_day.qtl.clust.mat<-pvclust(t(wr.all_day.qtl.clust.mat), method.dist=c("euclidean"), method.hclust=c("ward.D"))
#plot(PV_wr.all_day.qtl.clust.mat)

plot(fit) # display dendogram
rect.hclust(fit, k=4, border="red")

## Here we are making a plot of the euclidean distance between QTL
nHalf <- 100

Min <- min(r.all_day.qtl.clust.mat)
Max <- max(r.all_day.qtl.clust.mat)
Thresh <- 0

## Make vector of colors for values below threshold
rc1 <- colorRampPalette(colors = c("blue", "white"), space="Lab")(nHalf)    
## Make vector of colors for values above threshold
rc2 <- colorRampPalette(colors = c("white", "red"), space="Lab")(nHalf)
rampcols <- c(rc1, rc2)

rb1 <- seq(Min, Thresh, length.out=nHalf+1)
rb2 <- seq(Thresh, Max, length.out=nHalf+1)[-1]
rampbreaks <- c(rb1, rb2)

#heatmap.2(wr.all_day.qtl.clust.mat,  col=my_palette, Colv="NA", Rowv="NA", trace="none", density.info="none", key=F, margins=c(5,5))
#pdf("all_day_both.qtl.clust.pdf")
heatmap.2(as.matrix(r.all_day.qtl.clust.mat), distfun = function(x) dist(r.all_day.qtl.clust.mat,method = 'euclidean'), hclustfun = function(x) hclust(x,method = 'ward.D'), col=rampcols, breaks=rampbreaks, dendrogram=c("row"), Colv="NA", trace="none", density.info="none", key=F,margins=c(10,5))
#dev.off()


########## Lets try signed correlation distance now

## First calculate correlation between genes, then calculate a distance matrix based upon these

## By keeping values signed, this will identify general genetic trends that influence trait

cor.r.all_day.qtl.clust.mat<- cor(t(r.all_day.qtl.clust.mat))
cd <- dist(cor.r.all_day.qtl.clust.mat, method = "euclidean") # distance matrix

## Use wards method to cluster
c.fit <- hclust(cd, method="ward.D") 
plot(c.fit) # display dendogram

## Lets see how many clusters to expect using scree plot 
## wss is within sum of squares
c.wss <- (nrow(cor.r.all_day.qtl.clust.mat)-1)*sum(apply(cor.r.all_day.qtl.clust.mat,2,var))
for (i in 2:10) c.wss[i] <- sum(kmeans(cor.r.all_day.qtl.clust.mat,
                                       centers=i)$withinss)
plot(1:10, c.wss, type="b", xlab="Number of Clusters",
     ylab="Within groups sum of squares")

## Lets try using the pvclust algorithm
#PV_cor.wr.all_day.qtl.clust.mat<-pvclust(cor.wr.all_day.qtl.clust.mat, method.hclust=c("ward.D"))
#plot(PV_cor.wr.all_day.qtl.clust.mat)

groups <- cutree(fit, k=6) # cut tree into 6 clusters
# draw dendogram with red borders around the 6 clusters 
c.fit <- hclust(cd, method="ward.D") 
plot(c.fit) # display dendogram
rect.hclust(c.fit, k=6, border="red")


## Lets make heatmap of correlation distance between QTL 
## By keeping values signed, this will identify general genetic trends that influence trait
#pdf("all_day_both_signed.cor_qtl.clust.pdf")
heatmap.2(as.matrix(r.all_day.qtl.clust.mat),  col=rampcols, breaks=rampbreaks, distfun = function(x) dist(cor.r.all_day.qtl.clust.mat, method = 'euclidean'), hclustfun = function(x) hclust(x,method = 'ward.D'), dendrogram=c("row"), Colv="NA", trace="none", density.info="none", key=F,margins=c(10,5))
#dev.off()

########## Lets try unsigned correlation distance now

## Here we are looking for similar temporal trends in fx size regardless of directional effect
## Make negative values positive
p.r.all_day.qtl.clust.mat<-abs(r.all_day.qtl.clust.mat)
p.cor.r.all_day.qtl.clust.mat<- cor(t(p.r.all_day.qtl.clust.mat))
p.cd <- dist(p.cor.r.all_day.qtl.clust.mat, method = "euclidean") # distance matrix

## Use wards method to cluster
p.c.fit <- hclust(p.cd, method="ward.D") 
plot(p.c.fit) # display dendogram

## Lets see how many clusters to expect using scree plot 
## wss is within sum of squares
p.c.wss <- (nrow(p.cor.r.all_day.qtl.clust.mat)-1)*sum(apply(p.cor.r.all_day.qtl.clust.mat,2,var))
for (i in 2:10) p.c.wss[i] <- sum(kmeans(p.cor.r.all_day.qtl.clust.mat,
                                         centers=i)$withinss)
plot(1:10, p.c.wss, type="b", xlab="Number of Clusters",
     ylab="Within groups sum of squares")

## Lets try using the pvclust algorithm
#PV_cor.wr.all_day.qtl.clust.mat<-pvclust(cor.wr.all_day.qtl.clust.mat, method.hclust=c("ward.D"))
#plot(PV_cor.wr.all_day.qtl.clust.mat)

groups <- cutree(p.c.fit, k=6) # cut tree into 6 clusters
# draw dendogram with red borders around the 6 clusters 
p.c.fit <- hclust(cd, method="ward.D") 
plot(p.c.fit) # display dendogram
rect.hclust(p.c.fit, k=6, border="red")


## Lets make heatmap of correlation distance between QTL 
## Here we are looking for similar temporal trends in fx size regardless of directional effect
#pdf("all_day_both_unsigned.cor_qtl.clust.pdf")
heatmap.2(as.matrix(p.r.all_day.qtl.clust.mat),  col=rampcols, breaks=rampbreaks, distfun = function(x) dist(p.cor.r.all_day.qtl.clust.mat, method = 'euclidean'), hclustfun = function(x) hclust(x,method = 'ward.D'), dendrogram=c("row"), Colv="NA", trace="none", density.info="none", key=F,margins=c(10,5))
#dev.off()

####################################
## Lets split up into QTL of positive and negative effects and do the same thing
####################################
qtl.ave<-apply(r.all_day.qtl.clust.mat, 1, mean)
pos.fx<-names(qtl.ave[qtl.ave > 0])
neg.fx<-names(qtl.ave[qtl.ave < 0])

pos.r.all_day.qtl<-r.all_day.qtl.clust.mat[rownames(r.all_day.qtl.clust.mat) %in% pos.fx,]

cor.pos.r.all_day.qtl.clust.mat<- cor(t(pos.r.all_day.qtl))

cd <- dist(cor.pos.r.all_day.qtl.clust.mat, method = "euclidean") # distance matrix

## Use wards method to cluster
c.fit <- hclust(cd, method="ward.D") 
plot(c.fit) # display dendogram

## Lets see how many clusters to expect using scree plot 
## wss is within sum of squares
c.wss <- (nrow(cor.pos.r.all_day.qtl.clust.mat)-1)*sum(apply(cor.pos.r.all_day.qtl.clust.mat,2,var))
for (i in 2:6) c.wss[i] <- sum(kmeans(cor.pos.r.all_day.qtl.clust.mat,
                                      centers=i)$withinss)
plot(1:6, c.wss, type="b", xlab="Number of Clusters",
     ylab="Within groups sum of squares")

groups <- cutree(fit, k=4) # cut tree into 4 clusters
# draw dendogram with red borders around the 4 clusters 
c.fit <- hclust(cd, method="ward.D") 
plot(c.fit) # display dendogram
rect.hclust(c.fit, k=4, border="red")

my_palette <- colorRampPalette(c("white", "pink"))(n = 50)
## Lets make heatmap of correlation distance between QTL (Not as informative)
#pdf("all_day_both_pos.fx.cor_qtl.clust.pdf")
heatmap.2(as.matrix(pos.r.all_day.qtl), distfun = function(x) dist(cor.pos.r.all_day.qtl.clust.mat, method = 'euclidean'), hclustfun = function(x) hclust(x,method = 'ward.D'), col=my_palette, dendrogram=c("row"), Colv="NA", trace="none", density.info="none", key=F,margins=c(10,5))
#dev.off()

## Negative 
neg.r.all_day.qtl<-r.all_day.qtl.clust.mat[rownames(r.all_day.qtl.clust.mat) %in% neg.fx,]

cor.neg.r.all_day.qtl.clust.mat<- cor(t(neg.r.all_day.qtl))

cd <- dist(cor.neg.r.all_day.qtl.clust.mat, method = "euclidean") # distance matrix

## Use wards method to cluster
c.fit <- hclust(cd, method="ward.D") 
plot(c.fit) # display dendogram

## Lets see how many clusters to expect using scree plot 
## wss is within sum of squares
c.wss <- (nrow(cor.neg.r.all_day.qtl.clust.mat)-1)*sum(apply(cor.neg.r.all_day.qtl.clust.mat,2,var))
for (i in 2:6) c.wss[i] <- sum(kmeans(cor.neg.r.all_day.qtl.clust.mat,
                                      centers=i)$withinss)
plot(1:6, c.wss, type="b", xlab="Number of Clusters",
     ylab="Within groups sum of squares")

groups <- cutree(fit, k=3) # cut tree into 4 clusters
# draw dendogram with red borders around the 4 clusters 
c.fit <- hclust(cd, method="ward.D") 
plot(c.fit) # display dendogram
rect.hclust(c.fit, k=3, border="red")

my_palette <- colorRampPalette(c("blue", "white"))(n = 50)
## Lets make heatmap of correlation distance between QTL (Not as informative)
#pdf("all_day_both_neg.fx.cor_qtl.clust.pdf")
heatmap.2(as.matrix(neg.r.all_day.qtl), distfun = function(x) dist(cor.neg.r.all_day.qtl.clust.mat, method = 'euclidean'), hclustfun = function(x) hclust(x,method = 'ward.D'), col=my_palette, dendrogram=c("row"), Colv="NA", trace="none", density.info="none", key=F,margins=c(10,5))
#dev.off()



######################################################################################
## CUMULATIVE
######################################################################################

######################################################################################
## Clustering on wet traits  (cumulative)
######################################################################################

## Lets also subset into wet and dry
wr.all_total.qtl<-r.all_total.qtl[r.all_total.qtl$treatment == 'wet',]

traits<-unique(wr.all_total.qtl$trait)
wr.all_total.qtl.clust<-c()
for(t in 1:length(traits)){
  t.name<-traits[t]
  if (t == 1) {
    temp<-wr.all_total.qtl[wr.all_total.qtl$trait == t.name,]
    temp2<-temp[,c("qtl.name", "signed.prop.var")]
    colnames(temp2)[2]<-t.name
    wr.all_total.qtl.clust<-temp2
  }
  
  if (t > 1) {
    temp<-wr.all_total.qtl[wr.all_total.qtl$trait == t.name,]
    temp2<-temp[,c("qtl.name", "signed.prop.var")]
    colnames(temp2)[2]<-t.name
    wr.all_total.qtl.clust<-merge(wr.all_total.qtl.clust, temp2, by = c("qtl.name"))
    
  }
}

wr.all_total.qtl.clust.mat<-wr.all_total.qtl.clust
row.names(wr.all_total.qtl.clust.mat)<-wr.all_total.qtl.clust.mat$qtl.name
wr.all_total.qtl.clust.mat<-wr.all_total.qtl.clust.mat[,-c(1)]



## Start with euclidean distance between QTL
d <- dist(wr.all_total.qtl.clust.mat, method = "euclidean") # distance matrix
fit <- hclust(d, method="ward.D") 
plot(fit) # display dendogram

wss <- (nrow(wr.all_total.qtl.clust.mat)-1)*sum(apply(wr.all_total.qtl.clust.mat,2,var))
for (i in 2:10) wss[i] <- sum(kmeans(wr.all_total.qtl.clust.mat,
                                     centers=i)$withinss)
plot(1:10, wss, type="b", xlab="Number of Clusters",
     ylab="Within groups sum of squares")

## Lets try using the pvclust algorithm
#PV_wwr.all_total.qtl.clust.mat<-pvclust(t(wwr.all_total.qtl.clust.mat), method.dist=c("euclidean"), method.hclust=c("ward.D"))
#plot(PV_wwr.all_total.qtl.clust.mat)

plot(fit) # display dendogram
rect.hclust(fit, k=4, border="red")

## Here we are making a plot of the euclidean distance between QTL
nHalf <- 100

Min <- min(wr.all_total.qtl.clust.mat)
Max <- max(wr.all_total.qtl.clust.mat)
Thresh <- 0

## Make vector of colors for values below threshold
rc1 <- colorRampPalette(colors = c("blue", "white"), space="Lab")(nHalf)    
## Make vector of colors for values above threshold
rc2 <- colorRampPalette(colors = c("white", "red"), space="Lab")(nHalf)
rampcols <- c(rc1, rc2)

rb1 <- seq(Min, Thresh, length.out=nHalf+1)
rb2 <- seq(Thresh, Max, length.out=nHalf+1)[-1]
rampbreaks <- c(rb1, rb2)

#pdf("all_total.wet.qtl.clust.pdf")
heatmap.2(as.matrix(wr.all_total.qtl.clust.mat), distfun = function(x) dist(wr.all_total.qtl.clust.mat,method = 'euclidean'), hclustfun = function(x) hclust(x,method = 'ward.D'), col=rampcols, breaks=rampbreaks, dendrogram=c("row"), Colv="NA", trace="none", density.info="none", key=F,margins=c(10,5))
#dev.off()


########## Lets try signed correlation distance now

## First calculate correlation between genes, then calculate a distance matrix based upon these

## By keeping values signed, this will identify general genetic trends that influence trait

cor.wr.all_total.qtl.clust.mat<- cor(t(wr.all_total.qtl.clust.mat))
cd <- dist(cor.wr.all_total.qtl.clust.mat, method = "euclidean") # distance matrix

## Use wards method to cluster
c.fit <- hclust(cd, method="ward.D") 
plot(c.fit) # display dendogram

## Lets see how many clusters to expect using scree plot 
## wss is within sum of squares
c.wss <- (nrow(cor.wr.all_total.qtl.clust.mat)-1)*sum(apply(cor.wr.all_total.qtl.clust.mat,2,var))
for (i in 2:10) c.wss[i] <- sum(kmeans(cor.wr.all_total.qtl.clust.mat,
                                       centers=i)$withinss)
plot(1:10, c.wss, type="b", xlab="Number of Clusters",
     ylab="Within groups sum of squares")

## Lets try using the pvclust algorithm
#PV_cor.wwr.all_total.qtl.clust.mat<-pvclust(cor.wwr.all_total.qtl.clust.mat, method.hclust=c("ward.D"))
#plot(PV_cor.wwr.all_total.qtl.clust.mat)

groups <- cutree(fit, k=6) # cut tree into 6 clusters
# draw dendogram with red borders around the 6 clusters 
c.fit <- hclust(cd, method="ward.D") 
plot(c.fit) # display dendogram
rect.hclust(c.fit, k=6, border="red")


## Lets make heatmap of correlation distance between QTL 
## By keeping values signed, this will identify general genetic trends that influence trait
#pdf("all_total_wet_signed.cor_qtl.clust.pdf")
heatmap.2(as.matrix(wr.all_total.qtl.clust.mat),  col=rampcols, breaks=rampbreaks, distfun = function(x) dist(cor.wr.all_total.qtl.clust.mat, method = 'euclidean'), hclustfun = function(x) hclust(x,method = 'ward.D'), dendrogram=c("row"), Colv="NA", trace="none", density.info="none", key=F,margins=c(10,5))
#dev.off()

########## Lets try unsigned correlation distance now

## Here we are looking for similar temporal trends in fx size regardless of directional effect
## Make negative values positive
p.wr.all_total.qtl.clust.mat<-abs(wr.all_total.qtl.clust.mat)
p.cor.wr.all_total.qtl.clust.mat<- cor(t(p.wr.all_total.qtl.clust.mat))
p.cd <- dist(p.cor.wr.all_total.qtl.clust.mat, method = "euclidean") # distance matrix

## Use wards method to cluster
p.c.fit <- hclust(p.cd, method="ward.D") 
plot(p.c.fit) # display dendogram

## Lets see how many clusters to expect using scree plot 
## wss is within sum of squares
p.c.wss <- (nrow(p.cor.wr.all_total.qtl.clust.mat)-1)*sum(apply(p.cor.wr.all_total.qtl.clust.mat,2,var))
for (i in 2:10) p.c.wss[i] <- sum(kmeans(p.cor.wr.all_total.qtl.clust.mat,
                                         centers=i)$withinss)
plot(1:10, p.c.wss, type="b", xlab="Number of Clusters",
     ylab="Within groups sum of squares")

## Lets try using the pvclust algorithm
#PV_cor.wwr.all_total.qtl.clust.mat<-pvclust(cor.wwr.all_total.qtl.clust.mat, method.hclust=c("ward.D"))
#plot(PV_cor.wwr.all_total.qtl.clust.mat)

groups <- cutree(p.c.fit, k=6) # cut tree into 6 clusters
# draw dendogram with red borders around the 6 clusters 
p.c.fit <- hclust(cd, method="ward.D") 
plot(p.c.fit) # display dendogram
rect.hclust(p.c.fit, k=6, border="red")


## Lets make heatmap of correlation distance between QTL 
## Here we are looking for similar temporal trends in fx size regardless of directional effect
#pdf("all_total_wet_unsigned.cor_qtl.clust.pdf")
heatmap.2(as.matrix(p.wr.all_total.qtl.clust.mat),  col=rampcols, breaks=rampbreaks, distfun = function(x) dist(p.cor.wr.all_total.qtl.clust.mat, method = 'euclidean'), hclustfun = function(x) hclust(x,method = 'ward.D'), dendrogram=c("row"), Colv="NA", trace="none", density.info="none", key=F,margins=c(10,5))
#dev.off()

####################################
## Lets split up into QTL of positive and negative effects and do the same thing
####################################
qtl.ave<-apply(wr.all_total.qtl.clust.mat, 1, mean)
pos.fx<-names(qtl.ave[qtl.ave > 0])
neg.fx<-names(qtl.ave[qtl.ave < 0])

pos.wr.all_total.qtl<-wr.all_total.qtl.clust.mat[rownames(wr.all_total.qtl.clust.mat) %in% pos.fx,]

cor.pos.wr.all_total.qtl.clust.mat<- cor(t(pos.wr.all_total.qtl))

cd <- dist(cor.pos.wr.all_total.qtl.clust.mat, method = "euclidean") # distance matrix

## Use wards method to cluster
c.fit <- hclust(cd, method="ward.D") 
plot(c.fit) # display dendogram

## Lets see how many clusters to expect using scree plot 
## wss is within sum of squares
c.wss <- (nrow(cor.pos.wr.all_total.qtl.clust.mat)-1)*sum(apply(cor.pos.wr.all_total.qtl.clust.mat,2,var))
for (i in 2:6) c.wss[i] <- sum(kmeans(cor.pos.wr.all_total.qtl.clust.mat,
                                      centers=i)$withinss)
plot(1:6, c.wss, type="b", xlab="Number of Clusters",
     ylab="Within groups sum of squares")

groups <- cutree(fit, k=4) # cut tree into 4 clusters
# draw dendogram with red borders around the 4 clusters 
c.fit <- hclust(cd, method="ward.D") 
plot(c.fit) # display dendogram
rect.hclust(c.fit, k=4, border="red")

my_palette <- colorRampPalette(c("white", "pink"))(n = 50)
## Lets make heatmap of correlation distance between QTL (Not as informative)
#pdf("all_total_wet_pos.fx.cor_qtl.clust.pdf")
heatmap.2(as.matrix(pos.wr.all_total.qtl), distfun = function(x) dist(cor.pos.wr.all_total.qtl.clust.mat, method = 'euclidean'), hclustfun = function(x) hclust(x,method = 'ward.D'), col=my_palette, dendrogram=c("row"), Colv="NA", trace="none", density.info="none", key=F,margins=c(10,5))
#dev.off()

## Negative 
neg.wr.all_total.qtl<-wr.all_total.qtl.clust.mat[rownames(wr.all_total.qtl.clust.mat) %in% neg.fx,]

cor.neg.wr.all_total.qtl.clust.mat<- cor(t(neg.wr.all_total.qtl))

cd <- dist(cor.neg.wr.all_total.qtl.clust.mat, method = "euclidean") # distance matrix

## Use wards method to cluster
c.fit <- hclust(cd, method="ward.D") 
plot(c.fit) # display dendogram

## Lets see how many clusters to expect using scree plot 
## wss is within sum of squares
c.wss <- (nrow(cor.neg.wr.all_total.qtl.clust.mat)-1)*sum(apply(cor.neg.wr.all_total.qtl.clust.mat,2,var))
for (i in 2:6) c.wss[i] <- sum(kmeans(cor.neg.wr.all_total.qtl.clust.mat,
                                      centers=i)$withinss)
plot(1:6, c.wss, type="b", xlab="Number of Clusters",
     ylab="Within groups sum of squares")

groups <- cutree(fit, k=3) # cut tree into 4 clusters
# draw dendogram with red borders around the 4 clusters 
c.fit <- hclust(cd, method="ward.D") 
plot(c.fit) # display dendogram
rect.hclust(c.fit, k=3, border="red")

my_palette <- colorRampPalette(c("blue", "white"))(n = 50)
## Lets make heatmap of correlation distance between QTL (Not as informative)
#pdf("all_total_wet_neg.fx.cor_qtl.clust.pdf")
heatmap.2(as.matrix(neg.wr.all_total.qtl), distfun = function(x) dist(cor.neg.wr.all_total.qtl.clust.mat, method = 'euclidean'), hclustfun = function(x) hclust(x,method = 'ward.D'), col=my_palette, dendrogram=c("row"), Colv="NA", trace="none", density.info="none", key=F,margins=c(10,5))
#dev.off()



######################################################################################
## Clustering on dry traits  (cumulative)
######################################################################################



## Lets also subset into wet and dry
dr.all_total.qtl<-r.all_total.qtl[r.all_total.qtl$treatment == 'dry',]

traits<-unique(dr.all_total.qtl$trait)
dr.all_total.qtl.clust<-c()
for(t in 1:length(traits)){
  t.name<-traits[t]
  if (t == 1) {
    temp<-dr.all_total.qtl[dr.all_total.qtl$trait == t.name,]
    temp2<-temp[,c("qtl.name", "signed.prop.var")]
    colnames(temp2)[2]<-t.name
    dr.all_total.qtl.clust<-temp2
  }
  
  if (t > 1) {
    temp<-dr.all_total.qtl[dr.all_total.qtl$trait == t.name,]
    temp2<-temp[,c("qtl.name", "signed.prop.var")]
    colnames(temp2)[2]<-t.name
    dr.all_total.qtl.clust<-merge(dr.all_total.qtl.clust, temp2, by = c("qtl.name"))
    
  }
}

dr.all_total.qtl.clust.mat<-dr.all_total.qtl.clust
row.names(dr.all_total.qtl.clust.mat)<-dr.all_total.qtl.clust.mat$qtl.name
dr.all_total.qtl.clust.mat<-dr.all_total.qtl.clust.mat[,-c(1)]



## Start with euclidean distance between QTL
d <- dist(dr.all_total.qtl.clust.mat, method = "euclidean") # distance matrix
fit <- hclust(d, method="ward.D") 
plot(fit) # display dendogram

wss <- (nrow(dr.all_total.qtl.clust.mat)-1)*sum(apply(dr.all_total.qtl.clust.mat,2,var))
for (i in 2:10) wss[i] <- sum(kmeans(dr.all_total.qtl.clust.mat,
                                     centers=i)$withinss)
plot(1:10, wss, type="b", xlab="Number of Clusters",
     ylab="Within groups sum of squares")

## Lets try using the pvclust algorithm
#PV_wdr.all_total.qtl.clust.mat<-pvclust(t(wdr.all_total.qtl.clust.mat), method.dist=c("euclidean"), method.hclust=c("ward.D"))
#plot(PV_wdr.all_total.qtl.clust.mat)

plot(fit) # display dendogram
rect.hclust(fit, k=4, border="red")

## Here we are making a plot of the euclidean distance between QTL
nHalf <- 100

Min <- min(dr.all_total.qtl.clust.mat)
Max <- max(dr.all_total.qtl.clust.mat)
Thresh <- 0

## Make vector of colors for values below threshold
rc1 <- colorRampPalette(colors = c("blue", "white"), space="Lab")(nHalf)    
## Make vector of colors for values above threshold
rc2 <- colorRampPalette(colors = c("white", "red"), space="Lab")(nHalf)
rampcols <- c(rc1, rc2)

rb1 <- seq(Min, Thresh, length.out=nHalf+1)
rb2 <- seq(Thresh, Max, length.out=nHalf+1)[-1]
rampbreaks <- c(rb1, rb2)

#pdf("all_total.dry.qtl.clust.pdf")
heatmap.2(as.matrix(dr.all_total.qtl.clust.mat), distfun = function(x) dist(dr.all_total.qtl.clust.mat,method = 'euclidean'), hclustfun = function(x) hclust(x,method = 'ward.D'), col=rampcols, breaks=rampbreaks, dendrogram=c("row"), Colv="NA", trace="none", density.info="none", key=F,margins=c(10,5))
#dev.off()


########## Lets try signed correlation distance now

## First calculate correlation between genes, then calculate a distance matrix based upon these

## By keeping values signed, this will identify general genetic trends that influence trait

cor.dr.all_total.qtl.clust.mat<- cor(t(dr.all_total.qtl.clust.mat))
cd <- dist(cor.dr.all_total.qtl.clust.mat, method = "euclidean") # distance matrix

## Use wards method to cluster
c.fit <- hclust(cd, method="ward.D") 
plot(c.fit) # display dendogram

## Lets see how many clusters to expect using scree plot 
## wss is within sum of squares
c.wss <- (nrow(cor.dr.all_total.qtl.clust.mat)-1)*sum(apply(cor.dr.all_total.qtl.clust.mat,2,var))
for (i in 2:10) c.wss[i] <- sum(kmeans(cor.dr.all_total.qtl.clust.mat,
                                       centers=i)$withinss)
plot(1:10, c.wss, type="b", xlab="Number of Clusters",
     ylab="Within groups sum of squares")

## Lets try using the pvclust algorithm
#PV_cor.wdr.all_total.qtl.clust.mat<-pvclust(cor.wdr.all_total.qtl.clust.mat, method.hclust=c("ward.D"))
#plot(PV_cor.wdr.all_total.qtl.clust.mat)

groups <- cutree(fit, k=6) # cut tree into 6 clusters
# draw dendogram with red borders around the 6 clusters 
c.fit <- hclust(cd, method="ward.D") 
plot(c.fit) # display dendogram
rect.hclust(c.fit, k=6, border="red")


## Lets make heatmap of correlation distance between QTL 
## By keeping values signed, this will identify general genetic trends that influence trait
#pdf("all_total_dry_signed.cor_qtl.clust.pdf")
heatmap.2(as.matrix(dr.all_total.qtl.clust.mat),  col=rampcols, breaks=rampbreaks, distfun = function(x) dist(cor.dr.all_total.qtl.clust.mat, method = 'euclidean'), hclustfun = function(x) hclust(x,method = 'ward.D'), dendrogram=c("row"), Colv="NA", trace="none", density.info="none", key=F,margins=c(10,5))
#dev.off()

########## Lets try unsigned correlation distance now

## Here we are looking for similar temporal trends in fx size regardless of directional effect
## Make negative values positive
p.dr.all_total.qtl.clust.mat<-abs(dr.all_total.qtl.clust.mat)
p.cor.dr.all_total.qtl.clust.mat<- cor(t(p.dr.all_total.qtl.clust.mat))
p.cd <- dist(p.cor.dr.all_total.qtl.clust.mat, method = "euclidean") # distance matrix

## Use wards method to cluster
p.c.fit <- hclust(p.cd, method="ward.D") 
plot(p.c.fit) # display dendogram

## Lets see how many clusters to expect using scree plot 
## wss is within sum of squares
p.c.wss <- (nrow(p.cor.dr.all_total.qtl.clust.mat)-1)*sum(apply(p.cor.dr.all_total.qtl.clust.mat,2,var))
for (i in 2:10) p.c.wss[i] <- sum(kmeans(p.cor.dr.all_total.qtl.clust.mat,
                                         centers=i)$withinss)
plot(1:10, p.c.wss, type="b", xlab="Number of Clusters",
     ylab="Within groups sum of squares")

## Lets try using the pvclust algorithm
#PV_cor.wdr.all_total.qtl.clust.mat<-pvclust(cor.wdr.all_total.qtl.clust.mat, method.hclust=c("ward.D"))
#plot(PV_cor.wdr.all_total.qtl.clust.mat)

groups <- cutree(p.c.fit, k=6) # cut tree into 6 clusters
# draw dendogram with red borders around the 6 clusters 
p.c.fit <- hclust(cd, method="ward.D") 
plot(p.c.fit) # display dendogram
rect.hclust(p.c.fit, k=6, border="red")


## Lets make heatmap of correlation distance between QTL 
## Here we are looking for similar temporal trends in fx size regardless of directional effect
#pdf("all_total_dry_unsigned.cor_qtl.clust.pdf")
heatmap.2(as.matrix(p.dr.all_total.qtl.clust.mat),  col=rampcols, breaks=rampbreaks, distfun = function(x) dist(p.cor.dr.all_total.qtl.clust.mat, method = 'euclidean'), hclustfun = function(x) hclust(x,method = 'ward.D'), dendrogram=c("row"), Colv="NA", trace="none", density.info="none", key=F,margins=c(10,5))
#dev.off()

####################################
## Lets split up into QTL of positive and negative effects and do the same thing
####################################
qtl.ave<-apply(dr.all_total.qtl.clust.mat, 1, mean)
pos.fx<-names(qtl.ave[qtl.ave > 0])
neg.fx<-names(qtl.ave[qtl.ave < 0])

pos.dr.all_total.qtl<-dr.all_total.qtl.clust.mat[rownames(dr.all_total.qtl.clust.mat) %in% pos.fx,]

cor.pos.dr.all_total.qtl.clust.mat<- cor(t(pos.dr.all_total.qtl))

cd <- dist(cor.pos.dr.all_total.qtl.clust.mat, method = "euclidean") # distance matrix

## Use wards method to cluster
c.fit <- hclust(cd, method="ward.D") 
plot(c.fit) # display dendogram

## Lets see how many clusters to expect using scree plot 
## wss is within sum of squares
c.wss <- (nrow(cor.pos.dr.all_total.qtl.clust.mat)-1)*sum(apply(cor.pos.dr.all_total.qtl.clust.mat,2,var))
for (i in 2:6) c.wss[i] <- sum(kmeans(cor.pos.dr.all_total.qtl.clust.mat,
                                      centers=i)$withinss)
plot(1:6, c.wss, type="b", xlab="Number of Clusters",
     ylab="Within groups sum of squares")

groups <- cutree(fit, k=4) # cut tree into 4 clusters
# draw dendogram with red borders around the 4 clusters 
c.fit <- hclust(cd, method="ward.D") 
plot(c.fit) # display dendogram
rect.hclust(c.fit, k=4, border="red")

my_palette <- colorRampPalette(c("white", "pink"))(n = 50)
## Lets make heatmap of correlation distance between QTL (Not as informative)
#pdf("all_total_dry_pos.fx.cor_qtl.clust.pdf")
heatmap.2(as.matrix(pos.dr.all_total.qtl), distfun = function(x) dist(cor.pos.dr.all_total.qtl.clust.mat, method = 'euclidean'), hclustfun = function(x) hclust(x,method = 'ward.D'), col=my_palette, dendrogram=c("row"), Colv="NA", trace="none", density.info="none", key=F,margins=c(10,5))
#dev.off()

## Negative 
neg.dr.all_total.qtl<-dr.all_total.qtl.clust.mat[rownames(dr.all_total.qtl.clust.mat) %in% neg.fx,]

cor.neg.dr.all_total.qtl.clust.mat<- cor(t(neg.dr.all_total.qtl))

cd <- dist(cor.neg.dr.all_total.qtl.clust.mat, method = "euclidean") # distance matrix

## Use wards method to cluster
c.fit <- hclust(cd, method="ward.D") 
plot(c.fit) # display dendogram

## Lets see how many clusters to expect using scree plot 
## wss is within sum of squares
c.wss <- (nrow(cor.neg.dr.all_total.qtl.clust.mat)-1)*sum(apply(cor.neg.dr.all_total.qtl.clust.mat,2,var))
for (i in 2:6) c.wss[i] <- sum(kmeans(cor.neg.dr.all_total.qtl.clust.mat,
                                      centers=i)$withinss)
plot(1:6, c.wss, type="b", xlab="Number of Clusters",
     ylab="Within groups sum of squares")

groups <- cutree(fit, k=3) # cut tree into 4 clusters
# draw dendogram with red borders around the 4 clusters 
c.fit <- hclust(cd, method="ward.D") 
plot(c.fit) # display dendogram
rect.hclust(c.fit, k=3, border="red")

my_palette <- colorRampPalette(c("blue", "white"))(n = 50)
## Lets make heatmap of correlation distance between QTL (Not as informative)
#pdf("all_total_dry_neg.fx.cor_qtl.clust.pdf")
heatmap.2(as.matrix(neg.dr.all_total.qtl), distfun = function(x) dist(cor.neg.dr.all_total.qtl.clust.mat, method = 'euclidean'), hclustfun = function(x) hclust(x,method = 'ward.D'), col=my_palette, dendrogram=c("row"), Colv="NA", trace="none", density.info="none", key=F,margins=c(10,5))
#dev.off()


######################################################################################
## RATE
######################################################################################

######################################################################################
## Clustering on wet traits  (rate)
######################################################################################


## Lets also subset into wet and dry
wr.all_day.qtl<-r.all_day.qtl[r.all_day.qtl$treatment == 'wet',]

traits<-unique(wr.all_day.qtl$trait)
wr.all_day.qtl.clust<-c()
for(t in 1:length(traits)){
  t.name<-traits[t]
  if (t == 1) {
    temp<-wr.all_day.qtl[wr.all_day.qtl$trait == t.name,]
    temp2<-temp[,c("qtl.name", "signed.prop.var")]
    colnames(temp2)[2]<-t.name
    wr.all_day.qtl.clust<-temp2
  }
  
  if (t > 1) {
    temp<-wr.all_day.qtl[wr.all_day.qtl$trait == t.name,]
    temp2<-temp[,c("qtl.name", "signed.prop.var")]
    colnames(temp2)[2]<-t.name
    wr.all_day.qtl.clust<-merge(wr.all_day.qtl.clust, temp2, by = c("qtl.name"))
    
  }
}

wr.all_day.qtl.clust.mat<-wr.all_day.qtl.clust
row.names(wr.all_day.qtl.clust.mat)<-wr.all_day.qtl.clust.mat$qtl.name
wr.all_day.qtl.clust.mat<-wr.all_day.qtl.clust.mat[,-c(1)]



## Start with euclidean distance between QTL
d <- dist(wr.all_day.qtl.clust.mat, method = "euclidean") # distance matrix
fit <- hclust(d, method="ward.D") 
plot(fit) # display dendogram

wss <- (nrow(wr.all_day.qtl.clust.mat)-1)*sum(apply(wr.all_day.qtl.clust.mat,2,var))
for (i in 2:10) wss[i] <- sum(kmeans(wr.all_day.qtl.clust.mat,
                                     centers=i)$withinss)
plot(1:10, wss, type="b", xlab="Number of Clusters",
     ylab="Within groups sum of squares")

## Lets try using the pvclust algorithm
#PV_wwr.all_day.qtl.clust.mat<-pvclust(t(wwr.all_day.qtl.clust.mat), method.dist=c("euclidean"), method.hclust=c("ward.D"))
#plot(PV_wwr.all_day.qtl.clust.mat)

plot(fit) # display dendogram
rect.hclust(fit, k=4, border="red")

## Here we are making a plot of the euclidean distance between QTL
nHalf <- 100

Min <- min(wr.all_day.qtl.clust.mat)
Max <- max(wr.all_day.qtl.clust.mat)
Thresh <- 0

## Make vector of colors for values below threshold
rc1 <- colorRampPalette(colors = c("blue", "white"), space="Lab")(nHalf)    
## Make vector of colors for values above threshold
rc2 <- colorRampPalette(colors = c("white", "red"), space="Lab")(nHalf)
rampcols <- c(rc1, rc2)

rb1 <- seq(Min, Thresh, length.out=nHalf+1)
rb2 <- seq(Thresh, Max, length.out=nHalf+1)[-1]
rampbreaks <- c(rb1, rb2)

#pdf("all_day.wet.qtl.clust.pdf")
heatmap.2(as.matrix(wr.all_day.qtl.clust.mat), distfun = function(x) dist(wr.all_day.qtl.clust.mat,method = 'euclidean'), hclustfun = function(x) hclust(x,method = 'ward.D'), col=rampcols, breaks=rampbreaks, dendrogram=c("row"), Colv="NA", trace="none", density.info="none", key=F,margins=c(10,5))
#dev.off()


########## Lets try signed correlation distance now

## First calculate correlation between genes, then calculate a distance matrix based upon these

## By keeping values signed, this will identify general genetic trends that influence trait

cor.wr.all_day.qtl.clust.mat<- cor(t(wr.all_day.qtl.clust.mat))
cd <- dist(cor.wr.all_day.qtl.clust.mat, method = "euclidean") # distance matrix

## Use wards method to cluster
c.fit <- hclust(cd, method="ward.D") 
plot(c.fit) # display dendogram

## Lets see how many clusters to expect using scree plot 
## wss is within sum of squares
c.wss <- (nrow(cor.wr.all_day.qtl.clust.mat)-1)*sum(apply(cor.wr.all_day.qtl.clust.mat,2,var))
for (i in 2:10) c.wss[i] <- sum(kmeans(cor.wr.all_day.qtl.clust.mat,
                                       centers=i)$withinss)
plot(1:10, c.wss, type="b", xlab="Number of Clusters",
     ylab="Within groups sum of squares")

## Lets try using the pvclust algorithm
#PV_cor.wwr.all_day.qtl.clust.mat<-pvclust(cor.wwr.all_day.qtl.clust.mat, method.hclust=c("ward.D"))
#plot(PV_cor.wwr.all_day.qtl.clust.mat)

groups <- cutree(fit, k=6) # cut tree into 6 clusters
# draw dendogram with red borders around the 6 clusters 
c.fit <- hclust(cd, method="ward.D") 
plot(c.fit) # display dendogram
rect.hclust(c.fit, k=6, border="red")


## Lets make heatmap of correlation distance between QTL 
## By keeping values signed, this will identify general genetic trends that influence trait
#pdf("all_day_wet_signed.cor_qtl.clust.pdf")
heatmap.2(as.matrix(wr.all_day.qtl.clust.mat),  col=rampcols, breaks=rampbreaks, distfun = function(x) dist(cor.wr.all_day.qtl.clust.mat, method = 'euclidean'), hclustfun = function(x) hclust(x,method = 'ward.D'), dendrogram=c("row"), Colv="NA", trace="none", density.info="none", key=F,margins=c(10,5))
#dev.off()

########## Lets try unsigned correlation distance now

## Here we are looking for similar temporal trends in fx size regardless of directional effect
## Make negative values positive
p.wr.all_day.qtl.clust.mat<-abs(wr.all_day.qtl.clust.mat)
p.cor.wr.all_day.qtl.clust.mat<- cor(t(p.wr.all_day.qtl.clust.mat))
p.cd <- dist(p.cor.wr.all_day.qtl.clust.mat, method = "euclidean") # distance matrix

## Use wards method to cluster
p.c.fit <- hclust(p.cd, method="ward.D") 
plot(p.c.fit) # display dendogram

## Lets see how many clusters to expect using scree plot 
## wss is within sum of squares
p.c.wss <- (nrow(p.cor.wr.all_day.qtl.clust.mat)-1)*sum(apply(p.cor.wr.all_day.qtl.clust.mat,2,var))
for (i in 2:10) p.c.wss[i] <- sum(kmeans(p.cor.wr.all_day.qtl.clust.mat,
                                         centers=i)$withinss)
plot(1:10, p.c.wss, type="b", xlab="Number of Clusters",
     ylab="Within groups sum of squares")

## Lets try using the pvclust algorithm
#PV_cor.wwr.all_day.qtl.clust.mat<-pvclust(cor.wwr.all_day.qtl.clust.mat, method.hclust=c("ward.D"))
#plot(PV_cor.wwr.all_day.qtl.clust.mat)

groups <- cutree(p.c.fit, k=6) # cut tree into 6 clusters
# draw dendogram with red borders around the 6 clusters 
p.c.fit <- hclust(cd, method="ward.D") 
plot(p.c.fit) # display dendogram
rect.hclust(p.c.fit, k=6, border="red")


## Lets make heatmap of correlation distance between QTL 
## Here we are looking for similar temporal trends in fx size regardless of directional effect
#pdf("all_day_wet_unsigned.cor_qtl.clust.pdf")
heatmap.2(as.matrix(p.wr.all_day.qtl.clust.mat),  col=rampcols, breaks=rampbreaks, distfun = function(x) dist(p.cor.wr.all_day.qtl.clust.mat, method = 'euclidean'), hclustfun = function(x) hclust(x,method = 'ward.D'), dendrogram=c("row"), Colv="NA", trace="none", density.info="none", key=F,margins=c(10,5))
#dev.off()

####################################
## Lets split up into QTL of positive and negative effects and do the same thing
####################################
qtl.ave<-apply(wr.all_day.qtl.clust.mat, 1, mean)
pos.fx<-names(qtl.ave[qtl.ave > 0])
neg.fx<-names(qtl.ave[qtl.ave < 0])

pos.wr.all_day.qtl<-wr.all_day.qtl.clust.mat[rownames(wr.all_day.qtl.clust.mat) %in% pos.fx,]

cor.pos.wr.all_day.qtl.clust.mat<- cor(t(pos.wr.all_day.qtl))

cd <- dist(cor.pos.wr.all_day.qtl.clust.mat, method = "euclidean") # distance matrix

## Use wards method to cluster
c.fit <- hclust(cd, method="ward.D") 
plot(c.fit) # display dendogram

## Lets see how many clusters to expect using scree plot 
## wss is within sum of squares
c.wss <- (nrow(cor.pos.wr.all_day.qtl.clust.mat)-1)*sum(apply(cor.pos.wr.all_day.qtl.clust.mat,2,var))
for (i in 2:6) c.wss[i] <- sum(kmeans(cor.pos.wr.all_day.qtl.clust.mat,
                                      centers=i)$withinss)
plot(1:6, c.wss, type="b", xlab="Number of Clusters",
     ylab="Within groups sum of squares")

groups <- cutree(fit, k=4) # cut tree into 4 clusters
# draw dendogram with red borders around the 4 clusters 
c.fit <- hclust(cd, method="ward.D") 
plot(c.fit) # display dendogram
rect.hclust(c.fit, k=4, border="red")

my_palette <- colorRampPalette(c("white", "pink"))(n = 50)
## Lets make heatmap of correlation distance between QTL (Not as informative)
#pdf("all_day_wet_pos.fx.cor_qtl.clust.pdf")
heatmap.2(as.matrix(pos.wr.all_day.qtl), distfun = function(x) dist(cor.pos.wr.all_day.qtl.clust.mat, method = 'euclidean'), hclustfun = function(x) hclust(x,method = 'ward.D'), col=my_palette, dendrogram=c("row"), Colv="NA", trace="none", density.info="none", key=F,margins=c(10,5))
#dev.off()

## Negative 
neg.wr.all_day.qtl<-wr.all_day.qtl.clust.mat[rownames(wr.all_day.qtl.clust.mat) %in% neg.fx,]

cor.neg.wr.all_day.qtl.clust.mat<- cor(t(neg.wr.all_day.qtl))

cd <- dist(cor.neg.wr.all_day.qtl.clust.mat, method = "euclidean") # distance matrix

## Use wards method to cluster
c.fit <- hclust(cd, method="ward.D") 
plot(c.fit) # display dendogram

## Lets see how many clusters to expect using scree plot 
## wss is within sum of squares
c.wss <- (nrow(cor.neg.wr.all_day.qtl.clust.mat)-1)*sum(apply(cor.neg.wr.all_day.qtl.clust.mat,2,var))
for (i in 2:6) c.wss[i] <- sum(kmeans(cor.neg.wr.all_day.qtl.clust.mat,
                                      centers=i)$withinss)
plot(1:6, c.wss, type="b", xlab="Number of Clusters",
     ylab="Within groups sum of squares")

groups <- cutree(fit, k=3) # cut tree into 4 clusters
# draw dendogram with red borders around the 4 clusters 
c.fit <- hclust(cd, method="ward.D") 
plot(c.fit) # display dendogram
rect.hclust(c.fit, k=3, border="red")

my_palette <- colorRampPalette(c("blue", "white"))(n = 50)
## Lets make heatmap of correlation distance between QTL (Not as informative)
#pdf("all_day_wet_neg.fx.cor_qtl.clust.pdf")
heatmap.2(as.matrix(neg.wr.all_day.qtl), distfun = function(x) dist(cor.neg.wr.all_day.qtl.clust.mat, method = 'euclidean'), hclustfun = function(x) hclust(x,method = 'ward.D'), col=my_palette, dendrogram=c("row"), Colv="NA", trace="none", density.info="none", key=F,margins=c(10,5))
#dev.off()



######################################################################################
## Clustering on dry traits  (rate)
######################################################################################



## Lets also subset into wet and dry
dr.all_day.qtl<-r.all_day.qtl[r.all_day.qtl$treatment == 'dry',]

traits<-unique(dr.all_day.qtl$trait)
dr.all_day.qtl.clust<-c()
for(t in 1:length(traits)){
  t.name<-traits[t]
  if (t == 1) {
    temp<-dr.all_day.qtl[dr.all_day.qtl$trait == t.name,]
    temp2<-temp[,c("qtl.name", "signed.prop.var")]
    colnames(temp2)[2]<-t.name
    dr.all_day.qtl.clust<-temp2
  }
  
  if (t > 1) {
    temp<-dr.all_day.qtl[dr.all_day.qtl$trait == t.name,]
    temp2<-temp[,c("qtl.name", "signed.prop.var")]
    colnames(temp2)[2]<-t.name
    dr.all_day.qtl.clust<-merge(dr.all_day.qtl.clust, temp2, by = c("qtl.name"))
    
  }
}

dr.all_day.qtl.clust.mat<-dr.all_day.qtl.clust
row.names(dr.all_day.qtl.clust.mat)<-dr.all_day.qtl.clust.mat$qtl.name
dr.all_day.qtl.clust.mat<-dr.all_day.qtl.clust.mat[,-c(1)]



## Start with euclidean distance between QTL
d <- dist(dr.all_day.qtl.clust.mat, method = "euclidean") # distance matrix
fit <- hclust(d, method="ward.D") 
plot(fit) # display dendogram

wss <- (nrow(dr.all_day.qtl.clust.mat)-1)*sum(apply(dr.all_day.qtl.clust.mat,2,var))
for (i in 2:10) wss[i] <- sum(kmeans(dr.all_day.qtl.clust.mat,
                                     centers=i)$withinss)
plot(1:10, wss, type="b", xlab="Number of Clusters",
     ylab="Within groups sum of squares")

## Lets try using the pvclust algorithm
#PV_wdr.all_day.qtl.clust.mat<-pvclust(t(wdr.all_day.qtl.clust.mat), method.dist=c("euclidean"), method.hclust=c("ward.D"))
#plot(PV_wdr.all_day.qtl.clust.mat)

plot(fit) # display dendogram
rect.hclust(fit, k=4, border="red")

## Here we are making a plot of the euclidean distance between QTL
nHalf <- 100

Min <- min(dr.all_day.qtl.clust.mat)
Max <- max(dr.all_day.qtl.clust.mat)
Thresh <- 0

## Make vector of colors for values below threshold
rc1 <- colorRampPalette(colors = c("blue", "white"), space="Lab")(nHalf)    
## Make vector of colors for values above threshold
rc2 <- colorRampPalette(colors = c("white", "red"), space="Lab")(nHalf)
rampcols <- c(rc1, rc2)

rb1 <- seq(Min, Thresh, length.out=nHalf+1)
rb2 <- seq(Thresh, Max, length.out=nHalf+1)[-1]
rampbreaks <- c(rb1, rb2)

#pdf("all_day.dry.qtl.clust.pdf")
heatmap.2(as.matrix(dr.all_day.qtl.clust.mat), distfun = function(x) dist(dr.all_day.qtl.clust.mat,method = 'euclidean'), hclustfun = function(x) hclust(x,method = 'ward.D'), col=rampcols, breaks=rampbreaks, dendrogram=c("row"), Colv="NA", trace="none", density.info="none", key=F,margins=c(10,5))
#dev.off()


########## Lets try signed correlation distance now

## First calculate correlation between genes, then calculate a distance matrix based upon these

## By keeping values signed, this will identify general genetic trends that influence trait

cor.dr.all_day.qtl.clust.mat<- cor(t(dr.all_day.qtl.clust.mat))
cd <- dist(cor.dr.all_day.qtl.clust.mat, method = "euclidean") # distance matrix

## Use wards method to cluster
c.fit <- hclust(cd, method="ward.D") 
plot(c.fit) # display dendogram

## Lets see how many clusters to expect using scree plot 
## wss is within sum of squares
c.wss <- (nrow(cor.dr.all_day.qtl.clust.mat)-1)*sum(apply(cor.dr.all_day.qtl.clust.mat,2,var))
for (i in 2:10) c.wss[i] <- sum(kmeans(cor.dr.all_day.qtl.clust.mat,
                                       centers=i)$withinss)
plot(1:10, c.wss, type="b", xlab="Number of Clusters",
     ylab="Within groups sum of squares")

## Lets try using the pvclust algorithm
#PV_cor.wdr.all_day.qtl.clust.mat<-pvclust(cor.wdr.all_day.qtl.clust.mat, method.hclust=c("ward.D"))
#plot(PV_cor.wdr.all_day.qtl.clust.mat)

groups <- cutree(fit, k=6) # cut tree into 6 clusters
# draw dendogram with red borders around the 6 clusters 
c.fit <- hclust(cd, method="ward.D") 
plot(c.fit) # display dendogram
rect.hclust(c.fit, k=6, border="red")


## Lets make heatmap of correlation distance between QTL 
## By keeping values signed, this will identify general genetic trends that influence trait
#pdf("all_day_dry_signed.cor_qtl.clust.pdf")
heatmap.2(as.matrix(dr.all_day.qtl.clust.mat),  col=rampcols, breaks=rampbreaks, distfun = function(x) dist(cor.dr.all_day.qtl.clust.mat, method = 'euclidean'), hclustfun = function(x) hclust(x,method = 'ward.D'), dendrogram=c("row"), Colv="NA", trace="none", density.info="none", key=F,margins=c(10,5))
#dev.off()

########## Lets try unsigned correlation distance now

## Here we are looking for similar temporal trends in fx size regardless of directional effect
## Make negative values positive
p.dr.all_day.qtl.clust.mat<-abs(dr.all_day.qtl.clust.mat)
p.cor.dr.all_day.qtl.clust.mat<- cor(t(p.dr.all_day.qtl.clust.mat))
p.cd <- dist(p.cor.dr.all_day.qtl.clust.mat, method = "euclidean") # distance matrix

## Use wards method to cluster
p.c.fit <- hclust(p.cd, method="ward.D") 
plot(p.c.fit) # display dendogram

## Lets see how many clusters to expect using scree plot 
## wss is within sum of squares
p.c.wss <- (nrow(p.cor.dr.all_day.qtl.clust.mat)-1)*sum(apply(p.cor.dr.all_day.qtl.clust.mat,2,var))
for (i in 2:10) p.c.wss[i] <- sum(kmeans(p.cor.dr.all_day.qtl.clust.mat,
                                         centers=i)$withinss)
plot(1:10, p.c.wss, type="b", xlab="Number of Clusters",
     ylab="Within groups sum of squares")

## Lets try using the pvclust algorithm
#PV_cor.wdr.all_day.qtl.clust.mat<-pvclust(cor.wdr.all_day.qtl.clust.mat, method.hclust=c("ward.D"))
#plot(PV_cor.wdr.all_day.qtl.clust.mat)

groups <- cutree(p.c.fit, k=6) # cut tree into 6 clusters
# draw dendogram with red borders around the 6 clusters 
p.c.fit <- hclust(cd, method="ward.D") 
plot(p.c.fit) # display dendogram
rect.hclust(p.c.fit, k=6, border="red")


## Lets make heatmap of correlation distance between QTL 
## Here we are looking for similar temporal trends in fx size regardless of directional effect
#pdf("all_day_dry_unsigned.cor_qtl.clust.pdf")
heatmap.2(as.matrix(p.dr.all_day.qtl.clust.mat),  col=rampcols, breaks=rampbreaks, distfun = function(x) dist(p.cor.dr.all_day.qtl.clust.mat, method = 'euclidean'), hclustfun = function(x) hclust(x,method = 'ward.D'), dendrogram=c("row"), Colv="NA", trace="none", density.info="none", key=F,margins=c(10,5))
#dev.off()

####################################
## Lets split up into QTL of positive and negative effects and do the same thing
####################################
qtl.ave<-apply(dr.all_day.qtl.clust.mat, 1, mean)
pos.fx<-names(qtl.ave[qtl.ave > 0])
neg.fx<-names(qtl.ave[qtl.ave < 0])

pos.dr.all_day.qtl<-dr.all_day.qtl.clust.mat[rownames(dr.all_day.qtl.clust.mat) %in% pos.fx,]

cor.pos.dr.all_day.qtl.clust.mat<- cor(t(pos.dr.all_day.qtl))

cd <- dist(cor.pos.dr.all_day.qtl.clust.mat, method = "euclidean") # distance matrix

## Use wards method to cluster
c.fit <- hclust(cd, method="ward.D") 
plot(c.fit) # display dendogram

## Lets see how many clusters to expect using scree plot 
## wss is within sum of squares
c.wss <- (nrow(cor.pos.dr.all_day.qtl.clust.mat)-1)*sum(apply(cor.pos.dr.all_day.qtl.clust.mat,2,var))
for (i in 2:6) c.wss[i] <- sum(kmeans(cor.pos.dr.all_day.qtl.clust.mat,
                                      centers=i)$withinss)
plot(1:6, c.wss, type="b", xlab="Number of Clusters",
     ylab="Within groups sum of squares")

groups <- cutree(fit, k=4) # cut tree into 4 clusters
# draw dendogram with red borders around the 4 clusters 
c.fit <- hclust(cd, method="ward.D") 
plot(c.fit) # display dendogram
rect.hclust(c.fit, k=4, border="red")

my_palette <- colorRampPalette(c("white", "pink"))(n = 50)
## Lets make heatmap of correlation distance between QTL (Not as informative)
#pdf("all_day_dry_pos.fx.cor_qtl.clust.pdf")
heatmap.2(as.matrix(pos.dr.all_day.qtl), distfun = function(x) dist(cor.pos.dr.all_day.qtl.clust.mat, method = 'euclidean'), hclustfun = function(x) hclust(x,method = 'ward.D'), col=my_palette, dendrogram=c("row"), Colv="NA", trace="none", density.info="none", key=F,margins=c(10,5))
#dev.off()

## Negative 
neg.dr.all_day.qtl<-dr.all_day.qtl.clust.mat[rownames(dr.all_day.qtl.clust.mat) %in% neg.fx,]

cor.neg.dr.all_day.qtl.clust.mat<- cor(t(neg.dr.all_day.qtl))

cd <- dist(cor.neg.dr.all_day.qtl.clust.mat, method = "euclidean") # distance matrix

## Use wards method to cluster
c.fit <- hclust(cd, method="ward.D") 
plot(c.fit) # display dendogram

## Lets see how many clusters to expect using scree plot 
## wss is within sum of squares
c.wss <- (nrow(cor.neg.dr.all_day.qtl.clust.mat)-1)*sum(apply(cor.neg.dr.all_day.qtl.clust.mat,2,var))
for (i in 2:6) c.wss[i] <- sum(kmeans(cor.neg.dr.all_day.qtl.clust.mat,
                                      centers=i)$withinss)
plot(1:6, c.wss, type="b", xlab="Number of Clusters",
     ylab="Within groups sum of squares")

groups <- cutree(fit, k=3) # cut tree into 4 clusters
# draw dendogram with red borders around the 4 clusters 
c.fit <- hclust(cd, method="ward.D") 
plot(c.fit) # display dendogram
rect.hclust(c.fit, k=3, border="red")

my_palette <- colorRampPalette(c("blue", "white"))(n = 50)
## Lets make heatmap of correlation distance between QTL (Not as informative)
#pdf("all_day_dry_neg.fx.cor_qtl.clust.pdf")
heatmap.2(as.matrix(neg.dr.all_day.qtl), distfun = function(x) dist(cor.neg.dr.all_day.qtl.clust.mat, method = 'euclidean'), hclustfun = function(x) hclust(x,method = 'ward.D'), col=my_palette, dendrogram=c("row"), Colv="NA", trace="none", density.info="none", key=F,margins=c(10,5))
#dev.off()

setwd(wue_results.qtl.clustering.dir)
save.image("qtl_fixed_fx_size_clustering.Rdata")

rm(list=ls())
