# PtH2O2_heatmap_sigclust.R

# Created 5/1/15 by Bethanie Edwards, bedwards@whoi.edu
# Mod'fd 09/29/15 by BRE to prepare plot for JLR manuscript
# Mod'fd 1/17/16 by JRC to work with data from new (LOBSTAHS) pipeline

# Purpose: Generate heatmap, dendrogram, and some basic stats for PtH2O2 data using
# t-tests, hierarchical clustering and simprof analysis.

# Important: Assumes user has already run "PtH2O2_mz_rt_plots.R," and all objects generated using that script are still in current session.

################# setup and data prep #############

# load necessary packages

library(gplots)
library(RColorBrewer)
library(cluster)
library(vegan)
library("clustsig")

# get only "high confidence" assignments from ptdata.QA.norm.all (created by PtH2O2_mz-rt_plots.R)

hcPtassign = subset(ptdata.QA.norm.all, confcode==1)

# create single field to use as identifier

hcPtassign$plot.label = paste0(hcPtassign$compound_name, ", RT ", round(hcPtassign$peakgroup.rt/60,2), " min")

# subset data to only 0 and 150 uM treatments at 24 h

hcPtassign = hcPtassign[,-c(12,13,14,17,18,19,22:30)]

# restrict to compounds present in both replicates, or absent from both replicates

hcPtassign = hcPtassign[(apply(hcPtassign[,12:13], 1, function(x) sum(x>0))==2) | (apply(hcPtassign[,12:13], 1, function(x) sum(x>0))==0),]

hcPtassign = hcPtassign[(apply(hcPtassign[,14:15], 1, function(x) sum(x>0))==2) | (apply(hcPtassign[,14:15], 1, function(x) sum(x>0))==0),]

# also, don't want any compounds that aren't present at all!

hcPtassign = hcPtassign[(apply(hcPtassign[,12:15], 1, function(x) sum(x==0))!=4),]

# level scale according to van den Berg et al. (2006), BMC Bioinformatics

hcPtassign$X24h_150uM_levavg = apply(hcPtassign[,14:15]/apply(hcPtassign[,12:15], 1, mean), 1, mean)
hcPtassign$X24h_0uM_levavg = apply(hcPtassign[,12:13]/apply(hcPtassign[,12:15], 1, mean), 1, mean)

# set any 0 values = 1e-6 since log2 tranformation will generate an NA if given input of 0, and many hierarchical clustering functions will not handle NA values

hcPtassign$X24h_150uM_levavg[hcPtassign$X24h_150uM_levavg<1e-6] = 1e-6
hcPtassign$X24h_0uM_levavg[hcPtassign$X24h_0uM_levavg<1e-6] = 1e-6

################# obtain stats on molecular characteristics, by lipid class #############

# define groups for which we want to calculate stats

ttest.classes = list("DGCC","DGDG","DGTS_DGTA","MGDG","PC","PE","PG","SQDG","TAG",
                  c("DGDG","SQDG","MGDG","PG"),
                  c("PC","PE","PG"),
                  c("PE","PG"),
                  c("DGCC","DGTS_DGTA"),
                  unique(hcPtassign$species))

# preallocate structure for results

lipidclass.stats = matrix(data = NA, nrow = length(ttest.classes), ncol = 17)

# cycle through classes and calculate results

for (i in 1:length(ttest.classes)) {
  
  data.thisclass = hcPtassign[hcPtassign$species %in% ttest.classes[[i]],]
  
  # subset by state of regulation
  
  data.upreg = data.thisclass[data.thisclass$X24h_150uM_levavg>data.thisclass$X24h_0uM_levavg,]
  data.downreg = data.thisclass[data.thisclass$X24h_150uM_levavg<data.thisclass$X24h_0uM_levavg,]
  
  lipidclass.stats[i,1:2] = c(nrow(data.upreg), nrow(data.downreg)) 
  
  # run t-tests
  
  if (nrow(data.upreg)>1 & nrow(data.downreg)>1) {
    
  ttest.FA_C = t.test(data.upreg$FA_total_no_C,data.downreg$FA_total_no_C)
  ttest.FA_DB = t.test(data.upreg$FA_total_no_DB,data.downreg$FA_total_no_DB)
  ttest.degox = t.test(data.upreg$degree_oxidation,data.downreg$degree_oxidation)
  
  lipidclass.stats[i,3:5] = c(ttest.FA_C$estimate[1],
                              ttest.FA_DB$estimate[1],
                              ttest.degox$estimate[1])
  
  lipidclass.stats[i,6:8] = c(ttest.FA_C$estimate[2],
                              ttest.FA_DB$estimate[2],
                              ttest.degox$estimate[2])
  
  lipidclass.stats[i,9:11] = c(ttest.FA_C$p.value,
                              ttest.FA_DB$p.value,
                              ttest.degox$p.value)
  
  lipidclass.stats[i,12:14] = c(ttest.FA_C$conf.int[1],
                                ttest.FA_DB$conf.int[1],
                                ttest.degox$conf.int[1])
  
  lipidclass.stats[i,15:17] = c(ttest.FA_C$conf.int[2],
                                ttest.FA_DB$conf.int[2],
                                ttest.degox$conf.int[2])
  
  } else {
    
    lipidclass.stats[i,3:5] = c(mean(data.upreg$FA_total_no_C),
                                mean(data.upreg$FA_total_no_DB),
                                mean(data.upreg$degree_oxidation))
    
    lipidclass.stats[i,6:8] = c(mean(data.downreg$FA_total_no_C),
                                mean(data.downreg$FA_total_no_DB),
                                mean(data.downreg$degree_oxidation))
    
    lipidclass.stats[i,12:14] = c(sd(data.upreg$FA_total_no_C),
                                  sd(data.upreg$FA_total_no_DB),
    sd(data.upreg$degree_oxidation))

    lipidclass.stats[i,15:17] = c(sd(data.downreg$FA_total_no_C),
                              sd(data.downreg$FA_total_no_DB),
                              sd(data.downreg$degree_oxidation))

    
  }
  
}

# export raw results

write.csv(lipidclass.stats, file = "PtH202_ttest_stats_raw.csv")

# reshape data, export

for (i in 1:nrow(lipidclass.stats)) {
  
  if (i==1) {
    
    stats.forexport = round(rbind(lipidclass.stats[i,3:5],lipidclass.stats[i,6:8]),1)
    
  } else {
    
    stats.forexport = rbind(stats.forexport,
                            c("","",""),
                            round(rbind(lipidclass.stats[i,3:5],lipidclass.stats[i,6:8]),1))
    
  }
  
}

write.csv(stats.forexport, file = "PtH202_ttest_stats_fortable.csv")

################# look at distribution of peak area between various classes/timepoints #############

alltreat.tp.hc = subset(ptdata.QA.norm.all, confcode==1) # get high-confidence data from all treatments & timepoints
# 
# # restrict to compounds present in both replicates, or absent from both replicates
# 
# alltreat.tp.hc = alltreat.tp.hc[(apply(alltreat.tp.hc[,20:21], 1, function(x) sum(x>0))==2) | (apply(alltreat.tp.hc[,20:21], 1, function(x) sum(x>0))==0),]
# 
# alltreat.tp.hc = alltreat.tp.hc[(apply(alltreat.tp.hc[,15:16], 1, function(x) sum(x>0))==2) | (apply(alltreat.tp.hc[,15:16], 1, function(x) sum(x>0))==0),]

# increase in oxidized moieties from 4 to 24 h

# 0 uM H2O2

frac.ox.4h.0uM = sum(alltreat.tp.hc[alltreat.tp.hc$degree_oxidation!=0,c("X0uM_4h_Orbi_0476")])/sum(alltreat.tp.hc$X0uM_4h_Orbi_0476)

frac.ox.8h.0uM = mean(c(sum(alltreat.tp.hc[alltreat.tp.hc$degree_oxidation!=0,c("X0uM_8h_Orbi_0472")])/sum(alltreat.tp.hc$X0uM_8h_Orbi_0472),
                         sum(alltreat.tp.hc[alltreat.tp.hc$degree_oxidation!=0,c("X0uM_8h_Orbi_0477")])/sum(alltreat.tp.hc$X0uM_8h_Orbi_0477)))

frac.ox.24h.0uM = mean(c(sum(alltreat.tp.hc[alltreat.tp.hc$degree_oxidation!=0,c("X0uM_24h_Orbi_0473")])/sum(alltreat.tp.hc$X0uM_24h_Orbi_0473),
                     sum(alltreat.tp.hc[alltreat.tp.hc$degree_oxidation!=0,c("X0uM_24h_Orbi_0468")])/sum(alltreat.tp.hc$X0uM_24h_Orbi_0468)))

paste(frac.ox.4h.0uM,frac.ox.8h.0uM,frac.ox.24h.0uM)

# 150 uM H2O2

frac.ox.4h.150uM = sum(alltreat.tp.hc[alltreat.tp.hc$degree_oxidation!=0,c("X150uM_4h_Orbi_0478")])/sum(alltreat.tp.hc$X150uM_4h_Orbi_0478)

frac.ox.8h.150uM = mean(c(sum(alltreat.tp.hc[alltreat.tp.hc$degree_oxidation!=0,c("X150uM_8h_Orbi_0471")])/sum(alltreat.tp.hc$X150uM_8h_Orbi_0471),
                           sum(alltreat.tp.hc[alltreat.tp.hc$degree_oxidation!=0,c("X150uM_8h_Orbi_0467")])/sum(alltreat.tp.hc$X150uM_8h_Orbi_0467)))

frac.ox.24h.150uM = mean(c(sum(alltreat.tp.hc[alltreat.tp.hc$degree_oxidation!=0,c("X150uM_24h_Orbi_0469")])/sum(alltreat.tp.hc$X150uM_24h_Orbi_0469),
                     sum(alltreat.tp.hc[alltreat.tp.hc$degree_oxidation!=0,c("X150uM_24h_Orbi_0466")])/sum(alltreat.tp.hc$X150uM_24h_Orbi_0466)))

paste(frac.ox.4h.150uM,frac.ox.8h.150uM,frac.ox.24h.150uM)

# increase in peak area as TAG from 4 to 24 h

# 0 uM H2O2

frac.TAG.4h.0uM = sum(alltreat.tp.hc[alltreat.tp.hc$species=="TAG",c("X0uM_4h_Orbi_0476")])/sum(alltreat.tp.hc$X0uM_4h_Orbi_0476)

frac.TAG.8h.0uM = mean(c(sum(alltreat.tp.hc[alltreat.tp.hc$species=="TAG",c("X0uM_8h_Orbi_0472")])/sum(alltreat.tp.hc$X0uM_8h_Orbi_0472),
                        sum(alltreat.tp.hc[alltreat.tp.hc$species=="TAG",c("X0uM_8h_Orbi_0477")])/sum(alltreat.tp.hc$X0uM_8h_Orbi_0477)))

frac.TAG.24h.0uM = mean(c(sum(alltreat.tp.hc[alltreat.tp.hc$species=="TAG",c("X0uM_24h_Orbi_0473")])/sum(alltreat.tp.hc$X0uM_24h_Orbi_0473),
                         sum(alltreat.tp.hc[alltreat.tp.hc$species=="TAG",c("X0uM_24h_Orbi_0468")])/sum(alltreat.tp.hc$X0uM_24h_Orbi_0468)))

frac.TAG.8h.0uM.sd = sd(c(sum(alltreat.tp.hc[alltreat.tp.hc$species=="TAG",c("X0uM_8h_Orbi_0472")])/sum(alltreat.tp.hc$X0uM_8h_Orbi_0472),
                         sum(alltreat.tp.hc[alltreat.tp.hc$species=="TAG",c("X0uM_8h_Orbi_0477")])/sum(alltreat.tp.hc$X0uM_8h_Orbi_0477)))

frac.TAG.24h.0uM.sd = sd(c(sum(alltreat.tp.hc[alltreat.tp.hc$species=="TAG",c("X0uM_24h_Orbi_0473")])/sum(alltreat.tp.hc$X0uM_24h_Orbi_0473),
                          sum(alltreat.tp.hc[alltreat.tp.hc$species=="TAG",c("X0uM_24h_Orbi_0468")])/sum(alltreat.tp.hc$X0uM_24h_Orbi_0468)))

paste(frac.TAG.4h.0uM,frac.TAG.8h.0uM,frac.TAG.24h.0uM)
paste(frac.TAG.8h.0uM.sd,frac.TAG.24h.0uM.sd)

# 150 uM H2O2

frac.TAG.4h.150uM = sum(alltreat.tp.hc[alltreat.tp.hc$species=="TAG",c("X150uM_4h_Orbi_0478")])/sum(alltreat.tp.hc$X150uM_4h_Orbi_0478)

frac.TAG.8h.150uM = mean(c(sum(alltreat.tp.hc[alltreat.tp.hc$species=="TAG",c("X150uM_8h_Orbi_0471")])/sum(alltreat.tp.hc$X150uM_8h_Orbi_0471),
                          sum(alltreat.tp.hc[alltreat.tp.hc$species=="TAG",c("X150uM_8h_Orbi_0467")])/sum(alltreat.tp.hc$X150uM_8h_Orbi_0467)))

frac.TAG.24h.150uM = mean(c(sum(alltreat.tp.hc[alltreat.tp.hc$species=="TAG",c("X150uM_24h_Orbi_0469")])/sum(alltreat.tp.hc$X150uM_24h_Orbi_0469),
                           sum(alltreat.tp.hc[alltreat.tp.hc$species=="TAG",c("X150uM_24h_Orbi_0466")])/sum(alltreat.tp.hc$X150uM_24h_Orbi_0466)))

frac.TAG.8h.150uM.sd = sd(c(sum(alltreat.tp.hc[alltreat.tp.hc$species=="TAG",c("X150uM_8h_Orbi_0471")])/sum(alltreat.tp.hc$X150uM_8h_Orbi_0471),
                           sum(alltreat.tp.hc[alltreat.tp.hc$species=="TAG",c("X150uM_8h_Orbi_0467")])/sum(alltreat.tp.hc$X150uM_8h_Orbi_0467)))

frac.TAG.24h.150uM.sd = sd(c(sum(alltreat.tp.hc[alltreat.tp.hc$species=="TAG",c("X150uM_24h_Orbi_0469")])/sum(alltreat.tp.hc$X150uM_24h_Orbi_0469),
                            sum(alltreat.tp.hc[alltreat.tp.hc$species=="TAG",c("X150uM_24h_Orbi_0466")])/sum(alltreat.tp.hc$X150uM_24h_Orbi_0466)))

paste(frac.TAG.4h.150uM,frac.TAG.8h.150uM,frac.TAG.24h.150uM)
paste(frac.TAG.8h.150uM.sd,frac.TAG.24h.150uM.sd)

################# create heatmap #############

# create a matrix for the heatmap input

# create subsets of compounds for the heatmap, if desired
# different methods of subsetting follow

hcPtassign.plotsub = hcPtassign # use all

# hcPtassign.plotsub = hcPtassign[apply(hcPtassign[,12:15], 1, function (x) sum(x>0.015)>=2),] # take those where the peak area > 0.015 in >= 2 samples
# hcPtassign.plotsub = hcPtassign[abs(hcPtassign$X24h_150uM_levavg-hcPtassign$X24h_0uM_levavg)>1.5,] # take only those where the absolute difference between the two treatments exceeds > 1.5

# take only the top 25% of values by peak area
# quant.75 = quantile(unlist(hcPtassign[,12:15]))[4]
# hcPtassign.plotsub = hcPtassign[apply(hcPtassign[,12:15], 1, function (x) sum(x>quant.75)>=2),] 

# put in order from most upregulated in 150 uM treatment to most downregulated
hcPtassign.plotsub = hcPtassign.plotsub[order(hcPtassign.plotsub$X24h_150uM_levavg, decreasing = TRUE),]

hmdata <- hcPtassign.plotsub[,c("plot.label","X24h_0uM_levavg","X24h_150uM_levavg")]
rnames <- hmdata[,1]  
hmdata<-hmdata[,-1]
#create matrix
hmdata_matrix <-data.matrix(hmdata) 
#log 2 transform data (i.e. fold change)
hmdata_matrix <- log2(hmdata_matrix) 
# assign row names 
rownames(hmdata_matrix) <- rnames      

## creates a  color palette 
my_palette <- colorRampPalette(c("grey", "skyblue", "darkorchid3", "red")) (n = 9999)

## (optional) define the color breaks manually for a "skewed" color transition
## I run min(hmdata_matrix) and max(hmdata_matrix) to determine where to set my color breaks 

col_breaks = c(seq(-20, -5,length=2500),  #NA values are effectively grey         
               seq(-5, 0,length=2500),   #values between -5 and 0 are sky blue
               seq(0, 0.5, length= 2500), #between 0 and 0.5 is dark orchid 3
               seq(0.5, 1, length =2500))  #between 0.5 and 1 are red

# creates a x by y inch pdf
pdf("Assaf H2O2 24hr Control and 150uM_heatmap.pdf",            
    width = 10,        
    height = 22,
    pointsize = 4)       

#creates heat map
heatmap.2(hmdata_matrix,
          #main = "Grazing Experiment WT Normalized 6hr", # heat map title
          notecol="black",      # change font color of cell labels to black
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(12,9),     # widens margins around plot
          col=my_palette,       # use on color palette defined earlier 
          breaks=col_breaks,    # enable color transition at specified limits
          dendrogram="none",    # specify which dendrograms to create
          key = TRUE,           # color key
          Colv=FALSE,           #column dendrogram
          Rowv = FALSE)         #row dendrogram            

          dev.off()               # close the pdf device

#Creating custom color keys using heatmap.2 can be problematic. This snippet of code will throw up the following error:
# Error in seq.default(min.raw, max.raw, by = min(diff(breaks)/4)) : 
# invalid (to - from)/by in seq(.) 
#A color key can be generated by running the code with the color break off. Then in Adobe Illustrator paste the color key into the correct heatmap. Lastly edit the numbers to reflect the color breaks. 

####################                                              


###################
##Similarity profile analysis to determine the significant clusters of molecules
#run simprof and print the significant clusters to determine the number of groups so that the leaf colors can be choosen 
          
#Simprof is interative and permutational. Therefore it may give slightly different resuls each time the function is run.  
#Resulting dendrogram may need to rotated such that the heat map and dendrogram can be displayed side by side. For the Assaf H2O2 dataset the heatmap input data was ordered by most upregulated in the 150uM treatment to most downregulated.  Alternatively the order of the heat map input data can be altered to reflect the order of the significant clusters. 
          
sig.clus <- simprof(hmdata_matrix, num.expected=100, num.simulated=99,
                    method.cluster="complete", method.distance="euclidean",
                    method.transform="identity", alpha=0.01,
                    sample.orientation="row", const=0,
                    silent=FALSE, increment=1,
                    undef.zero=TRUE, warn.braycurtis=TRUE)

print(sig.clus$significantclusters)

# reorder the hclust dendrogram object created by sig.clus to align with our desired ordering (most to least upregulated in the 150 uM treatment)
# note that this will only reorder the leaf objects; R will restructure the dendrogram to preserve the clustering

sig.clus.plot = sig.clus

sig.clus.plot$hclust = reorder(sig.clus.plot$hclust, 1:896, FUN = mean)

#From the above output determine how many leaf colors are necessary and assign that number of colors to the value "colors"

# using a repeating pattern

colpat <- c("red3", "sienna1", "yellow", "springgreen3", "blue", "cyan", "darkviolet", "red3")
colors=c(rep(colpat,trunc(sig.clus.plot$numgroups/length(colpat))),colpat[1:(sig.clus.plot$numgroups-trunc(sig.clus.plot$numgroups/length(colpat))*length(colpat))])

# with a rainbow ramp
 
clust.ramp <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan","#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000")) # create a color ramp for shading the points according to lipid class
# clust.ramp <- colorRampPalette(c("black", "white")) # create a color ramp for shading the points according to lipid class
 
colors_alt = clust.ramp(sig.clus.plot$numgroups)
colors_alt.rand = sample(colors_alt)

# with a b&w ramp

colors_alt.bw = gray.colors(n = sig.clus.plot$numgroups, start = 0, end = 1)
colors_alt.bw.rand = sample(colors_alt.bw)

#create a x by y inch pdf
pdf("sigclus_compounds Assaf H2O2 control and 150uM_hierclust.dendrogram.pdf",    # create PNG for the heat map        
    width = 22,        
    height = 10,            
    pointsize = 4)       
#plot the similarity profile using the leafcolors designated above

simprof.plot(sig.clus.plot, leafcolors= colors_alt.rand, plot=TRUE, fill=TRUE,
             leaflab="perpendicular", siglinetype=1)
dev.off()

# finally, we need to re-generate the heatmap with appropriately spaced labels (and color annotation) for each group created during hierarchial clustering (this will have to be pasted in manually in Illustrator alongside the "good" heatmap and the dendrogram)

sig.clus.plotgrouplabels = sig.clus.plot

# will try four different coloring conventions

RowSideColors_1 = as.data.frame(matrix(ncol = length(sig.clus.plotgrouplabels$hclust$labels), nrow = 1))
RowSideColors_2 = as.data.frame(matrix(ncol = length(sig.clus.plotgrouplabels$hclust$labels), nrow = 1))
RowSideColors_3 = as.data.frame(matrix(ncol = length(sig.clus.plotgrouplabels$hclust$labels), nrow = 1))
RowSideColors_4 = as.data.frame(matrix(ncol = length(sig.clus.plotgrouplabels$hclust$labels), nrow = 1))

# set a counter
count = 1

for (i in 1:length(sig.clus.plotgrouplabels$significantclusters)) {
  
  if (count==1) {
    
    bgw.col = "black"
    
  } else if (count==2) {
    
    bgw.col = "grey"
    
  } else if (count==3) {
    
    bgw.col = "white"
    
  }
  
  center.element = sig.clus.plotgrouplabels$significantclusters[[i]][round(median(1:length(sig.clus.plotgrouplabels$significantclusters[[i]])))]
  
  RowSideColors_1[sig.clus.plotgrouplabels$hclust$labels==center.element] = colors[i]
  RowSideColors_2[sig.clus.plotgrouplabels$hclust$labels==center.element] = bgw.col
  RowSideColors_3[sig.clus.plotgrouplabels$hclust$labels==center.element] = colors_alt.rand[i]
  RowSideColors_4[sig.clus.plotgrouplabels$hclust$labels==center.element] = colors_alt.bw.rand[i]

    sig.clus.plotgrouplabels$hclust$labels[sig.clus.plotgrouplabels$hclust$labels==center.element] = i
  
  other.elements = sig.clus.plotgrouplabels$significantclusters[[i]][sig.clus.plotgrouplabels$significantclusters[[i]]!=center.element]
  
  RowSideColors_1[sig.clus.plotgrouplabels$hclust$labels %in% other.elements] = colors[i]
  RowSideColors_2[sig.clus.plotgrouplabels$hclust$labels %in% other.elements] = bgw.col
  RowSideColors_3[sig.clus.plotgrouplabels$hclust$labels %in% other.elements] = colors_alt.rand[i]
  RowSideColors_4[sig.clus.plotgrouplabels$hclust$labels %in% other.elements] = colors_alt.bw.rand[i]
  
  sig.clus.plotgrouplabels$hclust$labels[sig.clus.plotgrouplabels$hclust$labels %in% other.elements] = i
  
  # advance counter
  
  if (count<3) {
    
  count = count + 1
  
  } else {
    
    count = 1
    
  }
  
}

hmdata_matrix.plotgrouplabels = hmdata_matrix
rownames(hmdata_matrix.plotgrouplabels) <- sig.clus.plotgrouplabels$hclust$labels      

# creates a x by y inch pdf
pdf("Assaf H2O2 24hr Control and 150uM_grouplabels.pdf",            
    width = 10,        
    height = 22,
    pointsize = 4)       

#creates heat map
heatmap.2(hmdata_matrix.plotgrouplabels,
          #main = "Grazing Experiment WT Normalized 6hr", # heat map title
          notecol="black",      # change font color of cell labels to black
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(12,9),     # widens margins around plot
          col=my_palette,       # use on color palette defined earlier 
          breaks=col_breaks,    # enable color transition at specified limits
          dendrogram="none",    # specify which dendrograms to create
          key = TRUE,           # color key
          Colv=FALSE,           #column dendrogram
          Rowv = FALSE,        #row dendrogram      
          RowSideColors = unlist(c(RowSideColors_3)))

dev.off()               # close the pdf device


################# create heatmap using only a subset of the data #############

# create a matrix for the heatmap input

# create subset of compounds for the heatmap

hcPtassign.plotsub.sub = hcPtassign[hcPtassign$species %in% c("MGDG"),] # specific species
hcPtassign.plotsub.sub = hcPtassign.plotsub.sub[apply(hcPtassign.plotsub.sub[,12:15], 1, function (x) sum(x>0.015)>=2),] # take those where the peak area > 0.015 in >= 2 samples
#hcPtassign.plotsub.sub = hcPtassign.plotsub.sub[abs(hcPtassign.plotsub.sub$X24h_150uM_levavg-hcPtassign.plotsub.sub$X24h_0uM_levavg)>1.75,] # take only those where the absolute difference between the two treatments exceeds some number

# put in order from most upregulated in 150 uM treatment to most downregulated
hcPtassign.plotsub.sub = hcPtassign.plotsub.sub[order(hcPtassign.plotsub.sub$X24h_150uM_levavg, decreasing = TRUE),]

hmdata <- hcPtassign.plotsub.sub[,c("plot.label","X24h_0uM_levavg","X24h_150uM_levavg")]
rnames <- hmdata[,1]  
hmdata<-hmdata[,-1]
#create matrix
hmdata_matrix <-data.matrix(hmdata) 
#log 2 transform data (i.e. fold change)
hmdata_matrix <- log2(hmdata_matrix) 
# assign row names 
rownames(hmdata_matrix) <- rnames      

## creates a  color palette 
my_palette <- colorRampPalette(c("grey", "skyblue", "darkorchid3", "red")) (n = 9999)

## (optional) define the color breaks manually for a "skewed" color transition
## I run min(hmdata_matrix) and max(hmdata_matrix) to determine where to set my color breaks 

col_breaks = c(seq(-20, -5,length=2500),  #NA values are effectively grey         
               seq(-5, 0,length=2500),   #values between -5 and 0 are sky blue
               seq(0, 0.5, length= 2500), #between 0 and 0.5 is dark orchid 3
               seq(0.5, 1, length =2500))  #between 0.5 and 1 are red

# creates a x by y inch pdf
pdf("Assaf H2O2 24hr Control and 150uM_heatmap_subset.pdf",            
    width = 10,        
    height = 22,
    pointsize = 4)       

#creates heat map
heatmap.2(hmdata_matrix,
          #main = "Grazing Experiment WT Normalized 6hr", # heat map title
          notecol="black",      # change font color of cell labels to black
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(12,9),     # widens margins around plot
          col=my_palette,       # use on color palette defined earlier 
          breaks=col_breaks,    # enable color transition at specified limits
          dendrogram="none",    # specify which dendrograms to create
          key = TRUE,           # color key
          Colv=FALSE,           #column dendrogram
          Rowv = FALSE)         #row dendrogram            

dev.off()               # close the pdf device

################# export elements of clusters, and some stats #############

# compute some stats
# calculate whether elements of a particular cluster differ significantly between treatments

sigdiff = matrix(data = NA, nrow = sig.clus$numgroups, ncol = 1)
  
for (i in 1:sig.clus$numgroups) {
  
  elem.thisclust = sig.clus$significantclusters[[i]] # get elements in this cluster
  
  #   t.res = t.test(
  #     hcPtassign.plotsub[hcPtassign.plotsub$plot.label %in% elem.thisclust,14:15],
  #          hcPtassign.plotsub[hcPtassign.plotsub$plot.label %in% elem.thisclust,12:13])
  
  levavgs150 = hcPtassign.plotsub[hcPtassign.plotsub$plot.label %in% elem.thisclust,c("X24h_150uM_levavg")]
  levavgs0 = hcPtassign.plotsub[hcPtassign.plotsub$plot.label %in% elem.thisclust,c("X24h_0uM_levavg")]
  
  if (length(unique(levavgs150))>1 & length(unique(levavgs0))>1) { # can perform t-test
    
   t.res = t.test(levavgs150,levavgs0)
  
  if (t.res$p.value<0.00001) {
    
    if (t.res$estimate[1]>t.res$estimate[2]) {
      
      sigdiff[i] = "More abundant"
      
    } else if (t.res$estimate[2]>t.res$estimate[1]) {
      
      sigdiff[i] = "Less abundant"
      
    }
    
    
  } else {
    
    sigdiff[i] = "No significant difference"
    
  }
   
  } else {
    
    if (mean(levavgs150)>mean(levavgs0)) {
      
      sigdiff[i] = "More abundant"
      
    } else if (mean(levavgs0)>mean(levavgs150)) {
      
      sigdiff[i] = "Less abundant"
      
    }
    
  }
  
}

# generate table for export

for (i in 1:sig.clus$numgroups) {
  
  elem.thisclust = sig.clus$significantclusters[[i]] # get elements in this cluster
  
  remainder.3cols = (length(elem.thisclust)/3-trunc(length(elem.thisclust)/3))*3 # get remainder for number of elements/3
  
  if (remainder.3cols==0) { # number of elements is divisible by 3
    
    elem.exp = matrix(elem.thisclust, ncol=3)
    
  } else {
    
    elem.exp.a = matrix(elem.thisclust[1:(3*floor(length(elem.thisclust)/3))], ncol=3)
    elem.exp.b = matrix(c(elem.thisclust[(3*floor(length(elem.thisclust)/3)+1):length(elem.thisclust)],rep("",(3-round(remainder.3cols)))), ncol = 3)
    
    elem.exp = rbind(elem.exp.a,elem.exp.b)
    
  }
  
  # add in ancillary information
  
  anc.info = matrix(c(i,sigdiff[i],length(elem.thisclust)), ncol = 3)
  anc.cols = rbind(anc.info, matrix(rep(c("","",""),round(nrow(elem.exp)-1)), ncol = 3))
  
  expdata.thisclust = cbind(anc.cols,elem.exp)
  
  # concatenate table rows for export
  
  if (i==1) {
    
    elem.table = expdata.thisclust
    
  } else {
    
    elem.table = rbind(elem.table,expdata.thisclust)
    
  }
  
}
  
write.csv(elem.table, file = "sig.clusters_forexport.csv")