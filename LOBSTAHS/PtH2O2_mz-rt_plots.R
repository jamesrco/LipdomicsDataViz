# PtH2O2_mz-rt_plots.R

# Purpose: Analyze and generate some plots for the Graff Van Creveld et al. Pt H2O2 lipid data, after processing with LOBSTAHS pipeline

# Created 2/13/15 by J.R.C.
# Mod'fd 2/14/15, 5/6/15, 5/15/15 by J.R.C.
# Mod'fd 5/26/15 by J.R.C.
# Mod'fd 8/2/15 by J.R.C. to prepare plot for Lipids manuscript
# Mod'fd 12/30/15 by J.R.C. after retooling of LOBSTAHS pipeline (adaptation of "VML_lipidomics_stats_analysis_Assaf_H2O2_for_manuscript.R")

################# prep #############

library(TeachingDemos)
library(PtH2O2lipids)

# get data

ptdata = ptH2O2lipids$LOBSet
ptdata = getLOBpeaklist(ptdata)

setwd("/Users/jrcollins/Dropbox/Cruises & projects/High-Lat Lipid Peroxidation/Lipidomics methods paper/R analysis/New analysis with LOBSTAHS/")

################# QA/QC #############

# discard results with rt < 1 or > 26 mins

ptdata.QA <- ptdata[ptdata$peakgroup.rt>(1*60) & ptdata$peakgroup.rt<(26*60),]

################# normalize peak areas in screenedpeaks_bygroup.sum to DNPPE, if desired #############

ptdata.QA.norm <- ptdata.QA

DNPPE.data <- ptdata.QA.norm[ptdata.QA.norm$lipid_class=="DNPPE",] # get DNPPE data

# # if more than one DNPPE group, have to choose the right one and normalize to it
# #
# # ptdata.QA.norm <- ptdata.QA.norm[-171,] # get rid of bad one, if necessary 

# normalize all peak areas to the good DNPPE data

ptdata.QA.norm[,15:30] <- sweep(ptdata.QA.norm[,15:30], 2, as.numeric(DNPPE.data[,15:30]), "/")

# extract chl a data, will need this later to normalize

ptdata.QA.norm.chl.a <- ptdata.QA.norm[ptdata.QA.norm$compound_name =="Chl_a",]

# remove pigments, DNPPE

ptdata.QA.norm <- ptdata.QA.norm[grep("pigment", ptdata.QA.norm$lipid_class,invert=TRUE),]
ptdata.QA.norm <- ptdata.QA.norm[grep("DNPPE", ptdata.QA.norm$lipid_class,invert=TRUE),]

#   # ################# load in sample metadata for your experiment #############
#   
#   setwd("/Users/jrcollins/Dropbox/High-Lat Lipid Peroxidation/Lipidomics methods paper/R analysis/") # set working directory
#   sample_metadata_file <- "Assaf H2O2 culture experiments GBMF - metadata.csv" # location of sample metadata .csv format
#   sample_metadata_all <- read.csv(sample_metadata_file, skip = 0, header = T) # load file
#   sample_metadata_all <- sample_metadata_all[1:20,]
#   # colnames(sample_metadata_all) <- c("Sample.ID","Coll.DT","Extract.DT","Expt.ID","Treatment.ID","Vol.DNPPE.ul","Vol.samp_extracted.mL","CTD_station.ID","Depth.m","Notes","Old.style.vial","Samp_run.DT","Orbi_seq.ID","Instrument.method","HPLC_inject_vol.uL")
#   
#   # get number of unique treatments by timepoint
#   
#   treats_tp <- unique(sample_metadata_all[c("Treatment_uM_H2O2","Time_hrs")])
#   
#   # assign index numbers to treatments
#   treats_tp$Treat.ind <- seq(1,nrow(treats_tp))
#   # 
#   # # convert first two fields to characters
#   # treats_tp$Treatment.ID <- as.character(treats_tp$Treatment.ID)
#   # treats_tp$Coll.DT <- as.character(treats_tp$Coll.DT)

################# prep before plotting #############

# ppm.match.cut <- 2.5 # cutoff for ppm match, if desired
# 
# ptdata.QA.norm <- ptdata.QA.norm[abs(ptdata.QA.norm$mean_match_delta_ppm)<=ppm.match.cut,] # only ppm match <= ppm.match.cut

# coding by case assigned by the rules

# create a single code for "confidence" in the particular sample, based on the cases assigned by the rules

# set defaults

ptdata.QA.norm$confcode <- 0
ptdata.QA.norm$confcode.C3fandC3c <- 0
ptdata.QA.norm$confcode.C3r <- 0

ptdata.QA.norm[(rowSums(ptdata.QA.norm[,c("C1x","C4","C5","C3f","C3c","C6a","C6b","C2b")])==0),c("confcode")] <- 1 # select all which could be C1 or C2a (full satisfaction of rules) with no competing assignments

ptdata.QA.norm[((rowSums(ptdata.QA.norm[,c("C1x","C4","C5","C3f","C3c","C6a","C6b")])==0) & (ptdata.QA.norm[,c("C2b")]>=1)),c("confcode")] <- 2 # select all which could be 2b (partial satisfaction of rules) with no competing assignments

ptdata.QA.norm[((rowSums(ptdata.QA.norm[,c("C1x","C4","C5","C6a","C6b","C3c")])==0) & (ptdata.QA.norm[,c("C3f")]>=1)),c("confcode")] <- 3 # case 3fs that aren't case 3cs

ptdata.QA.norm[((rowSums(ptdata.QA.norm[,c("C1x","C4","C5","C6a","C6b","C3f")])==0) & (ptdata.QA.norm[,c("C3c")]>=1)),c("confcode")] <- 4 # case 3cs that aren't case 3fs

ptdata.QA.norm[((rowSums(ptdata.QA.norm[,c("C1x","C4","C5","C6a","C6b")])==0) & (ptdata.QA.norm[,c("C3c")]>=1) & (ptdata.QA.norm[,c("C3f")]>=1)),c("confcode.C3fandC3c")] <- 1 # features that are both case 3cs and case 3fs

ptdata.QA.norm[(ptdata.QA.norm[,c("C3r")]>=1),c("confcode.C3r")] <- 1 # case 3rs

################## assign colors and point symbols depending on properties #############

# set pch depending on oxylipin status
ptdata.QA.norm[ptdata.QA.norm$degree_oxidation==0,c("plot.pch")] <- 21
ptdata.QA.norm[ptdata.QA.norm$degree_oxidation==1,c("plot.pch")] <- 22
ptdata.QA.norm[ptdata.QA.norm$degree_oxidation==2,c("plot.pch")] <- 23
ptdata.QA.norm[ptdata.QA.norm$degree_oxidation==3,c("plot.pch")] <- 24
ptdata.QA.norm[ptdata.QA.norm$degree_oxidation==4,c("plot.pch")] <- 25

# basic prep for symbol fill color and point size

lipid.classes.plot <- c("DGCC","MGDG","TAG","PC","PE","SQDG","DGTS_DGTA","PG","DGDG") # lipid classes to plot

# set default distances for label proximity exclusion

distcuts = c(2.16, # DGCC
             2.13, # MGDG
             2.05, # TAG
             2.12, # PC
             2.12, # PE
             2.1, # SQDG
             2.15, # DGTS_DGTA
             2.1, # PG
             2.08) # DGDG

color_ramp.lipid.classes <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan","#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000")) # create a color ramp for shading the points according to lipid class

ptdata.QA.norm$plot.cex <- 0.85 # pch cex set to default uniform size, main plots
ptdata.QA.norm$plot.individual.lipid.classes.cex <- 1.5 # pch cex set to default uniform size, supplement plots

# basic prep for symbol border thickness and color

ptdata.QA.norm$plot.bord.col <- c("grey") # border color set to default uniform color, main plots
ptdata.QA.norm$plot.individual.lipid.classes.bord.col <- c("grey")# border color set to default uniform color, supplement plots

ptdata.QA.norm$plot.lwd <- 0.5 # border size set to default, main plots
ptdata.QA.norm$plot.individual.lipid.classes.lwd <- 0.5 # border size set to default, supplement plots

# specific adjustments for "big picture" plots (main text figure)

# set defaults to darkest color

ptdata.QA.norm$plot.lipid.classes.colors <- rev(color_ramp.lipid.classes(length(lipid.classes.plot)))[as.factor(as.character(ptdata.QA.norm$species))] # applies darkest shade of each color as a default

ptdata.QA.norm$plot.individual.lipid.classes.colors <- rev(color_ramp.lipid.classes(length(lipid.classes.plot)))[as.factor(as.character(ptdata.QA.norm$species))] # applies darkest shade of each color as a default

ptdata.QA.norm$plot.individual.lipid.classes.alt.colors <- rev(color_ramp.lipid.classes(length(lipid.classes.plot)))[as.factor(as.character(ptdata.QA.norm$species))] # applies darkest shade of each color as a default

# now, adjust to a lighter tone if it's a case 2b

for (i in 1:nrow(ptdata.QA.norm)) {
  
  this.darkestcolor <- ptdata.QA.norm$plot.lipid.classes.colors[i] # capture the current (dark) color for this lipid type
  
  new.plotramp <- colorRampPalette(c(this.darkestcolor,"white")) # generate a ramp to get the new colors below
  
  if (ptdata.QA.norm$confcode[i]==2) { # color should be slightly lighter
    
    lighter.plotcolor <- new.plotramp(20)[18]
    # generate a lighter shade for this point
    
    ptdata.QA.norm$plot.lipid.classes.colors[i] <- lighter.plotcolor # store the new color
    
  }
  
}

# specific adjustments to point color, border color, border thickness for individual species plots (to go in supplement)

# adjust tones for other than default cases

for (i in 1:nrow(ptdata.QA.norm)) {
  
  this.darkestcolor <- ptdata.QA.norm$plot.individual.lipid.classes.colors[i] # capture the current (dark) color for this lipid type
  
  new.plotramp <- colorRampPalette(c(this.darkestcolor,"white")) # generate a ramp to get the new colors below
  
  # generate colors for this point
  
  conf2.plotcolor <- new.plotramp(20)[6]
  conf3.plotcolor <- new.plotramp(20)[12]
  conf4.plotcolor <- new.plotramp(20)[18]
  
  if (ptdata.QA.norm$confcode[i]==2) { # color should be slightly lighter
    
    ptdata.QA.norm$plot.individual.lipid.classes.colors[i] <- conf2.plotcolor # store the new color
    
  } else if (ptdata.QA.norm$confcode[i]==3) { # color should be slightly lighter
    
    ptdata.QA.norm$plot.individual.lipid.classes.colors[i] <- conf3.plotcolor # store the new color
    
  } else if (ptdata.QA.norm$confcode[i]==4) { # color should be slightly lighter
    
    ptdata.QA.norm$plot.individual.lipid.classes.colors[i] <- conf4.plotcolor # store the new color
    
  }
  
  if (ptdata.QA.norm$confcode.C3fandC3c[i]==1) { # point shading should indicate both case C3f and C3c, need to adjust border color, fill color, and border thickness 
    
    ptdata.QA.norm$plot.individual.lipid.classes.colors[i] <- conf3.plotcolor # adjust fill color
    ptdata.QA.norm$plot.individual.lipid.classes.alt.colors[i] <- conf4.plotcolor  # adjust border color
    
  }
  
}

################## create final subsets before plotting #############

ptdata.QA.norm.all <- ptdata.QA.norm # create a subset including all data from this experiment

ptdata.QA.norm.sub.no_oxy <- ptdata.QA.norm[((ptdata.QA.norm$X0uM_24h_Orbi_0468>0) & (ptdata.QA.norm$X0uM_24h_Orbi_0473>0)),] # subset to only samples from 0 uM H2O2 treatment at 24 h

ptdata.QA.norm.sub.max_oxy <- ptdata.QA.norm[((ptdata.QA.norm$X30uM_24h_Orbi_0470>0) & (ptdata.QA.norm$X30uM_24h_Orbi_0479>0)),] # subset to only samples from 150 uM H2O2 treatment at 24 h ("max" oxidation)

################## universal plot settings #############

min.mz <- 575
max.mz <- 1250
min.rt <- 5
max.rt <- 22

# text offset for annotation on supplemental plots

text.offset.in <- 0.07 # specify only this value, absolute text offset in inches
text.offset.y <- 0 # ((max.mz-min.mz)/4.65)*text.offset.in
text.offset.x <- ((max.rt-min.rt)/6.93)*text.offset.in*1.5

# text size for annotation on supplemental plots

CtoDB_text.size <- 0.5 # in cex

################## single plots, all lipid classes #############

################# all lipid classes, just "best" matches #############

# ********* all compounds in screened & QA'd dataset

# # subset to just confcode = 1 or 2 (i.e., "best" or moderate matches... no competing assignments)
# ptdata.QA.norm.all.1 <- ptdata.QA.norm.all[ptdata.QA.norm.all$confcode==1 | ptdata.QA.norm.all$confcode==2,]

# # subset to just confcode = 1 (i.e., "best" matches... no competing assignments)
ptdata.QA.norm.all.1 <- ptdata.QA.norm.all[ptdata.QA.norm.all$confcode==1,] # alternate subset for figures in revised manuscript

graphics.off()
par(oma=c(0,0,0,0)) # set margins; large dataset seems to require this

# save image as pdf
pdf.filenm <- paste0("Rplot_all_classes_just_best_matches_all_compounds.pdf")
pdf(file=pdf.filenm,width=6.93,height=4.65,paper="special")

nf<-layout(matrix(c(1,1,1,1),2),c(1,1),c(1,1))

par(mar=c(5,5,1,1))
plot(ptdata.QA.norm.all.1$peakgroup.rt/60,ptdata.QA.norm.all.1$LOBdbase.mz,col=ptdata.QA.norm.all.1$plot.bord.col,bg=ptdata.QA.norm.all.1$plot.lipid.classes.colors,pch=ptdata.QA.norm.all.1$plot.pch,cex=ptdata.QA.norm.all.1$plot.cex,lwd=ptdata.QA.norm.all.1$plot.lwd,ylab=expression(italic(m/z)),xlab="Corrected retention time (min)",xlim=c(min.rt,max.rt),ylim=c(min.mz,max.mz),xaxs="i") # scatterplot, m/z vs r/t

# # turning this off, for now
# text(ptdata.QA.norm.all.1$peakgroup.rt/60, ptdata.QA.norm.all.1$LOBdbase.mz,paste0(ptdata.QA.norm.all.1$FA_total_no_C,":",ptdata.QA.norm.all.1$FA_total_no_DB),cex=0.6) # scatterplot, C and DB info

# generate some labels for the legend

IPL.labels <- unique(ptdata.QA.norm.all.1[,c("species","confcode")])[1]
IPL.labels <- as.character(IPL.labels$species)

confcode.labels <- unique(ptdata.QA.norm.all.1[,c("species","confcode")])[2]
confcode.labels <- as.character(confcode.labels$confcode)

IPL.confcode.labels <- paste(IPL.labels, confcode.labels)

#    # add plot legend
#    legend("topleft",c("Unoxidized","+1O","+2O","+3O","+4O"," "," "," "," "," ",IPL.confcode.labels),pch=c(21,22,23,24,25,rep(15,5),rep(15,nrow(unique(ptdata.QA.norm.all.1[,c("species","confcode")])))),col=c(rep("black",5),rep("white",5),unique(ptdata.QA.norm.all.1$plot.lipid.classes.colors)),pt.lwd=0.5,pt.bg=c(rep("white",5)),y.intersp=1.5,x.intersp=1.5,cex=1,pt.cex=c(rep(1.5,10),rep(3,nrow(unique(ptdata.QA.norm.all.1[,c("species","confcode")])))),xjust=1,yjust=1,bty="n",ncol=3,bg="white")

#   # add some descriptive text, upper left
#   text(min.rt,max.mz,paste0("All positive mode assignments, best matches only, all compounds in screened & QA'd dataset"),adj=c(-.025,.6),cex=2)
#   
#   # add some descriptive text, lower left
#   text(min.rt,min.mz,paste0("Database match at 2 ppm\nOnly assignments meeting r/t window criteria\nOnly assignments with even no. fatty acid C atoms\n"),pos=4,cex=0.75)

dev.off()

# ********* "no oxy" case 

# # subset to just confcode = 1 or 2 (i.e., "best" or moderate matches... no competing assignments)
# ptdata.QA.norm.sub.no_oxy.1 <- ptdata.QA.norm.sub.no_oxy[ptdata.QA.norm.sub.no_oxy$confcode==1 | ptdata.QA.norm.sub.no_oxy$confcode==2,]

# subset to just confcode = 1 (i.e., "best" matches... no competing assignments)
ptdata.QA.norm.sub.no_oxy.1 <- ptdata.QA.norm.sub.no_oxy[ptdata.QA.norm.sub.no_oxy$confcode==1,]

graphics.off()
par(oma=c(0,0,0,0)) # set margins; large dataset seems to require this

# save image as pdf
pdf.filenm <- paste0("Rplot_all_classes_just_best_matches_nooxy_0_uM_H2O2_24h.pdf")
pdf(file=pdf.filenm,width=6.93,height=4.65,paper="special")

nf<-layout(matrix(c(1,1,1,1),2),c(1,1),c(1,1))

par(mar=c(5,5,1,1))
plot(ptdata.QA.norm.sub.no_oxy.1$peakgroup.rt/60,ptdata.QA.norm.sub.no_oxy.1$LOBdbase.mz,col=ptdata.QA.norm.sub.no_oxy.1$plot.bord.col,bg=ptdata.QA.norm.sub.no_oxy.1$plot.lipid.classes.colors,pch=ptdata.QA.norm.sub.no_oxy.1$plot.pch,cex=ptdata.QA.norm.sub.no_oxy.1$plot.cex,lwd=ptdata.QA.norm.sub.no_oxy.1$plot.lwd,ylab=expression(italic(m/z)),xlab="Corrected retention time (min)",xlim=c(min.rt,max.rt),ylim=c(min.mz,max.mz),xaxs="i") # scatterplot, m/z vs r/t

# # turning this off, for now
# text(ptdata.QA.norm.sub.no_oxy.1$peakgroup.rt/60,ptdata.QA.norm.sub.no_oxy.1$LOBdbase.mz,paste0(ptdata.QA.norm.sub.no_oxy.1$FA_total_no_C,":",ptdata.QA.norm.sub.no_oxy.1$FA_total_no_DB),cex=0.6) # scatterplot, C and DB info

# generate some labels for the legend

IPL.labels <- unique(ptdata.QA.norm.sub.no_oxy.1[,c("species","confcode")])[1]
IPL.labels <- as.character(IPL.labels$species)

confcode.labels <- unique(ptdata.QA.norm.sub.no_oxy.1[,c("species","confcode")])[2]
confcode.labels <- as.character(confcode.labels$confcode)

IPL.confcode.labels <- paste(IPL.labels, confcode.labels)

#   # add plot legend
#   legend("topleft",c("Unoxidized","+1O","+2O","+3O","+4O"," "," "," "," "," ",IPL.confcode.labels),pch=c(21,22,23,24,25,rep(15,5),rep(15,nrow(unique(ptdata.QA.norm.sub.no_oxy.1[,c("species","confcode")])))),col=c(rep("black",5),rep("white",5),unique(ptdata.QA.norm.sub.no_oxy.1$plot.lipid.classes.colors)),pt.lwd=0.5,pt.bg=c(rep("white",5)),y.intersp=1.5,x.intersp=1.5,cex=1,pt.cex=c(rep(1.5,10),rep(3,nrow(unique(ptdata.QA.norm.sub.no_oxy.1[,c("species","confcode")])))),xjust=1,yjust=1,bty="n",ncol=3,bg="white")
#   
#   # add some descriptive text, upper left
#   text(min.rt,max.mz,paste0("All positive mode assignments, best matches only, from no oxy case"),adj=c(-.025,.6),cex=2)
#   
#   # add some descriptive text, lower left
#   text(min.rt,min.mz,paste0("Database match at 2 ppm\nOnly assignments meeting r/t window criteria\nOnly assignments with even no. fatty acid C atoms\n"),pos=4,cex=0.75)

dev.off()

# ********* "max oxy" case 

# # subset to just confcode = 1 or 2 (i.e., "best" or moderate matches... no competing assignments)
# ptdata.QA.norm.sub.max_oxy.1 <- ptdata.QA.norm.sub.max_oxy[ptdata.QA.norm.sub.max_oxy$confcode==1 | ptdata.QA.norm.sub.max_oxy$confcode==2,]

# subset to just confcode = 1 (i.e., "best" matches... no competing assignments)
ptdata.QA.norm.sub.max_oxy.1 <- ptdata.QA.norm.sub.max_oxy[ptdata.QA.norm.sub.max_oxy$confcode==1,]

graphics.off()
par(oma=c(0,0,0,0)) # set margins; large dataset seems to require this

# save image as pdf
pdf.filenm <- paste0("Rplot_all_classes_just_best_matches_maxoxy_150_uM_H2O2_24h.pdf")
pdf(file=pdf.filenm,width=6.93,height=4.65,paper="special")

nf<-layout(matrix(c(1,1,1,1),2),c(1,1),c(1,1))

par(mar=c(5,5,1,1))
plot(ptdata.QA.norm.sub.max_oxy.1$peakgroup.rt/60,ptdata.QA.norm.sub.max_oxy.1$LOBdbase.mz,col=ptdata.QA.norm.sub.max_oxy.1$plot.bord.col,bg=ptdata.QA.norm.sub.max_oxy.1$plot.lipid.classes.colors,pch=ptdata.QA.norm.sub.max_oxy.1$plot.pch,cex=ptdata.QA.norm.sub.max_oxy.1$plot.cex,lwd=ptdata.QA.norm.sub.max_oxy.1$plot.lwd,ylab=expression(italic(m/z)),xlab="Corrected retention time (min)",xlim=c(min.rt,max.rt),ylim=c(min.mz,max.mz),xaxs="i") # scatterplot, m/z vs r/t

# # turning this off for now
# text(ptdata.QA.norm.sub.max_oxy.1$peakgroup.rt/60,ptdata.QA.norm.sub.max_oxy.1$LOBdbase.mz,paste0(ptdata.QA.norm.sub.max_oxy.1$FA_total_no_C,":",ptdata.QA.norm.sub.max_oxy.1$FA_total_no_DB),cex=0.6) # scatterplot, C and DB info

#   # add plot legend
#   legend("topleft",c("Unoxidized","+1O","+2O","+3O","+4O"," "," "," "," "," ",IPL.confcode.labels),pch=c(21,22,23,24,25,rep(15,5),rep(15,nrow(unique(ptdata.QA.norm.sub.no_oxy.1[,c("species","confcode")])))),col=c(rep("black",5),rep("white",5),unique(ptdata.QA.norm.sub.no_oxy.1$plot.lipid.classes.colors)),pt.lwd=0.5,pt.bg=c(rep("white",5)),y.intersp=1.5,x.intersp=1.5,cex=1,pt.cex=c(rep(1.5,10),rep(3,nrow(unique(ptdata.QA.norm.sub.no_oxy.1[,c("species","confcode")])))),xjust=1,yjust=1,bty="n",ncol=3,bg="white")
#   
#   # add some descriptive text, upper left
#   text(min.rt,max.mz,paste0("All positive mode assignments, best matches only, from max oxy case"),adj=c(-.025,.6),cex=2)
#   
#   # add some descriptive text, lower left
#   text(min.rt,min.mz,paste0("Database match at 2 ppm\nOnly assignments meeting r/t window criteria\nOnly assignments with even no. fatty acid C atoms\n"),pos=4,cex=0.75)

dev.off()

################# use the subsets from the no and max oxy conditions to look at the data in a different way #############

# venn diagrams

library("VennDiagram")

# generate stats for individual diagrams, calculate relative size for figures

All.no = nrow(ptdata.QA.norm.sub.no_oxy.1)
All.max = nrow(ptdata.QA.norm.sub.max_oxy.1)
All.union = length(intersect(ptdata.QA.norm.sub.no_oxy.1$compound_name,ptdata.QA.norm.sub.max_oxy.1$compound_name))
All.scalefac = 1

DGCC.no = nrow(ptdata.QA.norm.sub.no_oxy.1[ptdata.QA.norm.sub.no_oxy.1$species=="DGCC",])
DGCC.max = nrow(ptdata.QA.norm.sub.max_oxy.1[ptdata.QA.norm.sub.max_oxy.1$species=="DGCC",])
DGCC.union = length(intersect(ptdata.QA.norm.sub.no_oxy.1[ptdata.QA.norm.sub.no_oxy.1$species=="DGCC",c("compound_name")],ptdata.QA.norm.sub.max_oxy.1[ptdata.QA.norm.sub.max_oxy.1$species=="DGCC",c("compound_name")]))
DGCC.scalefac = DGCC.no/All.no

DGDG.no = nrow(ptdata.QA.norm.sub.no_oxy.1[ptdata.QA.norm.sub.no_oxy.1$species=="DGDG",])
DGDG.max = nrow(ptdata.QA.norm.sub.max_oxy.1[ptdata.QA.norm.sub.max_oxy.1$species=="DGDG",])
DGDG.union = length(intersect(ptdata.QA.norm.sub.no_oxy.1[ptdata.QA.norm.sub.no_oxy.1$species=="DGDG",c("compound_name")],ptdata.QA.norm.sub.max_oxy.1[ptdata.QA.norm.sub.max_oxy.1$species=="DGDG",c("compound_name")]))
DGDG.scalefac = DGDG.no/All.no

DGTS.no = nrow(ptdata.QA.norm.sub.no_oxy.1[ptdata.QA.norm.sub.no_oxy.1$species=="DGTS_DGTA",])
DGTS.max = nrow(ptdata.QA.norm.sub.max_oxy.1[ptdata.QA.norm.sub.max_oxy.1$species=="DGTS_DGTA",])
DGTS.union = length(intersect(ptdata.QA.norm.sub.no_oxy.1[ptdata.QA.norm.sub.no_oxy.1$species=="DGTS_DGTA",c("compound_name")],ptdata.QA.norm.sub.max_oxy.1[ptdata.QA.norm.sub.max_oxy.1$species=="DGTS_DGTA",c("compound_name")]))
DGTS.scalefac = DGTS.no/All.no

MGDG.no = nrow(ptdata.QA.norm.sub.no_oxy.1[ptdata.QA.norm.sub.no_oxy.1$species=="MGDG",])
MGDG.max = nrow(ptdata.QA.norm.sub.max_oxy.1[ptdata.QA.norm.sub.max_oxy.1$species=="MGDG",])
MGDG.union = length(intersect(ptdata.QA.norm.sub.no_oxy.1[ptdata.QA.norm.sub.no_oxy.1$species=="MGDG",c("compound_name")],ptdata.QA.norm.sub.max_oxy.1[ptdata.QA.norm.sub.max_oxy.1$species=="MGDG",c("compound_name")]))
MGDG.scalefac = MGDG.no/All.no

PC.no = nrow(ptdata.QA.norm.sub.no_oxy.1[ptdata.QA.norm.sub.no_oxy.1$species=="PC",])
PC.max = nrow(ptdata.QA.norm.sub.max_oxy.1[ptdata.QA.norm.sub.max_oxy.1$species=="PC",])
PC.union = length(intersect(ptdata.QA.norm.sub.no_oxy.1[ptdata.QA.norm.sub.no_oxy.1$species=="PC",c("compound_name")],ptdata.QA.norm.sub.max_oxy.1[ptdata.QA.norm.sub.max_oxy.1$species=="PC",c("compound_name")]))
PC.scalefac = PC.no/All.no

PE.no = nrow(ptdata.QA.norm.sub.no_oxy.1[ptdata.QA.norm.sub.no_oxy.1$species=="PE",])
PE.max = nrow(ptdata.QA.norm.sub.max_oxy.1[ptdata.QA.norm.sub.max_oxy.1$species=="PE",])
PE.union = length(intersect(ptdata.QA.norm.sub.no_oxy.1[ptdata.QA.norm.sub.no_oxy.1$species=="PE",c("compound_name")],ptdata.QA.norm.sub.max_oxy.1[ptdata.QA.norm.sub.max_oxy.1$species=="PE",c("compound_name")]))
PE.scalefac = PE.no/All.no

PG.no = nrow(ptdata.QA.norm.sub.no_oxy.1[ptdata.QA.norm.sub.no_oxy.1$species=="PG",])
PG.max = nrow(ptdata.QA.norm.sub.max_oxy.1[ptdata.QA.norm.sub.max_oxy.1$species=="PG",])
PG.union = length(intersect(ptdata.QA.norm.sub.no_oxy.1[ptdata.QA.norm.sub.no_oxy.1$species=="PG",c("compound_name")],ptdata.QA.norm.sub.max_oxy.1[ptdata.QA.norm.sub.max_oxy.1$species=="PG",c("compound_name")]))
PG.scalefac = PG.no/All.no

SQDG.no = nrow(ptdata.QA.norm.sub.no_oxy.1[ptdata.QA.norm.sub.no_oxy.1$species=="SQDG",])
SQDG.max = nrow(ptdata.QA.norm.sub.max_oxy.1[ptdata.QA.norm.sub.max_oxy.1$species=="SQDG",])
SQDG.union = length(intersect(ptdata.QA.norm.sub.no_oxy.1[ptdata.QA.norm.sub.no_oxy.1$species=="SQDG",c("compound_name")],ptdata.QA.norm.sub.max_oxy.1[ptdata.QA.norm.sub.max_oxy.1$species=="SQDG",c("compound_name")]))
SQDG.scalefac = SQDG.no/All.no

TAG.no = nrow(ptdata.QA.norm.sub.no_oxy.1[ptdata.QA.norm.sub.no_oxy.1$species=="TAG",])
TAG.max = nrow(ptdata.QA.norm.sub.max_oxy.1[ptdata.QA.norm.sub.max_oxy.1$species=="TAG",])
TAG.union = length(intersect(ptdata.QA.norm.sub.no_oxy.1[ptdata.QA.norm.sub.no_oxy.1$species=="TAG",c("compound_name")],ptdata.QA.norm.sub.max_oxy.1[ptdata.QA.norm.sub.max_oxy.1$species=="TAG",c("compound_name")]))
TAG.scalefac = TAG.no/All.no

# generate diagrams

# all

graphics.off()
par(oma=c(0,0,0,0)) # set margins; large dataset seems to require this

# save image as pdf
pdf.filenm <- paste0("Venn_all.pdf")
pdf(file=pdf.filenm,width=6*All.scalefac,height=6*All.scalefac,paper="special")

grid.newpage()
draw.pairwise.venn(All.no, All.max, All.union, category = c("0 uM H2O2", "150 uM H2O2"),
                   lwd = rep(1,2), lty = rep(1,2), col = rep("grey",2),
                   fill = c("#FFFF00", "#f9d234"),
                   alpha = rep(0.5, 2), cat.pos = c(0,0),
                   cat.dist = rep(0.025, 2))

dev.off()

# DGCC

graphics.off()
par(oma=c(0,0,0,0)) # set margins; large dataset seems to require this

# save image as pdf
pdf.filenm <- paste0("Venn_DGCC.pdf")
pdf(file=pdf.filenm,width=6*DGCC.scalefac,height=6*DGCC.scalefac,paper="special")

grid.newpage()
draw.pairwise.venn(DGCC.no, DGCC.max, DGCC.union, category = c("0 uM H2O2", "150 uM H2O2"),
                   lwd = rep(1,2), lty = rep(1,2), col = rep("grey",2),
                   fill = c("#FFFF00", "#f9d234"),
                   alpha = rep(0.5, 2), cat.pos = c(0,0),
                   cat.dist = rep(0.025, 2))

dev.off()

# DGDG

graphics.off()
par(oma=c(0,0,0,0)) # set margins; large dataset seems to require this

# save image as pdf
pdf.filenm <- paste0("Venn_DGDG.pdf")
pdf(file=pdf.filenm,width=6*DGDG.scalefac,height=6*DGDG.scalefac,paper="special")

grid.newpage()
draw.pairwise.venn(DGDG.no, DGDG.max, DGDG.union, category = c("0 uM H2O2", "150 uM H2O2"),
                   lwd = rep(1,2), lty = rep(1,2), col = rep("grey",2),
                   fill = c("#FFFF00", "#f9d234"),
                   alpha = rep(0.5, 2), cat.pos = c(0,0),
                   cat.dist = rep(0.025, 2))

dev.off()

# DGTS_A

graphics.off()
par(oma=c(0,0,0,0)) # set margins; large dataset seems to require this

# save image as pdf
pdf.filenm <- paste0("Venn_DGTS.pdf")
pdf(file=pdf.filenm,width=6*DGTS.scalefac,height=6*DGTS.scalefac,paper="special")

grid.newpage()
draw.pairwise.venn(DGTS.no, DGTS.max, DGTS.union, category = c("0 uM H2O2", "150 uM H2O2"),
                   lwd = rep(1,2), lty = rep(1,2), col = rep("grey",2),
                   fill = c("#FFFF00", "#f9d234"),
                   alpha = rep(0.5, 2), cat.pos = c(0,0),
                   cat.dist = rep(0.025, 2))

dev.off()

# MGDG

graphics.off()
par(oma=c(0,0,0,0)) # set margins; large dataset seems to require this

# save image as pdf
pdf.filenm <- paste0("Venn_MGDG.pdf")
pdf(file=pdf.filenm,width=6*MGDG.scalefac,height=6*MGDG.scalefac,paper="special")

grid.newpage()
draw.pairwise.venn(MGDG.no, MGDG.max, MGDG.union, category = c("0 uM H2O2", "150 uM H2O2"),
                   lwd = rep(1,2), lty = rep(1,2), col = rep("grey",2),
                   fill = c("#FFFF00", "#f9d234"),
                   alpha = rep(0.5, 2), cat.pos = c(0,0),
                   cat.dist = rep(0.025, 2))

dev.off()

# PC

graphics.off()
par(oma=c(0,0,0,0)) # set margins; large dataset seems to require this

# save image as pdf
pdf.filenm <- paste0("Venn_PC.pdf")
pdf(file=pdf.filenm,width=6*PC.scalefac,height=6*PC.scalefac,paper="special")

grid.newpage()
draw.pairwise.venn(PC.no, PC.max, PC.union, category = c("0 uM H2O2", "150 uM H2O2"),
                   lwd = rep(1,2), lty = rep(1,2), col = rep("grey",2),
                   fill = c("#FFFF00", "#f9d234"),
                   alpha = rep(0.5, 2), cat.pos = c(0,0),
                   cat.dist = rep(0.025, 2))

dev.off()

# PE

graphics.off()
par(oma=c(0,0,0,0)) # set margins; large dataset seems to require this

# save image as pdf
pdf.filenm <- paste0("Venn_PE.pdf")
pdf(file=pdf.filenm,width=6*PE.scalefac,height=6*PE.scalefac,paper="special")

grid.newpage()
draw.pairwise.venn(PE.no, PE.max, PE.union, category = c("0 uM H2O2", "150 uM H2O2"),
                   lwd = rep(1,2), lty = rep(1,2), col = rep("grey",2),
                   fill = c("#FFFF00", "#f9d234"),
                   alpha = rep(0.5, 2), cat.pos = c(0,0),
                   cat.dist = rep(0.025, 2))

dev.off()

# PG

graphics.off()
par(oma=c(0,0,0,0)) # set margins; large dataset seems to require this

# save image as pdf
pdf.filenm <- paste0("Venn_PG.pdf")
pdf(file=pdf.filenm,width=6*PG.scalefac,height=6*PG.scalefac,paper="special")

grid.newpage()
draw.pairwise.venn(PG.no, PG.max, PG.union, category = c("0 uM H2O2", "150 uM H2O2"),
                   lwd = rep(1,2), lty = rep(1,2), col = rep("grey",2),
                   fill = c("#FFFF00", "#f9d234"),
                   alpha = rep(0.5, 2), cat.pos = c(0,0),
                   cat.dist = rep(0.025, 2))

dev.off()

# SQDG

graphics.off()
par(oma=c(0,0,0,0)) # set margins; large dataset seems to require this

# save image as pdf
pdf.filenm <- paste0("Venn_SQDG.pdf")
pdf(file=pdf.filenm,width=6*SQDG.scalefac,height=6*SQDG.scalefac,paper="special")

grid.newpage()
draw.pairwise.venn(SQDG.no, SQDG.max, SQDG.union, category = c("0 uM H2O2", "150 uM H2O2"),
                   lwd = rep(1,2), lty = rep(1,2), col = rep("grey",2),
                   fill = c("#FFFF00", "#f9d234"),
                   alpha = rep(0.5, 2), cat.pos = c(0,0),
                   cat.dist = rep(0.025, 2))

dev.off()

# TAG

graphics.off()
par(oma=c(0,0,0,0)) # set margins; large dataset seems to require this

# save image as pdf
pdf.filenm <- paste0("Venn_TAG.pdf")
pdf(file=pdf.filenm,width=6*TAG.scalefac,height=6*TAG.scalefac,paper="special")

grid.newpage()
draw.pairwise.venn(TAG.no, TAG.max, TAG.union, category = c("0 uM H2O2", "150 uM H2O2"),
                   lwd = rep(1,2), lty = rep(1,2), col = rep("grey",2),
                   fill = c("#FFFF00", "#f9d234"),
                   alpha = rep(0.5, 2), cat.pos = c(0,0),
                   cat.dist = rep(0.025, 2))

dev.off()


################# generate plots of mz vs rt by lipid type, with color coding for strength of match #############

for (i in 1:length(lipid.classes.plot)) {

  graphics.off()
  
  this.lipid.class <- lipid.classes.plot[i]
  
  ################# more subsetting of data #############
  
  IPL.subtype <- this.lipid.class # IPL subtype
  
  # create subsets for only one IPL type
  
  ptdata.QA.norm.sub.no_oxy.conf123.IPL_subtype <- ptdata.QA.norm.sub.no_oxy[ptdata.QA.norm.sub.no_oxy$species==IPL.subtype & ptdata.QA.norm.sub.no_oxy$confcode.C3fandC3c!=1,]
  
  ptdata.QA.norm.sub.no_oxy.conf123.regio.IPL_subtype <- ptdata.QA.norm.sub.no_oxy[ptdata.QA.norm.sub.no_oxy$species==IPL.subtype & ptdata.QA.norm.sub.no_oxy$confcode.C3r==1,]
  
  ptdata.QA.norm.sub.max_oxy.conf123.IPL_subtype <- ptdata.QA.norm.sub.max_oxy[ptdata.QA.norm.sub.max_oxy$species==IPL.subtype & ptdata.QA.norm.sub.max_oxy$confcode.C3fandC3c!=1,]
  
  ptdata.QA.norm.sub.max_oxy.conf123.regio.IPL_subtype <- ptdata.QA.norm.sub.max_oxy[ptdata.QA.norm.sub.max_oxy$species==IPL.subtype & ptdata.QA.norm.sub.max_oxy$confcode.C3r==1,]
  
  ptdata.QA.norm.sub.no_oxy.conf4.IPL_subtype <- ptdata.QA.norm.sub.no_oxy[ptdata.QA.norm.sub.no_oxy$species==IPL.subtype & ptdata.QA.norm.sub.no_oxy$confcode.C3fandC3c==1,]
  
  ptdata.QA.norm.sub.no_oxy.conf4.regio.IPL_subtype <- ptdata.QA.norm.sub.no_oxy[ptdata.QA.norm.sub.no_oxy$species==IPL.subtype & ptdata.QA.norm.sub.no_oxy$confcode.C3r==1,]
  
  ptdata.QA.norm.sub.max_oxy.conf4.IPL_subtype <- ptdata.QA.norm.sub.max_oxy[ptdata.QA.norm.sub.max_oxy$species==IPL.subtype & ptdata.QA.norm.sub.max_oxy$confcode.C3fandC3c==1,]
  
  ptdata.QA.norm.sub.max_oxy.conf4.regio.IPL_subtype <- ptdata.QA.norm.sub.max_oxy[ptdata.QA.norm.sub.max_oxy$species==IPL.subtype & ptdata.QA.norm.sub.max_oxy$confcode.C3r==1,]
  
  ptdata.QA.norm.all.IPL_subtype <- ptdata.QA.norm.all[ptdata.QA.norm.all$species==IPL.subtype,]
  
  # ********* "no oxy" case
  
  graphics.off()
  par(oma=c(0,0,0,0)) # set margins; large dataset seems to require this
  
  # save image as pdf
  pdf.filenm <- paste0("Rplot_",IPL.subtype,"_0_uM_H2O2_24h.pdf")
  pdf(file=pdf.filenm,width=6.93,height=4.65,paper="special")
  
  nf<-layout(matrix(c(1,1,1,1),2),c(1,1),c(1,1))
  
  par(mar=c(5,5,1,1))
  
  # first, plot points with confidence 1-3
  
  plot(ptdata.QA.norm.sub.no_oxy.conf123.IPL_subtype$peakgroup.rt/60,ptdata.QA.norm.sub.no_oxy.conf123.IPL_subtype$LOBdbase.mz,col=ptdata.QA.norm.sub.no_oxy.conf123.IPL_subtype$plot.individual.lipid.classes.bord.col,bg=ptdata.QA.norm.sub.no_oxy.conf123.IPL_subtype$plot.individual.lipid.classes.colors,pch=ptdata.QA.norm.sub.no_oxy.conf123.IPL_subtype$plot.pch,cex=ptdata.QA.norm.sub.no_oxy.conf123.IPL_subtype$plot.individual.lipid.classes.cex,lwd=ptdata.QA.norm.sub.no_oxy.conf123.IPL_subtype$plot.individual.lipid.classes.lwd,ylab=expression(italic(m/z)),xlab="Corrected retention time (min)",xlim=c(min.rt,max.rt),ylim=c(min.mz,max.mz),xaxs="i") # scatterplot, m/z vs r/t
  
  # overlay dots for regioisomers
  points(ptdata.QA.norm.sub.no_oxy.conf123.regio.IPL_subtype$peakgroup.rt/60,ptdata.QA.norm.sub.no_oxy.conf123.regio.IPL_subtype$LOBdbase.mz,bg="black",pch=20,cex=0.3)
  
  # overlay points with confidence code 4, if there are any
  
  if (nrow(ptdata.QA.norm.sub.no_oxy.conf4.IPL_subtype)>0) {
    
    points(ptdata.QA.norm.sub.no_oxy.conf4.IPL_subtype$peakgroup.rt/60,ptdata.QA.norm.sub.no_oxy.conf4.IPL_subtype$LOBdbase.mz,col=ptdata.QA.norm.sub.no_oxy.conf4.IPL_subtype$plot.individual.lipid.classes.bord.col,bg=ptdata.QA.norm.sub.no_oxy.conf4.IPL_subtype$plot.individual.lipid.classes.colors,pch=ptdata.QA.norm.sub.no_oxy.conf4.IPL_subtype$plot.pch,cex=ptdata.QA.norm.sub.no_oxy.conf4.IPL_subtype$plot.individual.lipid.classes.cex,lwd=ptdata.QA.norm.sub.no_oxy.conf4.IPL_subtype$plot.individual.lipid.classes.lwd)
    
    # overlay internal points for confidence code 4
    
    points(ptdata.QA.norm.sub.no_oxy.conf4.IPL_subtype$peakgroup.rt/60,ptdata.QA.norm.sub.no_oxy.conf4.IPL_subtype$LOBdbase.mz,col=ptdata.QA.norm.sub.no_oxy.conf4.IPL_subtype$plot.individual.lipid.classes.bord.col,bg=ptdata.QA.norm.sub.no_oxy.conf4.IPL_subtype$plot.individual.lipid.classes.alt.colors,pch=ptdata.QA.norm.sub.no_oxy.conf4.IPL_subtype$plot.pch,cex=0.9,lwd=ptdata.QA.norm.sub.no_oxy.conf4.IPL_subtype$plot.individual.lipid.classes.lwd)
    
    # overlay dots for regioisomers
    points(ptdata.QA.norm.sub.no_oxy.conf4.regio.IPL_subtype$peakgroup.rt/60,ptdata.QA.norm.sub.no_oxy.conf4.regio.IPL_subtype$LOBdbase.mz,bg="black",pch=20,cex=0.3)
    
  }
  
  # overlay labels on features using some imperfect code to that tries to limit overlap 
  
  if (nrow(ptdata.QA.norm.sub.no_oxy.conf4.IPL_subtype)>0){
    
    forlabeling = rbind(ptdata.QA.norm.sub.no_oxy.conf4.IPL_subtype,ptdata.QA.norm.sub.no_oxy.conf123.IPL_subtype)
    
  } else {
    
    forlabeling = ptdata.QA.norm.sub.no_oxy.conf123.IPL_subtype
    
  } 
  
  distmat = as.matrix(dist(cbind(forlabeling$peakgroup.rt/60+text.offset.x, forlabeling$LOBdbase.mz-text.offset.y), diag=TRUE, upper=TRUE, method = "euclidean")) # get euclidean distances
  
  dists <- apply(distmat<distcuts[i], 2, sum)  # count up the number of neighbors within distcut distance of each point
  
  labeledpoints = forlabeling[dists==1,] # subset to only points that should be labeled 
  
  #      text(labeledpoints$peakgroup.rt/60+text.offset.x,labeledpoints$LOBdbase.mz-text.offset.y,paste0(labeledpoints$FA_total_no_C,":",labeledpoints$FA_total_no_DB),cex=CtoDB_text.size) # scatterplot, C and DB info
  
  shadowtext(labeledpoints$peakgroup.rt/60+text.offset.x,labeledpoints$LOBdbase.mz-text.offset.y,paste0(labeledpoints$FA_total_no_C,":",labeledpoints$FA_total_no_DB),cex=CtoDB_text.size,col="black",bg="white",r=0.1,theta = seq(pi/4, 2 * pi, length.out = 30)) # scatterplot, C and DB info
  
  # generate some labels for the legend
  
  IPL.labels <- unique(ptdata.QA.norm.sub.no_oxy.conf123.IPL_subtype[,c("species","confcode")])[1]
  IPL.labels <- as.character(IPL.labels$species)
  
  confcode.labels <- unique(ptdata.QA.norm.sub.no_oxy.conf123.IPL_subtype[,c("species","confcode")])[2]
  confcode.labels <- as.character(confcode.labels$confcode)
  
  IPL.confcode.labels <- paste(IPL.labels, confcode.labels)
  
    # add plot legend
    legend("topleft",c("Unoxidized","+1O","+2O","+3O","+4O"," "," "," "," "," ",IPL.confcode.labels),pch=c(21,22,23,24,25,rep(15,5),rep(15,nrow(unique(ptdata.QA.norm.sub.no_oxy.conf123.IPL_subtype[,c("species","confcode")])))),col=c(rep("black",5),rep("white",5),unique(ptdata.QA.norm.sub.no_oxy.conf123.IPL_subtype$plot.individual.lipid.classes.colors)),pt.lwd=0.5,pt.bg=c(rep("white",5)),y.intersp=1.5,x.intersp=1.5,cex=1,pt.cex=c(rep(1.5,10),rep(3,nrow(unique(ptdata.QA.norm.sub.no_oxy.conf123.IPL_subtype[,c("species","confcode")])))),xjust=1,yjust=1,bty="n",ncol=3,bg="white")
  
  #   # add some descriptive text, upper left
  #   text(min.rt,max.mz,paste0(IPL.subtype, ", all positive mode assignments, 0 uM H2O2 at 24h"),adj=c(-.025,.6),cex=2)
  #   
  #   # add some descriptive text, lower left
  #   text(min.rt,min.mz,paste0("Database match at 2 ppm\nOnly assignments meeting r/t window criteria\nOnly assignments with even no. fatty acid C atoms\n"),pos=4,cex=0.75)
  
  dev.off()
  
  # ********* "max oxy" case
  
  graphics.off()
  par(oma=c(0,0,0,0)) # set margins; large dataset seems to require this
  
  # save image as pdf
  pdf.filenm <- paste0("Rplot_",IPL.subtype,"_150_uM_H2O2_24h.pdf")
  pdf(file=pdf.filenm,width=6.93,height=4.65,paper="special")
  
  nf<-layout(matrix(c(1,1,1,1),2),c(1,1),c(1,1))
  
  par(mar=c(5,5,1,1))
  
  # first, plot points with confidence 1-3
  
  plot(ptdata.QA.norm.sub.max_oxy.conf123.IPL_subtype$peakgroup.rt/60,ptdata.QA.norm.sub.max_oxy.conf123.IPL_subtype$LOBdbase.mz,col=ptdata.QA.norm.sub.max_oxy.conf123.IPL_subtype$plot.individual.lipid.classes.bord.col,bg=ptdata.QA.norm.sub.max_oxy.conf123.IPL_subtype$plot.individual.lipid.classes.colors,pch=ptdata.QA.norm.sub.max_oxy.conf123.IPL_subtype$plot.pch,cex=ptdata.QA.norm.sub.max_oxy.conf123.IPL_subtype$plot.individual.lipid.classes.cex,lwd=ptdata.QA.norm.sub.max_oxy.conf123.IPL_subtype$plot.individual.lipid.classes.lwd,ylab=expression(italic(m/z)),xlab="Corrected retention time (min)",xlim=c(min.rt,max.rt),ylim=c(min.mz,max.mz),xaxs="i") # scatterplot, m/z vs r/t
  
  # overlay dots for regioisomers
  points(ptdata.QA.norm.sub.max_oxy.conf123.regio.IPL_subtype$peakgroup.rt/60,ptdata.QA.norm.sub.max_oxy.conf123.regio.IPL_subtype$LOBdbase.mz,bg="black",pch=20,cex=0.3)
  
  # overlay points with confidence code 4, if there ase any
  
  if (nrow(ptdata.QA.norm.sub.max_oxy.conf4.IPL_subtype)>0) {
    
    points(ptdata.QA.norm.sub.max_oxy.conf4.IPL_subtype$peakgroup.rt/60,ptdata.QA.norm.sub.max_oxy.conf4.IPL_subtype$LOBdbase.mz,col=ptdata.QA.norm.sub.max_oxy.conf4.IPL_subtype$plot.individual.lipid.classes.bord.col,bg=ptdata.QA.norm.sub.max_oxy.conf4.IPL_subtype$plot.individual.lipid.classes.colors,pch=ptdata.QA.norm.sub.max_oxy.conf4.IPL_subtype$plot.pch,cex=ptdata.QA.norm.sub.max_oxy.conf4.IPL_subtype$plot.individual.lipid.classes.cex,lwd=ptdata.QA.norm.sub.max_oxy.conf4.IPL_subtype$plot.individual.lipid.classes.lwd)
    
    # overlay internal points for confidence code 4
    
    points(ptdata.QA.norm.sub.max_oxy.conf4.IPL_subtype$peakgroup.rt/60,ptdata.QA.norm.sub.max_oxy.conf4.IPL_subtype$LOBdbase.mz,col=ptdata.QA.norm.sub.max_oxy.conf4.IPL_subtype$plot.individual.lipid.classes.bord.col,bg=ptdata.QA.norm.sub.max_oxy.conf4.IPL_subtype$plot.individual.lipid.classes.alt.colors,pch=ptdata.QA.norm.sub.max_oxy.conf4.IPL_subtype$plot.pch,cex=0.9,lwd=ptdata.QA.norm.sub.max_oxy.conf4.IPL_subtype$plot.individual.lipid.classes.lwd)
    
    # overlay dots for regioisomers
    points(ptdata.QA.norm.sub.max_oxy.conf4.regio.IPL_subtype$peakgroup.rt/60,ptdata.QA.norm.sub.max_oxy.conf4.regio.IPL_subtype$LOBdbase.mz,bg="black",pch=20,cex=0.3)
    
  }
  
  # overlay labels on features using some imperfect code to that tries to limit overlap 
  
  if (nrow(ptdata.QA.norm.sub.max_oxy.conf4.IPL_subtype)>0){
    
    forlabeling = rbind(ptdata.QA.norm.sub.max_oxy.conf4.IPL_subtype,ptdata.QA.norm.sub.max_oxy.conf123.IPL_subtype)
    
  } else {
    
    forlabeling = ptdata.QA.norm.sub.max_oxy.conf123.IPL_subtype
    
  } 
  
  distmat = as.matrix(dist(cbind(forlabeling$peakgroup.rt/60+text.offset.x, forlabeling$LOBdbase.mz-text.offset.y), diag=TRUE, upper=TRUE, method = "euclidean")) # get euclidean distances
  
  dists <- apply(distmat<distcuts[i], 2, sum)  # count up the number of neighbors within distcut distance of each point
  
  labeledpoints = forlabeling[dists==1,] # subset to only points that should be labeled 
  
  #      text(labeledpoints$peakgroup.rt/60+text.offset.x,labeledpoints$LOBdbase.mz-text.offset.y,paste0(labeledpoints$FA_total_no_C,":",labeledpoints$FA_total_no_DB),cex=CtoDB_text.size) # scatterplot, C and DB info
  
  shadowtext(labeledpoints$peakgroup.rt/60+text.offset.x,labeledpoints$LOBdbase.mz-text.offset.y,paste0(labeledpoints$FA_total_no_C,":",labeledpoints$FA_total_no_DB),cex=CtoDB_text.size,col="black",bg="white",r=0.1,theta = seq(pi/4, 2 * pi, length.out = 30)) # scatterplot, C and DB info
  
  # generate some labels for the legend
  
  IPL.labels <- unique(ptdata.QA.norm.sub.max_oxy.conf123.IPL_subtype[,c("species","confcode")])[1]
  IPL.labels <- as.character(IPL.labels$species)
  
  confcode.labels <- unique(ptdata.QA.norm.sub.max_oxy.conf123.IPL_subtype[,c("species","confcode")])[2]
  confcode.labels <- as.character(confcode.labels$confcode)
  
  IPL.confcode.labels <- paste(IPL.labels, confcode.labels)
  
  #   # add plot legend
  #   legend("topleft",c("Unoxidized","+1O","+2O","+3O","+4O"," "," "," "," "," ",IPL.confcode.labels),pch=c(21,22,23,24,25,rep(15,5),rep(15,nrow(unique(ptdata.QA.norm.sub.max_oxy.conf123.IPL_subtype[,c("species","confcode")])))),col=c(rep("black",5),rep("white",5),unique(ptdata.QA.norm.sub.max_oxy.conf123.IPL_subtype$plot.individual.lipid.classes.colors)),pt.lwd=0.5,pt.bg=c(rep("white",5)),y.intersp=1.5,x.intersp=1.5,cex=1,pt.cex=c(rep(1.5,10),rep(3,nrow(unique(ptdata.QA.norm.sub.max_oxy.conf123.IPL_subtype[,c("species","confcode")])))),xjust=1,yjust=1,bty="n",ncol=3,bg="white")
  
  #   # add some descriptive text, upper left
  #   text(min.rt,max.mz,paste0(IPL.subtype, ", all positive mode assignments, 150 uM H2O2 at 24h"),adj=c(-.025,.6),cex=2)
  #   
  #   # add some descriptive text, lower left
  #   text(min.rt,min.mz,paste0("Database match at 2 ppm\nOnly assignments meeting r/t window criteria\nOnly assignments with even no. fatty acid C atoms\n"),pos=4,cex=0.75)
  
  dev.off()
  
}




