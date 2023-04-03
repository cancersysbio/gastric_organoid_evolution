#############################################################################################################################
##                                                                                                                      
##  CNA BASED FISHPLOT
##  FOR FIGURES 2C AND S8A
##                                                                                                                      
##  Date: 27 DECEMBER 2021                                                                                                                   
##  
##  Author: Kasper Karlsson
##
##                                                                                                                      
############################################################################################################################


library(devtools)
install_github("chrisamiller/fishplot")
library(fishplot)
library(RColorBrewer)


###### FOR D3C1 - FIGURE 2C ######


### SET PATH
setwd(file.path("/Users/kasper.karlsson/_Stanford/Papers/Evolution_paper/GitHub_Reproductions/Main/Figure_2/PLOTS/Figure_2C_Fishplot"))

### SET SEQUENCING TIME POINTS
timepoints=c(20,40,60,80,100) # EVEN TIMEPOINTS TO HARMONIZE WITH OTHER PANEL FIGURES


### SET PARENTS
parents = c(0,1,1,3,3,4,5,7) 

### DEFINE CLONE MATRIX BASED ON TABLE S6
frac.table = matrix(
  c(88,  50, 1,   0,  0,   0,  0,  0,
    100, 45, 52,  1,  1,   0,  0,  0,
    100, 0,  98,  78, 17,  1,  1,  0, 
    100, 0,  100, 68, 30,  40, 30, 1,
    100, 0,  100, 0,  100, 0,  100, 97),
  ncol=length(timepoints))

### CREATE FISH OBJECT
fish = createFishObject(frac.table,parents,timepoints=timepoints)
fish = layoutClones(fish)
fishPlot(fish)

### PLOT
pdf("Fish_D3C1_12.5p_100q.pdf",width=10,height=4,pointsize=0.1)
fishPlot(fish,shape="spline",bg.type="solid", bg.col="#FFFFFF")
dev.off()


###### FOR D2C2 - FIGURE ED4D ######

### SET PATH
setwd(file.path("/Users/kasperkarlsson/_Stanford/Papers/Evolution_paper/Script_and_figures_final/Figure_scripts_commented/Figure_2/PLOTS/Figure_S8A_Fishplot"))

### SET SEQUENCING TIME POINTS
timepoints=c(20,40,60,80,100) # EVEN TIMEPOINTS TO HARMONIZE WITH OTHER PANEL FIGURES

### SET PARENTS
parents = c(0,1,2,2,4,5) 

### DEFINE CLONE MATRIX BASED ON TABLE S6
frac.table = matrix(
  c(97, 51, 37, 1, 0, 0,
    100, 85, 45, 40, 1, 0,
    100, 93, 48, 45, 40, 0,
    100, 100, 0, 100, 99, 1,
    100, 100, 0, 100, 100, 62),
  ncol=length(timepoints))  

### CREATE FISH OBJECT
fish = createFishObject(frac.table,parents,timepoints=timepoints)
fish = layoutClones(fish)
fishPlot(fish)

pdf("Fish_D2C2_12.5p_100q.pdf",width=10,height=4,pointsize=0.1)
fishPlot(fish,shape="spline",bg.type="solid", bg.col="#FFFFFF")
dev.off()

svg("Fish_D2C2_12.5p_100q.svg",width=10,height=4)
fishPlot(fish,shape="spline",bg.type="solid", bg.col="#FFFFFF")
dev.off()
