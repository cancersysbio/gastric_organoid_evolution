#############################################################################################################################
##                                                                                                                      
##  SCRIPT TO CREATE FISHPLOT FOR FIGURE 6C
##                                                                                                                      
##  
##  Author: Kasper Karlsson
##
##                                                                                                                      
############################################################################################################################


library(devtools)
#install_github("chrisamiller/fishplot")
library(fishplot)

############### D2C2R2 - FIGURE 6C #####################
mypalette <- c("#D3D3D3","#b3cde0","#6497b1","#005b96","#011f4b",                                                                # CNV CLONE COLOR
               "#8B8378","#FF00FF","#8B008B","#8B5F65","#BBFFFF","#00EE76","#FFFACD","#FFA500","#FF3030","#D8292C","#9D1E20")    # ECB SUBCLONE COLOR
               
setwd(file.path("/Users/kasper.karlsson/_Stanford/Papers/Evolution_paper/GitHub_Reproductions/Main/Figure_6/Figure_6C/OUTPUT"))


timepoints = c(5,37, 48, 62,129,143,157,173,189,200,231,245,258,273,315) # Time points for D2C2R2
parents = c(0,1,2,3,4,1,2,4,4,4,5,5,5,5,14,15) # Vector with clone parent child relationship


### Note the subclone fraction table is slightly adjusted from the output "Figure_6C.1_fishplot_barcode_frequencies.py" to make plotting possible, and the bottleneck has been introduced
frac.table = matrix(
  c(100,	8,	      5,	     1,	 1,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,
    100,	60,	      25,	     8,	 5,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,
    100,	80,	      50,	     35, 20,	1,	1,	1,	1,	1,	1,	1,	1,	1,	0,	0,
    1,	0.8888145,	0.777629,	0.777462,	0.44424,	0.1110185,	0.1110185,	0.1110185,	0.1110185,	0.1110185,	0.1110185,	0.1110185,	0.1110185,	0.1110185,	0,	0,
    100,	97.27553,	82.60923,	82.57748,	54.71608,	2.69269,	14.63454,	10.4557,	9.9103,	  7.46364,	15.28721,	7.18589,	29.21122,	3,	      0,	      0,
    100,	89.88054,	79.43404,	79.40781,	53.36757,	10.09321,	10.42027,	12.6365,	7.36361,	6.0139,	  17.33474,	3.75862,	25.40514,	6.84284,	0,	      0,
    100,	88.43869,	78.79017,	78.7656,	55.39123,	11.53671,	9.62395,	11.53671,	6.89849,	4.91459,	18.6237,	5.00212,	23.30504,	8.4358,	  1,	      0,
    100,	87.52756,	81.01132,	80.98728,	58.2093,	12.44837,	6.4922,	  9.55085,	6.98853,	6.21457,	19.76341,	5.7843,  	22.21258,	10.42496,	4.89389,	0,
    100,	85.91841,	79.2347,	79.21179,	57.06394,	14.05865,	6.66079,	7.8967,	  9.90319,	4.32505,	18.80539,	6.23822,	20.90031,	11.09711,	7.61181,	0,
    100,	85.7131,	79.91336,	79.88879,	61.45822,	14.26231,	5.77516,	3.52077,	10.54733,	4.3379,	  21.18973,	5.00151,	21.81255,	13.42986,	8.16196,	0,
    100,	87.59914,	76.69633,	76.67257,	62.01651,	12.37708,	10.87906,	2.20007,	12.43222,	0,	      18.67264,	4.66499,	21.34801,	17.30712,	17.30712,	1,
    100,	90.11587,	75.78505,	75.75969,	64.23094,	9.85875,	14.30546,	0,	      11.50339,	0,	      18.07197,	2.88469,	22.69374,	20.55518,	20.55518,	5,
    100,	94.2322,	77.01223,	76.98353,	67.22044,	5.73908,	17.19127,	0,	      9.7344, 	0,	      17.0185,	0,	      24.72026,	25.45299,	25.45299,	24.59156,
    100,	95.90659,	78.86753,	78.83538,	73.26469,	4.06122,	17.0069,  0,	      5.53854,	0,	      17.21568,	0,	      26.04907,	29.96778,	29.96778,	29.27128,
    100,	99.96015,	79.14958,	79.10977,	79.06996,	0,	      20.77076,	0,	      0,	      0,	      14.32008,	0,	      25.90872,	38.80134,	38.80134,	38.60543
  ),          
  ncol=length(timepoints)) 

fish = createFishObject(frac.table,parents,timepoints=timepoints,fix.missing.clones=TRUE) # Create fishplot object
fish2 <- setCol(fish, col = mypalette) # Change colors according to mypalette
fish2 = layoutClones(fish2) # Set fishplot layout

### PLOT FISH PLOT

pdf("Fish_D2C2R2_10_4.pdf",width=10,height=4,pointsize=0.1)
fishPlot(fish2,bg.type="solid", bg.col="#FFFFFF",shape="spline")
dev.off()

pdf("Fish_D2C2R2_10_4_lables.pdf",width=10,height=4,pointsize=0.1)
fishPlot(fish2,bg.type="solid", bg.col="#FFFFFF",shape="spline",vlines=timepoints,vlab=vlable)
dev.off()

svg("Fish_D2C2R2_10_4.svg",width=10,height=4)
fishPlot(fish2,bg.type="solid", bg.col="#FFFFFF",shape="spline")
dev.off()



############### D2C1R1 - FIGURE S19B #####################
mypalette <- c("#D3D3D3","#FF3030","#00EE76","#8B4726","#1E90FF","#FFA500","#D8292C") 

setwd(file.path("/Users/kasperkarlsson/_Stanford/Papers/Evolution_paper/Script_and_figures_final/Figure_scripts_commented/Figure_6/Fishplot/D2C1R1/PLOT"))


timepoints = c(55, 101,114,129,143,157,173,189,200,245,258,315,329,357,386,399) # Time points
parents = c(0,1,1,1,1,1,2) # Vector with clone parent child relationship


frac.table = matrix(
  c(100,	0,	0,	0,	0,	0,	0,
    100,	1,	1,	1,	1,	1,	0,
    20,	3.9988,	3.9988,	3.9988,	3.9988,	3.9988,	0,
    100,	17.69526,	11.27058,	18.52811,	17.2319,	35.23824,	0,
    100,	19.53726,	15.68535,	17.56796,	17.88412,	29.29476,	0,
    100,	23.33369,	17.63146,	17.00831,	15.10809,	26.88922,	0,
    100,	22.54324,	17.57472,	17.2057,	17.79519,	24.85334,	0,
    100,	22.4311,	17.49712,	18.14572,	18.82606,	23.07294,	0,
    100,	25.41621,	17.80404,	17.94318,	16.26777,	22.54067,	3,
    100,	29.13761,	16.68667,	21.53854,	19.64591,	12.96059,	16,
    100,	32.63168,	14.23255,	20.30253,	16.37551,	16.42423,	20,
    100,	55.80266,	7.10474,	25.99153,	0,	11.04482,	40,
    100,	53.57378,	7.74052,	13.83237,	0,	24.79931,	48,
    100,	71.74031,	0.99194,	0.99194,	0,	26.20366,	68,
    100,	90.79465,	0,	0,	0,	9.11421,	90.25782,
    100,	84.27842,	0,	0,	0,	15.63697,	84.15522),          
  ncol=length(timepoints))  # Create fish plot matrix from adjusted frequencies


fish = createFishObject(frac.table,parents,timepoints=timepoints,fix.missing.clones=TRUE) # Create fishplot object
fish2 <- setCol(fish, col = mypalette) # Change colors according to mypalette
fish2 = layoutClones(fish2) # Set fishplot layout

### PLOT FISH PLOT

pdf("Fish_D2C1R1_2_10_4.pdf",width=10,height=4,pointsize=0.1)
fishPlot(fish2,bg.type="solid", bg.col="#FFFFFF",shape="spline")
dev.off()

pdf("Fish_D2C1R1_2_10_4_lables.pdf",width=10,height=4,pointsize=0.1)
fishPlot(fish2,bg.type="solid", bg.col="#FFFFFF",shape="spline",vlines=timepoints,vlab=vlable)
dev.off()

svg("Fish_D2C1R1_2_10_4.svg",width=10,height=4)
fishPlot(fish2,bg.type="solid", bg.col="#FFFFFF",shape="spline")
dev.off()


############### D2C1R2 - FIGURE S20B #####################
mypalette <- c("#D3D3D3","#FF3030","#00EE76","#8B4726","#1E90FF","#FFA500","#D8292C","#378805") 

setwd(file.path("/Users/kasperkarlsson/_Stanford/Papers/Evolution_paper/Script_and_figures_final/Figure_scripts_commented/Figure_6/Fishplot/D2C1R2/PLOT"))


timepoints = c(55, 101,114,129,143,157,173,189,200,231,245,258,273,315,329,343,357,386,399)# Time points
parents = c(0,1,1,1,1,1,2,3) # Vector with clone parent child relationship


frac.table = matrix(
  c(100,0,0,0,0,0,0,0,
    100,1,1,1,1,1,0,0,
    20,3.9988,3.9988,3.9988,3.9988,3.9988,0,0,
    100,17.69059,11.27398,18.5277,17.22713,35.24466,0,0,
    100,19.53798,15.68783,17.56402,17.88486,29.29475,0,0,
    100,23.09456,17.77442,18.03854,14.45115,26.61229,0,0,
    100,21.89264,17.01992,18.46431,17.78915,24.80631,0,0,
    100,21.84885,17.18247,20.10835,18.3142,22.51931,0,0,
    100,24.60432,17.72551,18.32872,16.93968,22.3741,0,0,
    100,24.46386,22.36571,24.91947,15.57503,12.64691,3,3,
    100,22.89617,24.05371,24.60644,17.40074,11.01401,5,8,
    100,24.76634,24.01451,24.6543,11.29817,15.2372,12,15,
    100,24.03655,27.14446,26.06912,12.62742,10.09137,15,20,
    100,26.14673,39.68657,25.91759,0,8.20819,20,28,
    100,26.44022,41.20275,22.4114,0,9.90349,24,35,
    100,25.65143,41.9607,18.79111,0,13.55404,25,41.9607,
    100,29.11222,45.7099,15.24904,0,9.88234,28.40534,45.7099,
    100,30.25105,65.24711,0,0,4.43612,29.76472,65.24711,
    100,18.36148,69.7456,0,0,11.82283,18.36148,69.7456),          
  ncol=length(timepoints))    # Create fish plot matrix from adjusted frequencies

fish = createFishObject(frac.table,parents,timepoints=timepoints,fix.missing.clones=TRUE) # Create fishplot object
fish2 <- setCol(fish, col = mypalette) # Change colors according to mypalette
fish2 = layoutClones(fish2) # Set fishplot layout

### PLOT FISH PLOT

pdf("Fish_D2C1R2_2_10_4.pdf",width=10,height=4,pointsize=0.1)
fishPlot(fish2,bg.type="solid", bg.col="#FFFFFF",shape="spline")
dev.off()

pdf("Fish_D2C1R2_2_10_4_lables.pdf",width=10,height=4,pointsize=0.1)
fishPlot(fish2,bg.type="solid", bg.col="#FFFFFF",shape="spline",vlines=timepoints,vlab=timepoints)
dev.off()

svg("Fish_D2C1R2_2_10_4.svg",width=10,height=4)
fishPlot(fish2,bg.type="solid", bg.col="#FFFFFF",shape="spline")
dev.off()


############### D3C2R1 - FIGURE S21B #####################
mypalette <- c("#D3D3D3","#b3cde0","#6497b1",
               "#FF3030","#00EE76","#8B4726","#98F5FF","#1E90FF","#8B8682",
               "#8B8B00","#FFA500","#9D1E20","#D8292C","#ce8b54","#d2a56d") 

setwd(file.path("/Users/kasperkarlsson/_Stanford/Papers/Evolution_paper/Script_and_figures_final/Figure_scripts_commented/Figure_6/Fishplot/D3C2R1/PLOT"))


timepoints = c(55, 115, 128, 157, 175, 189, 200, 231, 245, 273, 287, 301, 315, 371,427) # Time points
parents = c(0,1,2,3,2,3,3,2,3,3,3,4,4,6,6) # Vector with clone parent child relationship


frac.table = matrix(
  c(100,	5,	1,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,
    100,	60,	20,	1,	1,	1,	1,	1,	1,	1,	1,	0,	0,	0,	0,
    20,	19.99624,	14.9969,	2.5,	2.50136,	2.50136,	2.4976,	2.4976,	2.4976,	2.4976,	2.4976,	1,	0,	0,	0,
    100,	99.97313,	83.61987,	9.50709,	6.62996,	14.98902,	11.64211,	9.72061,	12.11109,	9.10775,	26.23595,	4.74012,	4.74012,	0,	0,
    100,	99.97984,	75.65437,	16.93172,	9.97932,	13.94783,	9.03386,	14.34414,	8.71905,	8.38224,	18.61953,	10.02447,	6.88711,	1,	1,
    100,	99.98392,	83.76696,	21.98054,	8.42061,	26.38794,	7.58509,	7.79475,	9.15142,	5.71911,	12.9268,	12.23171,	9.73277,	7.37387,	6.27513,
    100,	99.98281,	83.00969,	26.25816,	11.55437,	25.27417,	7.68201,	5.41703,	8.73099,	2.80361,	12.24358,	14.8294,	11.41158,	6.90123,	5.73044,
    100,	99.95625,	77.88224,	43.12448,	22.06964,	0.0437,	  13.41825,	0,	0,	0,	21.25211,	43.08078,	0,	0,	0,
    100,	99.95523,	81.73795,	44.14719,	18.2128,	0.04473,	12.73892,	0,	0,	0,	24.76239,	44.10246,	0,	0,	0,
    100,	99.94054,	83.07872,	59.04181,	16.85589,	0.0594,	  0,	0,	0,	0,	23.91811,	58.98242,	0,	0,	0,
    100,	99.93672,	86.56334,	62.9031,	13.36705,	0.06322,	0,	0,	0,	0,	23.5338,	62.83988,	0,	0,	0,
    100,	99.92632,	89.69908,	73.34717,	10.21988,	0.0736,	  0,	0,	0,	0,	16.20471,	73.27356,	0,	0,	0,
    100,	99.92248,	99.83729,	77.16999,	0.07744,	0.07744,	0,	0,	0,	0,	22.51241,	77.09254,	0,	0,	0,
    100,	99.91955,	99.83114,	80.11588,	0.08037,	0.08037,	0,	0,	0,	0,	19.55451,	80.03551,	0,	0,	0,
    100,	99.90369,	99.79785,	95.96355,	0.09622,	0.09622,	0,	0,	0,	0,	3.64187,	95.86733,	0,	0,	0),          
  ncol=length(timepoints))    # Create fish plot matrix from adjusted frequencies

fish = createFishObject(frac.table,parents,timepoints=timepoints,fix.missing.clones=TRUE) # Create fishplot object
fish2 <- setCol(fish, col = mypalette) # Change colors according to mypalette
fish2 = layoutClones(fish2) # Set fishplot layout

### PLOT FISH PLOT

pdf("Fish_D2C1R2_10_4.pdf",width=10,height=4,pointsize=0.1)
fishPlot(fish2,bg.type="solid", bg.col="#FFFFFF",shape="spline")
dev.off()

pdf("Fish_D2C1R2_10_4_lables.pdf",width=10,height=4,pointsize=0.1)
fishPlot(fish2,bg.type="solid", bg.col="#FFFFFF",shape="spline",vlines=timepoints,vlab=timepoints)
dev.off()

svg("Fish_D2C1R2_10_4.svg",width=10,height=4)
fishPlot(fish2,bg.type="solid", bg.col="#FFFFFF",shape="spline")
dev.off()
