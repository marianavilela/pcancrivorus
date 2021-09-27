######################################################################
######################################################################

# PLemes / jan2021

######################################################################
######################################################################
# R template for standard sdm 
# Using the sdm R package (Naimi & Araujo, 2016)
#install.packages("sdm")
#install.packages("Rtools")
devtools::install_github('babaknaimi/sdm') # to install the latest version of the sdm package from github
library(sdm)
#installAll()

library(dismo)
library(rgbif)
library(plyr)
library(usdm)
library(spThin)
library(rgdal)

setwd("") # Your directory

# ------- CLIMATE Past
setwd("C:/Users/Priscila/Documents/PLemes/Orientacao/Projetos PIBIC/Mariana/Climate")
dir()
# Check for name consistency
# LGM  ------------------------
setwd("LGM")
setwd("CCSM4")
bioL.cc <- stack(list.files(pattern = ".tif"))

setwd("..")
# LGM  ------------------------
setwd("MIROC-ESM")
bioL.mr <- stack(list.files(pattern = ".tif"))

setwd("")

# ------- CLIMATE Future
setwd("")
dir()
# Check for name consistency
# Future  ------------------------
setwd("Future")
setwd("CCSM4")
setwd("rcp45")
biof.cc.45 <- stack(list.files(pattern = ".tif"))

setwd("..")
setwd("rcp85")
biof.cc.85 <- stack(list.files(pattern = ".tif"))

setwd("")
setwd("Future")
setwd("MIROC-ESM")
setwd("rcp45")
biof.mr.45 <- stack(list.files(pattern = ".tif"))

setwd("..")
setwd("rcp85")
biof.mr.85 <- stack(list.files(pattern = ".tif"))

# --------- CLIMATE PRESENT
setwd("Current")
bio <- stack(list.files(pattern = ".tif"))


#-------- OCCURRENCES
setwd("")
sp_list= list.files() # List of species (binomials separated by space: "Fulanus beltranus")
sp_list 

setwd("")
# ------ SDM

resu <- matrix(nrow = length(sp_list), ncol = 10)
colnames(resu) <- c("sp_name","records","AUC","COR","Deviance","TSS","iniDist" ,
                    "finalDist1","finalDist2","finalDist3")

vars <- NULL # to save which variables were used

# Create a column of species  
t <- as.character(gsub(".txt","",sp_list[i]))

sp= read.table(paste("",sp_list[i],sep=""),h=T)
#sp <- dp[dp$species==t, c('longitude','latitude')]
sp$species <- 1

# Transform to spatial data
coordinates(sp) <- ~ longitude + latitude

# Create "response variable"
sp$species <- 1

# Remove collinear variables
# Climate
spx <- extract(bio, sp, na.omit=T) # extract from file
spx <- data.frame(spx) #convert to dataframe

v <- vifstep(spx)      # check collinearity (variance inflation and correlation)
bio_i <- exclude(bio, v) # exclude collinear predictors
names(bio_i)


# Bounding box  - species-specific background extent
bb <- bbox(sp)                                             #sets a bounding box o coordinate extremes
bb.buf <- extent(bb[1]-35, bb[3]+15, bb[2]-15, bb[4]+15)   # here adjust as much as appropriate (squares or rectangles)
bio_i <- crop(bio_i, bb.buf)                               # species-specific backgroud

# Define directory to save tmp files
setwd("")

# generate sdmData
d <- sdmData(species~., train=sp, predictors= bio_i, bg=list(n=10000))
d
write.sdm(d,"your_dataset")

# generate sdm model
m <- sdm(species ~ . , d, methods=c("domain", "bioclim.dismo", "glm","gam","svm","maxlike","rf","brt"), 
         replication='sub',test.percent=30, n=30,
         parallelSettings = list(ncore=4, method='parallel'))

m
write.sdm(m,"your_model")

# Get variables importance
# 
z <- getVarImp(m) 

svg (paste("", 
		t, "_variable_importance.svg"))
plot(z)
dev.off()


# Get the relationship with each predictor
m1 <- rcurve(m) 
svg (paste("", 
          t, "_predictor_curve.svg"))
plot(m1)
dev.off()


# Ensembling Present
en <- ensemble(m, bio_i, 
               setting=list(method='weighted',stat='TSS'),
               parallelSettings = list(ncore=4, method='parallel'))
plot(en)

writeRaster(en,  paste0("", 
                        t,"_pres.tif"), format = "GTiff", overwrite=TRUE)

# Evaluation
e <- getEvaluation(m)

# Save which variables were used and evaluation results for all species
resu[i, "sp_name"] <- gsub(" ", "_", t)
resu[i, "records"] <- nrow(as.data.frame(d))-10000
resu[i, "AUC"] <- round(mean(e$AUC), 2)
resu[i, "COR"] <- round(mean(e$COR), 2)
resu[i, "Deviance"] <- round(mean(e$Deviance),2)
resu[i, "TSS"] <- round(mean(e$TSS),2)

vars <- c(t, names(bio_i), vars)

# Find binarization threshold
dp <- data.frame(as.data.frame(d),coordinates(d)) # presence points and predictors associated
coords =cbind(dp$longitude,dp$latitude)
pr <- extract(en, coords)

ev <- evaluates(dp$species, pr) # evaluate prediction (observed vs expected) 
th <- ev@threshold_based$threshold[[2]] # threshold that maximizes sensitiv + specificity

# Binary prediction
pa <- en              
pa[] <- ifelse(pa[] >= th, 1, 0) # convert from continuous to binary

plot(pa, main = t)


# present PA
writeRaster(pa,   format = "GTiff",
            paste0("",
                   t,"_pres_PA.tif"),
            overwrite = T)

resu[i, "iniDist"] <- length(pa[pa==1])

# Ensembling past
# LGM --------------------------------------------------------------

enL.cc <- ensemble(m, crop(subset((bioL.cc), names(bio_i)), bb.buf),
                   setting=list(method='weighted',stat='TSS'),
                   parallelSettings = list(ncore=12, method='parallel'))
enL.mr <- ensemble(m, crop(subset((bioL.mr), names(bio_i)), bb.buf),
                   setting=list(method='weighted',stat='TSS'),
                   parallelSettings = list(ncore=12, method='parallel'))
beep(sound=3)

en.LGM <- mean(enL.cc, enL.mr)
#en.LGM <- sdmvspecies::rescale(en.LGM)
plot(en.LGM, main = t)


writeRaster(en.LGM,  paste0("",
                         t,"_LGM.tif"), format = "GTiff", overwrite=TRUE)

# Binary prediction
pa <- en.LGM              
pa[] <- ifelse(pa[] >= th, 1, 0) # convert from continuous to binary
writeRaster(pa,  paste0("",
                         t,"_LGM_PA.tif"), format = "GTiff", overwrite=TRUE)

resu[i, "finalDist1"] <- length(pa[pa==1])

# Ensembling Future rcp45
# Future --------------------------------------------------------------

enf.cc.45 = ensemble(m, crop(subset((biof.cc.45), names(bio_i)), bb.buf),
                   setting=list(method='weighted',stat='TSS'),
                   parallelSettings = list(ncore=12, method='parallel'))
enf.mr.45 = ensemble(m, crop(subset((biof.mr.45), names(bio_i)), bb.buf),
                   setting=list(method='weighted',stat='TSS'),
                   parallelSettings = list(ncore=12, method='parallel'))
beep(sound=3)

en.FUT45 <- mean(enf.cc.45, enf.mr.45)
#en.FUT45 <- sdmvspecies::rescale(en.FUT45)
plot(en.FUT45, main = t)


writeRaster(en.FUT45,  paste0("",
                         t,"_FUT45.tif"), format = "GTiff", overwrite=TRUE)

# Binary prediction
pa <- en.FUT45              
pa[] <- ifelse(pa[] >= th, 1, 0) # convert from continuous to binary
writeRaster(pa,  paste0("",
                         t,"_FUT45_PA.tif"), format = "GTiff", overwrite=TRUE)

resu[i, "finalDist2"] <- length(pa[pa==1])

# Ensembling Future rcp85
# Future --------------------------------------------------------------

enf.cc.85 = ensemble(m, crop(subset((biof.cc.85), names(bio_i)), bb.buf),
                   setting=list(method='weighted',stat='TSS'),
                   parallelSettings = list(ncore=12, method='parallel'))
enf.mr.85 = ensemble(m, crop(subset((biof.mr.85), names(bio_i)), bb.buf),
                   setting=list(method='weighted',stat='TSS'),
                   parallelSettings = list(ncore=12, method='parallel'))
beep(sound=3)

en.FUT85 <- mean(enf.cc.85, enf.mr.85)
#en.FUT85 <- sdmvspecies::rescale(en.FUT85)
plot(en.FUT85, main = t)

writeRaster(en.FUT85,  paste0("",
                         t,"_FUT85.tif"), format = "GTiff", overwrite=TRUE)

# Binary prediction
pa <- en.FUT85              
pa[] <- ifelse(pa[] >= th, 1, 0) # convert from continuous to binary
writeRaster(pa,  paste0("",
                         t,"_FUT85_PA.tif"), format = "GTiff", overwrite=TRUE)

resu[i, "finalDist3"] <- length(pa[pa==1])


write.csv(resu, "~/sdm/resu.csv")
write.csv(vars, "~/sdm/vars.csv")














