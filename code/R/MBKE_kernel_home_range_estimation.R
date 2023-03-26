## Movement-Base Kernel Estimation (MBKE)

# Load Packages
library(adehabitatHR)
library(dplyr)
library(lubridate)

# Set Working Directory
setwd("F:/C_files/C_Spring_2021/Real_Call_Localization/Kernel_Estimation_Activity_Space")

# Import data
fish <- read.csv('xxyy_564_to_937_model_clean.csv')
#fish <- read.csv('xxyy_1109_1527_model_clean.csv')

# Convert Date Time
t <- data.frame(Time = fish[,3])
t <- t %>% mutate(Time = dmy_hms(Time))

# Create Trajectory
tr <- as.ltraj(data.frame(x = fish[,1],y = fish[,2]),date=t[,1],id="a")

# Diffusion Coefficient
# 564 to 937 Tmax = 12
# 1109 to 1527 Tmax = 11
vv <- BRB.D(tr,Tmax=12,Lmin=0.25)

# Utilisation Distribution
# 564 to 937 hmin = 3.6519 for combo, 3.2043 for model
# 1109 to 1526 hmin = 2.8333 for combo, 4.5868 for model
ud <- BRB(tr,D = vv,Tmax = 12,tau = 1,Lmin = 0.25,hmin = 3.2043,grid = 120,b = 0,same4all = FALSE,extent = 1.5)

# Home Range
home95 <- getverticeshr(ud, percent=95, unin = c("m"), unout = c("m2"), standardize = FALSE)
home50 <- getverticeshr(ud, percent=50, unin = c("m"), unout = c("m2"), standardize = FALSE)

# 95% = Total Area, 50% = Core Area
total_area <- home95$area
core_area <- home50$area

# Get Coordinates
coords1_95 <- home95@polygons[[1]]@Polygons[[1]]@coords
coords2_95 <- home95@polygons[[1]]@Polygons[[2]]@coords

coords1_50 <- home50@polygons[[1]]@Polygons[[1]]@coords
coords2_50 <- home50@polygons[[1]]@Polygons[[2]]@coords

# Write Coordinates
write.csv(coords1_95,"coords1_model_564_937_95_MBKE.csv", row.names = FALSE)
write.csv(coords2_95,"coords2_model_564_937_95_MBKE.csv", row.names = FALSE)

write.csv(coords1_50,"coords1_model_564_937_50_MBKE.csv", row.names = FALSE)
write.csv(coords2_50,"coords2_model_564_937_50_MBKE.csv", row.names = FALSE)
