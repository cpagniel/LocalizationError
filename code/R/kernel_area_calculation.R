# Load Packages
library(adehabitatHR)
library(sp)

# Set Working Directory
setwd("F:/C_files/C_Spring_2021/Real_Call_Localization/Kernel_Estimation_Activity_Space")

# Import data
fish <- read.csv('xxyy_1109_1527_combo_clean.csv')

# Convert to Spatial Point Data Type
fish.sp <- fish[,c("xx","yy")]
coordinates(fish.sp) <- c("xx","yy")

# Kernel Density Estimation
kud <- kernelUD(fish.sp,h=2.8333)

#combo 1109 to 1527 h = 0.9022
#model 1109 to 1527 h = 1.7596
#combo 564 to 937 h = 1.3958
#model 564 to 937 h = 1.2983

# Kernel Home Range
home95 <- getverticeshr(kud, percent=95, unin = c("m"), unout = c("m2"), standardize = FALSE)
home50 <- getverticeshr(kud, percent=50, unin = c("m"), unout = c("m2"), standardize = FALSE)

# 95% = Total Area, 50% = Core Area
total_area <- home95$area
core_area <- home50$area

# Get Coordinates
coords1_95 <- home95@polygons[[1]]@Polygons[[1]]@coords
coords2_95 <- home95@polygons[[1]]@Polygons[[2]]@coords
coords3_95 <- home95@polygons[[1]]@Polygons[[3]]@coords
coords4_95 <- home95@polygons[[1]]@Polygons[[4]]@coords

coords1_60 <- home60@polygons[[1]]@Polygons[[1]]@coords
coords2_60 <- home60@polygons[[1]]@Polygons[[2]]@coords
coords3_60 <- home60@polygons[[1]]@Polygons[[3]]@coords

# Write Coordinates
write.csv(coords1_95,"coords1_model_1109_1527_95.csv", row.names = FALSE)
write.csv(coords2_95,"coords2_model_1109_1527_95.csv", row.names = FALSE)
write.csv(coords3_95,"coords3_model_1109_1527_95.csv", row.names = FALSE)
write.csv(coords4_95,"coords4_model_1109_1527_95.csv", row.names = FALSE)

write.csv(coords1_60,"coords1_model_1109_1527_60.csv", row.names = FALSE)
write.csv(coords2_60,"coords2_model_1109_1527_60.csv", row.names = FALSE)
write.csv(coords3_60,"coords3_model_1109_1527_60.csv", row.names = FALSE)
