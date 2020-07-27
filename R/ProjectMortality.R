###############################################################################
# Project mortality data  
# Period
# Time annuities
###############################################################################

rm(list = ls())
library(tidyverse)
library(data.table)
library(splines)
library(lubridate)
library(MortalitySmooth)
library(StMoMo)

Dx <- read_table2("./Data/Mortality/Deaths_1x1.txt",  skip = 2)
Nx <- read_table2("./Data/Mortality/Exposures_1x1.txt", skip = 2)


names(Dx) <- c("Year", "Age", "Dx.f", "Dx.m", "Dx.t")
names(Nx) <- c("Year", "Age", "Nx.f", "Nx.m", "Nx.t")

Dx$country <- "England and Wales"
Nx$country <- "England and Wales"

#Merge both datasets to obtain Dx
Mx <- merge(Dx,Nx, by =c("country", "Year", "Age"))

Mx <- data.frame(lapply(Mx, as.character), stringsAsFactors=FALSE)
Mx$Age[Mx$Age == "110+"] <-110
Mx[, c(2:9)] <- sapply(Mx[, c(2:9)], as.numeric)
Mx           <- arrange(Mx, country, Year, Age)
Mx[is.na(Mx)]<-0

Mx$Mx.f <- Mx$Dx.f / Mx$Nx.f
Mx$Mx.m <- Mx$Dx.m / Mx$Nx.m
Mx$Mx.t <- Mx$Dx.t / Mx$Nx.t

# Smooth lifetables -------------------------------------------------------

# Select the country and ages
startYear <- 1841
endYear   <- 2018
startAge <- 0
endAge   <- 110
sMx <- subset(Mx, Age >=startAge & Age<=endAge & Year >= startYear & Year <= endYear)

# Great Britain

Dxf <- subset(sMx[,c(1,2,3,4)],country == "England and Wales")
Nxf <- subset(sMx[,c(1,2,3,7)],country == "England and Wales")

Dxm <- subset(sMx[,c(1,2,3,5)],country == "England and Wales")
Nxm <- subset(sMx[,c(1,2,3,8)],country == "England and Wales")

# Females
# Convert data frames into matrices
Dxf <- as.data.frame(spread(Dxf[,-1],Year, Dx.f))
rownames(Dxf) <- Dxf$Age
#Dxf[Dxf == 0] <- 0.00001
Dxf <- as.matrix(Dxf[,-1])

Nxf <- as.data.frame(spread(Nxf[,-1],Year,Nx.f))
rownames(Nxf) <- Nxf$Age
Nxf <- as.matrix(Nxf[,-1])

Ages  <- as.integer(rownames(Dxf))
Years <- as.integer(colnames(Dxf))

# Males
Dxm <- as.data.frame(spread(Dxm[,-1],Year, Dx.m))
rownames(Dxm) <- Dxm$Age
Dxm <- as.matrix(Dxm[,-1])

Nxm <- as.data.frame(spread(Nxm[,-1],Year,Nx.m))
rownames(Nxm) <- Nxm$Age
Nxm <- as.matrix(Nxm[,-1])



# Project using a simple APC model -------------------------------------------------------
ages.fit <- 60:110
years.fit <- 1960:2018
  
wxt <- genWeightMat(ages.fit, years.fit, clip = 3, zeroCohorts = seq(1850,1860))
APCModelf <- fit(apc(), Dxt = Dxf, Ext = Nxf, 
                ages = startAge:endAge, years = startYear:endYear,
                ages.fit = ages.fit, years.fit = years.fit, wxt = wxt)
APCModelm <- fit(apc(), Dxt = Dxm, Ext = Nxm, 
                 ages = startAge:endAge, years = startYear:endYear,
                 ages.fit = ages.fit, years.fit = years.fit, wxt = wxt)



APCModelforf <- forecast(APCModelf)
APCModelform <- forecast(APCModelm)

# plot(APCModelforf, parametricbx = FALSE, nCol = 3)
# plot(APCModelform, parametricbx = FALSE, nCol = 3)
# # APCModelform$kt.f$model$drift
# # APCModelforf$kt.f$model$drift
# 
# APCModelform$gc.f$model
# APCModelforf$gc.f$model
# 
# plot(APCModelforf$rates[ , "2030"], log = "y")
# lines(APCModelform$rates[ , "2030"])

write.csv(as.data.frame(APCModelforf$rates), "./Data/Mortality/projectedMortalityFemales.csv" ,row.names = TRUE)
write.csv(as.data.frame(APCModelform$rates), "./Data/Mortality/projectedMortalityMales.csv" ,row.names = TRUE)
