###############################################################################
# Smooth mortality and interest data 
# Period
# Time annuities
###############################################################################

rm(list = ls())
library(tidyverse)
library(data.table)
library(splines)
library(lubridate)
library(MortalitySmooth)

setwd("C:/Git/TimeAnnuities/Data/")

iRate <- read.csv("C:/Git/TimeAnnuities/Data/Interest rates/UKlongtermrate.csv", sep=";")


Dx <- read_table2("C:/Git/TimeAnnuities/Data/Mortality/Deaths_1x1.txt",  skip = 2)
Nx <- read_table2("C:/Git/TimeAnnuities/Data/Mortality/Exposures_1x1.txt", skip = 2)


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
endYear   <- 2021
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
Nxf[Nxf == 0] <- 0.00001
Nxf <- as.matrix(Nxf[,-1])

Ages  <- as.integer(rownames(Dxf))
Years <- as.integer(colnames(Dxf))

# Smooth deaths
fit.f <- Mort2Dsmooth(x = Ages, y = Years, 
                          Z = Dxf, offset = log(Nxf))

smx.f <- as.data.frame(fit.f$fitted.values/Nxf)
smx.f$Age <- as.integer(rownames(smx.f))
# Convert matrix to dataframe again
smx.f <- gather(smx.f, Year, Mx,-Age)
smx.f$Year <- as.integer(smx.f$Year)

# Males
Dxm <- as.data.frame(spread(Dxm[,-1],Year, Dx.m))
rownames(Dxm) <- Dxm$Age
#Dxm[Dxm == 0] <- 0.00001
Dxm <- as.matrix(Dxm[,-1])

Nxm <- as.data.frame(spread(Nxm[,-1],Year,Nx.m))
rownames(Nxm) <- Nxm$Age
Nxm[Nxm == 0] <- 0.00001
Nxm <- as.matrix(Nxm[,-1])

Ages  <- as.integer(rownames(Dxm))
Years <- as.integer(colnames(Dxm))

# Smooth deaths
fit.m <- Mort2Dsmooth(x = Ages, y = Years, 
                          Z = Dxm, offset = log(Nxm))

smx.m <- as.data.frame(fit.m$fitted.values/Nxm)
smx.m$Age <- as.integer(rownames(smx.m))
# Convert matrix to dataframe again
smx.m <- gather(smx.m, Year, Mx,-Age)
smx.m$Year <- as.integer(smx.m$Year)

# Merge all ---------------------------------------------------------------

smx.f$Sex <- "Females"
smx.m$Sex <- "Males"   

smx.f$country <- "England and Wales"
smx.m$country <- "England and Wales"

smx <- (rbind(smx.f,
              smx.m))



# Smooth Death counts -----------------------------------------------------

# Females
sdx.f <- as.data.frame(fit.f$fitted.values)
sdx.f$Age <- as.integer(rownames(sdx.f))
# Convert matrix to dataframe again
sdx.f <- gather(sdx.f, Year, Dx,-Age)
sdx.f$Year <- as.integer(sdx.f$Year)

#Males
sdx.m <- as.data.frame(fit.m$fitted.values)
sdx.m$Age <- as.integer(rownames(sdx.m))
# Convert matrix to dataframe again
sdx.m <- gather(sdx.m, Year, Dx,-Age)
sdx.m$Year <- as.integer(sdx.m$Year)

sdx.f$Sex <- "Females"
sdx.m$Sex <- "Males"   

sdx.f$country <- "England and Wales"
sdx.m$country <- "England and Wales"

sdx <-      rbind(sdx.f, sdx.m)

smx <- as.data.frame(merge(smx,sdx, by = c("country","Sex","Year", "Age")))
smx <- arrange(smx, Sex, country, Year, Age)
smx$Mx <- as.numeric(smx$Mx)
smx$Nx <- smx$Dx /smx$Mx
smx <- data.table(smx)
# Split death rates -------------------------------------------------------

iterp <- function(fx, Age, s.age = 0, e.age= 110){
  ages <- seq( s.age, e.age, 0.01 )
  app <- approx(Age,fx, xout = ages)
  x    <- app$x
  y    <- app$y
  out <- data.frame(Age= x, hx = y)
  return(out)}

smx.cont<- smx[,iterp(fx =  Mx, Age = Age, s.age = 0),  by = list(Year, country, Sex)]
sdx.cont<- smx[,iterp(fx =  Dx, Age = Age, s.age = 0),  by = list(Year, country, Sex)]
names(sdx.cont)[5] <- "dx"
rawMx <- sMx

# Interest rates ----------------------------------------------------------
names(iRate) <- c("Year", "Month", "delta")
iRate$delta <- iRate$delta/100
# Annual rates
aRate <- iRate %>% group_by(Year)  %>% summarize(delta = mean(delta))

aRate$sdelta <- smooth.spline(aRate$delta, spar = 0.3)$y

ggplot(aRate)+
  geom_line(aes(Year,delta))+
  geom_line(aes(Year,sdelta), colour = "red")

hmx <- smx[,c("country","Sex","Year" ,"Age","Mx" )]

save(rawMx,hmx, aRate, file="SmoothedData.RData") # Only Mx

