########################################################################################
# Code to tidy up cause of death data
# Cause of death data
# England and Wales, Both sexes
########################################################################################


rm(list = ls())
library(tidyverse)
library(data.table)
library(hrbrthemes)
library(msm)
library(mltools)
library(pracma)
library(splines)
library(hrbrthemes)
library(lemon)
library(stats)
library(lubridate)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("C:/Git/TimeAnnuities/R/Decomp_FUN.R")
load("C:/Git/TimeAnnuities/Data/SmoothedData.RData")

cod <- as.data.frame(read_csv("C:/Git/TimeAnnuities/Data/Mortality/counts_cod.csv"))

cod <- cod %>% mutate(Cause = case_when(Cause %in% c(2)  ~ "Dx1", # Cancer
                                        Cause %in% c(7)  ~ "Dx2", # Hear
                                        Cause %in% c(8)  ~ "Dx3", # Cerebrovascular
                                        Cause %in% c(10,11)  ~ "Dx4", # Respiratory diseases, Cause 17 is the code for Covid19
                                        Cause %in% c(1,3,4,5,6,9,12,13,14,15,16) ~ "Dx5")) 


cod <- cod %>% group_by(country, Sex, Year, Cause) %>% summarise(d65 = sum(d65),
                                                                 d70 = sum(d70),
                                                                 d75 = sum(d75),
                                                                 d80 = sum(d80),
                                                                 d85 = sum(d85),
                                                                 d90 = sum(d90),
                                                                 d95 = sum(d95))

cod <- cod %>% gather(Age, Dx, -country,-Sex,-Year,-Cause )%>% spread(Cause, Dx)

cod[is.na(cod)] <- 0

cod$Age <- as.integer(substr(cod$Age,2,3))

cod$Dx <- apply(cod[,5:9], 1, sum)

cod$p1 <- cod$Dx1 / cod$Dx
cod$p2 <- cod$Dx2 / cod$Dx
cod$p3 <- cod$Dx3 / cod$Dx
cod$p4 <- cod$Dx4 / cod$Dx
cod$p5 <- cod$Dx5 / cod$Dx

save(cod, file = "causeOfDeathCounts.RData")
