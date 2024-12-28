########################################################################################
# Code to calculate time decomposition of life annuity factors
# Single interest rate factor, long term development
#
# Data: 
# England and Wales
# Interest data goes from 1731 to 2024
# Mortality data goes from 1841 to 2021
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
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("Decomp_FUN.R")


# Load data ---------------------------------------------------------------
# Historical data
load("C:/Git/TimeAnnuities/Data/SmoothedData.RData")
rawRates <- aRate
aRate <- aRate[,c("Year", "sdelta")]
names(aRate)[2] <- "delta" # Smoothed delta
aRate <- subset(aRate, Year<=2021)
# Historical data
hmxr <- data.table(merge(hmx,aRate, by = "Year"))
hmxr <- hmxr[,c("country", "Sex", "Year","Age", "Mx", "delta")]

smxr <- rbind(hmxr) # This is the full dataset
smxr <- data.table(arrange(smxr, country, Sex, Year, Age))
# Selection of years and ages to be decomposed ----------------------------

x <- 65 # Age to calculate annuities
Syear <- 1860
Fyear <- 2024# Final year to plot graphs
smxr <- subset(smxr, Age>=x) # Selection of interval

# Actuarial functions, entropy and durations ------------------------------

a <- smxr[, get.ax(Age = Age, mu = Mx, delta = delta), by = list(country,Sex, Year)]

# Rates of change, rho and phi --------------------------------------------

rho <- a[, get.rho(mx = mu, Year= Year, Age= Age), by = list(country, Sex)]
phi <- get.phi(delta = aRate$delta, Year = aRate$Year)

a <- merge(a, rho, by = c("country", "Sex", "Year", "Age"))
a <- merge(a, phi, by = c("Year"))

# Decomposition of changes over time --------------------------------------

adot <- subset(a)[, get.adot(ax=ax,Year=Year), by = list(country, Sex, Age)]
names(adot)[5] <- "adotc"
adotD <- a[,get.decomp(Age=Age, ax=ax, sEx=sEx, mu=mu,
                         delta=delta, rho=rho, phi=phi, H=Hp,D=Dp,Dc=Dc),
             by = list(country,  Sex, Year)]

adotD <- as.data.frame(arrange(adotD, country, Sex, Year))

# Check constant vs proportional changes in interest rates ----------------

adotD <- merge(adotD, aRate, by = "Year")
adotD <- merge(adotD, phi, by = "Year")

adotD$intc <- adotD$phibar * adotD$delta *adotD$Dc
adotD$Dp2 <- adotD$delta *adotD$Dc # Check Duration Prop changes

# Figures -----------------------------------------------------------------
brks <- seq(1860,2020,by = 20)
labs <- c("1860", "'80", "1900","'20","'40","'60","'80","2000", "'20")

intr <- "#2cb889"
long <- "#bf6524"
ann <- "#768F57"

# Interest rates
ggplot(rawRates)+
  geom_line( aes( Year,delta),size = 0.3,  alpha = 0.4)+
  geom_line( aes( Year,sdelta),size = 0.4, colour = intr)+##f5a6ff
  theme_minimal()+
  ylab("Interest rates (%)")+
  scale_y_continuous(breaks = seq(0,1,by = 0.04), expand = c(0,0))+
  scale_x_continuous(breaks = brks,labels = labs, expand = c(0,0))+
  theme_minimal()+
   coord_cartesian(ylim = c(0,0.16), xlim=c(Syear,Fyear))+
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),#  element_rect(fill = "grey30"),
        axis.ticks = element_line(size = 0.1, colour = "gray30"),
        text = element_text(size = 7,font_rc),
        panel.border = element_blank(),
        aspect.ratio = .5,
        strip.text = element_text(size = 7),
        legend.position = "none",
        panel.spacing.x = unit(1,"lines"))

ggsave("Fig/interestRates.svg", width = 2.5, height = 2)

# Raw death rates

ggplot(subset(rawMx, Age %in% c(65, 70,  75)))+
  geom_line( aes( Year,Mx.m, group = as.factor(Age)),size = 0.3,  alpha = 0.3)+
  geom_line(data= subset(smxr, Age %in% c(65,70, 75) & Sex == "Males"),
                         aes( Year,Mx, group = as.factor(Age),
                              colour = as.factor(Age)),
                         size = 0.4)+
  scale_color_manual(values = c("#fcc69f",  "#eb975b","#c27640" ))+
  theme_minimal()+
  ylab("Mortality rates (log)")+
  scale_y_continuous(breaks = c(0.01,seq(0,.09,by = 0.03)), trans = "log", expand = c(0,0))+
  scale_x_continuous(breaks = brks,labels = labs, expand = c(0,0))+
  theme_minimal()+
  coord_cartesian(ylim = c(0.01,0.12), xlim=c(Syear,Fyear))+
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),#  element_rect(fill = "grey30"),
        axis.ticks = element_line(size = 0.1, colour = "gray30"),
        text = element_text(size = 7,font_rc),
        panel.border = element_blank(),
        aspect.ratio = .5,
        strip.text = element_text(size = 7),
        legend.position = "none",
        panel.spacing.x = unit(1,"lines"))

ggsave("Fig/mortalityMales.svg", width = 2.5, height = 2)


# Annuity
ggplot(subset(a, Age %in% c(65) & Sex == "Males" & Year <= Fyear))+
  geom_line( aes( Year,ax, colour = as.factor(Age), group = as.factor(Age)),size = 0.4, colour = "black")+
  # geom_line( aes( Year,Hc, colour = as.factor(Age), group = as.factor(Age)),
  #            size = 0.3, linetype = 1, alpha = 0.7)+
  #scale_color_manual(values = c("#a1fffc", "#b1ff8a"))+
  facet_grid(~Sex, scales="free")+
  theme_minimal()+
  ylab("Life annuity factors at age 65")+
  scale_y_continuous(breaks = seq(2,20,by = 2), expand = c(0,0))+
  scale_x_continuous(breaks = brks,labels = labs, expand = c(0,0))+
  theme_minimal()+
 coord_cartesian(ylim = c(5,16.5), xlim=c(Syear,Fyear))+
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.ticks = element_line(size = 0.1, colour = "gray30"),
        text = element_text(size = 7,font_rc),
        panel.border = element_blank(),
        aspect.ratio = .5,
        strip.text = element_blank(),
        legend.position = "none",
        panel.spacing.x = unit(1,"lines"))

ggsave("Fig/annuityMales.svg", width = 2.5, height = 2)


# Derivative of ax
ggplot(subset(adot, Age %in% c(65) & Sex == "Males" & Year <= Fyear))+
  geom_hline(yintercept = 0, colour = "black", size = 0.1, alpha = 0.5)+
  #geom_line(aes(Year, adot), colour = "blue")+ # decomposition
  geom_line(aes(Year, adotc*100), color = ann,size = 0.4)+ 
  #scale_color_manual(values = c("#e4ff8a","#a1fffc", "#b1ff8a"))+
  facet_wrap(~Sex)+
  ylab("Change in life annuity factors (%)")+
  scale_y_continuous(breaks = seq(-6,6,by = 2), expand = c(0,0))+
  scale_x_continuous(breaks = brks,labels = labs, expand = c(0,0))+
  theme_minimal()+
  coord_cartesian(ylim = c(-6,6), xlim=c(Syear,Fyear))+
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.ticks = element_line(size = 0.1, colour = "gray30"),
        text = element_text(size = 7,font_rc),
        panel.border = element_blank(),
        aspect.ratio = .5,
        strip.text = element_blank(),
        legend.position = "none",
        panel.spacing.x = unit(1,"lines"))


ggsave("Fig/changeAnnuityMales.svg", width = 2.5, height = 2)

# Entropy vs Duration
ggplot(subset(adotD, Sex == "Males" & Year <= 2021) )+
  geom_line(aes(Year, D), colour = intr, size = 0.4)+
  geom_line(aes(Year, H), colour = long, size = 0.4)+
  facet_wrap(~Sex)+
  theme_minimal()+
  scale_y_continuous(breaks = seq(0,1,by = 0.2), expand = c(0,0))+
  scale_x_continuous(breaks = brks,labels = labs, expand = c(0,0))+
  theme_minimal()+
  coord_cartesian(ylim = c(0,0.8), xlim=c(Syear,Fyear))+
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.ticks = element_line(size = 0.1, colour = "gray30"),
        text = element_text(size = 7,font_rc),
        panel.border = element_blank(),
        aspect.ratio = .5,
        strip.text = element_blank(),
        panel.spacing.x = unit(1,"lines"))+
  ylab("Sensitivity to changes in mortality and interest rates")

ggsave("Fig/sensitivities.svg", width = 2.5, height = 2)

# Rates of change over time
ggplot(subset(adotD, Sex == "Males" & Year <= 2021))+
  geom_hline(yintercept = 0, colour = "black", size = 0.1, alpha = 0.5)+
  geom_line(aes(Year, rhobar*100), colour = long, size = 0.4)+
  geom_line(aes(Year, -phibar*delta*1000), colour = intr, size = 0.4)+
  #geom_line(aes(Year, phibar), colour = "white")+
  facet_wrap(~Sex)+
  scale_color_viridis_c()+
  scale_y_continuous(breaks = seq(-9,9,by = 3), expand = c(0,0))+
  scale_x_continuous(breaks = brks,labels = labs, expand = c(0,0))+
  theme_minimal()+
  coord_cartesian(ylim = c(-6,12), xlim=c(Syear,Fyear))+
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.ticks = element_line(size = 0.1, colour = "gray30"),
        text = element_text(size = 7,font_rc),
        panel.border = element_blank(),
        aspect.ratio = .5,
        strip.text = element_blank(),
        panel.spacing.x = unit(1,"lines"))+
  ylab("Changes in mortality and interest rates (%)")

ggsave("Fig/changeMortalityInterest.svg", width = 2.5, height = 2)


# Decomposition
adotD[,c("Year","country","Sex","mort","int")] %>% filter(Sex == "Males") %>% 
  gather(Component, Percentage, -country, -Sex, -Year) %>% 
  ggplot()+
  geom_area( aes( Year, Percentage*100, fill = Component), position = "stack")+
  # geom_bar(aes(x= Year, y = Percentage, fill = Component),
  #          position = "stack", stat = "identity")+
  scale_fill_manual(values = c(intr, long))+
  #facet_wrap(~Sex)+
  scale_y_continuous(breaks = seq(-10,10,by = 2), expand = c(0,0))+
  scale_x_continuous(breaks = brks,labels = labs, expand = c(0,0))+
  theme_minimal()+
  coord_cartesian(ylim = c(-6,6), xlim=c(Syear,Fyear))+
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.ticks = element_line(size = 0.1, colour = "gray30"),
        text = element_text(size = 7,font_rc),
        panel.border = element_blank(),
        legend.position = "none",
        aspect.ratio = .5,
        strip.text = element_text(size = 7),
        panel.spacing.x = unit(1,"lines"))+
  ylab("Contribution to the \n change in life annuity factors (%)")

  
ggsave("Fig/decompositionSingleInterestRate.svg", width = 2.5, height = 2)
  
