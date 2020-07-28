########################################################################################
# Code to calculate time decomposition of life annuity factors
# Historical and projected data
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
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("Decomp_FUN.R")


# Load data ---------------------------------------------------------------
# Historical data
load("SmoothedData.RData")

aRate <- aRate[,c("Year", "sdelta")]
names(aRate)[2] <- "delta" # Smoothed delta
aRate <- subset(aRate, Year<2020)
# Historical data
hmxr <- data.table(merge(hmx,aRate, by = "Year"))
hmxr <- hmxr[,c("country", "Sex", "Year","Age", "Mx", "delta")]

# Projected data
pRt     <- read_csv("C:/Users/jmartinez/OneDrive - Syddansk Universitet/TimeAnnuities/Data/Interest rates/projectedInterestRates.csv")
pMxf   <- as.data.frame(read_csv("C:/Users/jmartinez/OneDrive - Syddansk Universitet/TimeAnnuities/Data/Mortality/projectedMortalityFemales.csv"))
pMxm   <- as.data.frame(read_csv("C:/Users/jmartinez/OneDrive - Syddansk Universitet/TimeAnnuities/Data/Mortality/projectedMortalityMales.csv"))

pRt$rate <- pRt$rate/100
pRt <- subset(pRt, Year>=2020) # Take just forecasted rates
pRate <- pRt %>% group_by(Year)  %>% summarize(delta = mean(rate))
pRate <- data.frame(rbind(aRate, pRate))

# Transpose databases to long format 
pMxf <- pMxf %>% gather(Year, Mx,-X1)
pMxm <- pMxm %>% gather(Year, Mx,-X1)

pMxf$Sex <- "Females"
pMxm$Sex <- "Males"

pmx <- data.frame(rbind(pMxf, pMxm))
names(pmx)[1] <- "Age"
pmx$Year <- as.integer(pmx$Year)

# Projected mortality and interest rates
pmxr <- data.table(merge(pmx,pRate, by = "Year"))

pmxr$country <- "England and Wales"

pmxr <- pmxr[,c("country", "Sex", "Year","Age", "Mx", "delta")]

smxr <- rbind(hmxr, pmxr) # This is the full dataset
smxr <- data.table(arrange(smxr, country, Sex, Year, Age))
# Selection of years and ages to be decomposed ----------------------------

x <- 65 # Age to calculate annuities
Syear <- 1860
Fyear <- 2070# Final year to plot graphs
smxr <- subset(smxr, Age>=x) # Selection of interval


# Actuarial functions, entropy and durations ------------------------------

a <- smxr[, get.ax(Age = Age, mu = Mx, delta = delta), by = list(country,Sex, Year)]

# Rates of change, rho and phi --------------------------------------------

rho <- a[, get.rho(mx = mu, Year= Year, Age= Age), by = list(country, Sex)]
phi <- get.phi(delta = pRate$delta, Year = pRate$Year)

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


adotD <- merge(adotD, pRate, by = "Year")
adotD <- merge(adotD, phi, by = "Year")

adotD$intc <- adotD$phibar * adotD$delta *adotD$Dc
adotD$Dp2 <- adotD$delta *adotD$Dc # Check Duration Prop changes



# Figures -----------------------------------------------------------------
brks <- seq(1860,2060,by = 20)
labs <- c("1860", "'80", "1900","'20","'40","'60","'80","2000", "'20", "'40", "'60")

# Annuity
ggplot(subset(a, Age %in% c(65,70,75)))+
  geom_line( aes( Year,ax, colour = as.factor(Age), group = as.factor(Age)),size = 0.3)+
  # geom_line( aes( Year,Hc, colour = as.factor(Age), group = as.factor(Age)),
  #            size = 0.3, linetype = 1, alpha = 0.7)+
  scale_color_manual(values = c("#e4ff8a","#a1fffc", "#b1ff8a"))+
  facet_grid(~Sex, scales="free")+
  theme_minimal()+
  ylab("Life annuity factors")+
  scale_y_continuous(breaks = seq(2,20,by = 4), expand = c(0,0))+
  scale_x_continuous(breaks = brks,labels = labs, expand = c(0,0))+
  theme_minimal()+
# coord_cartesian(ylim = c(2,16), xlim=c(Syear,Fyear))+
  theme(panel.grid = element_blank(),
        panel.background = element_rect(fill = "grey30"),
        axis.ticks = element_line(size = 0.1, colour = "gray30"),
        text = element_text(size = 7,font_rc),
        panel.border = element_blank(),
        aspect.ratio = .5,
        strip.text = element_text(size = 7),
        legend.position = "none",
        panel.spacing.x = unit(1,"lines"))

ggsave("Fig/Proyax.pdf", width = 4, height = 2, device = cairo_pdf)

# Modified duration
ggplot(subset(a, Age %in% c(65,70,75)))+
  geom_line( aes( Year,Hc, colour = as.factor(Age), group = as.factor(Age)),size = 0.3)+
  #geom_line( aes( Year,Hc, colour = as.factor(Age), group = as.factor(Age)),
  #           size = 0.3, linetype = 3, alpha = 0.7)+
  #geom_line( aes( Year,Hc), colour = "#f5a57f",size = 0.3)+
  scale_color_manual(values = c("#e4ff8a","#a1fffc", "#b1ff8a"))+
  facet_grid(~Sex, scales="free")+
  theme_minimal()+
  ylab("Modified duration")+
  scale_y_continuous(breaks = seq(2,20,by = 4), expand = c(0,0))+
  scale_x_continuous(breaks = brks,labels = labs, expand = c(0,0))+
  theme_minimal()+
  coord_cartesian(ylim = c(2,16), xlim=c(Syear,Fyear))+
  theme(panel.grid = element_blank(),
        panel.background = element_rect(fill = "grey30"),
        axis.ticks = element_line(size = 0.1, colour = "gray30"),
        text = element_text(size = 7,font_rc),
        panel.border = element_blank(),
        aspect.ratio = .5,
        strip.text = element_text(size = 7),
        legend.position = "none",
        panel.spacing.x = unit(1,"lines"))

#ggsave("Fig/Dur.pdf", width = 4, height = 2, device = cairo_pdf)


# Derivative of ax
ggplot(subset(adot, Age %in% c(65,70,75)))+
  geom_hline(yintercept = 0, colour = "white", alpha = 0.7, size = 0.1)+
  #geom_line(aes(Year, adot), colour = "blue")+ # decomposition
  geom_line(aes(Year, adotc*100, color = as.factor(Age)),size = 0.2)+ 
  scale_color_manual(values = c("#e4ff8a","#a1fffc", "#b1ff8a"))+
  facet_wrap(~Sex)+
  ylab("Change over time (%)")+
  scale_y_continuous(breaks = seq(-6,6,by = 2), expand = c(0,0))+
  scale_x_continuous(breaks = brks,labels = labs, expand = c(0,0))+
  theme_minimal()+
  coord_cartesian(ylim = c(-6,6), xlim=c(Syear,Fyear))+
  theme(panel.grid = element_blank(),
        panel.background = element_rect(fill = "grey30"),
        axis.ticks = element_line(size = 0.1, colour = "gray30"),
        text = element_text(size = 7,font_rc),
        panel.border = element_blank(),
        aspect.ratio = .5,
        strip.text = element_text(size = 7),
        legend.position = "none",
        panel.spacing.x = unit(1,"lines"))


ggsave("Fig/Proyadot.pdf", width = 4, height = 2, device = cairo_pdf)

# Entropy vs Duration
ggplot(adotD)+
  geom_line(aes(Year, D), colour = "#f57f9e")+
  geom_line(aes(Year, H), colour = "#e1f57f")+
  facet_wrap(~Sex)+
  theme_minimal()+
  ylab("Sensitivity of life annuity factors")+
  scale_y_continuous(breaks = seq(0,1,by = 0.2), expand = c(0,0))+
  scale_x_continuous(breaks = brks,labels = labs, expand = c(0,0))+
  theme_minimal()+
  coord_cartesian(ylim = c(0,0.8), xlim=c(Syear,Fyear))+
  theme(panel.grid = element_blank(),
        panel.background = element_rect(fill = "grey30"),
        axis.ticks = element_line(size = 0.1, colour = "gray30"),
        text = element_text(size = 7,font_rc),
        panel.border = element_blank(),
        aspect.ratio = .5,
        strip.text = element_text(size = 7),
        panel.spacing.x = unit(1,"lines"))


ggsave("Fig/ProySens65.pdf", width = 4, height = 2, device = cairo_pdf)




# Rates of change over time

#adotD <- merge(adotD, aRate, by = "Year")


ggplot(adotD)+
  geom_hline(yintercept = 0, colour = "white", alpha = 0.7, size = 0.1)+
  geom_line(aes(Year, rhobar*100), colour = "#8bf7ab", size = 0.3)+
  geom_line(aes(Year, -phibar*delta*1000), colour = "#f5a6ff", size = 0.3)+
  #geom_line(aes(Year, phibar), colour = "white")+
  facet_wrap(~Sex)+
  scale_color_viridis_c()+
  scale_y_continuous(breaks = seq(-10,10,by = 2), expand = c(0,0))+
  scale_x_continuous(breaks = brks,labels = labs, expand = c(0,0))+
  theme_minimal()+
  coord_cartesian(ylim = c(-6,10), xlim=c(Syear,Fyear))+
  theme(panel.grid = element_blank(),
        panel.background = element_rect(fill = "grey30"),
        axis.ticks = element_line(size = 0.1, colour = "gray30"),
        text = element_text(size = 7,font_rc),
        panel.border = element_blank(),
        aspect.ratio = .5,
        strip.text = element_text(size = 7),
        panel.spacing.x = unit(1,"lines"))+
  ylab("Average change (%)")


ggsave("Fig/ProyPhiRho.pdf", width = 4, height = 2, device = cairo_pdf)


# Decomposition
adotD[,c("Year","country","Sex","mort","int")] %>% filter() %>% 
  gather(Component, Percentage, -country, -Sex, -Year) %>% 
  ggplot()+
  geom_area( aes( Year, Percentage*100, fill = Component), position = "stack")+
  # geom_bar(aes(x= Year, y = Percentage, fill = Component),
  #          position = "stack", stat = "identity")+
  scale_fill_manual(values = c("#fc6579", "#defc65"))+
  facet_wrap(~Sex)+
  scale_y_continuous(breaks = seq(-10,10,by = 2), expand = c(0,0))+
  scale_x_continuous(breaks = brks,labels = labs, expand = c(0,0))+
  theme_minimal()+
  coord_cartesian(ylim = c(-4,4), xlim=c(Syear,Fyear))+
  theme(panel.grid = element_blank(),
        panel.background = element_rect(fill = "grey30"),
        axis.ticks = element_line(size = 0.1, colour = "gray30"),
        text = element_text(size = 7,font_rc),
        panel.border = element_blank(),
        legend.position = "none",
        aspect.ratio = .5,
        strip.text = element_text(size = 7),
        panel.spacing.x = unit(1,"lines"))+
  ylab("Contribution (%)")

  
ggsave("Fig/ProyDescomposition.pdf", width = 4, height = 2, device = cairo_pdf)
  
