rm(list = ls())
library(tidyverse)
library(data.table)
library(hrbrthemes)
library(msm)
library(mltools)
library(pracma)
library(splines)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

source("Decomp_FUN.R")
load("SmoothedData.RData")


# Subset mortality data and merge it with interest rates
x <- 65 # Age to calculate annuities
smxr <- data.table(merge(subset(smx, Age>=x),subset(aRate, Year>=1841 & Year<=2015, by = "Year")))

# Actuarial functions, entropy and durations ------------------------------


a <- smxr[, get.ax(Age = Age, mu = Mx, delta = sdelta), by = list(Sex, country, Year)]

# Rates of change, rho and phi --------------------------------------------

rho <- a[, get.rho(mx = mu, Year= Year, Age= Age ), by = list(country, Sex)]
phi <- get.phi(delta = aRate$sdelta, Year = aRate$Year)

a <- merge(a, rho, by = c("country", "Sex", "Year", "Age"))
a <- merge(a, phi, by = c("Year"))


# Changes over time -------------------------------------------------------

adot <- subset(a, Age == x)[, get.adot(ax=ax,Year=Year), by = list(country, Sex)]
names(adot)[4] <- "adotc"
adotD <- a[,get.decomp(Age=Age, ax=ax, sEx=sEx, mu=mu,
                         delta=delta, rho=rho, phi=phi, H=Hp,D=Dp),
             by = list(country,  Sex, Year)]

adotD <- as.data.frame(arrange(adotD, country, Sex, Year))

adotD <- merge(adotD, adot, by = c("country", "Sex", "Year"))
adotD$mortp <- adotD$mort / adotD$adot
adotD$intp <- adotD$int / adotD$adot


# Figures -----------------------------------------------------------------

# Annuity
ggplot(subset(a, Age == x))+
  geom_line(aes(Year, ax))+
  geom_line(aes(Year, Dc), colour = "blue")+
  facet_wrap(~Sex)+
  theme_minimal()+
  ylab("Life annuity and modified duration at age 65")

ggsave("Fig/Annuity.pdf", width = 6, height = 4, device = cairo_pdf)



# Entropy vs Duration
ggplot(adotD)+
  geom_line(aes(Year, D), colour = "red")+
  geom_line(aes(Year, H), colour = "forestgreen")+
  facet_wrap(~Sex)+
  theme_minimal()+
  ylab("Sensitivity of a65")


ggsave("Fig/EntvsDur.pdf", width = 6, height = 4, device = cairo_pdf)


# Derivative of ax
ggplot(adotD)+
  geom_line(aes(Year, adot), colour = "blue")+ # decomposition
  geom_line(aes(Year, adotc))+ # real one
  facet_wrap(~Sex)+
  theme_minimal()+
  # coord_cartesian( xlim = c(1900,2020))+
  ylab("Relative derivative of a65 over time")

ggsave("Fig/adot.pdf", width = 6, height = 4, device = cairo_pdf)

# Rates of change over time

ggplot(adotD)+
  geom_line(aes(Year, rhobar), colour = "red")+
  geom_line(aes(Year, phibar))+
  facet_wrap(~Sex)+
  theme_minimal()+
 # coord_cartesian( xlim = c(1900,2020))+
  ylab("Average phi and rho (rhobar in red and phibar in black)")


ggsave("Fig/RhobarPhibar.pdf", width = 6, height = 4, device = cairo_pdf)


# Decomposition
adotD[,c(1,2,3,8,9)] %>% filter() %>% 
  gather(Component, Percentage, -country, -Sex, -Year) %>% 
  ggplot()+
  geom_bar(aes(x= Year, y = Percentage*100, fill = Component),
           position = "stack", stat = "identity")+
  scale_y_continuous(expand = c(0,0))+
  scale_x_continuous(expand = c(0,0))+
  scale_fill_viridis_d()+
  facet_wrap(~Sex)+
  coord_cartesian( xlim = c(1950,2020))+
  theme_minimal()+
  ylab("Decomposition of changes in adot over time")


ggsave("Fig/Decomposition.pdf", width = 6, height = 4, device = cairo_pdf)#

