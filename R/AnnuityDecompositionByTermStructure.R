########################################################################################
# Code to calculate time decomposition of life annuity factors
# Term structure of interest rates
#
# Data: 
# England and Wales
# Interest data goes from 1970 to 2024
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
library(stats)
library(lubridate)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("Decomp_FUN.R")

# Load data ---------------------------------------------------------------
# Mortality data
load("C:/Git/TimeAnnuities/Data/SmoothedData.RData")

# Term structure of interest rates
rstr <- as.data.frame(read_csv("C:/Git/TimeAnnuities/Data/Interest rates/Term Structure/Processed_Term_Structure.csv"))

# Extract year of each term-struc
rstr$Year <- as.integer(year(rstr$date))

# Average by year
rstr <- rstr %>% group_by(Year,term)  %>% summarize(rawdelta = mean(forward_force))
rstr$term  <- rstr$term - 1 # To start at age 0, e.g. age 65

# smooth rates over the term and over time
years <- 1970:2024 
rsm <- data.table(rstr)
rsm <- rsm[,smooth.spline(rawdelta, spar = 0)$y,by = Year ]
rsm$term <- rstr$term
rsm1 <- as.data.frame(rsm[,smooth.spline(V1, spar = 0.52)$y,by = term ])
rsm1$Year <- years
rstr1 <- merge(rstr, rsm1, by = c("Year", "term"))
names(rstr1)[4] <- "delta"
rstr1 <- arrange(rstr1, Year, term)

rsm2 <- as.data.frame(rsm[,smooth.spline(V1, spar = 0.42)$y,by = term ])
rsm2$Year <- years
rstr2 <- merge(rstr, rsm2, by = c("Year", "term"))
names(rstr2)[4] <- "delta"
rstr2 <- arrange(rstr2, Year, term)
yearcutoff <- 2015
rstr <- as.data.frame(rbind(subset(rstr1, Year<=yearcutoff),subset(rstr2, Year>yearcutoff) ))
# Merge interest rates with mortality data
x <- 65 # Age to calculate annuities

rstr$Age <- rstr$term + x

hmxr <- merge(hmx,rstr, by = c("Year", "Age"))
hmxr <- hmxr[,c("country", "Sex", "Year","Age","term", "Mx", "delta")]
hmxr <- data.table(arrange(hmxr, country, Sex, Year, Age))

# Actuarial functions, entropy and durations ------------------------------

a <- hmxr[, get.ax(Age = Age, mu = Mx, delta = delta), by = list(country,Sex, Year)]

# Rates of change, rho and phi --------------------------------------------

rho <- a[, get.rho(mx = mu, Year= Year, Age= Age), by = list(country, Sex)]

# This function is the one that captures the term structure of interest rates
phi <- a[, get.phi.str(delta = delta, Year= Year, Age= Age), by = list(country, Sex)]

a <- merge(a, rho, by = c("country", "Sex", "Year", "Age"))
a <- merge(a, phi, by = c("country", "Sex", "Year", "Age"))

# Decomposition by age ----------------------------------------------------

a.age <- a

a.age$xt <- case_when(a.age$Age %in% c(65:74) ~1,
                      a.age$Age %in% c(75:84) ~2,
                      a.age$Age %in% c(85:94) ~3,
                      a.age$Age %in% c(95:104)~4)


adotD.a <- a.age[,get.age.decomp(Age=Age, ax=ax, sEx=sEx, mu=mu,
                       delta=delta, rho=rho, phi=phi),
                       by = list(country,  Sex, Year, xt)]


adotD.age <- adotD.a[,c("country","Sex","Year","mort","int", "xt")] %>% filter() %>% 
             gather(Component, Percentage, -country, -Sex, -Year, -xt)


adotD.age$Component <- paste(adotD.age$Component, adotD.age$xt, sep= "")

adotD.age <- arrange(adotD.age, Sex,  Year)

# Full decomposition ------------------------------------------------------

adot <- a[, get.adot(ax=ax,Year=Year), by = list(country, Sex, Age)]
names(adot)[5] <- "adotc"
adotD <- a[,get.decomp(Age=Age, ax=ax, sEx=sEx, mu=mu,
                         delta=delta, rho=rho, phi=phi, H=Hp,D=Dp,Dc=Dc),
             by = list(country,  Sex, Year)]

adotD <- as.data.frame(arrange(adotD, country, Sex, Year))

# Figures -----------------------------------------------------------------
brks <- seq(1970,2020,by = 10)
labs <- c("1970", "'80","'90","2000","'10","'20")

Syear <- 1970
Fyear <- 2024# Final year to plot graphs

# Colour palette
#ii <- viridis::magma(4,alpha = 1,begin = 0.2, end =0.7,direction = -1)
#kk <- viridis::viridis(4, alpha = 1, begin = 0.6, end = 1, direction = -1)

ii <- c( "#006140", "#00bd7d","#92f7ae", "#bdffeb" )
kk <- c( "#ba3c06", "#d4742f", "#ffa463", "#ffdec7")
jj<- c(ii,kk)

# Age decomposition

  ggplot(subset(adotD.age, Sex == "Males"))+
 #geom_area( aes( Year, Percentage*100, fill = Component),position = "stack", alpha =1)+
  geom_bar(aes(x= Year, y = Percentage*100, fill = Component),
            position = "stack", stat = "identity",width = 1, alpha = 1)+
  scale_fill_manual(values = jj)+
  facet_wrap(~Sex)+
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
        aspect.ratio = 1,
        strip.text = element_blank(),
        panel.spacing.x = unit(1,"lines"))+
  ylab("Contribution to the \n change in life annuity factors (%)")

ggsave("Fig/ageTermStructureDescomposition.svg", width = 2.5, height = 2)

  
# Attribution of Durations and entropies 

adotD.a[,c("Year","country","Sex","H","D", "xt")] %>% filter(Sex == "Males" & Year %in% c(1990,2000,2015,2021)) %>% 
  gather(Component, Percentage, -country, -Sex, -Year, -xt) %>% 
  ggplot()+
  # geom_area( aes( Year, Percentage*100, fill = Component), position = "stack")+
  geom_bar(aes(x= xt, y = Percentage, fill = Component),
           position = "dodge", stat = "identity", width = 0.5)+
  scale_fill_manual(values = c("#2cb889", "#bf6524"))+
  facet_wrap(~Year, ncol = 2,scales = "free")+
  scale_y_continuous(breaks = seq(0,0.6,by = 0.1), expand = c(0,0))+
  scale_x_continuous(labels = c("65-74", "75-84", "85-94", "95-104"), expand = c(0,0))+
  theme_minimal()+
  coord_cartesian(ylim = c(0,0.6))+
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.ticks = element_line(size = 0.1, colour = "gray30"),
        text = element_text(size = 7,font_rc),
        panel.border = element_blank(),
        legend.position = "none",
        aspect.ratio = 1,
        strip.text = element_text(size = 7),
        panel.spacing.x = unit(2,"lines"))+
  ylab("Sensitivity to changes in mortality and interest rates")+
  xlab("Ages")

ggsave("Fig/attributionEntropyDuration.svg", width = 3, height = 5)
