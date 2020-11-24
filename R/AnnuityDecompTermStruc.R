########################################################################################
# Code to calculate time decomposition of life annuity factors
# Historical data using the term structure of interest rates
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
source("Decomp_FUN.R")


# Load data ---------------------------------------------------------------
# Mortality data
load("SmoothedData.RData")


# Term structure of interest rates
rstr <- as.data.frame(read_csv("Processed_Term_Structure.csv"))

# Extract year of each term-struc
rstr$Year <- as.integer(year(rstr$date))

# Average by year
rstr <- rstr %>% group_by(Year,term)  %>% summarize(rawdelta = mean(forward_force))
rstr$term  <- rstr$term - 1 # To start at age 0, e.g. age 65

# smooth rates over the term and over time
# I have to smooth them to have a continuous estimation of the derivative

rsm <- data.table(rstr)
rsm <- rsm[,smooth.spline(rawdelta, spar = 0.6)$y,by = Year ]
rsm$term <- rstr$term
rsm <- as.data.frame(rsm[,smooth.spline(V1, spar = 0.6)$y,by = term ])
years <- 1970:2020 # I HAVE TO FIND THE WAY TO DO IT AUTOMATICALLY
rsm$Year <- years

rstr <- merge(rstr, rsm, by = c("Year", "term"))
names(rstr)[4] <- "delta"
rstr <- arrange(rstr, Year, term)



# Merge interest rates with mortality data
x <- 65 # Age to calculate annuities

rstr$Age <- rstr$term + x

hmxr <- merge(hmx,rstr, by = c("Year", "Age"))
hmxr <- hmxr[,c("country", "Sex", "Year","Age","term", "Mx", "delta")]
hmxr <- data.table(arrange(hmxr, country, Sex, Year, Age))


# Actuarial functions, entropy and durations ------------------------------

# These ones stay the same, the functions are done to cope with the term structure
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
Fyear <- 2020# Final year to plot graphs


# Colour pallette
ii <- viridis::magma(4,alpha = 1,begin = 0.2, end =0.7,direction = -1)
kk <- viridis::viridis(4, alpha = 1, begin = 0.6, end = 1, direction = -1)
jj<- c(ii,kk)



# Age decomposition

  ggplot(adotD.age)+
 #geom_area( aes( Year, Percentage*100, fill = Component), position = "stack")+
   geom_bar(aes(x= Year, y = Percentage*100, fill = Component),
            position = "stack", stat = "identity",width = 1)+
  scale_fill_manual(values = jj)+
  facet_wrap(~Sex)+
  scale_y_continuous(breaks = seq(-10,10,by = 2), expand = c(0,0))+
  scale_x_continuous(breaks = brks,labels = labs, expand = c(0,0))+
  theme_minimal()+
  coord_cartesian(ylim = c(-4,4), xlim=c(Syear,Fyear))+
  theme(panel.grid = element_blank(),
        panel.background = element_rect(fill = "gray30"),
        axis.ticks = element_line(size = 0.1, colour = "gray30"),
        text = element_text(size = 7,font_rc),
        panel.border = element_blank(),
        legend.position = "none",
        aspect.ratio = .5,
        strip.text = element_text(size = 7),
        panel.spacing.x = unit(1,"lines"))+
  ylab("Contribution (%)")

ggsave("Fig/AgeDescomposition.pdf", width = 4, height = 2, device = cairo_pdf)

  
 
  
  
#Decomposition
adotD[,c("Year","country","Sex","mort","int")] %>% filter() %>% 
  gather(Component, Percentage, -country, -Sex, -Year) %>% 
  ggplot()+
 # geom_area( aes( Year, Percentage*100, fill = Component), position = "stack")+
   geom_bar(aes(x= Year, y = Percentage*100, fill = Component),
            position = "stack", stat = "identity", width = 1)+
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


#ggsave("Fig/StrucDescomposition.pdf", width = 4, height = 2, device = cairo_pdf)

# Attribution of Durations and entropies 

adotD.a[,c("Year","country","Sex","H","D", "xt")] %>% filter(Sex == "Males" & Year %in% c(1975, 2005)) %>% 
  gather(Component, Percentage, -country, -Sex, -Year, -xt) %>% 
  ggplot()+
  # geom_area( aes( Year, Percentage*100, fill = Component), position = "stack")+
  geom_bar(aes(x= xt, y = Percentage, fill = Component),
           position = "dodge", stat = "identity", width = 0.5)+
  scale_fill_manual(values = c("#fc6579", "#defc65"))+
  facet_wrap(~Year)+
  scale_y_continuous(breaks = seq(0,0.6,by = 0.1), expand = c(0,0))+
  scale_x_continuous(labels = c("65-74", "75-84", "85-94", "95-104"), expand = c(0,0))+
  theme_minimal()+
  coord_cartesian(ylim = c(0,0.6))+
  theme(panel.grid = element_blank(),
        panel.background = element_rect(fill = "grey30"),
        axis.ticks = element_line(size = 0.1, colour = "gray30"),
        text = element_text(size = 7,font_rc),
        panel.border = element_blank(),
        legend.position = "none",
        aspect.ratio = .5,
        strip.text = element_text(size = 7),
        panel.spacing.x = unit(1,"lines"))+
  ylab("Sensitivity")+
  xlab("Ages")


#ggsave("Fig/AttributionDH.pdf", width = 4, height = 2, device = cairo_pdf)


# Interest rates
ggplot(subset(rstr , term ==0))+
  geom_line( aes( Year,rawdelta),size = 0.3,  alpha = 0.4)+
  geom_line( aes( Year,delta),size = 0.3, colour = "#8600f2")+##f5a6ff
  theme_minimal()+
  ylab("Long-term interest rates (%)")+
  scale_y_continuous(breaks = seq(0,1,by = 0.02), expand = c(0,0))+
  scale_x_continuous(breaks = brks,labels = labs, expand = c(0,0))+
  theme_minimal()+
   coord_cartesian(ylim = c(0,0.16), xlim=c(Syear,2020))+
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),#  element_rect(fill = "grey30"),
        axis.ticks = element_line(size = 0.1, colour = "gray30"),
        text = element_text(size = 7,font_rc),
        panel.border = element_blank(),
        aspect.ratio = .5,
        strip.text = element_text(size = 7),
        legend.position = "none",
        panel.spacing.x = unit(1,"lines"))

#ggsave("Fig/delta.pdf", width = 4, height = 2, device = cairo_pdf)

# Raw death rates

ggplot(subset(rawMx, Age %in% c(65, 70,  75)))+
  geom_line( aes( Year,Mx.f, group = as.factor(Age)),size = 0.3,  alpha = 0.3)+
  geom_line(data= subset(hmxr, Age %in% c(65,70, 75) & Sex == "Females"),
                         aes( Year,Mx, group = as.factor(Age),
                              colour = as.factor(Age)),
                         size = 0.3)+
  # geom_line( aes( Year,Mx.m, group = as.factor(Age)),size = 0.3,  alpha = 0.3)+
  # geom_line(data= subset(hmxr, Age %in% c(65,75) & Sex == "Males"),
  #           aes( Year,Mx, group = as.factor(Age),
  #                colour = as.factor(Age)),
  #           size = 0.3)+
  
  #geom_line( aes( Year,sdelta),size = 0.3, colour = "#8600f2")+##f5a6ff
  scale_color_manual(values = c("#05cc02","#00c2f7", "#ccae02"))+
  theme_minimal()+
  ylab("Death rates (log)")+
  scale_y_continuous(breaks = c(0.02,seq(0,.1,by = 0.02)), trans = "log", expand = c(0,0))+
  scale_x_continuous(breaks = brks,labels = labs, expand = c(0,0))+
  theme_minimal()+
  coord_cartesian(ylim = c(0.005,0.12), xlim=c(Syear,2020))+
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),#  element_rect(fill = "grey30"),
        axis.ticks = element_line(size = 0.1, colour = "gray30"),
        text = element_text(size = 7,font_rc),
        panel.border = element_blank(),
        aspect.ratio = .5,
        strip.text = element_text(size = 7),
        legend.position = "none",
        panel.spacing.x = unit(1,"lines"))

#ggsave("Fig/muf.pdf", width = 4, height = 2, device = cairo_pdf)


# Annuity
ggplot(subset(a, Age %in% c(65) & Sex == "Males" & Year <= 2018))+
  geom_line( aes( Year,ax, colour = as.factor(Age), group = as.factor(Age)),size = 0.7)+
  # geom_line( aes( Year,Hc, colour = as.factor(Age), group = as.factor(Age)),
  #            size = 0.3, linetype = 1, alpha = 0.7)+
  scale_color_manual(values = c("#a1fffc", "#b1ff8a"))+
  facet_grid(~Sex, scales="free")+
  theme_minimal()+
  ylab("Life annuity factors")+
  scale_y_continuous(breaks = seq(2,20,by = 2), expand = c(0,0))+
  scale_x_continuous(breaks = brks,labels = labs, expand = c(0,0))+
  theme_minimal()+
 coord_cartesian(ylim = c(6,16), xlim=c(Syear,2020))+
  theme(panel.grid = element_blank(),
        panel.background = element_rect(fill = "grey30"),
        axis.ticks = element_line(size = 0.1, colour = "gray30"),
        text = element_text(size = 7,font_rc),
        panel.border = element_blank(),
        aspect.ratio = .5,
        strip.text = element_text(size = 7),
        legend.position = "none",
        panel.spacing.x = unit(1,"lines"))

#ggsave("Fig/axm.pdf", width = 4, height = 2, device = cairo_pdf)

# Modified duration
ggplot(subset(a, Age %in% c(65)))+
  geom_line( aes( Year,rho, colour = as.factor(Age), group = as.factor(Age)),size = 0.3)+
  #geom_line( aes( Year,Hc, colour = as.factor(Age), group = as.factor(Age)),
  #           size = 0.3, linetype = 3, alpha = 0.7)+
  #geom_line( aes( Year,Hc), colour = "#f5a57f",size = 0.3)+
  scale_color_manual(values = c("#e4ff8a","#a1fffc", "#b1ff8a"))+
  facet_grid(~Sex, scales="free")+
  theme_minimal()+
  ylab("Modified duration")+
  #scale_y_continuous(breaks = seq(2,20,by = 4), expand = c(0,0))+
  #scale_x_continuous(breaks = brks,labels = labs, expand = c(0,0))+
  theme_minimal()+
  #coord_cartesian(ylim = c(2,16), xlim=c(Syear,Fyear))+
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
ggplot(subset(adot, Age %in% c(65) & Sex == "Males" & Year <= 2018))+
  geom_hline(yintercept = 0, colour = "white", alpha = 0.7, size = 0.4)+
  #geom_line(aes(Year, adot), colour = "blue")+ # decomposition
  geom_line(aes(Year, -adotc*100, color = as.factor(Age)),size = 0.7)+ 
  scale_color_manual(values = c("#e4ff8a","#a1fffc", "#b1ff8a"))+
  facet_wrap(~Sex)+
  ylab("Change over time (%)")+
  scale_y_continuous(breaks = seq(-6,6,by = 2), expand = c(0,0))+
  scale_x_continuous(breaks = brks,labels = labs, expand = c(0,0))+
  theme_minimal()+
  coord_cartesian(ylim = c(-4,5), xlim=c(Syear,2020))+
  theme(panel.grid = element_blank(),
        panel.background = element_rect(fill = "grey30"),
        axis.ticks = element_line(size = 0.1, colour = "gray30"),
        text = element_text(size = 7,font_rc),
        panel.border = element_blank(),
        aspect.ratio = .5,
        strip.text = element_text(size = 7),
        legend.position = "none",
        panel.spacing.x = unit(1,"lines"))


#ggsave("Fig/adotm.pdf", width = 4, height = 2, device = cairo_pdf)

# Entropy vs Duration
ggplot(subset(adotD, Sex == "Males" & Year <= 2018) )+
  geom_line(aes(Year, D), colour = "#f57f9e")+
  geom_line(aes(Year, H), colour = "#8bf7ab")+
  facet_wrap(~Sex)+
  theme_minimal()+
  ylab("Sensitivity of life annuity factors")+
  scale_y_continuous(breaks = seq(0,1,by = 0.2), expand = c(0,0))+
  scale_x_continuous(breaks = brks,labels = labs, expand = c(0,0))+
  theme_minimal()+
  coord_cartesian(ylim = c(0,0.8), xlim=c(Syear,2020))+
  theme(panel.grid = element_blank(),
        panel.background = element_rect(fill = "grey30"),
        axis.ticks = element_line(size = 0.1, colour = "gray30"),
        text = element_text(size = 7,font_rc),
        panel.border = element_blank(),
        aspect.ratio = .5,
        strip.text = element_text(size = 7),
        panel.spacing.x = unit(1,"lines"))


#ggsave("Fig/Sens65.pdf", width = 4, height = 2, device = cairo_pdf)




# Rates of change over time

#adotD <- merge(adotD, aRate, by = "Year")


ggplot(subset(adotD, Sex == "Males" & Year <= 2018))+
  geom_hline(yintercept = 0, colour = "white", alpha = 0.7, size = 0.1)+
  geom_line(aes(Year, rhobar*100), colour = "#8bf7ab", size = 0.7)+
  geom_line(aes(Year, -phibar*100), colour = "#f5a6ff", size = 0.7)+
  #geom_line(aes(Year, phibar), colour = "white")+
  facet_wrap(~Sex)+
  scale_color_viridis_c()+
  scale_y_continuous(breaks = seq(-10,10,by = 2), expand = c(0,0))+
  scale_x_continuous(breaks = brks,labels = labs, expand = c(0,0))+
  theme_minimal()+
  coord_cartesian(ylim = c(-6,10), xlim=c(Syear,2020))+
  theme(panel.grid = element_blank(),
        panel.background = element_rect(fill = "grey30"),
        axis.ticks = element_line(size = 0.1, colour = "gray30"),
        text = element_text(size = 7,font_rc),
        panel.border = element_blank(),
        aspect.ratio = .5,
        strip.text = element_text(size = 7),
        panel.spacing.x = unit(1,"lines"))+
  ylab("Average change (%)")


#ggsave("Fig/PhiRho.pdf", width = 4, height = 2, device = cairo_pdf)

