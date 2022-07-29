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
rawRates <- aRate
aRate <- aRate[,c("Year", "sdelta")]
names(aRate)[2] <- "delta" # Smoothed delta
#aRate <- subset(aRate, Year<2020)
# Historical data
hmxr <- data.table(merge(hmx,aRate, by = "Year"))
hmxr <- hmxr[,c("country", "Sex", "Year","Age", "Mx", "delta")]


# Causes of death ---------------------------------------------------------

load("cod.RData")
hmxr <- hmxr %>% mutate(i = case_when(Age %in% c(65:69) ~ "65",
                              Age %in% c(70:74) ~ "70",
                              Age %in% c(75:79) ~ "75",
                              Age %in% c(80:84) ~ "80",
                              Age %in% c(85:89) ~ "85",
                              Age %in% c(90:94) ~ "90",
                              Age %in% c(95:110)~ "95"))



cod <- cod %>% mutate(i = case_when(Age %in% c(65:69) ~ "65",
                                      Age %in% c(70:74) ~ "70",
                                      Age %in% c(75:79) ~ "75",
                                      Age %in% c(80:84) ~ "80",
                                      Age %in% c(85:89) ~ "85",
                                      Age %in% c(90:94) ~ "90",
                                      Age %in% c(95:110)~ "95"))


# Merge

hmxr <- merge(hmxr, cod[,c("country","Sex","Year", "i",
                           "p1", "p2", "p3", "p4", "p5")],
              by = c("country","Sex","Year", "i"))


# Long format

smxr <- hmxr[,c("country", "Sex", "Year", "Age", "delta", "Mx", "p1","p2", "p3", "p4", "p5")]

#smxr <- hmxr # This is the full dataset
smxr <- data.table(arrange(smxr, country, Sex, Year, Age))
# Selection of years and ages to be decomposed ----------------------------

x <- 65 # Age to calculate annuities
Syear <- 2001
Fyear <- 2016# Final year to plot graphs
smxr <- subset(smxr, Age>=x) # Selection of interval


# Actuarial functions, entropy and durations ------------------------------

a <- smxr[, get.ax.cod(Age = Age, mu = Mx, delta = delta,
                       p1 = p1, p2 = p2, p3 = p3, p4 = p4, p5 = p5),
          by = list(country,Sex, Year)]

# Rates of change, rho and phi --------------------------------------------

phi <- get.phi(delta = aRate$delta, Year = aRate$Year)
rho <- a[, get.rho.cod(mu = mu, mu1 = mu1, mu2 = mu2, mu3 = mu3, mu4 = mu4, mu5 = mu5,
                       Year= Year, Age= Age), by = list(country, Sex)]

a <- merge(a, rho, by = c("country", "Sex", "Year", "Age"))
a <- merge(a, phi, by = c("Year"))


# Decomposition of changes over time --------------------------------------

adot <- subset(a)[, get.adot(ax=ax,Year=Year), by = list(country, Sex, Age)]
names(adot)[5] <- "adotc"

# Decomposition
adotD <- a[,get.decomp.cod(Age=Age, ax=ax, sEx=sEx, mu=mu, 
                           mu1 = mu1, mu2 = mu2, mu3 = mu3, mu4 = mu4, mu5 = mu5,
                           delta=delta, rho=rho,
                           rho1=rho1, rho2=rho2, rho3=rho3, rho4=rho4, rho5=rho5,
                           phi=phi, H=Hp,
                           H1=Hp1, H2=Hp2, H3=Hp3, H4=Hp4, H5=Hp5,
                           D=Dp),
             by = list(country,  Sex, Year)]

adotD <- as.data.frame(arrange(adotD, country, Sex, Year))


adotD.cause <- adotD[,c("country","Sex","Year","mort1","mort2","mort3","mort4","mort5",
                        "int")] %>% filter() %>% 
  gather(Component, Percentage, -country, -Sex, -Year)


adotD.cause <- arrange(adotD.cause, Sex,  Year)


# Figures -----------------------------------------------------------------

brks <- seq(1970,2020,by = 10)
labs <- c("1970", "'80","'90","2000","'10","'20")

Syear <- 2001
Fyear <- 2016# Final year to plot graphs


# Colour pallette
ii <- viridis::magma(1,alpha = 1,begin = 0.2, end =0.7,direction = -1)
kk <- viridis::viridis(5, alpha = 1, begin = 0.6, end = 1, direction = -1)
jj<- c("#76b3b8", "#542437", "#C02942", "#ECD078", "#45ba7a")  #c(ii,kk)

colores <- c("#68B3AF","#C02942","#E6AC27", "#A37E58","#254540")



#ggsave("Fig/ProyDescomposition.pdf", width = 4, height = 2, device = cairo_pdf)
  ggplot(subset(adotD.cause, Sex == "Males" & Component != "int" & Year %in% seq(2005,2016, by = 1)))+
  # geom_area( aes( Year, Percentage*100, fill = Component), position = "stack")+
  geom_hline(aes(yintercept = 0), colour = "black", size = 0.1, alpha =  0.5)+
    geom_bar(aes(x= Year, y = Percentage*100, fill = Component),
           position = "stack", stat = "identity",width = 0.7)+
  scale_fill_manual(values = colores)+
  facet_wrap(~Sex, ncol=2)+
  scale_y_continuous(breaks = seq(-1,1,by = 0.5), expand = c(0,0))+
  scale_x_continuous(breaks = seq(2000,2016, by = 5),  expand = c(0,0))+
  theme_minimal()+
  coord_flip(ylim = c(-0.5,1))+
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),#element_rect(fill = "grey30"),
        axis.ticks = element_line(size = 0.1, colour = "gray30"),
        text = element_text(size = 7,font_rc),
        panel.border = element_blank(),
        legend.position = "none",
        aspect.ratio = 1,
        strip.text = element_text(size = 7),
        panel.spacing.x = unit(2,"lines"))+
  ylab("Annual change in the longevity component (%)")+
  xlab("Years")
  
  adotD.cause$Component <- factor(adotD.cause$Component, levels = c("mort5","mort4","mort3","mort2","mort1", "int"))
  
  #ggsave("Fig/ProyDescomposition.pdf", width = 4, height = 2, device = cairo_pdf)
  ggplot(subset(adotD.cause, Sex == "Males" & Component != "int" & Year %in% seq(2005,2016, by = 5)))+
    # geom_area( aes( Year, Percentage*100, fill = Component), position = "stack")+
    geom_hline(aes(yintercept = 0), colour = "black", size = 0.1, alpha =  0.5)+
    geom_bar(aes(x= as.factor(Component), y = Percentage*100, fill = Component),
             position = "stack", stat = "identity",width = 0.8)+
    scale_fill_manual(values = rev(colores))+
    facet_wrap(~Year, ncol=3)+
    #scale_y_continuous(breaks = seq(-1,1,by = 0.5), expand = c(0,0))+
    #scale_x_reverse()+
    theme_minimal()+
    coord_flip(ylim = c(-0.3,0.5))+
    theme(panel.grid = element_blank(),
          panel.background = element_blank(),#element_rect(fill = "grey30"),
          axis.ticks = element_line(size = 0.1, colour = "gray30"),
          text = element_text(size = 7,font_rc),
          panel.border = element_blank(),
          legend.position = "none",
          aspect.ratio = 1,
          strip.text = element_text(size = 7),
          panel.spacing.x = unit(1,"lines"))+
    ylab("Annual change in the longevity component (%)")+
    xlab("")
 
ggsave("Fig/DecCodlabel.pdf", width = 4, height = 2, device = cairo_pdf)

