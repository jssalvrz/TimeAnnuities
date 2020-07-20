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

adot <- subset(a)[, get.adot(ax=ax,Year=Year), by = list(country, Sex, Age)]
names(adot)[5] <- "adotc"
adotD <- a[,get.decomp(Age=Age, ax=ax, sEx=sEx, mu=mu,
                         delta=delta, rho=rho, phi=phi, H=Hp,D=Dp,Dc=Dc),
             by = list(country,  Sex, Year)]

adotD <- as.data.frame(arrange(adotD, country, Sex, Year))


adotD$mortp <- adotD$mort / adotD$adot
adotD$intp <- adotD$int / adotD$adot


# Check constant vs proportional changes in interest rates ----------------


adotD <- merge(adotD, aRate, by = "Year")
adotD <- merge(adotD, phi, by = "Year")


adotD$intc <- adotD$phibar * adotD$sdelta *adotD$Dc
adotD$Dp2 <- adotD$sdelta *adotD$Dc
# Figures -----------------------------------------------------------------
brks <- seq(1860,2020,by = 20)
labs <- c("1860", "'80", "1900","'20","'40","'60","'80","2000", "'20")



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
  coord_cartesian(ylim = c(2,16), xlim=c(1850,2020))+
  theme(panel.grid = element_blank(),
        panel.background = element_rect(fill = "grey30"),
        axis.ticks = element_line(size = 0.1, colour = "gray30"),
        text = element_text(size = 7,font_rc),
        panel.border = element_blank(),
        aspect.ratio = .5,
        strip.text = element_text(size = 7),
        legend.position = "none",
        panel.spacing.x = unit(1,"lines"))

#ggsave("Fig/Ann.pdf", width = 4, height = 2, device = cairo_pdf)

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
  coord_cartesian(ylim = c(2,16), xlim=c(1850,2020))+
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
  coord_cartesian(ylim = c(-6,6), xlim=c(1850,2020))+
  theme(panel.grid = element_blank(),
        panel.background = element_rect(fill = "grey30"),
        axis.ticks = element_line(size = 0.1, colour = "gray30"),
        text = element_text(size = 7,font_rc),
        panel.border = element_blank(),
        aspect.ratio = .5,
        strip.text = element_text(size = 7),
        legend.position = "none",
        panel.spacing.x = unit(1,"lines"))


ggsave("Fig/adot.pdf", width = 4, height = 2, device = cairo_pdf)

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
  coord_cartesian(ylim = c(0,0.8), xlim=c(1850,2020))+
  theme(panel.grid = element_blank(),
        panel.background = element_rect(fill = "grey30"),
        axis.ticks = element_line(size = 0.1, colour = "gray30"),
        text = element_text(size = 7,font_rc),
        panel.border = element_blank(),
        aspect.ratio = .5,
        strip.text = element_text(size = 7),
        panel.spacing.x = unit(1,"lines"))


#ggsave("Fig/Sens75.pdf", width = 4, height = 2, device = cairo_pdf)


# Derivative of ax
ggplot(subset(adot, Age %in% c(65,70,75)))+
  #geom_line(aes(Year, adot), colour = "blue")+ # decomposition
  geom_line(aes(Year, adotc*100, group = Age), colour = "white")+ # real one
  facet_wrap(~Sex)+
  ylab("Change over time in life annuity factors (%)")+
  scale_y_continuous(breaks = seq(-6,6,by = 2), expand = c(0,0))+
  scale_x_continuous(breaks = brks,labels = labs, expand = c(0,0))+
  theme_minimal()+
  coord_cartesian(ylim = c(-6,6), xlim=c(1850,2020))+
  theme(panel.grid = element_blank(),
        panel.background = element_rect(fill = "grey30"),
        axis.ticks = element_line(size = 0.1, colour = "gray30"),
        text = element_text(size = 7,font_rc),
        panel.border = element_blank(),
        aspect.ratio = .5,
        strip.text = element_text(size = 7),
        panel.spacing.x = unit(1,"lines"))



# Rates of change over time

adotD <- merge(adotD, aRate, by = "Year")


ggplot(adotD)+
  geom_hline(yintercept = 0, colour = "white", alpha = 0.7, size = 0.1)+
  geom_line(aes(Year, rhobar*100), colour = "#8bf7ab", size = 0.3)+
  geom_line(aes(Year, -phibar*sdelta*1000), colour = "#f5a6ff", size = 0.3)+
  #geom_line(aes(Year, phibar), colour = "white")+
  facet_wrap(~Sex)+
  scale_color_viridis_c()+
  scale_y_continuous(breaks = seq(-10,10,by = 2), expand = c(0,0))+
  scale_x_continuous(breaks = brks,labels = labs, expand = c(0,0))+
  theme_minimal()+
  coord_cartesian(ylim = c(-6,10), xlim=c(1842,2020))+
  theme(panel.grid = element_blank(),
        panel.background = element_rect(fill = "grey30"),
        axis.ticks = element_line(size = 0.1, colour = "gray30"),
        text = element_text(size = 7,font_rc),
        panel.border = element_blank(),
        aspect.ratio = .5,
        strip.text = element_text(size = 7),
        panel.spacing.x = unit(1,"lines"))+
  ylab("Average change (%)")


ggsave("Fig/Fig3.pdf", width = 4, height = 2, device = cairo_pdf)


# Decomposition
adotD[,c(1,2,3,8,9)] %>% filter() %>% 
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
  coord_cartesian(ylim = c(-4,4), xlim=c(1851,2020))+
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

  
ggsave("Fig/D70.pdf", width = 4, height = 2, device = cairo_pdf)
  
#ggsave("Fig/Decomposition.pdf", width = 6, height = 4, device = cairo_pdf)

