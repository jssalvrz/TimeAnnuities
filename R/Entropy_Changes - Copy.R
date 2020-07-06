# Calculation of Demographic and Actuarial entropy

library(msm)
library(tidyverse)
library(data.table)
library(readr)
library(ggridges)
library(ungroup)
library(mltools)
library(pracma)
library(zoo)
library(hrbrthemes)


setwd("C:/Users/jmartinez/OneDrive - Syddansk Universitet/Danes Pensions/R")
source("LifeTableFUN.R") 
source("Ages_FUN.R")
source("Var_sim_FUN.R")
source("Hdecomp_FUN.R")
##### 0. Import datasets #####
load("datsims.RData")
lab_years <- c("1985", "'90","'95", "2000", "'05", "'10", "'15")
brk_years <- c(1985, 1990, 1995, 2000, 2005, 2010, 2015)

lab_years10 <- c("1985", "'95", "2005", "'15")
brk_years10 <- c(1985, 1995,  2005, 2015)




# 0. Retrieve target retirement ages --------------------------------------

retirementAges <- data.frame(Year = unique(sims.t$Year),
                             target= c(unique(sims.t$TQx.t)[1:17],68.37,
                                       unique(sims.t$TQx.t)[18:31]),
                             current=c(rep(67,15),rep(65,17))) # The age for 2002 was repeated, that's why I inserted it manually

# 1. lx at different target ages ------------------------------------------

# At age 67 and then at age 65 and then I will merge them
lx.67.f <- sims.f[, calc.lx(lt = lt, Qx  = 67 ), by =list(Year, SES)] # At the target retirement age
lx.67.m <- sims.m[, calc.lx(lt = lt, Qx  = 67 ), by =list(Year, SES)] # At the target retirement age
lx.67.T <- Tsims.t[, calc.lx(lt = lt, Qx = 67 ), by =list(Year)] # At the target retirement age

# At age 65
lx.65.f <- sims.f[, calc.lx(lt = lt, Qx  = 65 ), by =list(Year, SES)] # At the target retirement age
lx.65.m <- sims.m[, calc.lx(lt = lt, Qx  = 65 ), by =list(Year, SES)] # At the target retirement age
lx.65.T <- Tsims.t[, calc.lx(lt = lt, Qx = 65 ), by =list(Year)] # At the target retirement age

# At the target retirement age
lx.target.f <- sims.f[, calc.lx(lt = lt, Qx  = TQx.t ), by =list(Year, SES)] # At the target retirement age
lx.target.m <- sims.m[, calc.lx(lt = lt, Qx  = TQx.t ), by =list(Year, SES)] # At the target retirement age
lx.target.T <- Tsims.t[, calc.lx(lt = lt, Qx = TQx.t ), by =list(Year)] # At the target retirement age

# Get sexes
lx.67.f$Sex <- "Females"
lx.67.m$Sex <- "Males"
lx.67.T$Sex <- "Total"
lx.67.T$SES <- "Total"
lx.67.T <- lx.67.T[,c(1,7,2,3,4,5,6)]
lx.67 <- rbind(lx.67.f, lx.67.m, lx.67.T)


lx.65.f$Sex <- "Females"
lx.65.m$Sex <- "Males"
lx.65.T$Sex <- "Total"
lx.65.T$SES <- "Total"
lx.65.T <- lx.65.T[,c(1,7,2,3,4,5,6)]
lx.65 <- rbind(lx.65.f, lx.65.m, lx.65.T)

lx.target.f$Sex <- "Females"
lx.target.m$Sex <- "Males"
lx.target.T$Sex <- "Total"
lx.target.T$SES <- "Total"
lx.target.T <- lx.target.T[,c(1,7,2,3,4,5,6)]
lx.target <- rbind(lx.target.f, lx.target.m, lx.target.T)


# 2. Demographic values ---------------------------------------------------

delta <- 0 # interest rate

H.67  <- lx.67[, calc.aH(Age = Age, lx = lx, delta = delta ), by =list(Year, SES, Sex)]
H.65  <- lx.65[, calc.aH(Age = Age, lx = lx, delta = delta ), by =list(Year, SES, Sex)]
H.target  <- lx.target[, calc.aH(Age = Age, lx = lx, delta = delta ), by =list(Year, SES, Sex)]
H.current <- rbind(subset(H.67, Year < 2000), subset(H.65, Year>=2000))

H.target.T  <- subset(H.target, SES == "Total")[,c(1,4,5,6)]
H.current.T <- subset(H.current, SES == "Total")[,c(1,4,5,6)]
names(H.target.T) <- c("Year","T.ax","T.a.dagger", "T.H.a")
names(H.current.T) <- c("Year","T.ax","T.a.dagger", "T.H.a")

H.target <- merge(H.target, H.target.T, by = "Year")
H.current <- merge(H.current, H.current.T, by = "Year")

H.target$Diff.ax      <-(H.target$ax - H.target$T.ax) / H.target$T.ax
H.target$Diff.a.dagger<-(H.target$a.dagger - H.target$T.a.dagger) / H.target$T.a.dagger
H.target$Diff.H.a     <- (H.target$H.a - H.target$T.H.a) / H.target$T.H.a
H.target <- arrange(H.target, Sex, SES, Year)

H.current$Diff.ax      <- (H.current$ax - H.current$T.ax) / H.current$T.ax
H.current$Diff.a.dagger<- (H.current$a.dagger - H.current$T.a.dagger) /  H.current$T.a.dagger
H.current$Diff.H.a     <- (H.current$H.a - H.current$T.H.a) / H.current$T.H.a
H.current <- arrange(H.current, Sex, SES, Year)



# 3. Actuarial values, varying delta --------------------------------------

# 1%

H.67.1  <- lx.67[, calc.aH(Age = Age, lx = lx, delta = 0.01 ), by =list(Year, SES, Sex)]
H.65.1  <- lx.65[, calc.aH(Age = Age, lx = lx, delta = 0.01 ), by =list(Year, SES, Sex)]
H.target.1  <- lx.target[, calc.aH(Age = Age, lx = lx, delta = 0.01 ), by =list(Year, SES, Sex)]
H.current.1 <- rbind(subset(H.67.1, Year < 2000), subset(H.65.1, Year>=2000))

H.target.T.1  <- subset(H.target.1, SES == "Total")[,c(1,4,5,6)]
H.current.T.1 <- subset(H.current.1, SES == "Total")[,c(1,4,5,6)]
names(H.target.T.1) <- c("Year","T.ax","T.a.dagger", "T.H.a")
names(H.current.T.1) <- c("Year","T.ax","T.a.dagger", "T.H.a")

H.target.1 <- merge(H.target.1, H.target.T.1, by = "Year")
H.current.1 <- merge(H.current.1, H.current.T.1, by = "Year")

H.target.1$Diff.ax      <-(H.target.1$ax - H.target.1$T.ax) / H.target.1$T.ax
H.target.1$Diff.a.dagger<-(H.target.1$a.dagger - H.target.1$T.a.dagger) / H.target.1$T.a.dagger
H.target.1$Diff.H.a     <- (H.target.1$H.a - H.target.1$T.H.a) / H.target.1$T.H.a
H.target.1 <- arrange(H.target.1, Sex, SES, Year)

H.current.1$Diff.ax      <- (H.current.1$ax - H.current.1$T.ax) / H.current.1$T.ax
H.current.1$Diff.a.dagger<- (H.current.1$a.dagger - H.current.1$T.a.dagger) /  H.current.1$T.a.dagger
H.current.1$Diff.H.a     <- (H.current.1$H.a - H.current.1$T.H.a) / H.current.1$T.H.a
H.current.1 <- arrange(H.current.1, Sex, SES, Year)

# 5%

H.67.5  <- lx.67[, calc.aH(Age = Age, lx = lx, delta = 0.05 ), by =list(Year, SES, Sex)]
H.65.5  <- lx.65[, calc.aH(Age = Age, lx = lx, delta = 0.05 ), by =list(Year, SES, Sex)]
H.target.5  <- lx.target[, calc.aH(Age = Age, lx = lx, delta = 0.05 ), by =list(Year, SES, Sex)]
H.current.5 <- rbind(subset(H.67.5, Year < 2000), subset(H.65.5, Year>=2000))

H.target.T.5  <- subset(H.target.5, SES == "Total")[,c(1,4,5,6)]
H.current.T.5 <- subset(H.current.5, SES == "Total")[,c(1,4,5,6)]
names(H.target.T.5) <- c("Year","T.ax","T.a.dagger", "T.H.a")
names(H.current.T.5) <- c("Year","T.ax","T.a.dagger", "T.H.a")

H.target.5 <- merge(H.target.5, H.target.T.5, by = "Year")
H.current.5 <- merge(H.current.5, H.current.T.5, by = "Year")

H.target.5$Diff.ax      <-(H.target.5$ax - H.target.5$T.ax) / H.target.5$T.ax
H.target.5$Diff.a.dagger<-(H.target.5$a.dagger - H.target.5$T.a.dagger) / H.target.5$T.a.dagger
H.target.5$Diff.H.a     <- (H.target.5$H.a - H.target.5$T.H.a) / H.target.5$T.H.a
H.target.5 <- arrange(H.target.5, Sex, SES, Year)

H.current.5$Diff.ax      <- (H.current.5$ax - H.current.5$T.ax) / H.current.5$T.ax
H.current.5$Diff.a.dagger<- (H.current.5$a.dagger - H.current.5$T.a.dagger) /  H.current.5$T.a.dagger
H.current.5$Diff.H.a     <- (H.current.5$H.a - H.current.5$T.H.a) / H.current.5$T.H.a
H.current.5 <- arrange(H.current.5, Sex, SES, Year)


# 10%

H.67.10  <- lx.67[, calc.aH(Age = Age, lx = lx, delta = 0.1 ), by =list(Year, SES, Sex)]
H.65.10  <- lx.65[, calc.aH(Age = Age, lx = lx, delta = 0.1 ), by =list(Year, SES, Sex)]
H.target.10  <- lx.target[, calc.aH(Age = Age, lx = lx, delta = 0.1 ), by =list(Year, SES, Sex)]
H.current.10 <- rbind(subset(H.67.10, Year < 2000), subset(H.65.10, Year>=2000))

H.target.T.10  <- subset(H.target.10, SES == "Total")[,c(1,4,5,6)]
H.current.T.10 <- subset(H.current.10, SES == "Total")[,c(1,4,5,6)]
names(H.target.T.10) <- c("Year","T.ax","T.a.dagger", "T.H.a")
names(H.current.T.10) <- c("Year","T.ax","T.a.dagger", "T.H.a")

H.target.10 <- merge(H.target.10, H.target.T.10, by = "Year")
H.current.10 <- merge(H.current.10, H.current.T.10, by = "Year")

H.target.10$Diff.ax      <-(H.target.10$ax - H.target.10$T.ax) / H.target.10$T.ax
H.target.10$Diff.a.dagger<-(H.target.10$a.dagger - H.target.10$T.a.dagger) / H.target.10$T.a.dagger
H.target.10$Diff.H.a     <- (H.target.10$H.a - H.target.10$T.H.a) / H.target.10$T.H.a
H.target.10 <- arrange(H.target.10, Sex, SES, Year)

H.current.10$Diff.ax      <- (H.current.10$ax - H.current.10$T.ax) / H.current.10$T.ax
H.current.10$Diff.a.dagger<- (H.current.10$a.dagger - H.current.10$T.a.dagger) /  H.current.10$T.a.dagger
H.current.10$Diff.H.a     <- (H.current.10$H.a - H.current.10$T.H.a) / H.current.10$T.H.a
H.current.10 <- arrange(H.current.10, Sex, SES, Year)


# 5. Decomposition of H ---------------------------------------------------

rH <- data.table(H.target.5)
rax <- rH[, Roller(Rx= ax, Year = Year, fn = mean,
                                  nrol = 10, skip = 1),
                          by = list(Sex, SES)]

ra.dag <- rH[, Roller(Rx= a.dagger, Year = Year, fn = mean,
                                  nrol = 10, skip = 1),
                         by = list(Sex, SES)]

rH.a <- rH[, Roller(Rx= H.a, Year = Year, fn = mean,
                                  nrol = 10, skip = 1),
                            by = list(Sex, SES)]

Sx <- "Females"
ax1    <- subset(rax,    Sex == Sx & SES == 0)
ax2    <- subset(rax,    Sex == Sx & SES == 4)
a.dag1 <- subset(ra.dag, Sex == Sx & SES == 0)
a.dag2 <- subset(ra.dag, Sex == Sx & SES == 4)
H1     <- subset(rH.a,   Sex == Sx & SES == 0)
H2     <- subset(rH.a,   Sex == Sx & SES == 4)

decompEntropy <- H.decomp(Year = ax1$Year.rol, ax1 =  ax1$rol,
                                               ax2 =  ax2$rol,
                                         a.dagger1 =  a.dag1$rol,
                                         a.dagger2 =  a.dag2$rol,
                                              H.a1 =  H1$rol,
                                              H.a2 =  H2$rol)


ggplot(decompEntropy, aes(x = Year, y = Percentage, fill = Component))+
 # geom_hline(yintercept = 0,size = 0.2)+
  geom_bar(stat= "identity", position = "stack")+
  #coord_cartesian()+
  coord_cartesian(ylim = c(-0.02,0.02))+
  scale_x_continuous(expand=c(0,0), breaks = c(1986,2007), labels = c(1985,2016))+
  scale_y_continuous(expand=c(0,0), breaks =c(-0.02,0,0.02), sec.axis = dup_axis())+
  scale_fill_manual(values = c("#f2843a","#5eb5b1"))+
  theme_bw()+
  theme(strip.background = element_rect(fill="none"))+
  theme(#legend.title=element_blank(), 
    strip.background = element_blank(),
    panel.background = element_rect(fill = "purple"), #e7ebab
    axis.text.x  = element_text( hjust = 0.5),
    panel.spacing = unit(2, "lines"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line = element_blank(),
    axis.title.y = element_text(vjust = 2),
    axis.title.x = element_text(vjust = -1),
    # axis.text.y.left = element_blank(),
    # axis.ticks.y.left = element_blank(),
     axis.text.y.right = element_blank(),
     axis.ticks.y.right = element_blank(),
    aspect.ratio = 1,
    axis.ticks = element_line(size = 0.3, colour = "gray70"),
    text = element_text(size = 5,font_rc),
    legend.position = "none",
    plot.background = element_rect(fill = NA))+
  xlab("") +
  ylab("")+
  geom_text(aes(label = "\u03b4=5%"), 
                     x = 2003, y = -0.015, colour = "white", size = 2) 





ggsave("fig/Tf5.pdf", width = 2, height = 1, device = cairo_pdf)


# 6. Graphs for e(x) and H(x) ---------------------------------------------


# Target and current retirement ages

ggplot(retirementAges)+
  geom_line(aes(Year, current), size = 1.5, colour = "#f1ffa3")+
  geom_line(aes(Year, target), size = 1.5, colour = "purple4")+
  coord_cartesian(ylim= c(64,72))+
  scale_x_continuous(expand=c(0,0), breaks = brk_years, labels = lab_years)+
  scale_y_continuous(expand=c(0,0), breaks =c(seq(64,75, by = 1)))+
  #theme_minimal()+
  #theme_ft_rc()+
  #theme(strip.background = element_rect(fill="none"))+
  theme(#legend.title=element_blank(), 
    #strip.background = element_blank(),
    panel.background = element_rect(fill = "gray60"), 
    #panel.background = element_blank(),
    panel.border = element_blank(),
    axis.text.x  = element_text( hjust = 0.5),
    panel.spacing = unit(2, "lines"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_blank(),
    axis.title.y = element_text(vjust = 2),
    axis.title.x = element_text(vjust = -1),
    aspect.ratio = 1,
    text = element_text(size = 10,font_rc),
    strip.text.x = element_text(size = 18, colour = "black"),
    legend.position = "right",
    plot.background = element_rect(fill = NA))+
  xlab("Year") +
  ylab("Retirement age (in years)")

#ggsave("fig/Retirement ages.pdf", width = 2, height = 2, device = cairo_pdf)

# Life expectancy for the total population
#"lightsalmon2" for e(x) and a(x) and "lightseagreen" for H(x)
ggplot(H.current.T)+
  geom_line(aes(Year, T.ax), size = 1.2, colour = "#f1ffa3")+
  geom_hline(yintercept = 14.5, size = 1.2, colour = "purple4")+
  coord_cartesian(ylim= c(14,20))+
  scale_x_continuous(expand=c(0,0), breaks = brk_years, labels = lab_years)+
  scale_y_continuous(expand=c(0,0), breaks =c(seq(13,30, by = 1)))+
  #theme_minimal()+
  #theme(strip.background = element_rect(fill="none"))+
  theme(#legend.title=element_blank(), 
    #strip.background = element_blank(),
    panel.background = element_rect(fill = "#ff9d7d"), #"lightsalmon3""darkkhaki"
    #panel.background = element_blank(),
    #panel.border = element_rect(colour = "gray90"),
    axis.text.x  = element_text( hjust = 0.5),
    panel.spacing = unit(2, "lines"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_blank(),
    axis.title.y = element_text(vjust = 2),
    axis.title.x = element_text(vjust = -1),
    aspect.ratio = 1,
    text = element_text(size = 10,font_rc),
    strip.text.x = element_text(size = 18, colour = "black"),
    legend.position = "right",
    plot.background = element_rect(fill = NA))+
  xlab("Year") +
  ylab("Life expectancy")

#ggsave("fig/Life expectancy.pdf", width = 2, height = 2, device = cairo_pdf)

ggplot(H.current.T)+
  geom_line(data = H.current.T, aes(Year, T.H.a), colour = "#f1ffa3",size = 1.2)+
  geom_line(data = H.target.T, aes(Year, T.H.a), colour = "purple4", size = 1.2)+
  coord_cartesian(ylim= c(0.38,0.52))+
  scale_x_continuous(expand=c(0,0), breaks = brk_years, labels = lab_years)+
  scale_y_continuous(expand=c(0,0), breaks =c(seq(0.30,1, by = 0.02)))+
  #theme_minimal()+
  #theme(strip.background = element_rect(fill="none"))+
  theme(#legend.title=element_blank(), 
    #strip.background = element_blank(),
    panel.background = element_rect(fill = "lightseagreen"), 
    #panel.background = element_blank(),
    #panel.border = element_rect(colour = "gray90"),
    axis.text.x  = element_text( hjust = 0.5),
    panel.spacing = unit(2, "lines"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_blank(),
    axis.title.y = element_text(vjust = 2),
    axis.title.x = element_text(vjust = -1),
    aspect.ratio = 1,
    text = element_text(size = 10,font_rc),
    strip.text.x = element_text(size = 18, colour = "black"),
    legend.position = "right",
    plot.background = element_rect(fill = NA))+
  xlab("Year") +
  ylab("Lifetable entropy")

#ggsave("fig/Total entropy.pdf", width = 2, height = 2, device = cairo_pdf)


# Life expectancy by SES current retirement age

ggplot(subset(H.current, SES != "Total"))+
  geom_line(aes(Year, ax, colour = as.character(SES)), size = 0.8)+
  coord_cartesian(ylim= c(10,24))+
  facet_wrap(~Sex)+
  scale_x_continuous(expand=c(0,0), breaks = brk_years10, labels = lab_years10)+
  scale_y_continuous(expand=c(0,0), breaks =c(seq(10,24, by = 2)))+
  # scale_color_viridis_d(direction = -1,
  #                       option = "D",labels = c("First", "Second", "Third", "Fourth", "Fifth"))+
  scale_colour_grey(start = 1, end = 0, na.value = "red",
                  labels = c("First", "Second", "Third", "Fourth", "Fifth"),
                  aesthetics = "colour")+
  theme_bw()+
   theme(#legend.title=element_blank(), 
     strip.background = element_blank(),
     panel.background = element_rect(fill = "#ff9d7d"), #"lightsalmon3""darkkhaki""darkcyan
     #panel.background = element_blank(),
     panel.border = element_blank(),
     axis.text.x  = element_text( hjust = 0.5),
     panel.spacing = unit(2, "lines"),
     panel.grid.major = element_blank(),
     panel.grid.minor = element_blank(),
     axis.line = element_blank(),
     axis.title.y = element_text(vjust = 2),
     axis.title.x = element_text(vjust = -1),
     aspect.ratio = 1,
     text = element_text(size = 10,font_rc),
     strip.text.x = element_text(size = 18, colour = "black"),
     legend.position = "none",
     plot.background = element_rect(fill = NA))+
  xlab("Year") +
  ylab("at age c" )+#H(r,\u03b4) with \u03b4 = 1%
  labs(colour = "SES group")+
  geom_text(data = . %>% distinct(Sex) %>% 
              mutate(over = c("A", "B")),
            aes(label = over), 
            x = 1988, y = 23, size = 6, colour = "black")


#ggsave("fig/Life expectancy SES current.pdf", width = 4, height = 2, device = cairo_pdf)


# Life expectancy target retirement age

ggplot(subset(H.target, SES != "Total"))+
  #geom_line(data= H.current.T,aes(Year, T.ax), colour = "firebrick", size = 2, alpha = 0.8)+
  geom_line(aes(Year, ax, colour = as.character(SES)), size = 0.8)+
  coord_cartesian(ylim= c(10,24))+
  facet_wrap(~Sex)+
  scale_x_continuous(expand=c(0,0), breaks = brk_years10, labels = lab_years10)+
  scale_y_continuous(expand=c(0,0), breaks =c(seq(10,24, by = 2)))+
  # scale_color_viridis_d(direction = -1,
  #                       option = "D",labels = c("First", "Second", "Third", "Fourth", "Fifth"))+
  scale_colour_grey(start = 1, end = 0, na.value = "red",
                    labels = c("First", "Second", "Third", "Fourth", "Fifth"),
                    aesthetics = "colour")+
  theme_bw()+
  theme(#legend.title=element_blank(), 
    strip.background = element_blank(),
    panel.background = element_rect(fill = "#ff9d7d"), #"lightsalmon3""darkkhaki""darkcyan
    #panel.background = element_blank(),
    panel.border = element_blank(),
    axis.text.x  = element_text( hjust = 0.5),
    panel.spacing = unit(2, "lines"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_blank(),
    axis.title.y = element_text(vjust = 2),
    axis.title.x = element_text(vjust = -1),
    aspect.ratio = 1,
    text = element_text(size = 10,font_rc),
    strip.text.x = element_text(size = 18, colour = "black"),
    legend.position = "right",
    plot.background = element_rect(fill = NA))+
  xlab("Year") +
  ylab("at age t" )+#H(r,\u03b4) with \u03b4 = 1%
  labs(colour = "SES group")+
  geom_text(data = . %>% distinct(Sex) %>% 
              mutate(over = c("C", "D")),
            aes(label = over), 
            x = 1988, y = 23, size = 6, colour = "black")

#ggsave("fig/SES.pdf", width = 4, height = 2, device = cairo_pdf)


# Entropy at the current retirement age

ggplot(subset(H.current, SES != "Total"))+
  geom_line(aes(Year, H.a, colour = as.character(SES)), size = 0.8)+
  coord_cartesian(ylim= c(0.3,0.65))+
  facet_wrap(~Sex)+
  scale_x_continuous(expand=c(0,0), breaks = brk_years10, labels = lab_years10)+
  scale_y_continuous(expand=c(0,0), breaks =c(seq(0,1, by = 0.1)))+
  scale_colour_grey(start = 1, end = 0, na.value = "red",
                    labels = c("First", "Second", "Third", "Fourth", "Fifth"),
                    aesthetics = "colour")+
  theme_bw()+
  theme(#legend.title=element_blank(), 
    strip.background = element_blank(),
    panel.background = element_rect(fill = "lightseagreen"), #"lightsalmon3""darkkhaki""darkcyan
    #panel.background = element_blank(),
    panel.border = element_blank(),
    axis.text.x  = element_text( hjust = 0.5),
    panel.spacing = unit(2, "lines"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_blank(),
    axis.title.y = element_text(vjust = 2),
    axis.title.x = element_text(vjust = -1),
    aspect.ratio = 1,
    text = element_text(size = 10,font_rc),
    strip.text.x = element_text(size = 18, colour = "black"),
    legend.position = "none",
    plot.background = element_rect(fill = NA))+
  xlab("Year") +
  ylab("at age t" )+#H(r,\u03b4) with \u03b4 = 1%
  labs(colour = "SES group")+
  geom_text(data = . %>% distinct(Sex) %>% 
              mutate(over = c("A", "B")),
            aes(label = over), 
            x = 1987, y = 0.61, size = 6, colour = "black")

#ggsave("fig/Lifespan inequality SES current.pdf", width = 4, height = 2, device = cairo_pdf)

# Entropy at the target retirement age

ggplot(subset(H.target, SES != "Total"))+
  geom_line(aes(Year, H.a, colour = as.character(SES)), size = 0.8)+
  coord_cartesian(ylim= c(0.3,0.65))+
  facet_wrap(~Sex)+
  scale_x_continuous(expand=c(0,0), breaks = brk_years10, labels = lab_years10)+
  scale_y_continuous(expand=c(0,0), breaks =c(seq(0,1, by = 0.1)))+
  scale_colour_grey(start = 1, end = 0, na.value = "red",
                    labels = c("First", "Second", "Third", "Fourth", "Fifth"),
                    aesthetics = "colour")+
  theme_bw()+
  theme(#legend.title=element_blank(), 
    strip.background = element_blank(),
    panel.background = element_rect(fill = "lightseagreen"), #"lightsalmon3""darkkhaki""darkcyan
    #panel.background = element_blank(),
    panel.border = element_blank(),
    axis.text.x  = element_text( hjust = 0.5),
    panel.spacing = unit(2, "lines"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_blank(),
    axis.title.y = element_text(vjust = 2),
    axis.title.x = element_text(vjust = -1),
    aspect.ratio = 1,
    text = element_text(size = 10,font_rc),
    strip.text.x = element_text(size = 18, colour = "black"),
    legend.position = "none",
    plot.background = element_rect(fill = NA))+
  xlab("Year") +
  ylab("at age t" )+#H(r,\u03b4) with \u03b4 = 1%
  labs(colour = "SES group")+
  geom_text(data = . %>% distinct(Sex) %>% 
              mutate(over = c("C", "D")),
            aes(label = over), 
            x = 1987, y = 0.61, size = 6, colour = "black")
#ggsave("fig/Lifespan inequality SES target.pdf", width = 4, height = 2, device = cairo_pdf)




# 4. Differences in a(x) and H(x) -----------------------------------------

# Difference in ax
Hh <- unique(subset(H.target, SES != "Total")[c("Sex")])


ggplot(subset(H.target.5, SES != "Total"))+
  # geom_rect(data = Hh,aes(fill = Sex),xmin = -Inf,xmax = Inf,
  #           ymin = -Inf,ymax = Inf,alpha = 0.2) +
  geom_hline(yintercept = 0,  size = 0.3, colour = "purple4", linetype = "dashed")+
  geom_line(aes(Year, Diff.ax*100, colour = as.character(SES)), size =0.3)+
  coord_cartesian(ylim= c(-40,45))+
  facet_wrap(~Sex)+
  scale_x_continuous(expand=c(0,0), breaks = c(1985,2016))+
  scale_y_continuous(expand=c(0,0), breaks =c(seq(-40,40, by = 20)), sec.axis = dup_axis())+
  scale_colour_grey(start = 1, end = 0, na.value = "red",
                    labels = c("First", "Second", "Third", "Fourth", "Fifth"),
                    aesthetics = "colour")+
  theme_bw()+
  theme(#legend.title=element_blank(), 
    strip.background = element_blank(),
    panel.background = element_rect(fill = "#ff9d7d"), #"lightsalmon3""darkkhaki""darkcyan
    #panel.background = element_blank(),
    panel.border = element_blank(),
    axis.text.x  = element_text( hjust = 0.5),
    axis.text.y.left = element_blank(),
    axis.ticks.y.left = element_blank(),
    panel.spacing = unit(0.8, "lines"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_blank(),
    axis.title.y = element_text(vjust = 2),
    axis.title.x = element_text(vjust = -1),
    aspect.ratio = 1,
    axis.ticks = element_line(size = 0.3, colour = "gray70"),
    text = element_text(size = 5,font_rc),
    strip.text.x = element_blank(),
    legend.position = "none",
    plot.background = element_rect(fill = NA))+
  xlab("") +
  ylab("" )+#H(r,\u03b4) with \u03b4 = 1%
  labs(colour = "SES group")+
  geom_text(data = . %>% distinct(Sex) %>% 
              mutate(over = c("", "")),
            aes(label = over), 
            x = 1989, y = 38, size = 2)+
  geom_text(data = . %>% distinct(Sex) %>% 
              mutate(over = c("\u03b4=5%", "\u03b4=5%")),
            aes(label = over), 
            x = 2009, y = 38, size = 2) 
  
ggsave("fig/at5.pdf", width = 2, height = 1, device = cairo_pdf)

# Difference in H(x)

ggplot(subset(H.target, SES != "Total"))+
  geom_hline(yintercept = 0,  size = 0.3, colour = "purple4")+
  geom_line(aes(Year, Diff.H.a*100, colour = as.character(SES)), size =0.3)+
  coord_cartesian(ylim= c(-40,45))+
  facet_wrap(~Sex)+
  scale_x_continuous(expand=c(0,0), breaks = c(1985,2016))+
  scale_y_continuous(expand=c(0,0), breaks =c(seq(-40,40, by = 20)), sec.axis = dup_axis())+
  scale_colour_grey(start = 1, end = 0, na.value = "red",
                    labels = c("First", "Second", "Third", "Fourth", "Fifth"),
                    aesthetics = "colour")+
  theme_bw()+
  theme(#legend.title=element_blank(), 
    strip.background = element_blank(),
    panel.background = element_rect(fill = "lightseagreen"), #"lightsalmon3""darkkhaki""darkcyan
    panel.border = element_blank(),
    axis.text.x  = element_text( hjust = 0.5),
    axis.text.y.left = element_blank(),
    axis.ticks.y.left = element_blank(),
    panel.spacing = unit(0.8, "lines"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_blank(),
    axis.title.y = element_text(vjust = 2),
    axis.title.x = element_text(vjust = -1),
    aspect.ratio = 1,
    axis.ticks = element_line(size = 0.3, colour = "gray70"),
    text = element_text(size = 5,font_rc),
    strip.text.x = element_blank(),
    legend.position = "none",
    plot.background = element_rect(fill = NA))+
  xlab("") +
  ylab("" )+#H(r,\u03b4) with \u03b4 = 1%
  labs(colour = "SES group")+
  geom_text(data = . %>% distinct(Sex) %>% 
              mutate(over = c("", "")),
            aes(label = over), 
            x = 1989, y = 38, size = 2)+
  geom_text(data = . %>% distinct(Sex) %>% 
              mutate(over = c("\u03b4=0%", "\u03b4=0%")),
            aes(label = over), 
            x = 2009, y = 38, size = 2) 


ggsave("fig/Ht0.pdf", width = 2, height = 1, device = cairo_pdf)



# 5. Values of a(x) and H(x) for the total population ---------------------

ggplot()+
  geom_line(data = H.target.T.1, aes(Year, T.ax),colour = "purple4", size = 0.8)+
  #geom_line(data = H.current.T.10,aes(Year, T.ax),colour = "firebrick1", size = 1.2, alpha = 1)+
  coord_cartesian(ylim= c(0,20))+
  scale_x_continuous(expand=c(0,0), breaks = c(1985,2016))+
  scale_y_continuous(expand=c(0,0), breaks =c(seq(0,20, by = 5)))+
  theme_bw()+
  theme(strip.background = element_rect(fill="none"))+
  theme(#legend.title=element_blank(), 
    strip.background = element_blank(),
    panel.background = element_rect(fill = "#ff9d7d"), #"lightsalmon3""darkkhaki""darkcyan
    axis.text.x  = element_text( hjust = 0.5),
    panel.spacing = unit(2, "lines"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line = element_blank(),
    axis.title.y = element_text(vjust = 2),
    axis.title.x = element_text(vjust = -1),
    aspect.ratio = 1,
    axis.ticks = element_line(size = 0.3, colour = "gray70"),
    text = element_text(size = 5,font_rc),
    legend.position = "right",
    plot.background = element_rect(fill = NA))+
  xlab("") +
  ylab("")+geom_text(aes(label = "\u03b4=1%"), 
                     x = 1992, y = 18, size = 2) 


ggsave("fig/Dac1.pdf", width = 1, height = 1, device = cairo_pdf)

ggplot()+
  geom_line(data = H.target.T, aes(Year, T.H.a),colour = "purple4", size = 0.8)+
  #geom_line(data = H.current.T.10,aes(Year, T.ax),colour = "firebrick1", size = 1.2, alpha = 1)+
  coord_cartesian(ylim= c(0.2,0.52))+
  scale_x_continuous(expand=c(0,0), breaks = c(1985,2016))+
  scale_y_continuous(expand=c(0,0), breaks =c(seq(0,1, by = 0.1)))+
  theme_bw()+
  theme(strip.background = element_rect(fill="none"))+
  theme(#legend.title=element_blank(), 
    strip.background = element_blank(),
    panel.background = element_rect(fill = "lightseagreen"), #"lightsalmon3""darkkhaki""darkcyan
    axis.text.x  = element_text( hjust = 0.5),
    panel.spacing = unit(2, "lines"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line = element_blank(),
    axis.title.y = element_text(vjust = 2),
    axis.title.x = element_text(vjust = -1),
    aspect.ratio = 1,
    axis.ticks = element_line(size = 0.3, colour = "gray70"),
    text = element_text(size = 5,font_rc),
    legend.position = "right",
    plot.background = element_rect(fill = NA))+
  xlab("") +
  ylab("")+geom_text(aes(label = "\u03b4=0%"), 
                     x = 2009, y = 0.5, size = 2) 


ggsave("fig/DHt0.pdf", width = 1, height = 1, device = cairo_pdf)
