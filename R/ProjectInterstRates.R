###############################################################################
# Project interest rate data using CIR model
# Period
# Time annuities
###############################################################################


#Load data
lngRate <- read_csv("./Data/Interest rates/UKlongtermrate.csv")
lngRate <- lngRate %>%
  mutate(date = dmy(paste("1", lngRate$Month, lngRate$Year, sep = "-")))


#Estimate parameter of CIR model
YearsStart <- 1753
YearsProj <- 50

dataEst <- lngRate %>% 
  filter(Year >= YearsStart) %>%
  mutate(t = Year + month(date)/12,
         rate = rate/100, 
         Y = (lead(rate) -rate)/sqrt(rate), 
         X1 = (lead(t) - t)/sqrt(rate),
         X2 = -(lead(t) - t)*sqrt(rate))

CIRreg <- lm(Y ~ -1  + X1 + X2, data = dataEst)

aCIR <-  CIRreg$coefficients[2]
bCIR <-  CIRreg$coefficients[1]/CIRreg$coefficients[2]


#Projection
r0 <- tail(lngRate$rate,1)/100


tProj <- seq(1/12,YearsProj, 1/12)
rProj <- r0*exp(-aCIR*tProj) + bCIR*(1 - exp(-aCIR*tProj))
forRate <- tibble(rate = rProj*100, 
                  date = add_with_rollback(tail(lngRate$date, 1), months(seq(1,YearsProj*12, 1))),
                  Year = year(date),
                  Month = month(date, label = TRUE, abbr = TRUE)) 

rateDF <- bind_rows(lngRate, forRate) %>% 
  select(Year, Month, rate)


write_csv(rateDF, path = "./Data/Interest rates/projectedInterestRates.csv")





