# This function is to calculate population rhos
get.rho <- function(mx, Year, Age){ 
  MX<- data.frame(mx, Year, Age)
  MX <- MX %>% spread(Year, mx)
  mx.2   <- MX[ ,-c(1)]
  mx.1   <- MX[ ,-c(ncol(MX))]
  rho <- cbind(Age=MX$Age,-log(mx.2/mx.1)) # rho populaton by age
  rho <- rho %>% gather(Year, rho, -Age )
  
  return(data.frame(rho))
} 


get.phi <- function(delta, Year){ 
  rx<- data.frame(delta, Year)
  rx <- rx %>% spread(Year, delta)
  rx.2   <- rx[ ,-c(1)]
  rx.1   <- rx[ ,-c(ncol(rx))]
  phi <- -log(rx.2/rx.1)
  phi <- phi %>% gather(Year, phi)
  phi$Year <- as.integer(phi$Year)
  return(data.frame(phi))
}



get.adot <- function(ax, Year){ 
  a <- data.frame(ax, Year)
  a <- a %>% spread(Year, ax)
  ax.2   <- a[ ,-c(1)]
  ax.1   <- a[ ,-c(ncol(a))]
  adot <- ax.2 - ax.1
  #adot <- -log(ax.2/ax.1) * ax[-1]
  adot <-adot %>% gather(Year,adot)
  adot$Year <- as.integer(adot$Year)
  return(data.frame(adot))
}







calc.fun <- function(Age, mu, delta){
# Mortality
 Mu  <- cumtrapz(Age,mu) # Cumulative hazard
 spx <- exp(-Mu)
 sqx <- 1-spx
 fx  <- mu * spx

# Interest
 Delta <- cumtrapz(Age,delta) # Cumulative interest
 vs <- exp(-Delta)
 dat <- data.frame(Age,mu, Mu, spx, sqx, fx, delta, Delta, vs)
 
 return(dat)}



# Calculate annuities, entropies and durations

ss <- subset(r, Year == 2010 & Sex == "Females")

Age <- ss$Age
vs <- ss$vs
spx <- ss$spx
mu <- ss$mu
fx <- ss$fx
delta <- ss$delta
sEx <- spx * vs 


#assume SOA example life table to be load
data(soaLt)
act=with(soaLt, new("actuarialtable",interest=delta,
                         x=Age,lx=spx,name="SOA2008"))
#evaluate and life-long annuity for an aged 65
axn(soa08Act, x=50)

 
 
# life expectancy
 
Lx  <- c(spx[-1],0) + 0.5*fx
Tx  <- rev(cumsum(rev(Lx)))

Tx/spx # Esta hace mas sentido, lo que tengo que hacer es todo en base a life tables

 rev(cumtrapz(Age, spx )) *spx

 rev(cumsum(rev(spx  )))
 
  
 # Entropy mort
 #THIS IS CORRECT
 -rev(cumtrapz(Age, spx * log(spx) ))[1]
 -rev(cumsum(rev(spx * log(spx) )))[1]

 

 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 calc.aH <- function(Age, spx, vs){
  

  # Annuities
  ax <- trapz(Age, spx * vs )
  cum.ax <-rev(cumtrapz(Age, spx * vs )) * (spx *vs)
  sum(spx *vs)
  
  cum.ax <- rev(cumsum(rev(spx * vs)))  * (spx *vs)
  
  # Numerators entropy
  
 - trapz(Age, spx * log(spx) * vs)
  
  trapz(Age, mu *cum.ax)
  
  
  a.dagger <- -trapz(Age, lx * log(lx) * exp(-delta * (Age))) #- min(Age) )))
  H.a <- a.dagger / ax
  out <- data.frame(ax, a.dagger, H.a)
  return(out)
}





