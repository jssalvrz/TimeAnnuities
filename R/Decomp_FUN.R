# Rates of mortality improvemen
get.rho <- function(mx, Year, Age){ 
  MX<- data.frame(mx, Year, Age)
  MX <- MX %>% spread(Year, mx)
  mx.2   <- MX[ ,-c(1,2)]
  mx.1   <- MX[ ,-c(1,ncol(MX))]
  rho <- cbind(Age=MX$Age,-log(mx.2/mx.1)) # rho populaton by age
  rho <- rho %>% gather(Year, rho, -Age )
  rho$Year <- as.integer(rho$Year)
  return(data.frame(rho))
} 

# Change in a single interest rate
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


# Phi using the terms structure
get.phi.str <- function(delta, Year, Age){ 
    Delta<- data.frame(delta, Year, Age)
    Delta <- Delta %>% spread(Year, delta)
    delta.2   <- Delta[ ,-c(1,2)]
    delta.1   <- Delta[ ,-c(1,ncol(Delta))]
    
    phi <- cbind(Age=Delta$Age,-log(delta.2/delta.1)) # phi by age
    phi <- phi %>% gather(Year, phi, -Age )
    phi$Year <- as.integer(phi$Year)
    return(data.frame(phi))
} 




get.adot <- function(ax, Year){ 
  a <- data.frame(ax, Year)
  a <- a %>% spread(Year, ax)
  ax.2   <- a[ ,-c(1)]
  ax.1   <- a[ ,-c(ncol(a))]
  adot <- (ax.2 - ax.1)/ax.1
  #adot <- -log(ax.2/ax.1) #* ax[-1]
  adot <-adot %>% gather(Year,adot)
  adot$Year <- as.integer(adot$Year)
  return(data.frame(adot))
}



# Annuity, entropy and durations
get.ax <- function(Age, mu, delta){
# Mortality
# I calculate spx from mu, as in a lifetable where 
# it is assumed that deaths occur at the half of the interval
  m <- length(Age)
  qx  <- pmin(mu / (1 + 0.5* mu),1) # Mx to qx, I capped it to one to avoid negative numbers
  qx[m] <- 1
  px  <- 1-qx
  lx  <- cumprod(c(1,px))
  dx  <- -diff(lx)
  #Lx  <- lx[-1] + 0.5*dx # I need all these functions to calculate ax
  lx <- lx[-(m+1)]
  #Lx[m] <- lx[m]/mu[m]
  #Lx[is.na(Lx)] <- 0 ## in case of NA values
  #Lx[is.infinite(Lx)] <- 0 ## in case of Inf values
  #Tx  <- rev(cumsum(rev(Lx)))
  #ex  <- Tx/lx
  
  # Interest
  Delta <- cumtrapz(Age,delta) # Cumulative interest
  vs <- exp(-Delta)
  
  # Actuarial PV
  sEx <- lx * vs 
  
  # Annuity
  ax <- rev(cumsum(rev(sEx)))/ sEx
  
  ax[is.na(ax)] <- 0 ## in case of NA values
  ax[is.infinite(ax)] <- 0 ## in case of Inf values
  
  # Entropy
  mEax <- mu * sEx * ax
  hp <- rev(cumsum(rev( mEax )))
  hc <- rev(cumsum(rev( (Age-Age[1]+1) * sEx )))
  Hp <- hp / ax
  Hc <- hc / ax
  
  Hp[is.na(Hp)] <- 0 ## in case of NA values
  Hp[is.infinite(Hp)] <- 0 ## in case of Inf values
  
  Hc[is.na(Hc)] <- 0 ## in case of NA values
  Hc[is.infinite(Hc)] <- 0 ## in case of Inf values
  
  
  # Duration
  dEax <- delta *sEx * ax
  dp <-   rev(cumsum(rev( dEax )))
  dc <- rev(cumsum(rev( (Age-Age[1]+1) * sEx )))
  Dp <- dp / ax
  Dc <- dc / ax
  
  Dp[is.na(Dp)] <- 0 ## in case of NA values
  Dp[is.infinite(Dp)] <- 0 ## in case of Inf values
  
  Dc[is.na(Dc)] <- 0 ## in case of NA values
  Dc[is.infinite(Dc)] <- 0 ## in case of Inf values
  
  
 dat <- data.frame(Age, mu, delta, fx=dx, spx=lx, vs, sEx, ax, hp, hc, dp, dc, Hp, Hc, Dp, Dc)
 return(dat)}


# Full decomposition
get.decomp <- function(Age, ax, sEx, mu, delta, rho, phi, H,D,Dc){

# Weights
sMx <- mu * sEx * ax
sWx <- delta * sEx * ax

# Average change over time 
rhobar.num   <- sum( rho * sMx )
rhobar.denom <- sum( sMx )
rhobar <- rhobar.num / rhobar.denom

phibar.num   <- sum( phi * sWx )
phibar.denom <- sum( sWx )
phibar <- phibar.num / phibar.denom

# Decomposition
mort <- rhobar * H[1]
int <- phibar * D[1]
adot <- mort+int

out <- data.frame(ax=ax[1],rhobar, phibar, H=H[1],D=D[1],Dc=Dc[1], mort, int, adot)
return(out)}



# Full decomposition
get.age.decomp <- function(Age, ax, sEx, mu, delta, rho, phi){
  
  # Weights
  sMx <- mu * sEx * ax
  sWx <- delta * sEx * ax
  
  # Average change over time 
  rhobar.num   <- sum( rho * sMx )
  rhobar.denom <- sum( sMx )
  rhobar <- rhobar.num / rhobar.denom
  
  phibar.num   <- sum( phi * sWx )
  phibar.denom <- sum( sWx )
  phibar <- phibar.num / phibar.denom
  
  
  # Entropy
  hp <- rev(cumsum(rev( sMx )))
  Hp <- hp / ax

  Hp[is.na(Hp)] <- 0 ## in case of NA values
  Hp[is.infinite(Hp)] <- 0 ## in case of Inf values
  
 
  # Duration
  dp <-   rev(cumsum(rev( sWx )))
  Dp <- dp / ax

  Dp[is.na(Dp)] <- 0 ## in case of NA values
  Dp[is.infinite(Dp)] <- 0 ## in case of Inf values
  
  H<-Hp[1]
  D<-Dp[1]
 
  # Decomposition
  mort <- rhobar * H
  int  <- phibar * D
  adot <- mort+int
  
  out <- data.frame(ax=ax[1],rhobar, phibar, H, D, mort, int, adot)
  return(out)}




