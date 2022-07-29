# Rates of mortality improvement
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





# Annuity, entropy and durations
get.ax.cod <- function(Age, mu, delta, p1, p2, p3, p4, p5){
  # Total mortality
  m <- length(Age)
  qx  <- pmin(mu / (1 + 0.5* mu),1) # Mx to qx, I capped it to one to avoid negative numbers
  qx[m] <- 1
  px  <- 1-qx
  lx  <- cumprod(c(1,px))
  dx  <- -diff(lx)
  lx <- lx[-(m+1)]
  
  # Interest, this is for all causes
  Delta <- cumtrapz(Age,delta) # Cumulative interest
  vs <- exp(-Delta)
  
  # Actuarial PV
  sEx <- lx * vs 
  
  # Annuity
  ax <- rev(cumsum(rev(sEx)))/ sEx
  ax[is.na(ax)] <- 0 ## in case of NA values
  ax[is.infinite(ax)] <- 0 ## in case of Inf values
  
  # Entropy
  mEax <- mu * sEx * ax # This is the important one
  hp <- rev(cumsum(rev( mEax )))
  hc <- rev(cumsum(rev( (Age-Age[1]+1) * sEx )))
  Hp <- hp / ax
  
  Hp[is.na(Hp)] <- 0 ## in case of NA values
  Hp[is.infinite(Hp)] <- 0 ## in case of Inf values
  
  
  # Duration, this is for the total population
  dEax <- delta *sEx * ax
  dp <-   rev(cumsum(rev( dEax )))
  Dp <- dp / ax
  
  Dp[is.na(Dp)] <- 0 ## in case of NA values
  Dp[is.infinite(Dp)] <- 0 ## in case of Inf values
  
  # Cause-specific entropies, as Preston et. al (2001)
  # First cause
  mu1 <- mu * p1
  mu2 <- mu * p2
  mu3 <- mu * p3
  mu4 <- mu * p4
  mu5 <- mu * p5
  
  mEax1 <- mu1 * sEx * ax # cause specific mu times differed annuity
  hp1 <- rev(cumsum(rev( mEax1 )))
  Hp1 <- hp1 / ax
  Hp1[is.na(Hp1)] <- 0 ## in case of NA values
  Hp1[is.infinite(Hp1)] <- 0 ## in case of Inf values
  
  mEax2 <- mu2 * sEx * ax # cause specific mu times differed annuity
  hp2 <- rev(cumsum(rev( mEax2 )))
  Hp2 <- hp2 / ax
  Hp2[is.na(Hp2)] <- 0 ## in case of NA values
  Hp2[is.infinite(Hp2)] <- 0 ## in case of Inf values

  mEax3 <- mu3 * sEx * ax # cause specific mu times differed annuity
  hp3 <- rev(cumsum(rev( mEax3 )))
  Hp3 <- hp3 / ax
  Hp3[is.na(Hp3)] <- 0 ## in case of NA values
  Hp3[is.infinite(Hp3)] <- 0 ## in case of Inf values
  
  mEax4 <- mu4 * sEx * ax # cause specific mu times differed annuity
  hp4 <- rev(cumsum(rev( mEax4 )))
  Hp4 <- hp4 / ax
  Hp4[is.na(Hp4)] <- 0 ## in case of NA values
  Hp4[is.infinite(Hp4)] <- 0 ## in case of Inf values
  
  mEax5 <- mu5 * sEx * ax # cause specific mu times differed annuity
  hp5 <- rev(cumsum(rev( mEax5 )))
  Hp5 <- hp5 / ax
  Hp5[is.na(Hp5)] <- 0 ## in case of NA values
  Hp5[is.infinite(Hp5)] <- 0 ## in case of Inf values
  
  dat <- data.frame(Age, mu, mu1, mu2, mu3, mu4, mu5, 
                    delta, fx=dx, spx=lx, vs, sEx, ax, dp, Dp,
                    hp, hp1, hp2, hp3, hp4, hp5,
                    Hp, Hp1, Hp2, Hp3, Hp4, Hp5)
  return(dat)}


# Rates of mortality improvement by cause of death
get.rho.cod <- function(mu,mu1,mu2,mu3,mu4,mu5, Year, Age){ 
  
  # Total
  MX<- data.frame(mu, Year, Age)
  MX <- MX %>% spread(Year, mu)
  mx.2   <- MX[ ,-c(1,2)]
  mx.1   <- MX[ ,-c(1,ncol(MX))]
  rho <- cbind(Age=MX$Age,-log(mx.2/mx.1)) # rho populaton by age
  rho <- rho %>% gather(Year, rho, -Age )
  rho$Year <- as.integer(rho$Year)
  
  # First cause
  MX1<- data.frame(mu1, Year, Age)
  MX1 <- MX1 %>% spread(Year, mu1)
  mx1.2   <- MX1[ ,-c(1,2)]
  mx1.1   <- MX1[ ,-c(1,ncol(MX1))]
  rho1 <- cbind(Age=MX1$Age,-log(mx1.2/mx1.1)) # rho populaton by age
  rho1 <- rho1 %>% gather(Year, rho1, -Age )
  rho1$Year <- as.integer(rho1$Year)
  
  # Second cause
  MX2<- data.frame(mu2, Year, Age)
  MX2 <- MX2 %>% spread(Year, mu2)
  mx2.2   <- MX2[ ,-c(1,2)]
  mx2.1   <- MX2[ ,-c(1,ncol(MX2))]
  rho2 <- cbind(Age=MX2$Age,-log(mx2.2/mx2.1)) # rho populaton by age
  rho2 <- rho2 %>% gather(Year, rho2, -Age )
  rho2$Year <- as.integer(rho2$Year)
  
  # Third cause
  MX3<- data.frame(mu3, Year, Age)
  MX3 <- MX3 %>% spread(Year, mu3)
  mx3.2   <- MX3[ ,-c(1,2)]
  mx3.1   <- MX3[ ,-c(1,ncol(MX3))]
  rho3 <- cbind(Age=MX3$Age,-log(mx3.2/mx3.1)) # rho populaton by age
  rho3 <- rho3 %>% gather(Year, rho3, -Age )
  rho3$Year <- as.integer(rho3$Year)
  
  # Fourth cause
  MX4<- data.frame(mu4, Year, Age)
  MX4 <- MX4 %>% spread(Year, mu4)
  mx4.2   <- MX4[ ,-c(1,2)]
  mx4.1   <- MX4[ ,-c(1,ncol(MX4))]
  rho4 <- cbind(Age=MX4$Age,-log(mx4.2/mx4.1)) # rho populaton by age
  rho4 <- rho4 %>% gather(Year, rho4, -Age )
  rho4$Year <- as.integer(rho4$Year)
  
  # Fifth cause
  MX5<- data.frame(mu5, Year, Age)
  MX5 <- MX5 %>% spread(Year, mu5)
  mx5.2   <- MX5[ ,-c(1,2)]
  mx5.1   <- MX5[ ,-c(1,ncol(MX5))]
  rho5 <- cbind(Age=MX5$Age,-log(mx5.2/mx5.1)) # rho populaton by age
  rho5 <- rho5 %>% gather(Year, rho5, -Age )
  rho5$Year <- as.integer(rho5$Year)
  
  
  out <- data.frame(merge(rho, rho1, by = c("Age", "Year")))
  out <- data.frame(merge(out, rho2, by = c("Age", "Year")))
  out <- data.frame(merge(out, rho3, by = c("Age", "Year")))
  out <- data.frame(merge(out, rho4, by = c("Age", "Year")))
  out <- data.frame(merge(out, rho5, by = c("Age", "Year")))
  out <- arrange(out, Year, Age)
  return(out)
} 



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


# Cause-specific decomposition
get.decomp.cod <- function(Age, ax, sEx, mu, mu1, mu2, mu3, mu4, mu5,
                           delta, rho, rho1, rho2, rho3, rho4, rho5,
                           phi, H, H1, H2, H3, H4, H5, D){
  
  # Weights
  sMx <- mu * sEx * ax
  sMx1 <- mu1 * sEx * ax
  sMx2 <- mu2 * sEx * ax
  sMx3 <- mu3 * sEx * ax
  sMx4 <- mu4 * sEx * ax
  sMx5 <- mu5 * sEx * ax
  
  sWx <- delta * sEx * ax
  
  # Average change over time 
  phibar.num   <- sum( phi * sWx )
  phibar.denom <- sum( sWx )
  phibar <- phibar.num / phibar.denom
  
  rhobar.num   <- sum( rho * sMx )
  rhobar.denom <- sum( sMx )
  rhobar <- rhobar.num / rhobar.denom
  
  # by causes of death
  
  rhobar1.num   <- sum( rho1 * sMx1 )
  rhobar1.denom <- sum( sMx1 )
  rhobar1 <- rhobar1.num / rhobar1.denom
  
  rhobar2.num   <- sum( rho2 * sMx2 )
  rhobar2.denom <- sum( sMx2 )
  rhobar2 <- rhobar2.num / rhobar2.denom
  
  rhobar3.num   <- sum( rho3 * sMx3 )
  rhobar3.denom <- sum( sMx3 )
  rhobar3 <- rhobar3.num / rhobar3.denom

  rhobar4.num   <- sum( rho4 * sMx4 )
  rhobar4.denom <- sum( sMx4 )
  rhobar4 <- rhobar4.num / rhobar4.denom
  
  rhobar5.num   <- sum( rho5 * sMx5 )
  rhobar5.denom <- sum( sMx5 )
  rhobar5 <- rhobar5.num / rhobar5.denom
  
  # Decomposition
  int  <- phibar * D[1]
  mort <- rhobar * H[1]
  mort1 <- rhobar1 * H1[1]
  mort2 <- rhobar2 * H2[1]
  mort3 <- rhobar3 * H3[1]
  mort4 <- rhobar4 * H4[1]
  mort5 <- rhobar5 * H5[1]
  
  
  adot <- mort+int
  adot.cod <- mort1+mort2+mort3+mort4+mort5+int
  
  out <- data.frame(ax=ax[1],rhobar,rhobar1, rhobar2,rhobar3,rhobar4,rhobar5,
                    phibar, H=H[1],D=D[1],
                    mort, mort1, mort2,mort3,mort4,mort5,
                    int, adot, adot.cod)
  return(out)}



# Age decomposition
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




