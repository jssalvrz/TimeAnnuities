

# Calculate lx
calc.lx <- function(lt, Qx){
  dat     <- data.frame(lt, Qx)
  dat$id  <- ifelse(dat$lt>= dat$Qx, 1,0)
  dat <- subset(dat, id == 1)
  dat.cdf <- empirical_cdf(dat$lt, 
                          ubounds=seq(round(min(dat$lt), 2),
                                      round(max(dat$lt), 2)-0.01,.01)) # Calculate CDF
  names(dat.cdf) <- c("Age", "dx.mil", "dx")
  dat.cdf$lx <- 1 - dat.cdf$dx # Survival function
  return(dat.cdf)}


# Calculate e-dagger, H and ex
calc.H <- function(Age, lx){
  ex <- trapz(Age, lx) # life expectancy; trapz is used to integrate
  e.dagger <- -trapz(Age, lx * log(lx))
  H <- e.dagger / ex
  out <- data.frame(ex, e.dagger, H)
  return(out)
  }


# Calculate annuities and the entropy of annuities

calc.aH <- function(Age, lx, delta){
ax <- trapz(Age, lx * exp(-delta * (Age))) #- min(Age) )) ) 
a.dagger <- -trapz(Age, lx * log(lx) * exp(-delta * (Age))) #- min(Age) )))
H.a <- a.dagger / ax
out <- data.frame(ax, a.dagger, H.a)
return(out)
}




