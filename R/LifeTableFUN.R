## function for constructing a classic (& rather general) lifetable
lifetable.mx <- function(x, mx, sex="M", ax=NULL){
  m <- length(x)
  mx  <- mx
  n <- c(diff(x), NA)
  if(is.null(ax)){
    ax <- rep(0,m)
    if(x[1]!=0 | x[2]!=1){
      ax <- n/2
      ax[m] <- 1 / mx[m]
    }else{    
      if(sex=="F"){
        if(mx[1]>=0.107){
          ax[1] <- 0.350
        }else{
          ax[1] <- 0.053 + 2.800*mx[1]
        }
      }
      if(sex=="M"){
        if(mx[1]>=0.107){
          ax[1] <- 0.330
        }else{
          ax[1] <- 0.045 + 2.684*mx[1]
        }
      }
      ax[-1] <- n[-1]/2
      ax[m] <- 1 / mx[m]
    }
  }
  qx  <- pmin((n*mx / (1 + (n - ax) * mx)), rep(1,m))
  qx[m] <- 1
  
  px  <- 1-qx
  lx  <- cumprod(c(1,px))*100000
  dx  <- -diff(lx)
  Lx  <- n*lx[-1] + ax*dx
  lx <- lx[-(m+1)]
  Lx[m] <- lx[m]/mx[m]
  Lx[is.na(Lx)] <- 0 ## in case of NA values
  Lx[is.infinite(Lx)] <- 0 ## in case of Inf values
  Tx  <- rev(cumsum(rev(Lx)))
  ex  <- Tx/lx
  return.df <- data.frame(x, n, mx, ax, qx, px, lx, dx, Lx, Tx, ex)
  return(return.df)
}


## function for constructing a classic (& rather general) lifetable
lifetable <- function(x, Nx, Dx, sex="M", ax=NULL){
  m <- length(x)
  mx  <- Dx/Nx
  n <- c(diff(x), NA)
  if(is.null(ax)){
    ax <- rep(0,m)
    if(x[1]!=0 | x[2]!=1){
      ax <- n/2
      ax[m] <- 1 / mx[m]
    }else{    
      if(sex=="F"){
        if(mx[1]>=0.107){
          ax[1] <- 0.350
        }else{
          ax[1] <- 0.053 + 2.800*mx[1]
        }
      }
      if(sex=="M"){
        if(mx[1]>=0.107){
          ax[1] <- 0.330
        }else{
          ax[1] <- 0.045 + 2.684*mx[1]
        }
      }
      ax[-1] <- n[-1]/2
      ax[m] <- 1 / mx[m]
    }
  }
  qx  <- n*mx / (1 + (n - ax) * mx)
  qx[m] <- 1
  px  <- 1-qx
  lx  <- cumprod(c(1,px))*100000
  dx  <- -diff(lx)
  Lx  <- n*lx[-1] + ax*dx
  lx <- lx[-(m+1)]
  Lx[m] <- lx[m]/mx[m]
  Lx[is.na(Lx)] <- 0 ## in case of NA values
  Lx[is.infinite(Lx)] <- 0 ## in case of Inf values
  Tx  <- rev(cumsum(rev(Lx)))
  ex  <- Tx/lx
  return.df <- data.frame(x, n, Nx, Dx, mx, ax, qx, px, lx, dx, Lx, Tx, ex)
  return(return.df)
}

## function for constructing a classic (& rather general) lifetable using IDs


lifetable.id <- function(id, datos, sex = "M"){
  datos <- datos
  ax<-NULL
  country    <- datos[[id]][ ,1]
  Year       <- as.numeric(datos[[id]][ ,2])
  
  x    <- as.numeric(datos[[id]][ ,3])
  if(sex =="M"){
    Dx <- as.numeric(datos[[id]][ ,5]) 
    Nx <- as.numeric(datos[[id]][ ,8]) 
  }else{
    Dx <- as.numeric(datos[[id]][ ,4])
    Nx <- as.numeric(datos[[id]][ ,8]) 
  }
  m <- length(x)
  mx  <- Dx/Nx
  n <- c(diff(x), NA)
  if(is.null(ax)){
    ax <- rep(0,m)
    if(x[1]!=0 | x[2]!=1){
      ax <- n/2
      ax[m] <- 1 / mx[m]
    }else{    
      if(sex=="F"){
        if(mx[1]>=0.107){
          ax[1] <- 0.350
        }else{
          ax[1] <- 0.053 + 2.800*mx[1]
        }
      }
      if(sex=="M"){
        if(mx[1]>=0.107){
          ax[1] <- 0.330
        }else{
          ax[1] <- 0.045 + 2.684*mx[1]
        }
      }
      ax[-1] <- n[-1]/2
      ax[m] <- 1 / mx[m]
    }
  }
  qx  <- n*mx / (1 + (n - ax) * mx)
  qx[m] <- 1
  px  <- 1-qx
  lx  <- cumprod(c(1,px))*100000
  dx  <- -diff(lx)
  Lx  <- n*lx[-1] + ax*dx
  lx <- lx[-(m+1)]
  Lx[m] <- lx[m]/mx[m]
  Lx[is.na(Lx)] <- 0 ## in case of NA values
  Lx[is.infinite(Lx)] <- 0 ## in case of Inf values
  Tx  <- rev(cumsum(rev(Lx)))
  ex  <- Tx/lx
  id <- paste(country, Year, sex, sep = "")
  return.df <- data.frame(x, n, Nx, Dx, mx, ax, qx, px, lx, dx, Lx, Tx, ex,country, Year, sex, id)
  return(return.df)
}

#Life tables using death rates from HMD

lifetable.id.r <- function(id, datos, sex = "M"){
  datos <- datos
  ax<-NULL
  country    <- datos[[id]][ ,6]
  Year       <- as.numeric(datos[[id]][ ,1])
  
  x    <- as.numeric(datos[[id]][ ,2])
  if(sex =="M"){
    mx <- as.numeric(datos[[id]][ ,4]) 
  }else{
    mx <- as.numeric(datos[[id]][ ,3])
  }
  m <- length(x)
  n <- c(diff(x), NA)
  if(is.null(ax)){
    ax <- rep(0,m)
    if(x[1]!=0 | x[2]!=1){
      ax <- n/2
      ax[m] <- 1 / mx[m]
    }else{    
      if(sex=="F"){
        if(mx[1]>=0.107){
          ax[1] <- 0.350
        }else{
          ax[1] <- 0.053 + 2.800*mx[1]
        }
      }
      if(sex=="M"){
        if(mx[1]>=0.107){
          ax[1] <- 0.330
        }else{
          ax[1] <- 0.045 + 2.684*mx[1]
        }
      }
      ax[-1] <- n[-1]/2
      ax[m] <- 1 / mx[m]
    }
  }
  qx  <- n*mx / (1 + (n - ax) * mx)
  qx[m] <- 1
  px  <- 1-qx
  lx  <- cumprod(c(1,px))*100000
  dx  <- -diff(lx)
  Lx  <- n*lx[-1] + ax*dx
  lx <- lx[-(m+1)]
  Lx[m] <- lx[m]/mx[m]
  Lx[is.na(Lx)] <- 0 ## in case of NA values
  Lx[is.infinite(Lx)] <- 0 ## in case of Inf values
  Tx  <- rev(cumsum(rev(Lx)))
  ex  <- Tx/lx
  id <- paste(country, Year, sex, sep = "")
  return.df <- data.frame(x, n, mx, ax, qx, px, lx, dx, Lx, Tx, ex,country, Year, sex, id)
  return(return.df)
}



## function for constructing a lifetable starting from probabilities
lifetable.qx <- function(x, qx, sex="M", ax=NULL, last.ax=5.5){
  m <- length(x)
  n <- c(diff(x), NA)
  qx[is.na(qx)] <- 0
  if(is.null(ax)){
    ax <- rep(0,m)
    if(x[1]!=0 | x[2]!=1){
      ax <- n/2
      ax[m] <- last.ax
    }else{    
      if(sex=="F"){
        if(qx[1]>=0.1){
          ax[1] <- 0.350
        }else{
          ax[1] <- 0.05 + 3*qx[1]
        }
      }
      if(sex=="M"){
        if(qx[1]>=0.1){
          ax[1] <- 0.33
        }else{
          ax[1] <- 0.0425 + 2.875*qx[1]
        }
      }
      ax[-1] <- n[-1]/2
      ax[m] <- last.ax
    }
  }
  px  <- 1-qx
  lx  <- cumprod(c(1,px))*100000
  dx  <- -diff(lx)
  Lx  <- n*lx[-1] + ax*dx
  lx <- lx[-(m+1)]
  Lx[m] <- lx[m]*last.ax
  Lx[is.na(Lx)] <- 0 ## in case of NA values
  Lx[is.infinite(Lx)] <- 0 ## in case of Inf values
  Tx  <- rev(cumsum(rev(Lx)))
  ex  <- Tx/lx
  return.df <- data.frame(x, n, ax, qx, px, lx, dx, Lx, Tx, ex)
  return(return.df)
}




## x=agesA
## Nx=myE
## Dx=myD
## ax=rep(0.5,mA)
## ns=100
## level=0.90
## which.x = 31
## sex = "M"
## last.ax=0.5
## general function for getting CI of life expectancy at x
CIex <- function(x, Nx, Dx, sex = "M", ax = NULL,
                 which.x=0, ns=1000, level=0.95){
  ## point-estimated lifetable
  LT <- lifetable(x, Nx, Dx, sex = sex, ax = ax)
  ## number of ages
  m <- nrow(LT)
  ## estimated probs
  qx <- LT$qx
  ## trials for binomial, rounded
  Ntil <- round(Dx/qx)
  ## ax for last age
  last.ax <- LT$ax[m]
  ## simulated death counts
  ## from Binomial distribution
  Y <- suppressWarnings(matrix(rbinom(m*ns,
                                      Ntil,
                                      qx),
                               m,ns))
  ## simulated probabilities
  QX <- Y/Ntil
  ## which age?
  wh <- which(x==which.x)
  ## for all replicates, compute life expectancy
  ## ## by a for-loop
  ## exsimA <- rep(0,ns)
  ## for(s in 1:ns){
  ##   exsimA[s] <-lifetable.qx(x, qx=QX[,s], sex,
  ##                        last.ax=last.ax)$ex[wh]
  ## }
  ## by apply command
  ## (slighly faster, less intuitive)
  fun.ex <- function(qx){
    return(lifetable.qx(x=x, qx, sex=sex,
                        last.ax=last.ax)$ex[wh])
  }
  exsim <- apply(QX, 2, fun.ex)

  ## confidence interval
  CI <- quantile(exsim,
                 probs = c((1-level)/2,
                           1 - (1-level)/2))
  ## output
  out <- list(ex=LT$ex[wh],
              meanex=mean(exsim),
              CIex=CI,
              exsim=exsim,
              which.x=which.x)
  return(out)
}
