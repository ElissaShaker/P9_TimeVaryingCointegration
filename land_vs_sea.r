start_time <- Sys.time()

library(tsDyn)
library(data.table)

# Read your CSV files
berkeley <- fread("data/Land vs sea/Berkeley_GlobTemp.csv")
setnames(
  berkeley,
  old = names(berkeley),
  new = c(
    "year","month",
    "anomaly_monthly","unc_monthly",
    "anomaly_annual","unc_annual",
    "anomaly_5yr","unc_5yr",
    "anomaly_10yr","unc_10yr",
    "anomaly_20yr","unc_20yr"
  )
)

hadsst <- fread("data/Land vs sea/HadSST.4.1.0.0_monthly_GLOBE.csv")

# Merge datasets by year and month
merged_data <- merge(berkeley, hadsst, by = c("year", "month"))

# Example: pick Berkeley and HadSST monthly anomalies
y1 <- merged_data$anomaly_monthly
y2 <- merged_data$anomaly

# Bind into a matrix
y <- cbind(y1, y2)
y <- as.data.table(y)
y <- na.omit(y) # Remove rows with missing values


n <- nrow(y)
# Number of endogenous variables in the analysis
k <- ncol(y)
# Helping functions
# Function to calculate the eigenvalues and eigenvectors with the QR algorithm
my.qr.eigen <- function(x){
  if(is.matrix(x) == FALSE){
    stop("\nPlease provide a quadratic matrix.\n")
  }
  x <- as.matrix(x)
  pQ <- diag(1, dim(x)[1]);
  # iterate 
  while(sum(abs(x[lower.tri(x, diag = FALSE)])) > 0.0000001){
    d <- qr(x);
    Q <- qr.Q(d);
    pQ <- pQ %*% Q;
    x <- qr.R(d) %*% Q; 
    
  }
  # Eigenvalues
  eig.values <- round(diag(x),5)
  # Eigenvectors
  eig.vec <- round(pQ,5)
  return(list(values = eig.values, vectors = eig.vec))
  
}

# Function to create lagged data.tables
# functions validated against proc varlags 
lag.data.table <- function(data, max.lag = 1){
  data[, shift(data, max.lag, give.names = TRUE)]
}
lags <- 1
# Output conformable to the the output of varlags
cbind(y[-(1:lags),] ,na.exclude(lag.data.table(y, max.lag = 1)))

# Function diffmatrix to calculate the first differences
diffmatrix <- function(x,max.diff = 1,max.lag = 1){
  #Add if condition to make it possible to differentiate between matrix and vector                  
  if(is.vector(x) == TRUE ){
    myx <- embed(c(rep(NA,max.lag),diff(x,max.lag,max.diff)),max.diff)
    colnames(myx) <- paste("v1.d",max.diff, sep=".")
    return(myx)
  }
  
  else if(is.matrix(x) == TRUE){
    myx <- rbind(matrix(rep(NA,max.lag), ncol = ncol(x)), matrix(diff(x,max.lag,max.diff), ncol = ncol(x)))
    mycolnames <- colnames(x)
    colnames(myx) <- paste(mycolnames,"d",max.diff, sep=".")
    return(myx)
  }
  #if matrix the result should be 0, if only a vector it should be 1
  else if(as.integer(is.null(ncol(x))) == 0 ){
    myx <- rbind(matrix(rep(NA,max.lag), ncol = ncol(x)), matrix(diff(x,max.lag,max.diff), ncol = ncol(x)))
    mycolnames <- colnames(x)
    colnames(myx) <- paste(mycolnames,"d",max.diff, sep=".")
    return(myx)
    
  }
}

as.data.table(diffmatrix(as.matrix(y), max.diff = 1, max.lag = 1))
# Removal of na.s
na.exclude(as.data.table(diffmatrix(as.matrix(y), max.diff = 1, max.lag = 1)))


# In gauss a preceding . means that the following calculations are meant to be carried out element-wise
# Function to calculate the time-varying cointegration
tvcoint <- function(y,p,m){
  # local ystar_1,ysub,y_1,dy,dylags,T,betau,resu,betav,resv,S00,S01,S10,S11,S00inv,S11inv,A,eval,evec,evs;
  y <- as.matrix(y)
  ystar_1 <- ycheb(y,varord = p,chebdim = m);		
  # dy is the same as the Z0 matrix in ca.jo
  # ylags would be the same as ZK
  # dylags is equivalent to Z1
  # checked the dimensions and content of all preliminary matrices against original gauss code
  dy <- na.exclude(as.data.table(diffmatrix(as.matrix(y), max.diff = 1, max.lag = 1)))
  dylags <-  na.exclude(lag.data.table(dy, max.lag = p))
  dy <- dy[-(1:lags)]
  myT <- nrow(dylags);             
  dylags <- cbind(rep(1,myT), dylags) 
  # a= b/c
  # a = (c'c)^(-1)c'b
  # betau <- dy/dylags
  # betau <- (dylags'dylags)^-1 dylags dy
  betau <- as.matrix(solve(crossprod(as.matrix(dylags))) %*% crossprod(as.matrix(dylags), as.matrix(dy)))
  resu <- as.matrix(dy)-(as.matrix(dylags) %*% betau)     
  # betav <- ystar_1/dylags;     
  # betav <- (dylags'dylags)^-1 dylags ystar_1
  betav <-  as.matrix(solve(crossprod(as.matrix(dylags))) %*% crossprod(as.matrix(dylags), as.matrix(ystar_1)))
  resv <- ystar_1-(as.matrix(dylags) %*% as.matrix(betav));
  
  S00 <- (crossprod(resu))/ myT;        
  S01 <- (crossprod(resu,resv))/myT;        
  S10 <- t(S01);                
  S11 <- (crossprod(resv))/myT;       
  S00inv <- solve(S00);       
  S11inv <- solve(S11);       
  A <- S11inv %*% S10 %*% S00inv %*% S01; 
  # calculation of eigenvalues and -vectors with QR algorithm
  valeigen <- my.qr.eigen(A);
  evs <- cbind(valeigen$values,valeigen$vectors)[order(valeigen$values, decreasing = TRUE), ]
  evec = t(evs[,2:ncol(evs)])
  detS00 <- det(S00)
  output <- list("eigenvalues" = evs[, 1], "eigenvectors" = evec, "determinant" = detS00)
  return(output)
}


# function to calculate the chebyshev polynomials            
# function validated against the original gauss code
ycheb <- function(mydata,varord,chebdim){
  #  local i,yst,k,n,nn,ind;
  k <- ncol(mydata);
  nn <- nrow(mydata)
  yst <- matrix(nrow = nn-varord-1, ncol = (chebdim+1)*k)
  yst[, 1:k] <- as.matrix(mydata[(varord+1):(nn-1), ])
  if(chebdim == 0){
    return(yst);
  } else {
    n <- length(yst[, 1]);
    ind <- seq(from = varord+2, by = 1, length.out = n);
    i = 1;
    for(i in 1:chebdim){
      yst[, (i*k+1):((i+1)*k)] <- sqrt(2)*cos(i*pi*(ind-0.5)/n)*yst[,1:k]
    }
    return(yst);
  }
}

mmax <- round(n/10)	# maximum dimension of Chebishev Polynomials
p <- 1;			# VAR order

# Setting up the matrices to hold the values

lrtvc <- matrix(nrow = mmax, ncol = k);	     # TVC Stat (row) m=1,...,mmax and (col) r=1,..,k
lrtvcpv <- matrix(nrow = mmax, ncol = k);	     # P-values using the Asymp Distr (chisquare(mrk))
betat <- matrix(nrow = n-p-1, ncol = k*mmax)  # b_t for r=1.  (row) i= observation; Cols 1 up to k= b1~...~bk (m=1); Cols k+1 up to 2k=b1~...~bk (m=2); ...
lnlikm <- matrix(nrow = mmax, ncol = k);	     # log likelihood for different m and r 
aic <- matrix(nrow = mmax, ncol = k);	     # AIC Akaike model selection criterion 
bic <- matrix(nrow = mmax, ncol = k);	     # BIC Schwarz msc 
hann <- matrix(nrow = mmax, ncol = k);	     # Hannan-Quinn msc 



# Baseline standard cointegration analysis without chebyshev polynomials
return.tvcoint.0 <- tvcoint(y,p,0)

# Log-likelihood of baseline cointegration
ll0=log(1 - return.tvcoint.0$eigenvalues, base = exp(1));

# Main Function call, using all defined functions
for(m in 1:mmax){
  result.tvcoint.m <- tvcoint(y,p,m); # OUT: Lambda; Eigenvectors q1...qr...q(m+1)k; det(S00)
  evect <- result.tvcoint.m$eigenvectors
  eval <- result.tvcoint.m$eigenvalues
  llm <- log(1-eval, base = exp(1));	
  ind <- seq(from = p+2, by = 1, length.out = n-p-1); 
  beta1sum <- matrix(nrow = m, ncol = length(ind));	
  beta2sum <- matrix(nrow = m, ncol = length(ind));	
  beta3sum <- matrix(nrow = m, ncol = length(ind));	
  
  betaSum <- matrix(0, nrow = length(ind), ncol = k)
  for (j in 1:k) {
    for (mm in 1:m) {
      betaSum[, j] <- betaSum[, j] + 
        evect[k * mm + j, 1] * sqrt(2) * cos(mm * pi * (ind - 0.5) / (n - p - 1))
    }
  }
  betat[, (k*(m-1)+1):(k*m)] <- evect[1:k, 1] + t(betaSum)
  
  for(r in 1:k){
    lrtvc[m,r] <- (n-p-1)*sum(ll0[1:r]) - (n-p-1)*sum(llm[1:r]);	
    lrtvcpv[m,r] <- 1-pchisq(lrtvc[m,r],m*r*k)
    lnlikm[m,r] <- (log(r, base = exp(1))-k-k*log((2*pi), base = exp(1)))*(n-p-1)/2 - sum(llm[1:r])*(n-p-1)/2 - (log(result.tvcoint.m$determinant, base = exp(1)))*(n-p-1)/2;				        
    npar <- (m+1)*k*r+r*k+k^2+(k+(p-1)*k^2);		
    aic[m,r] <-  -2*lnlikm[m,r]/(n-p-1)+npar*2/(n-p-1);					
    bic[m,r] <-  -2*lnlikm[m,r]/(n-p-1)+npar*(log(n-p-1, base= exp(1)))/(n-p-1);			
    hann[m,r] <-  -2*lnlikm[m,r]/(n-p-1)+npar*(log(log((n-p-1), base = exp(1)), base = exp(1)))*2/(n-p-1);	
  }
  
}

mminhan <- matrix(nrow = k, ncol = 1);
for(r in 1:k){
  mminhan[r,1] <- which.min(hann[,r])
}

m <- max(mminhan)

#save.image("myTVCData.RData", compress = "xz")

end_time <- Sys.time()
time <- end_time - start_time
cat("It took", time, "min to run the code")

