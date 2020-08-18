
#######################################################################
####  Initial look-see simulation for PO clustering methods.
####  Scott Aug 2020.
####
#######################################################################

library( raster)
library( RandomFields)

set.seed( 001)

####  Define survey area and covariates
lenny <- 250
xSeq <- seq( from=0, to=1, length=lenny)
ySeq <- seq( from=0, to=1, length=lenny)
X <- expand.grid( x=xSeq, y=ySeq)
#define two covariates within the study area
Mod1 <- RMgauss( var=1, scale=0.2) + RMnugget( var=0.01)
Mod2 <- RMgauss( var=1, scale=0.5) + RMnugget( var=0.01)
simmy1 <- RFsimulate( Mod1, x=xSeq, y=ySeq)
simmy2 <- RFsimulate( Mod2, x=xSeq, y=ySeq)
X <- cbind( X, as.numeric( as.matrix( simmy1)), as.numeric( as.matrix( simmy2)))
X[,-(1:2)] <- apply( X[,-(1:2)], 2, scale)
colnames( X) <- c("x","y","covar1","covar2")
#plot the covariates
Xbrick <- brick( rasterFromXYZ( X[,c(1,2,3)]), rasterFromXYZ( X[,c(1,2,4)]))
plot( Xbrick, zlim=range( values( Xbrick)))

####  Define the distribution of the groups (thinning process)
multinomialLogit <- function(linPreds){
  #multinomial logistic function.  Takes K-1 linear predictors and turns them into K probabilitites
  K <- length( linPreds) + 1
  pis <- rep( NA, K)
  pis[1:(K-1)] <- exp( linPreds[1:(K-1)])
  pis[1:(K-1)] <- pis[1:(K-1)] /  (1+sum( pis[1:(K-1)]))
  pis[K] <- 1 - sum( pis[1:(K-1)])
  
  return( pis)
}

####  Data set up
S <- 20 #number of spp
K <- 5  #number of RCPs
p <- 2

modMat <- model.matrix( ~1+covar1+covar2, data=X)
beta <- matrix( rnorm( (K-1)*(p+1)), nrow=p+1, ncol=K-1)
rownames( beta) <- colnames( modMat)
linPred <- modMat %*% beta
pis <- t( apply( linPred, 1, multinomialLogit))
colnames( pis) <- paste0( "pi",1:K)

piBrick <- rasterFromXYZ( cbind( X[,1:2], pis))
plot( piBrick, zlim=c(0,1))

####  Define the unthinned process
averageNumberPresent <- rnbinom( n=S, mu=10, size=0.95)
average.alpha <- log( averageNumberPresent + 1)

#the number of each species
alphas <- rnorm( n=S, mean=average.alpha)  #Note that the mean of back-transformed won't be average.mean.prev (but in the ball park)
#the number within each group
taus <- matrix( rnorm( n=K*S), ncol=K, nrow=S)
taus <- t( apply( taus, 1, function(x) x - mean( x)))
lambdas <- exp( kronecker( matrix( 1, ncol=K, nrow=1), alphas) + taus)

####  Define the OBSERVED K thinned processes.
n.of.each.spp <- matrix( rpois( K*S, lambdas), ncol=K, nrow=S)
tmp.combs <- expand.grid( spp=1:S, rcp=1:K)
loc.of.each.spp <- tmp.combs[rep( 1:(S*K), n.of.each.spp),]
loc.of.each.spp$ID <- NA
for( ss in 1:S){
  for( kk in 1:K){
    if( n.of.each.spp[ss,kk] > 0){
      pts.ids <- which( apply( loc.of.each.spp, 1, function(x) all( x[1:2]==c(ss,kk))))
      loc.of.each.spp[pts.ids,"ID"] <- sample( x=1:nrow( X), size=n.of.each.spp[ss,kk], prob=pis[,kk])
    }
  }
}
loc.of.each.spp <- cbind( loc.of.each.spp, X[loc.of.each.spp$ID,c("x","y")])

#the locations of the thinned values match the pis
par( mfrow=c(2,3))
for( kk in 1:K){
  plot( piBrick[[kk]])
  tmpDat <- loc.of.each.spp[loc.of.each.spp$rcp==kk,]
  text( tmpDat[,c("x","y")], labels=tmpDat$spp, cex=0.5, col=tmpDat$spp)
}

#The distribution of spp is spread between RCPs
par( mfrow=c(5,4))
for( ss in 1:S){
  plot( extent( piBrick[[1]]), ylab="", xlab="")
  tmpDat <- loc.of.each.spp[loc.of.each.spp$spp==ss,]
  text( tmpDat[,c("x","y")], labels=tmpDat$rcp, cex=0.5, col=tmpDat$rcp)
}

sppDat <- loc.of.each.spp

rm( list=ls()[! ls() %in% c("sppDat", "X", "Xbrick", "pis", "piBrick", "alphas", "taus", "lambdas", "betas")])


library( "qrbp")

bkpts <- ppmData( npoints = 10000,
                    presences = sppDat[,c("x","y","spp")],
                    window = Xbrick[[1]],
                    covariates = Xbrick,
                    resolution = 0.00001,
                    method = "quasirandom",
                    interpolation = "bilinear",
                    coord = c("X","Y"),
                    control = ppmData.control( cores=7))

