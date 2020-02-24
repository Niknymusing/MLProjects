library(foreign)

library(EMD)

library(foreach)

library(doParallel)

registerDoParallel(4)

# For parallelizing computations using the %dopar% operator bellow.

data <- read.dta("/Users/nilskakoseosnystrom/Downloads/spotmarket.dta")

sig <- data[,28]

sigemd <- emd(sig, tt=NULL, tol=sd(sig)*0.1^2, max.sift=20, stoprule="type1",

boundary="periodic", sm="none", smlevels=c(1), spar=NULL, alpha=NULL,

check=FALSE, max.imf=10, plot.imf=FALSE, interm=NULL, weight=NULL)

# Calculates the Hilbert-Huang transform of the signal, which gives a decomposition of the signal 
# into a number of (10 in this case) different "IMF" components, which are characteristic oscillatiory components present 
# in the signal) and a residue, as described in [N.Huang 2003, "Applications of Hilbert–Huang transform to non‐stationary 
# financial time series analysis"]. Summing the IMFs and the residue reconstruct the signal completely. 
# By decomposing the signal in this way the predictive analysis is more robust to noise, non-linear and non-stationary 
# characteristics of the signal.

j=960 
# length of forecast, look-ahead timewindow is j hours/pricepoints.

l=24 
# length of look-back timewindow, for computing the k nearest neighboors of consecutive length-j segments.

k=1060 
# number of nearest neigboors (number of weighted training examples to add together to construct the forecast)

forecast <- vector(mode = "numeric", length=j)

foreach(i in 1:10) %dopar% {

sigimf <- sigemd$imf[,i]
query_timewindow <- sigimf[(length(sigimf) - j - l):(length(sigimf) - j - 1)] 
Training_data <- matrix(ncol = l, nrow = 0)

for(i in 1:(length(sigimf) - j - l)) {        

Training_data <- rbind(Training_data, t(sigimf[(i):(i + l - 1)]))
}
# Constructs training data set (for iteratively selecting the K nearest neigbhoors based on distance 
# to query_timewindow) from l-segments of the given timeseries.

query_dist <- matrix(ncol = 0, nrow = 1)


for(i in 1:(length(sigimf) - j - 2*l)) {

D <- dist(rbind(Training_data[i,], t(query_timewindow)))        
query_dist <- cbind(query_dist,  D)         
 
}
 
# Computes the euclidean distance of [all l-segments in the training dataset] and [the last segment of 
# the timseries (query_timewindow above)], and puts the resulting distance-coefficients in the query_dist object. 

indices_of_k_nearest_neighboors <- which(query_dist[1,] %in% sort(query_dist[1,])[1:k])

# Indices in the training dataset of the k nearest neighboors to query_timewindow

f <- vector(mode = "numeric", length=j)

for(i in 1:k) {

q <- 1/(dist(rbind(Training_data[indices_of_k_nearest_neighboors[i],], t(query_timewindow))) + 0.01)
f <- f + q*sigimf[(indices_of_k_nearest_neighboors[i] + l):(indices_of_k_nearest_neighboors[i] + l + j - 1)]

}
# Computes the distance-weights for each of the k nearest neigboors and sums the weighted k nearest neigbhoors.
 
c = 0

for(i in 1:k) {

c <-  c + 1/(dist(rbind(Training_data[indices_of_k_nearest_neighboors[i]], t(query_timewindow))) + 0.01) 

}
# Computes normalisation constant to form a convex combination of weighted k-nearest neighboors.


forecast <- forecast + (1/c)*f


} 

# The lines bellow gives the same computations as above but applied to the residue component of the EMD, 
# which added to the forecast gives the resulting prediction of the model.

sigres <- sigemd$residue

query_timewindow <- sigres[(length(sigres) - j - l):(length(sigres) - j - 1)] 

Training_data <- matrix(ncol = l, nrow = 0)


for(i in 1:(length(sigres) - j - l)) {        

Training_data <- rbind(Training_data, t(sigres[(i):(i + l - 1)]))

}

L <- matrix(ncol = l, nrow = 0)

query_dist <- matrix(ncol = 0, nrow = 1)

for(i in 1:(length(sigimf) - j - 2*l)) {


 D <- dist(rbind(Training_data[i,], t(query_timewindow)))        
query_dist <- cbind(query_dist,  D)         

}

indices_of_k_nearest_neighboors <- which(query_dist[1,] %in% sort(query_dist[1,])[1:k])

f <- vector(mode = "numeric", length=j)

for(i in 1:k) {

q <- 1/(dist(rbind(Training_data[indices_of_k_nearest_neighboors[i],], t(query_timewindow))) + 0.01)
f <- f + q*sigres[(indices_of_k_nearest_neighboors[i]+l):(indices_of_k_nearest_neighboors[i] + l + j - 1)]

}
 
c = 0

for(i in 1:k) {

c <-  c + 1/(dist(rbind(Training_data[indices_of_k_nearest_neighboors[i]], t(query_timewindow))) + 0.01) 

} 
 
forecast <- forecast + (1/c)*f
