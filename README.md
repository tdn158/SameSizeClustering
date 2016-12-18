# SameSizeClustering

A function for creating fixed size clusters, based on a 'distance trade' algorithm. 
Used for geographic applications, in which distance determined by latitude and longitude 
is being minimized. However, can be modified to include different distance dimensions.

The base algorithm was developed by Wes Stevenson (http://statistical-research.com/page/3/), 
but was made scalable and more optimized by incorporating an initial kMeans split of the data
into subgroups. This reduces distance between points by significantly decreasing the "worst case matching"
in which the starting seeds are very far apart.

Additionally, this can allow for parallelization of the algorithm (although this isn't supported at this time).

# Example Usage

```R
# Generate sample data
sample_data_chicago_distribution <- function(Mu1, Sigma1, Mu2, Sigma2, n1, n2) {
  library(MASS)
  dist1 <- data.frame(mvrnorm(n = n1, Mu1, Sigma1))
  names(dist1) <- c("y", "x")
  dist2 <- data.frame(mvrnorm(n = n2, Mu2, Sigma2))
  names(dist2) <- c("y", "x")
  dist <- rbind(dist1, dist2)
  dist$UserID <- 1:length(dist$x)
  return(dist)
}

sample_data_users_chicago <- sample_data_chicago_users(c(41.855970, -87.68684), 
                              matrix(c(0.01125,-0.005,-.005,0.002975),2,2), 
                              c(41.855970, -87.76684),
                              matrix(c(0.003,0,-0.08,0.003),2,2),
                              ceiling(.7*2000),
                              ceiling(.3*2000))

distance_clustered_data <- SameSizeClustering(sample_data_users_chicago)
```
