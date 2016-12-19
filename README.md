# SameSizeClustering

A function for creating tight, same sized clusters, based on a 'distance trade' algorithm. 
Used for geographic applications, in which distance determined by latitude and longitude 
is being minimized. However, this can be modified to include different distance dimensions.

The base algorithm was developed by Wes Stevenson (http://statistical-research.com/page/3/), 
but was made scalable and more optimized by incorporating an initial kmeans split of the data
into subgroups. This reduces distance between points by significantly decreasing the "worst case matching"
in which the starting seeds are very far apart.

Additionally, this can allow for parallelization of the algorithm (although this isn't supported at this time).

It returns a dataframe with several new columns:
  1. BigCluster -- which initial kmeans group was it included in
  2. subCluster -- for each initial kmeans group, which subcluster was it located in after optimization
  3. clusterOriginal -- for each initial kmeans group, this was the initialized subcluster
  4. distance.from.center -- the distance of the point from the center of the final, optimized cluster
  5. distance.from.center.original -- the distance of the point from the initialized center
  6. cluster_final -- the name of the final cluster, leading with a letter to preserve data integrity

# Example Usage

```R
# Generate sample data
sample_data_distribution <- function(Mu1, Sigma1, Mu2, Sigma2, n1, n2) {
  library(MASS)
  dist1 <- data.frame(mvrnorm(n = n1, Mu1, Sigma1))
  names(dist1) <- c("y", "x")
  dist2 <- data.frame(mvrnorm(n = n2, Mu2, Sigma2))
  names(dist2) <- c("y", "x")
  dist <- rbind(dist1, dist2)
  dist$UserID <- 1:length(dist$x)
  return(dist)
}

# Generate sample data to replicate pop distribution of Chicago

sample_data_users_chicago <- sample_data_distribution(c(41.855970, -87.68684), 
                              matrix(c(0.01125,-0.005,-.005,0.002975),2,2), 
                              c(41.855970, -87.76684),
                              matrix(c(0.003,0,-0.08,0.003),2,2),
                              ceiling(.7*2000),
                              ceiling(.3*2000))

head(sample_data_users_chicago)
### y           x           UserID
### 41.75815    -87.64354     1
### 41.92560    -87.69849     2
### 41.73634    -87.61671     3


distance_clustered_data <- SameSizeClustering(sample_data_users_chicago)

head(distance_clustered_data)

###     y            x          UserID         BigCluster     subCluster      clusterOriginal    distance.from.center
### 1   41.75815      -87.64354      1              1           39               26                0.2726409
### 2   41.92560      -87.69849      2              8           51               77                0.2513240
### 3   41.73634      -87.61671      3              1           26               36                0.1123567
 
###   distance.from.center_original     cluster_final
### 1                     1.2549047              a39
### 2                     0.6990344              h51
### 3                     0.9942355              a26


```
