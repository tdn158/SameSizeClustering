# SameSizeClustering

A function for creating fixed size clusters, based on a 'distance trade' algorithm. 
Used for geographic applications, in which distance determined by latitude and longitude 
is being minimized. However, can be modified to include different distance dimensions.

The base algorithm was developed by Wes Stevenson (http://statistical-research.com/page/3/), 
but was made scalable and more optimized by incorporating an initial kMeans split of the data
into subgroups. This reduces distance between points by significantly decreasing the "worst case matching"
in which the starting seeds are very far apart.

Additionally, this can allow for parallelization of the algorithm (although this isn't supported at this time).


```R
x <- sample(nrow(data), 500)
```
