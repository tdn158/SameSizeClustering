SameSizeClusters <- function(orig.data, k=5, max.iter =50, tolerance = 1, plot.iter=FALSE) {
  
  #### This is a function that takes a dataframe with lat / lon 
  #### (must be named 'lat' and 'lon') coordinates and
  #### 1. Initializes groups of size at most k
  #### 2. Minimizes distance between the groups by
  ###### a. For every observation, the closest cluster center is found (k_take). 
  ######    This observation now belongs to this cluster (k_take).
  ###### b. The cluster that gave the observation (k_give) receives the observation
  ######    closest to its center from k_take
  ###### c. This is repeated until convergence or max.iter is reached
  
  fr = to = NULL
  
  #build initial clusters from sampled points from the data.
  r.k.start = sample(seq(1:k))
  n = nrow( orig.data )
  k.size = ceiling(n/k)
  initial.clusters = rep(r.k.start, k.size)
  
  if(n%%length(initial.clusters)!=0){
    exclude.k = length(initial.clusters) - n%%length(initial.clusters)
  } else {
    exclude.k = 0
  }
  orig.data$cluster = initial.clusters[1:(length(initial.clusters)-exclude.k)]
  orig.data$cluster_original = orig.data$cluster
  
  ## Calc centers and merge
  mu = cbind( by(orig.data$y, orig.data$cluster, mean), by(orig.data$x, orig.data$cluster, mean), seq(1:k) )
  tmp1 = matrix( match(orig.data$cluster, mu[,3]) )
  orig.data.centers = cbind(as.matrix(orig.data), mu[tmp1,])[,c(1:2,4:6)]
  
  ## Calc initial distance from centers
  fr$lat = orig.data.centers[,3]; fr$lon = orig.data.centers[,4]
  to$lat = orig.data.centers[,1]; to$lon = orig.data.centers[,2]
  orig.data$distance.from.center = calc_dist(fr, to)$distance_miles
  orig.data$distance.from.center_original = orig.data$distance.from.center
  
  ## Set some initial configuration values
  is.converged = FALSE
  iteration = 0
  error.old = Inf
  error.curr = Inf
  
  while ( !is.converged && iteration < max.iter ) { # Iterate until threshold or maximum iterations
    
    if(plot.iter==TRUE){
      plot(orig.data$x, orig.data$y, col=orig.data$cluster, pch=16, cex=.6,
           xlab="Longitude",ylab="Latitude")
    }
    iteration = iteration + 1
    start.time = as.numeric(Sys.time())
    cat("Iteration ", iteration,sep="")
    for( i in 1:n ) {
      # Iterate over each observation and measure the distance each observation' from its mean center
      # Produces an exchange. It takes the observation closest to it's mean and in return it gives the observation
      # closest to the giver, k, mean
      fr = to = distances = NULL
      for( j in 1:k ){
        # Determine the distance from each k group
        fr$lat = orig.data$y[i]; fr$lon = orig.data$x[i]
        to$lat = mu[j,1]; to$lon = mu[j,2]
        distances[j] = as.numeric( calc_dist(fr, to) )
      }
      
      # Which k cluster is the observation closest.
      which.min.distance = which(distances==min(distances), arr.ind=TRUE)
      previous.cluster = orig.data$cluster[i]
      orig.data$cluster[i] = which.min.distance # Replace cluster with closest cluster
      
      # Trade an observation that is closest to the giving cluster
      if(previous.cluster != which.min.distance){
        new.cluster.group = orig.data[orig.data$cluster==which.min.distance,]
        
        fr$lat = mu[previous.cluster,1]; fr$lon = mu[previous.cluster,2]
        to$lat = new.cluster.group$y; to$lon = new.cluster.group$x
        new.cluster.group$tmp.dist = calc_dist(fr, to)$distance_miles
        
        take.out.new.cluster.group = which(new.cluster.group$tmp.dist==min(new.cluster.group$tmp.dist), arr.ind=TRUE)
        UserID = new.cluster.group$UserID[take.out.new.cluster.group]
        orig.data$cluster[orig.data$UserID == UserID] = previous.cluster
      }
      
    }
    
    # Calculate new cluster means
    mu = cbind( by(orig.data$y, orig.data$cluster, mean), by(orig.data$x, orig.data$cluster, mean), seq(1:k) )
    tmp1 = matrix( match(orig.data$cluster, mu[,3]) )
    orig.data.centers = cbind(as.matrix(orig.data), mu[tmp1,])[,c(1:2,4:6)]
    mu = cbind( by(orig.data$y, orig.data$cluster, mean), by(orig.data$x, orig.data$cluster, mean), seq(1:k) )
    
    ## Calc initial distance from centers
    fr$lat = orig.data.centers[,3]; fr$lon = orig.data.centers[,4]
    to$lat = orig.data.centers[,1]; to$lon = orig.data.centers[,2]
    orig.data$distance.from.center = calc_dist(fr, to)$distance_miles
    
    # Test for convergence. Is the previous distance within the threshold of the current total distance from center
    error.curr = sum(orig.data$distance.from.center)
    
    error.diff = abs( error.old - error.curr )
    error.old = error.curr
    if( !is.nan( error.diff ) && error.diff < tolerance ) {
      is.converged = TRUE
    }
    
    # Set a time to see how long the process will take is going through all iterations
    stop.time = as.numeric(Sys.time())
    hour.diff = (((stop.time - start.time) * (max.iter - iteration))/60)/60
    cat("\n Error ",error.diff," Hours remain from iterations ",hour.diff,"\n")
    
    # Write out iterations. Can later be used as a starting point if iterations need to pause
    #write.table(orig.data, paste("C:\\optimize_iteration_",iteration,"_data.csv", sep=""), sep=",", row.names=F)
  }
  
  centers = data.frame(mu)
  ret.val = list("centers" = centers, "cluster" = factor(orig.data$cluster), "UserID" = orig.data$UserID,
                 "Latitude" = orig.data$y, "Longitude" = orig.data$x, "Distance" = orig.data$distance.from.center,
                 "k" = k, "iterations" = iteration, "error.diff" = error.diff)
  
  return(ret.val)
}