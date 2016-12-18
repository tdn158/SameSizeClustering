SameSizeClusters <- function(orig.data, group_size=5, max.iter =50, tolerance = 1, plot.iter=FALSE) {
  require(dplyr)
  
  #### Citation for base algorithm: http://statistical-research.com/page/3/
  #### My contribution adds scalability and further minimizes distance than baseline algorithm.
  #### This is a function that takes a dataframe with lat / lon 
  #### (must be named 'y' and 'x' respectively) coordinates and ID (must be named 'UserID')
  #### and initially splits into nrow(data)/250 clusters for scalability.
  #### Then, for every subcluster, the algorithm
  #### 1. Initializes n/k clusters with groups of size at most k
  #### 2. Minimizes distance between the groups by
  ###### a. For every observation, the closest cluster center is found (k_take). 
  ######    This observation now belongs to this cluster (k_take).
  ###### b. The cluster that gave the observation (k_give) receives the observation
  ######    closest to its center from k_take
  ###### c. This is repeated for every obs multiple times until convergence or max.iter is reached
  
  
  # Check conditions
  if (!is.data.frame(orig.data)) {
    stop('argument [orig.data] must be a data frame')
  }
  if (!require(dplyr)) {
    stop('packages [dplyr] are required for this function')
  }
  if (!("x" %in% names(orig.data) && "y" %in% names(orig.data) && "UserID" %in% names(orig.data))) {
    stop('requires latitude column named "y", longitude column named "x", and an ID column named "UserID"')
  }
      
  # Letters needed for cluster key, a measure of data integrity
  letters <- 'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ'
  
  # Set number of initial clusters
  if (nrow(orig.data) < 1000) {
    bigK <- 4
  }
  else {
    bigK <- ceiling(nrow(orig.data)/250)
  }
  
  kmeansFit <- kmeans(dplyr::select(orig.data, x, y), centers = bigK, nstart = 10)
  orig.data$BigClusters <- kmeansFit$cluster
  
  for (cluster in 1:bigK) {
    
    cat("\n\n\n ------------- Cluster ", cluster, 'of ', bigK, ' ------------- ')
    if (cluster == 1) {
      tmp.data.old <- NULL
    }
    else {
      tmp.data.old <- tmp.data.bind
    }
    
    tmp.data <- subset(orig.data, BigClusters == cluster)
    fr = to = NULL
    
    #build initial clusters from sampled points from the data.
    number_groups = ceiling(nrow(tmp.data)/group_size)
    r.k.start = sample(seq(1:number_groups))
    n = nrow( tmp.data )
    initial.clusters = rep(r.k.start, group_size)
    
    if(n%%length(initial.clusters)!=0){
      exclude.k = length(initial.clusters) - n%%length(initial.clusters)
    } else {
      exclude.k = 0
    }
    tmp.data$cluster = initial.clusters[1:(length(initial.clusters)-exclude.k)]
    tmp.data$cluster_original = tmp.data$cluster
    
    ## Calc centers and merge
    mu = cbind( by(tmp.data$y, tmp.data$cluster, mean), by(tmp.data$x, tmp.data$cluster, mean), seq(1:number_groups) )
    tmp1 = matrix( match(tmp.data$cluster, mu[,3]) )
    orig.data.centers = cbind(as.matrix(tmp.data), mu[tmp1,])
    
    ## Calc initial distance from centers
    fr$lat = as.numeric(orig.data.centers[,ncol(orig.data.centers)-2])
    fr$lon = as.numeric(orig.data.centers[,ncol(orig.data.centers)-1])
    fr <- as.data.frame(fr)
    to$lat = as.numeric(orig.data.centers[,1])
    to$lon = as.numeric(orig.data.centers[,2])
    to <- as.data.frame(to)
    tmp.data$distance.from.center = calc_dist(fr, to)$distance_miles
    tmp.data$distance.from.center_original = tmp.data$distance.from.center
    
    ## Set some initial configuration values
    is.converged = FALSE
    iteration = 0
    error.old = Inf
    error.curr = Inf
    
    while ( !is.converged && iteration < max.iter ) { # Iterate until threshold or maximum iterations
      
      if(plot.iter==TRUE){
        plot(tmp.data$x, tmp.data$y, col=tmp.data$cluster, pch=16, cex=.6,
             xlab="Longitude",ylab="Latitude")
      }
      
      iteration = iteration + 1
      start.time = as.numeric(Sys.time())
      cat("\nIteration ", iteration,sep="")
      for( i in 1:n ) {
        # Iterate over each observation and measure the distance each observation' from its mean center
        # Produces an exchange. It takes the observation closest to it's mean and in return it gives the observation
        # closest to the giver, k, mean
        fr = to = distances = NULL
        for( j in 1:number_groups ){
          # Determine the distance from each k group
          fr$lat = tmp.data$y[i]; fr$lon = tmp.data$x[i]
          to$lat = mu[j,1]; to$lon = mu[j,2]
          distances[j] = as.numeric( calc_dist(fr, to) )
        }
        
        # Which k cluster is the observation closest.
        which.min.distance = which(distances==min(distances), arr.ind=TRUE)
        previous.cluster = tmp.data$cluster[i]
        tmp.data$cluster[i] = which.min.distance # Replace cluster with closest cluster
        
        # Trade an observation that is closest to the giving cluster
        if(previous.cluster != which.min.distance){
          new.cluster.group = tmp.data[tmp.data$cluster==which.min.distance,]
          fr$lat = mu[previous.cluster,1]; fr$lon = mu[previous.cluster,2]
          to$lat = new.cluster.group$y; to$lon = new.cluster.group$x
          new.cluster.group$tmp.dist = calc_dist(fr, to)$distance_miles
          
          take.out.new.cluster.group = which(new.cluster.group$tmp.dist==min(new.cluster.group$tmp.dist), arr.ind=TRUE)
          UserID = new.cluster.group$UserID[take.out.new.cluster.group]
          tmp.data$cluster[tmp.data$UserID == UserID] = previous.cluster
        }
        
      }
      
      # Calculate new cluster means
      mu = cbind( by(tmp.data$y, tmp.data$cluster, mean), by(tmp.data$x, tmp.data$cluster, mean), seq(1:number_groups) )
      tmp1 = matrix( match(tmp.data$cluster, mu[,3]) )
      orig.data.centers = cbind(as.matrix(tmp.data), mu[tmp1,])
      mu = cbind( by(tmp.data$y, tmp.data$cluster, mean), by(tmp.data$x, tmp.data$cluster, mean), seq(1:number_groups) )
      
      ## Calc initial distance from centers
      fr$lat = as.numeric(orig.data.centers[,ncol(orig.data.centers)-2]) 
      fr$lon = as.numeric(orig.data.centers[,ncol(orig.data.centers)-1])
      fr <- as.data.frame(fr)
      to$lat = as.numeric(orig.data.centers[,1])
      to$lon = as.numeric(orig.data.centers[,2])
      to <- as.data.frame(to)
      tmp.data$distance.from.center = calc_dist(fr, to)$distance_miles
      
      # Test for convergence. Is the previous distance within the threshold of the current total distance from center
      error.curr = sum(tmp.data$distance.from.center)
      
      error.diff = abs( error.old - error.curr )
      error.old = error.curr
      if( !is.nan( error.diff ) && error.diff < tolerance ) {
        is.converged = TRUE
      }
      # Set a time to see how long the process will take is going through all iterations
      tmp.data$cluster_kmeans <- paste0(substr(letters, cluster, cluster), tmp.data$cluster, sep = '')
      stop.time = as.numeric(Sys.time())
      hour.diff = (((stop.time - start.time) * (max.iter - iteration))/60)/60
      cat("\n Error ",error.diff," Hours remain from iterations ",hour.diff,"\n")
      
    }
    tmp.data.bind <- rbind(tmp.data.old, tmp.data)

  }
  
    return(tmp.data.bind)
}
  




