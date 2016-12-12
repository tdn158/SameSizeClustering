calc_dist <- function(fr, to) {
  lat1 = fr$lat * (pi/180)
  lon1 = fr$lon * (pi/180)
  lat2 = to$lat * (pi/180)
  lon2 = to$lon * (pi/180)
  a = 3963.191;
  b = 3949.903;
  numerator = ( a^2 * cos(lat2) )^2 + ( b^2 * sin(lat2) ) ^2
  denominator = ( a * cos(lat2) )^2 + ( b * sin(lat2) )^2
  radiusofearth = sqrt(numerator/denominator) #Accounts for the ellipticity of the earth.
  d = radiusofearth * acos( sin(lat1) * sin(lat2) + cos(lat1)*cos(lat2)*cos(lon2 - lon1) )
  d.return = list(distance_miles=d)
  return(d.return)
}
