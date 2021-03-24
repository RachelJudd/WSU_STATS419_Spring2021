github.monte.raw = "https://raw.githubusercontent.com/MonteShaffer/";
include.setup = paste0(github.monte.raw, "humanVerse/main/include.setup.R");
source(include.setup);  # Maybe comment this out, if error 
 
github.monte.http = "https://github.com/MonteShaffer/";
######## we will parse this page to get a list of the .R functions to include ########
github.monte.humanVerse = paste0(github.monte.http, "humanVerse/tree/main/humanVerse/R/"); 

######## you can pass flag `force.download = TRUE` if you want to make certain it is not coming from cache ########

###### R::humanVerse #####
includeGithubFolder(github.monte.humanVerse, force.download = TRUE); 
 
  


folder.owell = paste0(github.monte.raw, "humanVerse/main/data/o-well/");

file.location = paste0(folder.owell, "wells-location-plus.txt");

wells.pipe = readFromPipe(file.location);  

wells.pipe = replaceFactorColumnWithIndicatorVariables(wells.pipe, "geology", use.boolean=FALSE);

wells.pipe;








file.location = paste0(folder.owell, "well-23.R");
source(file.location);

names(wells);

# "sodium" is Table 6

Wilcox = wells$sodium$Wilcox;
Wilcox = str_replace("C","", Wilcox);
Wilcox = str_replace("S","", Wilcox);
tmp = explodeMe("-",Wilcox,NULL);
Wilcox.C = as.numeric( getElementsInList(tmp, 1) ); 
Wilcox.S = as.numeric( getElementsInList(tmp, 2) ); 

wells$sodium$Wilcox.C = Wilcox.C;
wells$sodium$Wilcox.S = Wilcox.S;

wells;









wells.df = wells.pipe;
  wells.df = merge(wells.df, wells$metals, by="well");
  wells.df = merge(wells.df, wells$chem, by="well");
  wells.df = merge(wells.df, wells$sodium, by="well");

dim(wells.df); 
nrow(wells.df);
ncol(wells.df);

wells.df;  # do we have any variable collision?
  
library(geosphere);  
library(measurements); 
library(pracma);    
# http://www.csgnetwork.com/degreelenllavcalc.html 
# default "mi" is miles as in statute miles (5280 feet)
# GREAT CIRCLE ... http://en.wikipedia.org/wiki/Great_circle_distance
# ??distHaversine 
convertLatitudeLongitudeToDistance = function(lat, format="mi")
  {
	m1 = 111132.92; m2 = -559.82; m3 = 1.175; m4 = -0.0023;
	p1 = 111412.84; p2 = -93.5;   p3 = 0.118;
		
	lat = pracma::deg2rad(lat);
	# lon = pracma::deg2rad(lon);
	# note:  longitude is not required
	# below computes in meters
	latlen  = m1 + (m2 * cos(2 * lat)) + (m3 * cos(4 * lat)) +	(m4 * cos(6 * lat));
	lonlen = (p1 * cos(lat)) + (p2 * cos(3 * lat)) + (p3 * cos(5 * lat));
	
	res = list();
	  res$format     = format;
	  res$lat.factor = measurements::conv_unit( latlen,  "m", format);
	  res$lon.factor = measurements::conv_unit( lonlen,  "m", format);
	res;
  }

latlon = convertLatitudeLongitudeToDistance( mean(wells.df$latitude) );
wells.df$lat.mi = latlon$lat.factor * wells.df$latitude;
wells.df$lon.mi = latlon$lon.factor * wells.df$longitude;


# names(wells.df);
# d1 = conv_unit(  distm( wells.df[,3:2], fun=distHaversine) ,  "m", "mi");
# d2 = as.matrix ( dist( wells.df[,55:56], 
#               method="euclidean", diag=TRUE, upper=TRUE) );
# sum( d2 - d1 ) / (23 * 23);



library(measurements); # install.packages("measurements");
conv_unit(2.54, "cm", "inch");

wells.df$altitude.mi = conv_unit(wells.df$altitude.ft, "ft", "mi");
wells.df$sea.mi = conv_unit(wells.df$sea.m, "m", "mi");
wells.df$fault.mi = conv_unit(wells.df$fault.m, "m", "mi");
  