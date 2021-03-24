github.monte.raw = "https://raw.githubusercontent.com/MonteShaffer/";
include.setup = paste0(github.monte.raw, "humanVerse/main/include.setup.R");
source(include.setup);

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
  
  
  
  
  
# https://gis.stackexchange.com/questions/142326/calculating-longitude-length-in-miles#:~:text=Each%20degree%20of%20latitude%20is,111.699%20km)%20at%20the%20poles.
getMileMarkerFromLatitude = function(lat, factor=69)
  {
  # is there a formula that would figure this out? I could just linear extract
  # https://gis.stackexchange.com/questions/110730/mercator-scale-factor-is-changed-along-the-meridians-as-a-function-of-latitude
  # MERCATOR
  factor.equator  = 68.703;
  factor.poles    = 69.407;
  
  lat*factor;
  }
  

library(pracma);
getMileMarkerFromLongitude = function(long, lat, factor=69.172)
  {
  # https://gis.stackexchange.com/questions/142326/calculating-longitude-length-in-miles
  
  one = 1 / ( factor * cos(deg2rad(lat)) );  # this is one degree
  long * one;
  }

wells.df$lat.mi = getMileMarkerFromLatitude(wells.df$latitude);
wells.df$lon.mi = getMileMarkerFromLatitude(wells.df$longitude, wells.df$latitude);






library(measurements); # install.packages("measurements");
conv_unit(2.54, "cm", "inch");

wells.df$altitude.mi = conv_unit(wells.df$altitude.ft, "ft", "mi");
wells.df$sea.mi = conv_unit(wells.df$sea.m, "m", "mi");
wells.df$fault.mi = conv_unit(wells.df$fault.m, "m", "mi");
  