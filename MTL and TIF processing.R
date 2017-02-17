
# purpose: batch process GeoTIFF metadata and satellite band files
# these files contain the metadata for the GeoTIFF image files 
# used to calculate NDVI
# author: Andy Lathrop, BlueGranite, Inc.

library(tidyverse)
library(lubridate)
library(stringr)
library(parallel)
library(doParallel)
library(raster)
library(rgdal)

# list and check files
#####
data_path <- "E:/Data"
processed_path <- "E:/XDF/band-unfiltered"

MTL_IDs   <- str_sub(list.files(data_path, "*.txt"),    -29, -9)
TIF_4_IDs <- str_sub(list.files(data_path, "*_B4.TIF"), -29, -8)
TIF_5_IDs <- str_sub(list.files(data_path, "*_B5.TIF"), -29, -8)

# check if all the vectors of IDs are the same
all.equal(MTL_IDs, TIF_4_IDs, TIF_5_IDs)

## read and transform metadata files
#####
# example from http://serialmentor.com/blog/2016/6/13/reading-and-combining-many-tidy-data-files-in-R

MTL_data_path <- "E:Data"

# find all file names ending in .txt
files <- dir(MTL_data_path, pattern = "*.txt")

# process individual metadata files to data frame
start.time <- proc.time()
  # create a data frame holding the file names
  MTL_data <- data_frame(filename = files) %>%  
    
    # read file contents into a new data column (nested data frame)
    mutate(file_contents = map(filename,    
                               ~ read_delim(file.path(MTL_data_path, . ),
                                            delim = "=",
                                            col_names = c("myColName", "myValue"),
                                            trim_ws = TRUE))) %>% 
    
    # unnest to use in regular data frame format, then transpose with 'spread'
    unnest %>%
    
    # filter out columns used to identify column hierarchies
    filter(!(myColName %in% c("END_GROUP", "END", "GROUP"))) %>% 
    
    # transpose
    spread(key = myColName, value = myValue) %>% 
    
    # DATE_ACQUIRED and FILE_DATE as date column 
    mutate(DATE_ACQUIRED = as_date(DATE_ACQUIRED)) %>% 
    mutate(FILE_DATE     = as_date(FILE_DATE)) %>% 
    
    # remove leading and training quote marks
    mutate(LANDSAT_SCENE_ID = str_sub(LANDSAT_SCENE_ID,2, -2)) %>% 
    
    # just the columns we need
    dplyr::select(LANDSAT_SCENE_ID, DATE_ACQUIRED, FILE_DATE)
proc.time() - start.time
# end

############

# read GeoTIFF files to prepare for modeling
#####

TIF_data_path <- "E:/Data"

# account for corrupt and already processed IDs
processed_B4_IDs <- str_sub(unique(dir(processed_path, pattern = "*_B4.xdf")), 1, -8)
processed_B5_IDs <- str_sub(unique(dir(processed_path, pattern = "*_B5.xdf")), 1, -8)

# remove corrupt and processed IDs
B4_IDs_to_process <- MTL_IDs[!(MTL_IDs %in% processed_B4_IDs)]
B4_IDs_to_process <- B4_IDs_to_process[!(B4_IDs_to_process %in% corrupt_IDs)]

B5_IDs_to_process <- MTL_IDs[!(MTL_IDs %in% processed_B5_IDs)]
B5_IDs_to_process <- B5_IDs_to_process[!(B5_IDs_to_process %in% corrupt_IDs)]

# filenames for all band 4 and band 5 .TIF files needed to calculate NDVI
files <- c(paste0(B5_IDs_to_process, "_B5.TIF"),
           paste0(B4_IDs_to_process, "_B4.TIF"))  

##
# read and transform GeoTIFF band files to points, write to XDF
# make sure to check output folder on last line of function
# make sure to check filtering
# based on example at 
# http://neondataskills.org/R/Extract-NDVI-From-Rasters-In-R/
# and
# https://geoscripting-wur.github.io/IntroToRaster/

TIF_to_points <- function(myFile) {
  
  myFile %>% 
    map(stack) %>% 
    map(overlay, fun = min) %>% 
    map(rasterToPoints) %>% 
    data.frame %>% 
    as_tibble %>% 
    
    # filter to subset of pixels
    filter(x > 625000 & x < 625200 & y > 3650000 & y < 3650500) %>% 
    
    # keep only layer values > 0 (remove empty space in image) 
    filter(layer > 0) %>% 
    
    # create ID column
    mutate(LANDSAT_SCENE_ID = str_sub(myFile, -28, -8)) %>% 
    
    # create new cols for x and y since unite will remove original cols
    # mutate(pixel_x = x, pixel_y = y) %>% 
    # unite (x, y, col = "pixelID") %>% 
    
    # get band value from filename
    mutate(band = str_sub(myFile, -6, -5)) %>% 
    
    # mutate(write_time = Sys.time()) %>%
    
    # join TIF and meatadata to add image date column
    left_join(MTL_data) %>% 
    
    # write data to XDF
    rxDataStep(paste0("E:/XDF/band-unfiltered/", str_sub(myFile, -28, -5), ".xdf"),
               overwrite = TRUE)
  
} 

##
# iterate processing over all TIF files
# benchmark 186s for 1 TIFF on DSVM 56GB mem - filtered x,y
# benchmark 337s for 2 TIFFs on DSVM 56GB mem - filtered x,y
# benchmark 2964s (50min) for 16 TIFFs on DSVM 56GB mem - filtered x,y
# benchmark 3321s (55min) for 18 TIFFs on DSVM 56GB mem - filtered x,y
# benchmark 654s (10min) for 2 TIFFs on local 16GB - not filtered
# benchmark 4642s (77min) for 16 TIFFs on DSVM - not filtered
# benchmark 4878s (81min) for 18 TIFFs on DSVM - not filtered
# benchmark 312s (5.2min) for 4 TIFFs on DSVM *parallel* - not filtered
start.time <- proc.time()
  files %>%  
    map(~ TIF_to_points(paste0(TIF_data_path, "/", files)))
proc.time() - start.time 

# check if file was written correctly
rm(myDataFrame)
myDataFrame <- rxImport(inData = "E:/XDF/band-filtered/LC80390372017002LGN00_B4.xdf")

#####

# R Server parallel processing
#####
# how many cores are available
computeCores <- detectCores()

registerDoParallel(cores=computeCores)
rxSetComputeContext(RxForeachDoPar())

# TIFF to XDF
start.time <- proc.time()
  rxExec(TIF_to_points, 
         myFile = rxElemArg(paste0(TIF_data_path, "/", files)), 
         elemType = "cores",
         # taskChunkSize=length(i)/computeCores,
         packagesToLoad = c("tidyverse",
                            "dplyr",
                            "lubridate",
                            "stringr",
                            "parallel",
                            "doParallel",
                            "RevoScaleR",
                            "rgdal",
                            "raster"),
         execObjects = c("TIF_to_points","files", "MTL_data"))
proc.time() - start.time 

stopImplicitCluster()

# end