#current version: 20/08/2024

# --------------------------------------------------------------------------------

# Developed and owned by Flanders Marine Institute (VLIZ)
# Please do not change anything in this file
# For any questions, bugs or requests, please contact francisco.hernandez@vliz.be

# --------------------------------------------------------------------------------


#packages
library(foreach)
library(doParallel)
library(stars)
library(terra)
library(magrittr)
library(dplyr)
library(Rarr)
library(stringr)
library(lubridate)
library(arrow)


#function used to specify an s3 client based on the path, and assuming anonymous access, 
#required for the rarr functions because of a conflict with S3 credentails specified on the data lab virtual machines 
s3_client <- function(endpoint) {paws.storage::s3(
  config = list(
    credentials = list(anonymous = TRUE), 
    region = "auto",
    endpoint = endpoint)
)}


#retrieve STAC catalog in EDITOSTAC
acf <- S3FileSystem$create(scheme = "https", endpoint_override = "s3.waw3-1.cloudferro.com", anonymous = T)
EDITOSTAC = arrow::open_dataset(acf$path("emodnet/edito_stac_cache.parquet")) %>% collect()


#function to return the closest value in a list, needed to select the time and depth slice from the zarr file
closest<-function(xv,sv){
  which(abs(xv-sv)==min(abs(xv-sv))) }


#wrap a gdal request string
toGDAL<-function(fn,par=""){
  if(par!="") url=sprintf('ZARR:"/vsicurl?list_dir=no&retry_delay=60&max_retry=3&url=%s":/%s',fn,par)
  else   url= sprintf('ZARR:"/vsicurl?list_dir=no&retry_delay=60&max_retry=3&url=%s"',fn)
  return(url)
}


#derive info
gdalinfo<-function(fn,par=""){
  url=toGDAL(fn,par)
  jsonlite::fromJSON(gdal_utils('mdiminfo', url, quiet = T))
}


#track and report progress including from forked processes
options("outputdebug"=c('L','M'))
loglist=""
dbl<-function(...){
  m=paste(Sys.time(), ":", ..., sep=" ", collapse= "\t") 
  if(...length()==0) loglist<<-""
  if("L" %in% getOption("outputdebug")) if(class(loglist)=='list') loglist<<-append(loglist,c(m))   else loglist<<-c(m)
  if("P" %in% getOption("outputdebug")) print(m)
  if("M" %in% getOption("outputdebug")) message(m)
}


#used in parallel processing; drop results of incomplete parameter instead of losing all data
CombineC<-function(lx,ly)
{
  if(length(lx)==0) return(ly)
  if(is.null(lx$df) || nrow(lx$df)==0)
  { lx$log$progess=paste(lx$log$progess, "first table is empty")
  return(list('log'=append(lx$log,list(ly$log))  , 'df'= ly$df)) 
  }
  if(is.null(ly$df) || nrow(ly$df)==0)
  { 
    return(list('log'=append(lx$log,list(ly$log))  , 'df'= lx$df)) 
  }
  if(nrow(lx$df)!=nrow(ly$df) )
  { ly$log$progess=paste(ly$log$progress, paste("new dataframe has different rowcount ", nrow(lx$df), "><", nrow(ly$df) ), sep='\n')
  return(list('log'=append(lx$log,list(ly$log))  , 'df'= lx$df)) 
  }
  
  newdf=left_join(lx$df,ly$df,by=c('Latitude','Longitude','Time'))
  
  list('log'=append(lx$log,ly$log), 'df'= newdf )
}

# needs error trapping 
CombineR<-function(x,y)
{
  if(is.null(x)) return(y)
  
  dplyr::bind_rows(x,y)
}

#function to retry reading from aws s3, because sometime it just fails, without reason 
#reads json and makes a named list
retryJson<-function(url){

  resp=httr::RETRY(url=url,verb='GET',times=3, encode='json', httr::content_type('json'))
  if(resp$status_code!=200)
    return("")
  
  return(jsonlite::fromJSON(httr::content(resp,'text',encoding = 'UTF-8')))
}

#uses terra::rast
getParFromZarrwInfo<-function(usePar, coords, atTime, atDepth, zinfo, isCategory=F)
{ 
  param = usePar[[1]]
  params = usePar
  
  dsn=toGDAL(zinfo$href)
  
  if(is.na(atDepth)) atDepth=0
  
  #check if parameter is available in zarr file
  if(param %in% names(zinfo$meta$arrays)) 
    sdsn=sdsn_par=sprintf('%s:/%s',dsn,param)
  
  closest_time=NULL
  #is there a time component
  if('time' %in% zinfo$meta$dimensions$name & !is.null(zinfo$times))
  {
    if(length(zinfo$times) > 1)
    {  
      closest_time=closest(zinfo$times,atTime)
      
      sdsn=sprintf('%s:%s',sdsn,closest_time-1)
    }
  }
  
  #is the data in several depth or elevation levels
  closest_level=NULL
  if('elevation' %in% zinfo$meta$dimensions$name & !is.null(zinfo$levels) )
  { 
    if(length(zinfo$levels) > 1)
    {
      closest_level=closest(zinfo$levels,atDepth)
      sdsn=sprintf('%s:%s',sdsn,closest_level-1)
      closest_level = zinfo$levels[closest_level]
    }
  }
  
  r=rast(sdsn)
  
  if(isCategory)  {#r=terra::flip(r, "vertical") --> currently does not work
    #manually flip raster
    r <- read_stars(sdsn,driver = "zarr")
    values_matrix <- as.matrix(r[[1]])
    flipped_matrix <- values_matrix[, ncol(values_matrix):1]
    raster_flipped <- st_as_stars(flipped_matrix)
    st_dimensions(raster_flipped) <- st_dimensions(r)
    r <- st_as_stars(raster_flipped)
    r <- as(r, "SpatRaster")
  }
  
  if(crs(r)=="") {
    dbl("extent and coordinate system missing, assuming epsg:4326, extent from stac catalogue")
    crs(r)='epsg:4326'
    ext(r)=c(zinfo$lonmin[1],zinfo$lonmax[1],zinfo$latmin[1],zinfo$latmax[1])
  }
  
  parTabel=dplyr::tibble()
  
  try({

    #1. no buffer is provided
    if(!"buffer" %in% names(params)) {
      parTabel=terra::extract(x = r, y = dplyr::select(coords,x,y), ID= F, xy=T) 
    } 
    #2. buffer value is provided
    else {  
      bufferSize = as.numeric(params[["buffer"]])
      fun = params[["fun"]]
      bufferedY = terra::buffer(vect(cbind(coords$x, coords$y), crs="+proj=longlat"), bufferSize)
      
      # categorical variable
      # derive the most frequent category for categorical variables
      if(isCategory) {
        if(params[["fun"]] == "most_frequent") {
          dbl("using most frequent value in buffer to look up categorical data")
          
          parTabel = cbind(as.numeric(colnames(par2)[max.col(par2)]), dplyr::select(coords,x,y))}
        else if(params[["fun"]] == "exact") {
          if(!"category" %in% names(params)) {
            dbl("please provide a category to match to")
            stop()
          }
          dbl("calculating percentage area in buffer matching provided category")
          
          par2 = terra::extract(x=r, y=bufferedY, table , na.rm=T)
          row_sums <- rowSums(par2[, sapply(par2, is.numeric)])
          par2_percent <- par2[, sapply(par2, is.numeric)] / row_sums
          # select category that was asked for
          
          parTabel = cbind(par2_percent[,which(names(par2_percent) == params[["category"]])], 
                           dplyr::select(coords,x,y))
        }
      } 
      # numerical variable
      else {
        # regular extraction with specified fun
        if (!"convert_from_timestep" %in% names(params)) {
          dbl("using buffer to look up data")
          
          parTabel = cbind(terra::extract(x=r, y=bufferedY, fun, na.rm=T, ID= F),
                           dplyr::select(coords,x,y))
        } 
        # an averaging needs to be done towards a coarser time step before extraction
        else {
          for(d in 1:lubridate::days_in_month(atTime)) {
            closest_time=closest(zinfo$times,atTime + lubridate::days(d) - 1)
            sdsn_d=sprintf('%s:%s',sdsn_par,max(closest_time-1,1))
            r_d = rast(sdsn_d)
            if(d == 1) r = r_d
            else r = c(r, r_d)
          }
          parTabel = cbind(rowMeans(terra::extract(x=r, y=bufferedY, fun, na.rm=T, ID= F)),
                           dplyr::select(coords,x,y))
        }
      }
    }
  }, silent=T)
  
  cn=c('','_x','_y')
  
  if(nrow(parTabel) > 0) 
  { if(!is.null(closest_level))
  {  parTabel$DataEvelation= closest_level
  cn=c(cn,'_z')
  }
    if(!is.null(closest_time))
    { 
      #parTabel$DataTime = as.POSIXct(zinfo$times[closest_time], origin="1970-01-01", tz='UTC')
      parTabel$DataTime = zinfo$times[closest_time]
      cn=c(cn,'_t')
    }  
    colnames(parTabel)<-c(sprintf("%s%s",param,cn))
  }
  return(parTabel)
}


#uses rarr::read_array
getTimeSeriesFromZarr<-function(usePar, coords, from=NULL, till=NULL, zinfo){
  
  param = usePar[[1]]
  params = usePar
  
  buoyposlat =  (coords$y[1] - zinfo$latmin) / zinfo$latstep
  buoyposlon =  (coords$x[1] - zinfo$lonmin) / zinfo$lonstep
  endpoint = paste0("https://", str_extract(zinfo$href, "(?<=//)[^\\s/:]+"))
  
  ilist=list( buoyposlat   , buoyposlon)
  if( zinfo$meta$arrays[[param]]$dimensions[2] =='/elevation' ) # read last level, normally sea surface
  {
    levels=length(zinfo$levels)
    ilist=append(levels,ilist)
  }
  
  if(!is.null(zinfo$timeunit))
  {
    
    timeStart=closest(zinfo$times, as.POSIXct(from,origin="1970-01-01",  tz="Z"))
    timeEnd=closest(zinfo$times, as.POSIXct(till,origin="1970-01-01",  tz="Z"))
    times=zinfo$times[timeStart:timeEnd]
    ilist=append(list(timeStart:timeEnd), ilist)
  } 
  
  sc=1
  if(!is.null(zinfo$meta$arrays[[param]]$scale)) sc=zinfo$meta$arrays[[param]]$scale
  offs=0
  if(!is.null(zinfo$meta$arrays[[param]]$offset )) offs=zinfo$meta$arrays[[param]]$offset 
  
  parValues=NULL
  parTable=tibble()
  try({ 
    parValues=read_zarr_array(sprintf("%s/%s", zinfo$href, param ), index= ilist, s3_client = s3_client(endpoint)) *sc   + offs
  })      
  
  if(!is.null(zinfo$timeunit) )
  {
    if(is.null(nrow(parValues))) parValues=rep(NA,length(times))
    parTable= dplyr::tibble("Time" = times, par = parValues)
    colnames(parTable)<-c('Time', param)
  }
  else
  { parTable=dplyr::tibble(par=parValues)
  colnames(parTable)<-c(param)
  }
  
  return(parTable)
}


getLastInfoFromZarr<-function(href,ori=NULL)
{
  #info from zarr file specified in href
  meta=jsonlite::fromJSON(gdal_utils('mdiminfo', toGDAL(href), quiet = T))
  
  zi=list()
  endpoint = paste0("https://", str_extract(href, "(?<=//)[^\\s/:]+"))
  
  #info from stac catalog specified in ori
  if(!is.null(ori) && !is.na(ori)){
    
    zm = retryJson(ori)
    if(length(zm) > 1 ) {
      step3=NULL
      for( a in names(zm$assets)) if(zm$assets[[a]]$href == href  ) step3=zm$assets[[a]]
      
      zi$latmin=step3$viewDims$latitude$coords$min
      zi$latmax=step3$viewDims$latitude$coords$max
      zi$latstep=step3$viewDims$latitude$coords$step
      
      zi$lonmin=step3$viewDims$longitude$coords$min
      zi$lonmax=step3$viewDims$longitude$coords$max
      zi$lonstep=step3$viewDims$longitude$coords$step
      
      zi$timemin=step3$viewDims$time$coords$min
      zi$timemax=step3$viewDims$time$coords$max
      zi$timestep=step3$viewDims$time$coords$step
      zi$timetype=step3$viewDims$time$coords$type
      
      zi$timeunit=step3$viewDims$time$units
      
      zi$levels=step3$timeChunked$viewDims$elevation$len
    }
  }
  
  if(is.null(zi$latstep) | is.null(zi$latmin)) 
  {
    zi$latstep=meta$arrays$latitude$attributes$step
    zi$latmin=meta$arrays$latitude$attributes$valid_min
    zi$latmax=meta$arrays$latitude$attributes$valid_max
    zi$lonstep=meta$arrays$longitude$attributes$step
    zi$lonmin=meta$arrays$longitude$attributes$valid_min
    zi$lonmax=meta$arrays$longitude$attributes$valid_max
  }
  
  if(is.null(zi$latstep) | is.null(zi$latmin)) 
  {
    zi$latmin = read_zarr_array(sprintf("%s%s", href, '/latitude' ), index=list(1), s3_client = s3_client(endpoint))
    zi$latsize=(meta$dimensions %>% dplyr::filter(name=='latitude'))$size
    
    zi$latmax = read_zarr_array(sprintf("%s%s", href, '/latitude' ), index=list(zi$latsize), s3_client = s3_client(endpoint))
    zi$latstep= (zi$latmax-zi$latmin)/zi$latsize
    
    zi$lonmin = read_zarr_array(sprintf("%s%s", href, '/longitude' ), index=list(1), s3_client = s3_client(endpoint))
    zi$lonsize=(meta$dimensions %>% dplyr::filter(name=='longitude'))$size
    
    zi$lonmax = read_zarr_array(sprintf("%s%s", href, '/longitude' ), index=list(zi$lonsize), s3_client = s3_client(endpoint))
    zi$lonstep= (zi$lonmax-zi$lonmin)/zi$lonsize
  }
  if(is.null(zi$timestep))  
    if(!is.null(meta$arrays$time) )
    {
      t1 = read_zarr_array(sprintf("%s%s", href, '/time' ), index=list(1:2), s3_client = s3_client(endpoint))
      zi$timemin = t1[1]
      zi$timestep = t1[2] - t1[1]
      zi$timemax = zi$timemin + (meta$arrays$time$dimension_size -1) * zi$timestep
      zi$timeunit=meta$arrays$time$unit
    }
  
  #is there a time component
  if('time' %in% meta$dimensions$name)
  {
    if(!is.character(meta$arrays$time$unit)) unit='seconds since 1970-01-01'
    else unit=meta$arrays$time$unit 

    origin="1970-01-01"
    if(stringr::str_detect(unit,"2000-01-01")) origin="2000-01-01"
    else if(stringr::str_detect(unit,"1970-01-01")) origin="1970-01-01"
    else if(stringr::str_detect(unit,"1950-01-01")) origin="1950-01-01"
    
    m=1
    if(stringr::str_detect(unit,"milliseconds"))  m=1000
    if(stringr::str_detect(unit,"hours"))         m=1/3600
    if(stringr::str_detect(unit,"minutes"))       m=1/60
    
    zi$times = as.POSIXct(read_zarr_array(sprintf("%s%s", href, '/time'), index=list(1:meta$arrays$time$dimension_size), s3_client = s3_client(endpoint))/m,
                          origin=origin,
                          tz='Z')
  }
  
  if('elevation' %in% meta$dimensions$name)
  { 
    elevation = dplyr::filter(meta$dimensions,name=='elevation')
    zi$levels=read_zarr_array(sprintf("%s%s", href, '/elevation' ), index=list(1:meta$arrays$elevation$dimension_size), s3_client = s3_client(endpoint))
  }
  
  zi$meta=meta
  zi$href=href
  
  return(zi)
}


# main function, adds the requested parameter (usePar) to the data frame (pts) 
# parameters = dslist , usePar, pts
lookupParameter<-function(dslist=NULL, usePar, pts, atDepth=0)
{
  if(is.null(dslist)) return()
  
  param = usePar[[1]]
  params = usePar
  
  registerDoParallel(cores=16)
  mcoptions=list(preschedule=F, silent=F)

  #depending on the time resolution of the dataset, calculate a period to group points by timeslice
  #for some parameters we have climatology data that can be per month and other timesteps.. so check only the first one will not work..
  
  #TODELETE add NA again here once monthly = NA is solved
  if(! dslist$timestep[1] %in% c(0
                                 # , NA
  ) ) {
    
    #whatif the parameter has a time dimension but the pts not?
    if('Time' %in% colnames(pts))
    {  tslist=dplyr::filter(dslist,start_datetime <= min(pts$Time) & end_datetime >= max(pts$Time))  
    if(nrow(tslist)>1) dslist=tslist
    
    dslist=dplyr::arrange(dslist,asset,timestep,latstep)
    
    step=dslist$timestep[1]
    if(is.na(step)) units='months'
    else if(step== 900000 ) units='15 mins'
    else if(step== 1800000 ) units='30 mins'
    else if(step== 3600000 ) units='hours'
    else if(step== 10800000 ) units='3 hours'
    else if(step== 21600000 ) units='6 hours'
    else if(step== 86400000) units='days'
    else if(step== 604800000 ) units='weeks'
    else units='months' # most are defined as specific iso timestep ... extract should be changed
    
    pts=dplyr::arrange(pts,Time)
    pts$Period=lubridate::round_date(pts$Time,units)
    #difference between consecutive times.. needed to group in slices
    
    dbl("Product timestep:", units)
    }  
  }
  #parameter has no timeresolution eg bathymetry
  else
  {# use most recent value.. could be a different choice
    
    dslist=dplyr::filter(dslist, end_datetime == max(end_datetime) )
    pts$Period=round.POSIXt(dslist$end_datetime[1] )
    dslist=dplyr::arrange(dslist,latstep)
  }
  
  #group by timeslice
  periods=c(unique(pts$Period))
  
  dsng=dslist$href[1]
  
  zinfo=getLastInfoFromZarr(dslist$href[1], dslist$ori[1])

  if("Time" %in% colnames(pts))
  {
    pts$diff= as.numeric(difftime(pts$Time,pts$Time[1], 'secs'))
    pts$Slice=cumsum(c(T, diff(pts$diff) > 24*60*60)) # threshold of days
    locs=dplyr::group_by(pts,Longitude,Latitude,Slice) %>% dplyr::summarise(from=min(Period),till=max(Period), cnt=n(), .groups='drop')
  }
  else
    locs=dplyr::group_by(pts,Longitude,Latitude) %>% dplyr::summarise(from=min(Period),till=max(Period), cnt=n(), .groups='drop')
  
  
  lookup=c("Time","Longitude","Latitude","loop","method")
  names(lookup)=c(sprintf("vb_%s_t", param),sprintf("vb_%s_x", param),sprintf("vb_%s_y",param),sprintf("vb_%s_l",param),sprintf("vb_%s_m",param))
  
  if(nrow(locs) < length(periods)) {
    #loop locations
    resulting= foreach(p=1:nrow(locs),.options.multicore=mcoptions, .combine='CombineR', .packages = c("terra","stars","magrittr","dplyr","Rarr")) %dopar% 
      { 
        tble=getTimeSeriesFromZarr(param, tibble('x'= locs$Longitude[p],'y'=locs$Latitude[p]), from=locs$from[p], till=locs$till[p], zinfo=zinfo )
        
        thistble=dplyr::filter(pts, Latitude==locs$Latitude[p],Longitude==locs$Longitude[p] ) %>% dplyr::select(Time=Period,Latitude,Longitude)
        if("Time" %in% colnames(tble)) tble=dplyr::left_join(thistble, tble, by="Time" ) 
        else if(nrow(tble)==1) {v=tble[[param]][1]; thistble[[param]]=v; tble=thistble }
        else tble=bind_cols(thistble,tble)
        tble$loop=p
        tble$method="ts"
        tble=dplyr::rename(tble, any_of(lookup))
        
        tble
      }
    dbl("zarr extracted tble: ", nrow(resulting))
  }   
  else {
    
    #loop time slices   
    resulting= foreach(p=1:length(periods),.options.multicore=mcoptions, .combine='CombineR', .packages = c("terra","stars","magrittr","dplyr","Rarr")) %dopar% 
      { 
        tble=NULL
        
        thistble=dplyr::filter(pts,Period==periods[p])
        
        if(nrow(thistble)>0) 
        {
          if("Elevation" %in% colnames(pts)) atDepth = mean(pts$Elevation, na.rm=T) else
            if("Bathymetry" %in% colnames(pts)) atDepth= mean(pts$Bathymetry * -1, na.rm=T)
            if(is.na(atDepth)) atDepth=0
            if(atDepth >0 ) atDepth=0-atDepth
            tble=getParFromZarrwInfo(params, coords = dplyr::select(thistble, x=Longitude,y=Latitude) , atTime = periods[p] , atDepth = atDepth, zinfo = zinfo, isCategory =  (!dslist$categories[1] %in% c(NA,0)) )
            # veryvery slow, moved to sandbox
            # tble=getMultipleTimeSeriesFromZarr(params,coords = dplyr::select(thistble, x=Longitude,y=Latitude), from=periods[p], till=periods[p], zinfo=zinfo )
        }  
        
        #this prevents cbind errors when row counts don't match
        if(nrow(tble) != nrow(thistble) ) { tble = tibble( "par" = rep(NA , nrow(thistble)  ) ) ; colnames(tble)<-c(param) }
        rm(thistble)
        tble$loop=p
        tble$method="rast"
        tble=dplyr::rename(tble, any_of(lookup))
        tble
      }
  }  
  resulting = cbind(pts,  resulting  )
  
  newpiece=list(usePar=list("href"= dsng, "par"=param,"nr"= nrow(resulting),"points"=nrow(locs),"periods"=length(periods),"stacinfo"= dplyr::slice(dslist,1),"zarrmeta"=zinfo$meta ,"progress"=loglist ) )
  names(newpiece)=c(param)
  
  #for category data, add the description found in the zarr meta data
  if(!dslist$categories[1] %in% c(NA,0))
  {
    catdesc=zinfo$meta$attributes[[param]]
    resulting[[paste0(param,"_Description")]]=ifelse(resulting[[param]] %in% catdesc,  names(catdesc)[match(resulting[[param]], catdesc)] ,'NA')
  }
  
  resulting <- resulting %>% dplyr::select( -any_of(c("Period","Slice","diff") ))
  
  #done with this parameter, save memory
  rm(zinfo)
  list('log'=newpiece , 'df'=resulting)
}


enhanceDF<-function(inputPoints, requestedParameters, requestedTimeSteps, stacCatalogue, verbose="", atDepth = NA, select_layers = NULL)
{
  #we need to group the sampling points to reduce the data lookups
  #rounding to 3 deg decimals for lat/lon ot minutes for time is a very crude approximation
  #better is to use a snap to grid approach based on the spatiotemporal resolution of the parameter  
  # so this needs to move to the enhanceDF function
  # select timechunked or geochunked zarr files for the cmems / copernicus marine data store Zarr files to optimise lookup performance
  # emodnet zarr files are just chunked
  
  if("Time" %in% colnames(inputPoints)) {
    inputPoints$Time = as.POSIXct(inputPoints$Time, tz="UTC")
    inputPoints=dplyr::mutate(inputPoints,RTime=round.POSIXt(Time,units="mins"), RLatitude= round(Latitude, 3), RLongitude= round(Longitude,3) )
    optimalChunking=c('timeChunked','geoChunked','chunked')
  } else {
    inputPoints=dplyr::mutate(inputPoints,RLatitude= round(Latitude, 3), RLongitude= round(Longitude,3) )
    optimalChunking=c('timeChunked','chunked')
  }
  
  pts=dplyr::select(inputPoints,Latitude=RLatitude,Longitude=RLongitude,any_of(c("Time"="RTime")),any_of(c("Bathymetry","Elevation"))) %>% dplyr::distinct()
  
  extracted=c()
  
  sources=c()
  
  for(currentPar in 1:length(requestedParameters))
  { 
    parameter=requestedParameters[[currentPar]]
    
    param = ifelse(length(parameter)==1, parameter, parameter[[1]])
    
    dbl("Deriving",param)

    #check available zarr assets in the stack catalogue for the region and period in the data file
    #order by resolution and take the first record
    
    dslist=dplyr::filter(stacCatalogue, par==param & 
                           latmin <= min(pts$Latitude) & 
                           latmax >= max(pts$Latitude) & 
                           lonmin <= min(pts$Longitude) & 
                           lonmax >= max(pts$Longitude) &
                           chunktype %in% c('timeChunked','geoChunked','chunked'))
    #TODELETE
    if(param == "elevation") dslist=dplyr::filter(stacCatalogue, par==param &
                                                    chunktype %in% c('timeChunked','geoChunked','chunked'))
    
    #if the variable is categorical, do not filter based on start and end time (should be changed to dynamic vs static variable in future)
    if(dslist$categories[1] %in% c(NA,0))    dslist=dplyr::filter(dslist, start_datetime <= min(pts$Time) & end_datetime >= max(pts$Time)) 
    
    if(! "Time" %in% colnames(pts))
      dslist=dplyr::filter(dslist,timestep %in% c(NA,0))
    
    if(nrow(dslist)==0) {
      if(param %in% stacCatalogue$par)
        dbl("no data assets fit the criteria")
      else
        dbl("requested parameter ", param, " not found")
      next
    }
    
    # if a preferred timestep is specified, use it
    if(! is.null(requestedTimeSteps))
    {
      tlist=dplyr::filter(dslist, timestep %in% requestedTimeSteps)
      if("convert_from_timestep" %in% names(parameter)) {
        dbl(paste0("Converting original product with timestep ", as.numeric(parameter["convert_from_timestep"]), " towards a timestep of ", requestedTimeSteps))
        requestedTimeSteps2 = as.numeric(parameter["convert_from_timestep"])
        tlist=dplyr::filter(dslist, timestep %in% requestedTimeSteps2)
      }
      if(nrow(tlist)==0) {
        dbl("no data assets fit the criteria")
        return()
      }
      if(nrow(tlist)>0) dslist=tlist
    }
    
    # geochunked, timechunked or just chunked, to be optimized in future versions
    if(! is.null(optimalChunking))
    {
      tlist=dplyr::filter(dslist, asset %in% optimalChunking)
      if(nrow(tlist)>0) dslist=tlist
    }
    
    #add time filter here (see lookupParameter function)
    
    if(is.null(select_layers) & nrow(dslist) > 1) {
      dbl("Multiple options are available for your input parameters:")
      
      print(dslist %>% dplyr::select(title, latmin, latmax, lonmin, lonmax, latstep, lonstep, timestep, start_datetime, end_datetime, chunktype))
      
      ind <- readline(prompt = paste0("choose a product (number in the range 1 - ", nrow(dslist),"): "))
      dslist <- dslist[ind,]
    } else if(!is.null(select_layers)) dslist <- dslist[select_layers[currentPar],]
    
    # this object will have a dataframe $df and a logfile $log, the dataframe has several verbose columns pointing to the nearest lat,lon,time and depth found in the requested data set    
    extracted=CombineC(extracted, lookupParameter(dslist=dslist,usePar=parameter,pts=pts,atDepth=atDepth) )
    
    sources = c(sources,dslist$title)
  }
  
  pts=extracted$df
  
  #add everything to the user data frame
  #join the results, depending on available columns
  cl=c("RLatitude"="Latitude","RLongitude"="Longitude")
  if("Time" %in% colnames(pts)) cl=c(cl,"RTime"="Time")
  if("Bathymetry" %in% colnames(pts)) cl=c(cl,"Bathymetry"="Bathymetry")
  
  #remove verbose columns , unless specified
  if(verbose == "") {
    pts = pts %>% dplyr::select( -ends_with(c('_x','_y','_z','_t'  )))
    pts = pts %>% dplyr::select( -starts_with('vb_'))
  }
  #add to the input dataframe, remove the processing columns
  inputPoints=dplyr::left_join(inputPoints,pts,by =cl) %>% dplyr::select(-any_of(c("RLatitude","RLongitude","RTime"))) 
  
  #clean up
  rm(pts)
  
  message(paste0("Extraction done. The following datasets were used:\n - ", paste(sources, collapse = "\n - ")))
  return(inputPoints)
  
}

getRasterSlice <- function(requestedParameter='thetao', lon_min=-10, lon_max=10, lat_min=50, lat_max=60, 
                           requestedTimeSteps = NULL, date = NULL, atDepth = NULL, stacCatalogue, 
                           select_layers = NULL) {
  
  if(is.null(atDepth)) atDepth=0
  if(! is.null(date))
  { date <- as.POSIXct(date, format = "%Y-%m-%d")
  if(is.na(date)) {
    dbl("failed to recognise date format, please use format %Y-%m-%d")
    date=NULL
  }}
  
  if(! is.null(date))  
    dslist=dplyr::filter(stacCatalogue, par==requestedParameter & 
                           latmin <= lat_min & 
                           latmax >= lat_max & 
                           lonmin <= lon_min & 
                           lonmax >= lon_max &
                           start_datetime <= date & end_datetime >= date &
                           chunktype %in% c('timeChunked','chunked'))
  
  else  
    dslist=dplyr::filter(stacCatalogue, par==requestedParameter & 
                           latmin <= lat_min & 
                           latmax >= lat_max & 
                           lonmin <= lon_min & 
                           lonmax >= lon_max &
                           chunktype %in% c('timeChunked','chunked'))
  
  
  # if a prefered timestep specified, use it
  if(! is.null(requestedTimeSteps))
  {
    tlist=dplyr::filter(dslist, timestep %in% requestedTimeSteps)
    if(nrow(tlist)>0) dslist=tlist
  }
  
  if(is.null(select_layers) & nrow(dslist) > 1) {
    dbl("Multiple options are available for your input parameters:")
    
    print(dslist %>% dplyr::select(title, latmin, latmax, lonmin, lonmax, latstep, lonstep, timestep, start_datetime, end_datetime))
    
    ind <- readline(prompt = paste0("choose a product (number in the range 1 - ", nrow(dslist),"): "))
    dslist <- dslist[ind,]
  }
  else if(!is.null(select_layers))  dslist <- dslist[select_layers,]
  
  zinfo=getLastInfoFromZarr(dslist$href[1], dslist$ori[1])
  isCategory =  (!dslist$categories[1] %in% c(NA,0)) 
  
  dsn=toGDAL(zinfo$href)
  
  #check if parameter is available in zarr file
  if(requestedParameter %in% names(zinfo$meta$arrays)) 
    sdsn=sprintf('%s:/%s',dsn,requestedParameter)
  
  closest_time=NULL
  #is there a time component
  if('time' %in% zinfo$meta$dimensions$name & !is.null(zinfo$times))
  {
    if(length(zinfo$times) > 1)
    { 
      if(is.null(date)) closest_time=1 #take first available data
      else closest_time=closest(zinfo$times,date)
      
      sdsn=sprintf('%s:%s',sdsn,closest_time-1)
    }
  }
  
  #is the data in several depth or elevation levels
  closest_level=NULL
  if('elevation' %in% zinfo$meta$dimensions$name & !is.null(zinfo$levels) )
  { 
    if(length(zinfo$levels) > 1)
    {
      closest_level=closest(zinfo$levels,atDepth)
      sdsn=sprintf('%s:%s',sdsn,closest_level-1)
      closest_level = zinfo$levels[closest_level]
    }
  }
  
  r=rast(sdsn)

  if(isCategory)  {#r=terra::flip(r, "vertical") --> currently does not work
    #manually flip raster
    r <- read_stars(sdsn,driver = "zarr")
    values_matrix <- as.matrix(r[[1]])
    flipped_matrix <- values_matrix[, ncol(values_matrix):1]
    raster_flipped <- st_as_stars(flipped_matrix)
    st_dimensions(raster_flipped) <- st_dimensions(r)
    r <- st_as_stars(raster_flipped)
    r <- as(r, "SpatRaster")
  }
  
  dbl("extent and coordinate system missing, assuming epsg:4326, extent from stac catalogue")
  crs(r)='epsg:4326'
  ext(r)=c(zinfo$lonmin[1],zinfo$lonmax[1],
           min(zinfo$latmin[1],zinfo$latmax[1]),
           max(zinfo$latmin[1],zinfo$latmax[1]))
  
  r <- crop(r, ext(lon_min, lon_max, lat_min, lat_max))
  
  return(r)
}

