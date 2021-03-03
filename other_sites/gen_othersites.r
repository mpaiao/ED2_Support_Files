#==========================================================================================#
#==========================================================================================#
#     Reset the session.                                                                   #
#------------------------------------------------------------------------------------------#
rm(list=ls())
graphics.off()
options(warn=0)
#==========================================================================================#
#==========================================================================================#




#----- Define some paths. -----------------------------------------------------------------#
main    = dirname(getwd())           # Main directory
here    = getwd()                    # Current path 
srcdir  = file.path(main,"Rsc")      # Path with additional R functions
outroot = file.path(here,"figures")  # Output directory
#------------------------------------------------------------------------------------------#




#------------------------------------------------------------------------------------------#
#     Session information, for saving/loading data already read/analysed.                  #
#------------------------------------------------------------------------------------------#
reload.clean     = TRUE                             # Reload the cleaned data? 
save.clean       = TRUE                             # Save cleaned data if read?
clean.rdata      = "raw_othersites-2008-2016.RData" # R object with the cleaned data
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#     Define the site-specific information for reading the dataset.                        #
#------------------------------------------------------------------------------------------#
yeara              = 2008                         # First year to use
yearz              = 2016                         # Last year to use
lon.ref            = -52.376                      # Tanguro Longitude (degrees east)
lat.ref            = -13.082                      # Tanguro Latitude (degrees north)
dtdat              = 3600.                        # Minimum time step between observation
imetavg            = 1                            # How do the averages relate to time?
                                                  # 1. Ending at the time stamp
                                                  # 2. Beginning at time stamp
                                                  # 3. Centred around time stamp
discard.suspicious = TRUE                         # Remove data that are suspicious based
                                                  #   on visual inspection.
#------------------------------------------------------------------------------------------#




#------------------------------------------------------------------------------------------#
#     Plot options.                                                                        #
#------------------------------------------------------------------------------------------#
outform        = "pdf"          # Formats for output file (case insensitive):
                                #   - "X11" - for printing on screen
                                #   - "eps" - for postscript printing
                                #   - "png" - for PNG printing
byeold         = TRUE           # Remove old files of the given format?
depth          = 96             # PNG resolution, in pixels per inch
paper          = "square"       # Paper size, to define the plot shape
ptsz           = 20             # Font size.
plotgrid       = TRUE           # Should I plot the grid in the background?
legwhere       = "topleft"      # Where should I place the legend?
inset          = 0.01           # inset distance between legend and edge of plot region.
legbg          = "white"        # Legend background colour.
plotsd         = TRUE          # Plot the standard deviation (T/F)
angle          = 45             # Angle for shaded 1-SD
dens           = 40             # Density of lines
llwd           = 2.5            # Line width
shwd           = 1.0            # Width of shaded lines
fracexp        = 0.2            # Expansion factor
f.leg          = 1/6            # Fraction of the plotting area for legend.
#------------------------------------------------------------------------------------------#




#------------------------------------------------------------------------------------------#
#     List of sites.  Variables within the list elements should be:                        #
#  - iata  -- 3-letter code identifying the site.                                          #
#  - short -- longer name handle to identify site                                          #
#  - desc  -- nice-looking name to appear in plot titles                                   #
#  - csv   -- csv file with data for this site                                             #
#  - lon   -- site longitude                                                               #
#  - lat   -- site altitude                                                                #
#  - alt   -- site altitude (above sea level)                                              #
#  - off.utc -- offset for UTC.  If data are in local time, the time stamp corresponds to  #
#               the end of the averaging period, and the time zone is UTC-4, set           #
#               off.utc=-4.  If the data are already in UTC, set off.utc=0.  If the time   #
#               is UTC but the time stamp is the beginning of the hourly average period,   #
#               set off.utc=-1                                                             #
#  - dtime   -- average period, in seconds (e.g., if data are hourly averages, set it to   #
#               3600.                                                                      #
#  - dtype   -- data type.  Currently only "inmet" is provided, but feel free to define    #
#               others, and look for inmet occurrences throughout this script for the      #
#               appropriate changes.                                                       #
#------------------------------------------------------------------------------------------#
n          = 0
sites      = list()
n          = n + 1
sites[[n]] = list( iata    = "que"
                 , short   = "querencia"
                 , desc    = "Querencia, MT"
                 , csv     = "inmet_a916_merged_2008-2016.csv"
                 , lon     =  -52.221
                 , lat     = -12.627
                 , alt     = 361.
                 , off.utc = 0
                 , dtime   = 3600.
                 , dtype   = "inmet"
                 )#end list
n          = n + 1
sites[[n]] = list( iata    = "aub"
                 , short   = "agua_boa"
                 , desc    = "Agua Boa, MT"
                 , csv     = "inmet_a908_merged_2007-2016.csv"
                 , lon     =  -52.212
                 , lat     = -14.016
                 , alt     = 440.
                 , off.utc = 0
                 , dtime   = 3600.
                 , dtype   = "inmet"
                 )#end list
#------------------------------------------------------------------------------------------#


#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#
#            CHANGES BEYOND THIS POINT ARE FOR ADJUSTING THE INPUT FILE ONLY.              #
#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#


#------------------------------------------------------------------------------------------#
#     Load some useful packages and functions.                                             #
#------------------------------------------------------------------------------------------#
srcdir = (srcdir[file.exists(srcdir)])[1]
source(file.path(srcdir,"load.everything.r"))
#------------------------------------------------------------------------------------------#



#----- Create output directory in case it doesn't exist. ----------------------------------#
if (! file.exists(outroot)) dir.create(outroot)
#------------------------------------------------------------------------------------------#



#----- Define plot window size ------------------------------------------------------------#
f.ext   = f.leg / (1. - f.leg)
sq.size = plotsize(proje=FALSE,paper=paper)
ex.size = plotsize(proje=FALSE,paper=paper,extendfc="lat",extfactor=f.ext)
#------------------------------------------------------------------------------------------#



#----- Determine how many formats we must output. -----------------------------------------#
outform = tolower(outform)
nout    = length(outform)
#------------------------------------------------------------------------------------------#


#------------------------------------------------------------------------------------------#
#     List of variables from INMET.  This allows variable transformation.                  #
#------------------------------------------------------------------------------------------#
n          = 0
inmet      = list()
n          = n + 1
inmet[[n]] = list( vname = "year"
                 , type  = "integer"
                 , mult  = 1.
                 , add   = 0.
                 )#end list
n          = n + 1
inmet[[n]] = list( vname = "month"
                 , type  = "integer"
                 , mult  = 1.
                 , add   = 0.
                 )#end list
n          = n + 1
inmet[[n]] = list( vname = "day"
                 , type  = "integer"
                 , mult  = 1.
                 , add   = 0.
                 )#end list
n          = n + 1
inmet[[n]] = list( vname = "hour"
                 , type  = "numeric"
                 , mult  =  1./day.hr
                 , add   = 0.
                 )#end list
n          = n + 1
inmet[[n]] = list( vname = "atm.tmp"
                 , type  = "numeric"
                 , mult  = 1.
                 , add   = t00
                 )#end list
n          = n + 1
inmet[[n]] = list( vname = "atm.rhv"
                 , type  = "numeric"
                 , mult  = 0.01
                 , add   = 0.
                 )#end list
n          = n + 1
inmet[[n]] = list( vname = "atm.tdew"
                 , type  = "numeric"
                 , mult  = 1.
                 , add   = t00
                 )#end list
n          = n + 1
inmet[[n]] = list( vname = "atm.prss"
                 , type  = "numeric"
                 , mult  = 100.
                 , add   = 0.
                 )#end list
n          = n + 1
inmet[[n]] = list( vname = "atm.vdir"
                 , type  = "numeric"
                 , mult  = 1.
                 , add   = 0.
                 )#end list
n          = n + 1
inmet[[n]] = list( vname = "atm.vels"
                 , type  = "numeric"
                 , mult  = 1.
                 , add   = 0.
                 )#end list
n          = n + 1
inmet[[n]] = list( vname = "rain"
                 , type  = "numeric"
                 , mult  = 1./hr.sec
                 , add   = 0.
                 )#end list
n          = n + 1
inmet[[n]] = list( vname = "rshort.in"
                 , type  = "numeric"
                 , mult  = 1000. / hr.sec
                 , add   = 0.
                 )#end list
#------------------------------------------------------------------------------------------#


#------------------------------------------------------------------------------------------#
#      List of variables to be plotted.                                                    #
#------------------------------------------------------------------------------------------#
n         = 0
tvar      = list()
n         = n + 1
tvar[[n]] = list( vnam      = "rshort.in"
                , desc      = "Incident SW"
                , unit      = "wom2"
                , add       = 0
                , mult      = 1.
                , colmean   = "darkorange3"
                , colerr    = "gold"
                , coledge   = "gray31"
                , plt       = TRUE
                , edge      = TRUE
                , out.hour  = TRUE
                , out.all   = TRUE
                )#end list
n         = n + 1
tvar[[n]] = list( vnam      = "rlong.in"
                , desc      = "Incident LW"
                , unit      = "wom2"
                , add       = 0
                , mult      = 1.
                , colmean   = "darkorange3"
                , colerr    = "gold"
                , coledge   = "gray31"
                , plt       = TRUE
                , edge      = TRUE
                , out.hour  = TRUE
                , out.all   = TRUE
                )#end list
n         = n + 1
tvar[[n]] = list( vnam      = "atm.tmp"
                , desc      = "Air temperature"
                , unit      = "degC"
                , add       = -t00
                , mult      = 1.
                , colmean   = "darkorange3"
                , colerr    = "gold"
                , coledge   = "gray31"
                , plt       = TRUE
                , edge      = TRUE
                , out.hour  = TRUE
                , out.all   = TRUE
                )#end list
n         = n + 1
tvar[[n]] = list( vnam      = "atm.rhv"
                , desc      = "Relative humidity"
                , unit      = "pc"
                , add       = 0.
                , mult      = 100.
                , colmean   = "midnightblue"
                , colerr    = "steelblue"
                , coledge   = "gray31"
                , plt       = TRUE
                , edge      = TRUE
                , out.hour  = TRUE
                , out.all   = TRUE
                )#end list
n         = n + 1
tvar[[n]] = list( vnam      = "atm.pvap"
                , desc      = "Vapour pressure"
                , unit      = "hpa"
                , add       = 0.
                , mult      = 0.01
                , colmean   = "midnightblue"
                , colerr    = "steelblue"
                , coledge   = "gray31"
                , plt       = TRUE
                , edge      = TRUE
                , out.hour  = TRUE
                , out.all   = TRUE
                )#end list
n         = n + 1
tvar[[n]] = list( vnam      = "rain"
                , desc      = "Precipitation rate"
                , unit      = "kgwom2oday"
                , add       = 0.
                , mult      = day.sec
                , colmean   = "midnightblue"
                , colerr    = "steelblue"
                , coledge   = "gray31"
                , plt       = TRUE
                , edge      = FALSE
                , out.hour  = FALSE
                , out.all   = FALSE
                )#end list
n         = n + 1
tvar[[n]] = list( vnam      = "atm.vels"
                , desc      = "Wind speed"
                , unit      = "mos"
                , add       = 0.
                , mult      = 1.
                , colmean   = "purple4"
                , colerr    = "orchid"
                , coledge   = "gray31"
                , plt       = TRUE
                , edge      = TRUE
                , out.hour  = TRUE
                , out.all   = TRUE
                )#end list
n         = n + 1
tvar[[n]] = list( vnam      = "atm.prss"
                , desc      = "Atmospheric Pressure"
                , unit      = "hpa"
                , add       = 0.
                , mult      = 0.01
                , colmean   = "purple4"
                , colerr    = "orchid"
                , coledge   = "gray31"
                , plt       = TRUE
                , edge      = TRUE
                , out.hour  = TRUE
                , out.all   = TRUE
                )#end list
#------------------------------------------------------------------------------------------#


#----- Transform lists into data frames. --------------------------------------------------#
sites   = list.2.data.frame(sites)
inmet   = list.2.data.frame(inmet)
tvar    = list.2.data.frame(tvar )
#------------------------------------------------------------------------------------------#


#----- Determine how many places we must compile and how many variables we must plot. -----#
nsites = nrow(sites)
nepro  = nrow(epro )
ninmet = nrow(inmet)
ntvar  = nrow(tvar )
#------------------------------------------------------------------------------------------#


#------------------------------------------------------------------------------------------#
#    Limits for the combined dataset.                                                      #
#------------------------------------------------------------------------------------------#
whena = chron(dates=paste( 1, 1,yeara,sep="/"),times=paste( 0, 0, 0,sep=":"))
whenz = chron(dates=paste(12,31,yearz,sep="/"),times=paste(23, 0, 0,sep=":"))
#------------------------------------------------------------------------------------------#


#----- List of variables that MUST exist in the data frame (even if it's all NA). ---------#
must.exist = c("atm.tmp","atm.prss","atm.pvap","rshort.in","rlong.in","atm.vels"
              ,"atm.vdir","rain")
#------------------------------------------------------------------------------------------#


#----- Define the R Data file name with cleaned data. -------------------------------------#
session.clean    = file.path(here,clean.rdata)
#------------------------------------------------------------------------------------------#


#------------------------------------------------------------------------------------------#
#      Function that helps visualising data in a short window.                             #
#------------------------------------------------------------------------------------------#
inspect <<- function( aws
                    , vnam  = stop("Provide a variable name")
                    , whena = min(aws$when)
                    , whenz = max(aws$when)
                    , ...
                    ){

   #----- Find time bounds. ---------------------------------------------------------------#
   whena   = chron(whena)
   whenz   = chron(whenz)
   sel     = aws$when %wr% c(whena,whenz)
   #---------------------------------------------------------------------------------------#


   #----- Find neat x and y limits. -------------------------------------------------------#
   xlim    = c(whena,whenz)
   xstuff  = pretty.time(when=c(whena,whenz),n=8)
   ylim    = range(aws[[vnam]][sel],finite=TRUE)
   yat     = pretty(ylim)
   ylabels = sprintf("%g",yat)
   #---------------------------------------------------------------------------------------#


   #----- Open device. --------------------------------------------------------------------#
#   pdf("inspect.pdf",width=9,height=9,pointsize=20,onefile=FALSE,family="Helvetica")
   par(par.user)
   par(mar=c(4.6,4.6,1.1,1.1))
   plot.new()
   plot.window(xlim=xlim,ylim=ylim)
   axis(side=1,las=1,at=xstuff$levels,labels=xstuff$labels)
   axis(side=2,las=1,at=yat,labels=ylabels)
   abline(h=yat,v=xstuff$levels,col="grey80",lty="solid")
   lines(x=aws$when[sel],y=aws[[vnam]][sel],...)
   title(ylab=vnam,xlab="Time [UTC]")
   box()
#   graphics.off()
   #---------------------------------------------------------------------------------------#

   invisible()
}#end function inspect
#------------------------------------------------------------------------------------------#


#------------------------------------------------------------------------------------------#
#    Check whether we should read the data or reload from the background.                  #
#------------------------------------------------------------------------------------------#
if (reload.clean && file.exists(session.clean)){
   cat0(" + Load cleaned data from file ",basename(session.clean),".")
   load(session.clean)
}else{
   #---------------------------------------------------------------------------------------#
   #     Make the complete list of David Fitzjarrald's AWS.                                #
   #---------------------------------------------------------------------------------------#
   fitz.aws = list()
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #    Determine the bounds for the times we are interested in.                           #
   #---------------------------------------------------------------------------------------#
   whena = chron(dates=paste( 1, 1,yeara,sep="/"),times=paste( 0, 0, 0,sep=":"))
   whenz = chron(dates=paste(12,31,yearz,sep="/"),times=paste(23, 0, 0,sep=":"))
   zero  = chron(paste(12,31,yeara-1,sep="/"))
   #---------------------------------------------------------------------------------------#




   #----- Define the time vector that will be used by the full time series. ---------------#
   cat0(" + Define the time vector.")
   when = NULL
   for (yr in seq(from=yeara,to=yearz,by=1)){
      for (mon in sequence(12)){
         dmax     = daymax(mon,yr)
         when.now = chron(dates = rep(paste(  mon, 1:dmax, yr,sep="/") ,each  = 24  )
                         ,times = rep(paste( 0:23,      0,  0,sep=":") ,times = dmax) )
         when     = c(when,when.now)
      }#end for
   }#end for
   when     = chron(when)
   dt.scale = hr.sec / min(sites$dtime)
   elapsed  = round(day.hr * dt.scale * (when - zero)) / dt.scale
   #---------------------------------------------------------------------------------------#


   #----- Find auxiliary variables that will help to make the combined dataset. -----------#
   ntimes = length(when)
   empty  = rep(NA,times=ntimes)
   #---------------------------------------------------------------------------------------#

   #---------------------------------------------------------------------------------------#
   #    Make a template data frame.                                                        #
   #---------------------------------------------------------------------------------------#
   template            = list()
   template$when       = when
   #----- Initialise the other variables with NA. -----------------------------------------#
   template$atm.tmp    = empty
   template$atm.prss   = empty
   template$atm.pvap   = empty
   template$rshort.in  = empty
   template$rlong.in   = empty
   template$atm.vels   = empty
   template$atm.vdir   = empty
   template$rain       = empty
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #     Loop over all sites.                                                              #
   #---------------------------------------------------------------------------------------#
   for (s in sequence(nsites)){
      #----- Aliases for some useful variables. -------------------------------------------#
      s.iata    = sites$iata   [s]
      s.short   = sites$short  [s]
      s.desc    = sites$desc   [s]
      s.csv     = sites$csv    [s]
      s.lon     = sites$lon    [s]
      s.lat     = sites$lat    [s]
      s.alt     = sites$alt    [s]
      s.off.utc = sites$off.utc[s]
      s.dtime   = sites$dtime  [s]
      s.dtype   = sites$dtype  [s]
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #    "aws" is the current dataset we are working, which will eventually become       #
      # part of fitz.aws.                                                                  #
      #------------------------------------------------------------------------------------#
      aws              = list()
      aws$short        = s.short
      aws$longname     = s.desc
      aws$iata         = s.iata
      aws$lon          = s.lon
      aws$lat          = s.lat
      aws$alt          = s.alt
      aws$dist.ref     = 1000. * rdist.earth( x1    = rbind( c(aws$lon,aws$lat)
                                                           , c(lon.ref,lat.ref)
                                                           )#end rbind
                                            , miles = FALSE
                                            )[1,2]
      cat0(" + Read data from ",aws$longname,".")
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #     Create three structures:                                                       #
      # - aws is going to be the eventual data, but initially it will contain the hourly   #
      #   average ending at the time stamp;                                                #
      # - aws.p30 will hold the data with the half-hour averages centred at the half       #
      #   past the hour                                                                    #
      # - aws.60 will hold the data with the 30-minute averages centred at the full        #
      #   hour before the hour of interest.                                                #
      #------------------------------------------------------------------------------------#
      cat0("   - Initialise data frames.")
      #----- Append the template to the list. ---------------------------------------------#
      aws = c(aws,template)
      aws = alltimes( datin     = aws
                    , lon       = aws$lon
                    , lat       = aws$lat
                    , ed21      = TRUE
                    , zeronight = FALSE
                    , meanval   = TRUE
                    , imetavg   = 1
                    , nmean     = 120
                    , na.rm     = TRUE
                    )#end alltimes
      #------------------------------------------------------------------------------------#



      #----- Read the raw file. -----------------------------------------------------------#
      rawfile = file.path(here,s.csv)
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #       Define how to read the data set based on the file type.                      #
      #------------------------------------------------------------------------------------#
      if (s.dtype %in% "epro"){
         #----- Read in the data set. -----------------------------------------------------#
         cat0("   - Read EddyPRO data from ",s.csv,".")
         nlines = as.integer( system( command = paste0("cat ",rawfile," | wc -l")
                                    , intern  = TRUE
                                    )#end system
                            )#end as.integer
         raw    = read.csv( file             = rawfile
                          , header           = TRUE
                          , sep              = ","
                          , comment.char     = ""
                          , nrows            = nlines
                          , colClasses       = c("character",rep("numeric",194),"character"
                                                ,rep("numeric",7),"character","numeric")
                          , stringsAsFactors = FALSE
                          )#end read.table
         #---------------------------------------------------------------------------------#


         #----- Throw away columns that we won't use. -------------------------------------#
         raw  = raw[,epro$vorig,drop=FALSE]
         names(raw) = epro$vname
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #      Standardise scale of variables that must be scaled.                        #
         #---------------------------------------------------------------------------------#
         cat0("   - Standardise variable scales:")
         for (e in which(epro$type %in% "numeric")){
            e.vnam        = epro$vname[e]
            e.mult        = epro$mult [e]
            e.add         = epro$add  [e]
            raw[[e.vnam]] = e.add + e.mult * raw[[e.vnam]]
            cat0("     * ",e.vnam,".")
            
         }#end for (e in which(epro$type %in% "numeric"))
         #---------------------------------------------------------------------------------#
         
         #---------------------------------------------------------------------------------#
         #       Create time variable.                                                     #
         #---------------------------------------------------------------------------------#
         raw$when    = chron( chron(paste(12,31,raw$year-1,sep="/"))
                            + raw$doy + raw$hfrac - s.off.utc / day.hr
                            )#end chron
         raw$elapsed = round(day.hr * dt.scale * (raw$when - zero)) / dt.scale
         #---------------------------------------------------------------------------------#
      }else{
         #----- Read in the data set. -----------------------------------------------------#
         cat0("   - Read INMET data from ",s.csv,".")
         nlines = as.integer( system( command = paste0("cat ",rawfile," | wc -l")
                                    , intern  = TRUE
                                    )#end system
                            )#end as.integer
         raw    = read.csv( file             = rawfile
                          , header           = TRUE
                          , sep              = ","
                          , comment.char     = ""
                          , nrows            = nlines
                          , colClasses       = rep("numeric",12)
                          , stringsAsFactors = FALSE
                          )#end read.table
         #---------------------------------------------------------------------------------#


         #----- Throw away columns that we won't use. -------------------------------------#
         raw  = raw[,inmet$vname,drop=FALSE]
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #      Standardise scale of variables that must be scaled.                        #
         #---------------------------------------------------------------------------------#
         cat0("   - Standardise variable scales:")
         for (i in which(inmet$type %in% "numeric")){
            i.vnam        = inmet$vname[i]
            i.mult        = inmet$mult [i]
            i.add         = inmet$add  [i]
            cat0("     * ",i.vnam,".")
            raw[[i.vnam]] = i.add + i.mult * raw[[i.vnam]]
         }#end for (i in which(inmet$type %in% "numeric"))
         #---------------------------------------------------------------------------------#
         
         #---------------------------------------------------------------------------------#
         #       Create time variable.                                                     #
         #---------------------------------------------------------------------------------#
         raw$when    = chron( chron(paste(raw$month,raw$day,raw$year,sep="/"))
                            + raw$hour - s.off.utc / day.hr
                            )#end chron
         raw$elapsed = round(day.hr * dt.scale * (raw$when - zero)) / dt.scale
         #---------------------------------------------------------------------------------#
      }#end if (s.dtype)
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #      Fill in data with all relevant time information.                              #
      #------------------------------------------------------------------------------------#
      cat0("   - Obtain additional information on time.")
      raw = alltimes( datin     = raw
                    , lon       = aws$lon
                    , lat       = aws$lat
                    , ed21      = TRUE
                    , zeronight = FALSE
                    , meanval   = TRUE
                    , imetavg   = 1
                    , nmean     = 120
                    , na.rm     = TRUE
                    )#end alltimes
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #     Check that all variables that must exist are assigned in 'raw'.  In case any   #
      # variable is missing, make a dummy vector.                                          #
      #------------------------------------------------------------------------------------#
      cat0("   - Ensure that all critical variables exist in 'raw'.")
      for (m in which(! must.exist %in% names(raw))){
         v.must = must.exist[m]
         raw[[v.must]] = rep(NA,times=length(raw$when))
      }#end for (m in which(! must.exist %in% names(raw))
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #    Discard obviously wrong or quite suspicious data.  This is based on visual      #
      # inspection plus looking at the log for the sites, and it is usually done after     #
      # running a first iteration keeping everything.  Right now, this is has to be hard-  #
      # -coded unfortunately.                                                              #
      #------------------------------------------------------------------------------------#
      if (discard.suspicious){
         cat0("   - Discard suspicious data.")
         if (aws$iata %in% "que"){
            #------------------------------------------------------------------------------#
            #     Querencia.                                                               #
            #------------------------------------------------------------------------------#


            #------------------------------------------------------------------------------#
            #     Pressure shifted after January 2014.                                     #
            #------------------------------------------------------------------------------#
            del = raw$when %>=% chron("01/01/2014")
            raw$atm.prss[del] = NA
            #------------------------------------------------------------------------------#


            #------------------------------------------------------------------------------#
            #     Radiation in 2014 has lots of missing data and some strange remaining    #
            # values.  Remove everything.                                                  #
            #------------------------------------------------------------------------------#
            del = raw$when %wr% chron(c("01/01/2014","01/01/2015"))
            raw$rshort.in[del] = NA
            #------------------------------------------------------------------------------#


            #------------------------------------------------------------------------------#
            #     Rainfall late 2009/early 2010 is extremely low between 2009 and 2010,    #
            # whereas it looks alright in Agua Boa.  Discard data.                         #
            #------------------------------------------------------------------------------#
            del = raw$when %wr% chron(c("09/15/2009","11/15/2010"))
            raw$rain[del] = NA
            #------------------------------------------------------------------------------#


            #------------------------------------------------------------------------------#
            #     Rainfall is pretty much zero for all 2014.                               #
            #------------------------------------------------------------------------------#
            del = raw$when %wr% chron(c("01/22/2014","01/01/2015"))
            raw$rain[del] = NA
            #------------------------------------------------------------------------------#

         }else if (aws$iata %in% "aub"){
            #------------------------------------------------------------------------------#
            #     Agua Boa.                                                                #
            #------------------------------------------------------------------------------#


            #------------------------------------------------------------------------------#
            #     Radiation looks totally off from 15 Mar 2013 onwards.                    #
            #------------------------------------------------------------------------------#
            del = raw$when %>=% chron("03/15/2013")
            raw$rshort.in[del] = NA
            #------------------------------------------------------------------------------#


            #------------------------------------------------------------------------------#
            #     Rainfall is pretty much zero for all 2014.                               #
            #------------------------------------------------------------------------------#
            del = raw$when %wr% chron(c("01/22/2014","01/01/2015"))
            raw$rain[del] = NA
            #------------------------------------------------------------------------------#
         }#end if (aws$iata %in% "que")
         #---------------------------------------------------------------------------------#
      }#end if (discard.suspicious)
      #------------------------------------------------------------------------------------#



      #----- Calculate partial pressure of water vapour. ----------------------------------#
      cat0("   - Find vapour pressure.")
      fill         = is.na(raw$atm.pvap)
      if ("atm.tdew" %in% names(raw)){
         raw$atm.pvap[fill] = eslif(raw$atm.tdew[fill])
         fill               = is.na(raw$atm.pvap)
      }#end if("atm.tdew" %in% names(raw))
      if ("atm.rhv" %in% names(raw)){
         raw$atm.pvap[fill] = raw$atm.rhv[fill] * eslif(raw$atm.tmp[fill])
         fill               = is.na(raw$atm.pvap)
      }#end if("atm.rhv" %in% names(raw))
      if ("atm.shv" %in% names(raw)){
         raw$atm.pvap[fill] = ( raw$atm.prss[fill] * raw$atm.shv[fill]
                              / ( ep + (1. - ep)*raw$atm.shv[fill] )
                              )#end atm.pvap
      }#end if ("atm.shv" %in% names(raw))
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #      Find wind speed components (so we can interpolate direction if needed be.     #
      #------------------------------------------------------------------------------------#
      raw$atm.uspd = - raw$atm.vels * sin(raw$atm.vdir * pio180)
      raw$atm.vspd = - raw$atm.vels * cos(raw$atm.vdir * pio180)
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #     Make sure data are hourly, keeping time stamp for the end of averaging period. #
      #------------------------------------------------------------------------------------#
      if (s.dtime < hr.sec){
         cat0("   - Obtain hourly averages.")
         #----- This index will be the same for times that can be combined. ---------------#
         elapidx = ceiling(raw$elapsed)
         #---------------------------------------------------------------------------------#

         #---------------------------------------------------------------------------------#
         #       The way we will represent the hourly average is by making an average of   #
         # the staggered hourly averages.                                                  #
         #---------------------------------------------------------------------------------#
         h.elapsed   = tapply(X=elapidx      ,INDEX=elapidx,FUN=max ,na.rm=FALSE)
         h.atm.tmp   = tapply(X=raw$atm.tmp  ,INDEX=elapidx,FUN=mean,na.rm=FALSE)
         h.atm.prss  = tapply(X=raw$atm.prss ,INDEX=elapidx,FUN=mean,na.rm=FALSE)
         h.atm.pvap  = tapply(X=raw$atm.pvap ,INDEX=elapidx,FUN=mean,na.rm=FALSE)
         h.rshort.in = tapply(X=raw$rshort.in,INDEX=elapidx,FUN=mean,na.rm=FALSE)
         h.rlong.in  = tapply(X=raw$rlong.in ,INDEX=elapidx,FUN=mean,na.rm=FALSE)
         h.atm.vels  = tapply(X=raw$atm.vels ,INDEX=elapidx,FUN=mean,na.rm=FALSE)
         h.atm.uspd  = tapply(X=raw$atm.uspd ,INDEX=elapidx,FUN=mean,na.rm=FALSE)
         h.atm.vspd  = tapply(X=raw$atm.vspd ,INDEX=elapidx,FUN=mean,na.rm=FALSE)
         h.atm.vdir  = ( 270. - atan2(x=h.atm.uspd,y=h.atm.vspd) / pio180 ) %% 360.
         h.rain      = tapply(X=raw$rain     ,INDEX=elapidx,FUN=mean,na.rm=FALSE)
         #---------------------------------------------------------------------------------#


         #---------------------------------------------------------------------------------#
         #       Standardise missing values to NA.                                         #
         #---------------------------------------------------------------------------------#
         h.elapsed   = ifelse(test=is.finite(h.elapsed  ),yes=h.elapsed  ,no=NA)
         h.atm.tmp   = ifelse(test=is.finite(h.atm.tmp  ),yes=h.atm.tmp  ,no=NA)
         h.atm.prss  = ifelse(test=is.finite(h.atm.prss ),yes=h.atm.prss ,no=NA)
         h.atm.pvap  = ifelse(test=is.finite(h.atm.pvap ),yes=h.atm.pvap ,no=NA)
         h.rshort.in = ifelse(test=is.finite(h.rshort.in),yes=h.rshort.in,no=NA)
         h.rlong.in  = ifelse(test=is.finite(h.rlong.in ),yes=h.rlong.in ,no=NA)
         h.atm.vels  = ifelse(test=is.finite(h.atm.vels ),yes=h.atm.vels ,no=NA)
         h.atm.uspd  = ifelse(test=is.finite(h.atm.uspd ),yes=h.atm.uspd ,no=NA)
         h.atm.vspd  = ifelse(test=is.finite(h.atm.vspd ),yes=h.atm.vspd ,no=NA)
         h.atm.vdir  = ifelse(test=is.finite(h.atm.vdir ),yes=h.atm.vdir ,no=NA)
         h.rain      = ifelse(test=is.finite(h.rain     ),yes=h.rain     ,no=NA)
         #---------------------------------------------------------------------------------#
      }else{
         #---------------------------------------------------------------------------------#
         #       Data are already hourly. Just standardise missing values.                 #
         #---------------------------------------------------------------------------------#
         h.elapsed   = ifelse(test=is.finite(raw$elapsed  ),yes=raw$elapsed  ,no=NA)
         h.atm.tmp   = ifelse(test=is.finite(raw$atm.tmp  ),yes=raw$atm.tmp  ,no=NA)
         h.atm.prss  = ifelse(test=is.finite(raw$atm.prss ),yes=raw$atm.prss ,no=NA)
         h.atm.pvap  = ifelse(test=is.finite(raw$atm.pvap ),yes=raw$atm.pvap ,no=NA)
         h.rshort.in = ifelse(test=is.finite(raw$rshort.in),yes=raw$rshort.in,no=NA)
         h.rlong.in  = ifelse(test=is.finite(raw$rlong.in ),yes=raw$rlong.in ,no=NA)
         h.atm.vels  = ifelse(test=is.finite(raw$atm.vels ),yes=raw$atm.vels ,no=NA)
         h.atm.vdir  = ifelse(test=is.finite(raw$atm.vdir ),yes=raw$atm.vdir ,no=NA)
         h.atm.uspd  = ifelse(test=is.finite(raw$atm.uspd ),yes=raw$atm.uspd ,no=NA)
         h.atm.vspd  = ifelse(test=is.finite(raw$atm.vspd ),yes=raw$atm.vspd ,no=NA)
         h.rain      = ifelse(test=is.finite(raw$rain     ),yes=raw$rain     ,no=NA)
         #---------------------------------------------------------------------------------#
      }#end if
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #     Map raw data onto the standardised output.                                     #
      #------------------------------------------------------------------------------------#
      cat0("   - Map raw data set onto standardised output.")
      idx = match(h.elapsed,elapsed)
      sel = is.finite(idx)
      aws$atm.tmp   [idx[sel]] =  h.atm.tmp  [sel]
      aws$atm.prss  [idx[sel]] =  h.atm.prss [sel]
      aws$atm.pvap  [idx[sel]] =  h.atm.pvap [sel]
      aws$rshort.in [idx[sel]] =  h.rshort.in[sel]
      aws$rlong.in  [idx[sel]] =  h.rlong.in [sel]
      aws$atm.vels  [idx[sel]] =  h.atm.vels [sel]
      aws$atm.vdir  [idx[sel]] =  h.atm.vdir [sel]
      aws$atm.uspd  [idx[sel]] =  h.atm.uspd [sel]
      aws$atm.vspd  [idx[sel]] =  h.atm.vspd [sel]
      aws$rain      [idx[sel]] =  h.rain     [sel]
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #      Find derived variables.                                                       #
      #------------------------------------------------------------------------------------#
      aws$atm.rhv      =   aws$atm.pvap / eslif(aws$atm.tmp)
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #     Find components of incoming radiation.                                         #
      #------------------------------------------------------------------------------------#
      cat0("   - Find shortwave radiation components.")
      rsbreak                     = rshort.bdown( rad.in   = aws$rshort.in
                                                , atm.prss = aws$atm.prss
                                                , cosz     = aws$cosz
                                                , rad.type = "rshort")
      aws$par.diff                = rsbreak$par.diff
      aws$nir.diff                = rsbreak$nir.diff
      aws$par.beam                = rsbreak$par.beam
      aws$nir.beam                = rsbreak$nir.beam
      #------------------------------------------------------------------------------------#




      #------------------------------------------------------------------------------------#
      #    Find the maximum possible radiation.  We will use the maximum measured pressure #
      # as the reference value.                                                            #
      #------------------------------------------------------------------------------------#
      pmean = range(aws$atm.prss,na.rm=TRUE)
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #    Find the lowest pressure amongst all measurements, and use it as the reference  #
      # pressure, to be on the conservative side.                                          #
      #------------------------------------------------------------------------------------#
      cat0("   - Find the maximum radiation.")
      press.min      = min(aws$atm.prss,na.rm=TRUE) + 0. * aws$cosz
      sunny          = rshort.bdown( rad.in   = aws$rshort.in
                                   , atm.prss = press.min
                                   , cosz     = aws$cosz
                                   , rad.type = "rshort"
                                   )
      aws$rshort.pot = sunny$rshort.max
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #    Force night-time data to be NA before filtering the data (they will become 0    #
      # afterwards).                                                                       #
      #------------------------------------------------------------------------------------#
      cat0("   - Temporarily discard nighttime radiation.")
      aws$rshort.in[aws$nighttime] = NA
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #      Loop over all variables.                                                      #
      #------------------------------------------------------------------------------------#
      cat0("   - Remove suspicious data from variables:")
      for (v in sequence(ntvar)){
         #----- Handy aliases. ------------------------------------------------------------#
         this.vnam      = tvar$vnam    [v]
         this.desc      = tvar$desc    [v]
         this.unit      = tvar$unit    [v]
         this.add       = tvar$add     [v]
         this.mult      = tvar$mult    [v]
         this.colmean   = tvar$colmean [v]
         this.colerr    = tvar$colerr  [v]
         this.coledge   = tvar$coledge [v]
         this.plt       = tvar$plt     [v]
         this.edge      = tvar$edge    [v]
         this.out.hour  = tvar$out.hour[v]
         this.out.all   = tvar$out.all [v]
         #---------------------------------------------------------------------------------#


         #----- Pick the variable and call the outlier remover function. ------------------#
         cat0("     > ",this.desc," (",this.vnam,").")
         thisvar          = aws[[this.vnam]]
         aws[[this.vnam]] = del.outliers( x    = thisvar
                                        , when = aws$when
                                        , out.hour = this.out.hour
                                        , out.all  = this.out.all
                                        )#end del.outliers
         #---------------------------------------------------------------------------------#
      }#end for (v in sequence(ntvar))
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #      There were two checks for humidity, so make sure to delete both if one of     #
      # them was considered weird.  In addition, delete relative humidity values in case   #
      # temperature was considered weird.                                                  #
      #------------------------------------------------------------------------------------#
      dl1               = ( is.na(aws$atm.pvap) 
                          | ( is.na(aws$atm.rhv) & is.finite(aws$atm.tmp) ) )
      dl2               = is.finite(aws$atm.rhv) & is.na(aws$atm.tmp)
      aws$atm.pvap[dl1] = NA
      aws$atm.rhv [dl1] = NA
      aws$atm.rhv [dl2] = NA
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #      Here we will remove data points and/or periods that we don't believe.         #
      #------------------------------------------------------------------------------------#
      #----- Remove suspicious incoming shortwave radiation. ------------------------------#
      cat0("   - Remove days with suspicious incoming shortwave radiation.")
      suspect             = as.integer(aws$daytime & aws$rshort.in  < -10.)
      testdays            = dates(sort(unique(aws$today)))
      nsuspects           = tapply(X=suspect,INDEX=aws$today,FUN=sum,na.rm=TRUE)
      weird               = testdays[nsuspects > 0]
      del                 = aws$today %in% weird
      aws$rshort.in [del] = NA
      cat0("   - Remove suspicious data from variables:")
      aurora                      = aws$nighttime & aws$rshort.in > 0.
      aurora[is.na(aurora)]       = FALSE
      aws$rshort.in [aurora]      = NA
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #    Make nighttime radiation zero.                                                  #
      #------------------------------------------------------------------------------------#
      aws$rshort.in[aws$nighttime] = 0.
      aws$par.beam [aws$nighttime] = 0.
      aws$nir.beam [aws$nighttime] = 0.
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #      Find the aggregated variables.                                                #
      #------------------------------------------------------------------------------------#
      cat0("   - Find daily and monthly time vectors.")
      aws$today   = chron(unique(aws$today))
      aws$tomonth = chron(unique(dates(aws$when[aws$day == 1])))
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #     Find several auxiliary variables for daily and monthly statistics...           #
      #------------------------------------------------------------------------------------#
      cat0("   - Find time-related variables for the daily and monthly statistics.")
      mm    = aws$month
      hh    = aws$hour
      hhmm  = paste(sprintf("%2.2i",hh),sprintf("%2.2i",mm),sep="")
      yy    = aws$year
      hrs   = sort(unique(hh))
      mons  = mon2mmm(sort(unique(mm)))
      yrs   = unique(aws$year)
      nhrs  = length(hrs)
      nmons = length(mons)
      nyrs  = length(yrs)
      #----- Find the maximum number of days for normal and leap years. -------------------#
      emean.dmax = daymax(nummonths(aws$tomonth),numyears(aws$tomonth))
      mmean.dmax = tapply(X=emean.dmax,INDEX=nummonths(aws$tomonth),FUN=mean,na.rm=TRUE)
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #      Find the mean diurnal cycle.                                                  #
      #------------------------------------------------------------------------------------#
      cat0("   - Find the mean diurnal cycles:")
      for (v in sequence(ntvar)){
         #----- Handy aliases. ------------------------------------------------------------#
         varname = tvar$vnam[v]
         vardesc = tvar$desc[v]
         varval  = aws[[varname]]
         varval2 = varval * varval
         cat0("    * ",vardesc,".")
         #---------------------------------------------------------------------------------#


         #---------------------------------------------------------------------------------#
         #     Find the statistics of the diurnal cycle.  These are for a skewed normal    #
         # distribution, so the variables are location, scale, and shape.  For rainfall,   #
         # we bypass the calculation because the distribution is extremely skewed.         #
         #---------------------------------------------------------------------------------#
         if (varname %in% "rain"){
            qlocation = rep(NA,times=length(unique(hhmm)))
            qscale    = rep(NA,times=length(unique(hhmm)))
            qshape    = rep(NA,times=length(unique(hhmm)))
         }else{
            qlocation = tapply(X=thisvar,INDEX=hhmm,FUN=sn.location ,na.rm=TRUE)
            qscale    = tapply(X=thisvar,INDEX=hhmm,FUN=sn.scale    ,na.rm=TRUE)
            qshape    = tapply(X=thisvar,INDEX=hhmm,FUN=sn.shape    ,na.rm=TRUE)
         }#end if
         #---------------------------------------------------------------------------------#



         okvar    = is.finite(varval)
         #----- Mean diurnal cycle. -------------------------------------------------------#
         qmean    = aggregate(varval ,by=list(mm,hh),FUN=mean,na.rm=TRUE,simplify=TRUE)
         qmean    = matrix   (qmean[,3],nrow=nmons,ncol=nhrs,byrow=FALSE)
         dimnames(qmean) = list(mons,hrs)
         #----- Mean sum of squares (diurnal cycle). --------------------------------------#
         qmsqu    = aggregate(varval2,by=list(mm,hh),FUN=mean,na.rm=TRUE,simplify=TRUE)
         qmsqu    = matrix   (qmsqu[,3],nrow=nmons,ncol=nhrs,byrow=FALSE)
         dimnames(qmsqu) = list(mons,hrs)
         #----- Standard deviation. -------------------------------------------------------#
         qsdev    = aggregate(varval  ,by=list(mm,hh),FUN=sd  ,na.rm=TRUE,simplify=TRUE)
         qsdev    = matrix   (qsdev[,3],nrow=nmons,ncol=nhrs,byrow=FALSE)
         dimnames(qsdev) = list(mons,hrs)
         #----- Valid data. ------------------------------------------------------#
         qccnt    = aggregate(okvar    ,by=list(mm,hh),FUN=sum  ,na.rm=TRUE,simplify=TRUE)
         qccnt    = matrix   (qccnt[,3],nrow=nmons,ncol=nhrs,byrow=FALSE)
         dimnames(qccnt) = list(mons,hrs)
         #---------------------------------------------------------------------------------#


         #----- Time series of monthly means. ---------------------------------------------#
         a.mean = aggregate(x=varval ,by=list(hh,mm,yy),FUN=mean,na.rm=TRUE,simplify=TRUE)
         a.msqu = aggregate(x=varval2,by=list(hh,mm,yy),FUN=mean,na.rm=TRUE,simplify=TRUE)
         a.ccnt = aggregate(x=okvar  ,by=list(hh,mm,yy),FUN=sum ,na.rm=TRUE,simplify=TRUE)
         nagg   = length(a.mean[,4])
         ind    = paste(sprintf("%4.4i",a.mean[,3]),sprintf("%2.2i",a.mean[,2])
                       ,sep="",concatenate="")
         emean   = tapply(X=a.mean[,4],INDEX=ind,FUN=mean,na.rm=FALSE)
         emsqu   = tapply(X=a.msqu[,4],INDEX=ind,FUN=mean,na.rm=FALSE)
         eccnt   = tapply(X=a.ccnt[,4],INDEX=ind,FUN=sum ,na.rm=FALSE)
         srnorm1 = sqrt(1./(1. - 1./eccnt))
         esdev   = sqrt((emsqu - emean*emean)*srnorm1)
         #---------------------------------------------------------------------------------#



         #----- Daily mean and standard deviation. ----------------------------------------#
         today = dates(aws$when)
         dmean = tapply(X=varval ,INDEX=today,FUN=mean,na.rm=FALSE)
         dmsqu = tapply(X=varval2,INDEX=today,FUN=mean,na.rm=FALSE)
         dsdev = tapply(X=varval ,INDEX=today,FUN=sd  ,na.rm=FALSE)
         dccnt = tapply(X=okvar  ,INDEX=today,FUN=sum ,na.rm=FALSE)
         #---------------------------------------------------------------------------------#



         #----- Monthly mean and standard deviation. --------------------------------------#
         mmean   = rowMeans(qmean,na.rm=TRUE)
         mmsqu   = rowMeans(qmsqu,na.rm=TRUE)
         mccnt   = rowMeans(qccnt,na.rm=TRUE)
         srnorm1 = sqrt(1./(1. - 1./qccnt))
         msdev   = rowMeans(sqrt(qmsqu-qmean*qmean)*srnorm1,na.rm=TRUE)
         #---------------------------------------------------------------------------------#



         #----- Change units for precipitation (mm/month) for monthly means. --------------#
         if (varname %in% "rain"){
            mmean = mmean * mmean.dmax
            mmsqu = mmsqu * mmean.dmax
            msdev = msdev * mmean.dmax
            emean = emean * emean.dmax
            emsqu = emsqu * emean.dmax
            esdev = esdev * emean.dmax
         }#end if
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #     Make unique names for the statistics.                                       #
         #---------------------------------------------------------------------------------#
         qlocation.name = paste("qlocation",varname,sep=".")
         qscale.name    = paste("qscale"   ,varname,sep=".")
         qshape.name    = paste("qshape"   ,varname,sep=".")
         qmean.name     = paste("qmean"    ,varname,sep=".")
         qmsqu.name     = paste("qmsqu"    ,varname,sep=".")
         qsdev.name     = paste("qsdev"    ,varname,sep=".")
         qccnt.name     = paste("qccnt"    ,varname,sep=".")
         dmean.name     = paste("dmean"    ,varname,sep=".")
         dmsqu.name     = paste("dmsqu"    ,varname,sep=".")
         dsdev.name     = paste("dsdev"    ,varname,sep=".")
         dccnt.name     = paste("dccnt"    ,varname,sep=".")
         mmean.name     = paste("mmean"    ,varname,sep=".")
         mmsqu.name     = paste("mmsqu"    ,varname,sep=".")
         msdev.name     = paste("msdev"    ,varname,sep=".")
         mccnt.name     = paste("mccnt"    ,varname,sep=".")
         emean.name     = paste("emean"    ,varname,sep=".")
         emsqu.name     = paste("emsqu"    ,varname,sep=".")
         esdev.name     = paste("esdev"    ,varname,sep=".")
         eccnt.name     = paste("eccnt"    ,varname,sep=".")
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #     Append the statistics to the list.                                          #
         #---------------------------------------------------------------------------------#
         aws[[qlocation.name]] = qlocation
         aws[[qscale.name   ]] = qscale
         aws[[qshape.name   ]] = qshape
         aws[[qmean.name    ]] = qmean
         aws[[qmsqu.name    ]] = qmsqu
         aws[[qsdev.name    ]] = qsdev
         aws[[qccnt.name    ]] = qccnt
         aws[[dmean.name    ]] = dmean
         aws[[dmsqu.name    ]] = dmsqu
         aws[[dsdev.name    ]] = dsdev
         aws[[dccnt.name    ]] = dccnt
         aws[[mmean.name    ]] = mmean
         aws[[mmsqu.name    ]] = mmsqu
         aws[[msdev.name    ]] = msdev
         aws[[mccnt.name    ]] = mccnt
         aws[[emean.name    ]] = emean
         aws[[emsqu.name    ]] = emsqu
         aws[[esdev.name    ]] = esdev
         aws[[eccnt.name    ]] = eccnt
         #---------------------------------------------------------------------------------#
      }#end for (v in sequence(ntvar))
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #     Save the dataset at the full structure.                                        #
      #------------------------------------------------------------------------------------#
      fitz.aws[[aws$iata]] = aws
      #------------------------------------------------------------------------------------#
   }#end for (s in sequence(nsites))
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Save the compounded dataset so we don't load them again.                          #
   #---------------------------------------------------------------------------------------#
   if (save.clean){
      cat0(" + Save AWS data to file ",basename(session.clean),".")
      save(fitz.aws,file=session.clean)
   }#end if
   #---------------------------------------------------------------------------------------#
}#end if (reload.clean && file.exists(session.clean)){
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#     Plotting output.                                                                     #
#------------------------------------------------------------------------------------------#
cat0(" + Plot averages:")
for (s in sequence(nsites)){
   #----- Handy aliases. ------------------------------------------------------------------#
   s.iata = sites$iata[s]
   aws    = fitz.aws[[s.iata]]
   #---------------------------------------------------------------------------------------#

   #----- Copy structure to local variables. ----------------------------------------------#
   outplace  = file.path(outroot,aws$short)
   placename = aws$long
   thisiata  = aws$iata
   thislon   = aws$lon
   thislat   = aws$lat
   yeara     = min(aws$year)
   yearz     = max(aws$year)

   outmdcyc   = file.path(outplace,"qmean")
   outmmean   = file.path(outplace,"mmean")
   outdmean   = file.path(outplace,"dmean")
   outemean   = file.path(outplace,"emean")
   outhours   = file.path(outplace,"hours")

   if (! file.exists(outplace)) dir.create(outplace)
   if (! file.exists(outmdcyc)) dir.create(outmdcyc)
   if (! file.exists(outmmean)) dir.create(outmmean)
   if (! file.exists(outdmean)) dir.create(outdmean)
   if (! file.exists(outemean)) dir.create(outemean)
   if (! file.exists(outhours)) dir.create(outhours)

   #---------------------------------------------------------------------------------------#
   #    Plot the diurnal cycle by month.                                                   #
   #---------------------------------------------------------------------------------------#
   cat0("   - Plot mean diurnal cycle for ",placename,".")
   for (tt in sequence(ntvar)){
      #----- Handy aliases. ---------------------------------------------------------------#
      vname   = tvar$vnam   [tt]
      desc    = tvar$desc   [tt]
      unit    = tvar$unit   [tt]
      unit    = untab[[unit]]
      add0    = tvar$add    [tt]
      mult0   = tvar$mult   [tt]
      colmean = tvar$colmean[tt]
      colerr  = tvar$colerr [tt]
      coledge = tvar$coledge[tt]
      edge    = tvar$edge   [tt]
      plotit  = tvar$plt    [tt]
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #     Check whether to plot it.                                                      #
      #------------------------------------------------------------------------------------#
      if (plotit){
         cat0("     > ",desc,".")


         #----- Create directory for this variable. ---------------------------------------#
         outvar = file.path(outmdcyc,vname)
         if (! file.exists(outvar)) dir.create(outvar)
         #---------------------------------------------------------------------------------#


         #----- Find average and standard deviation. --------------------------------------#
         qmean.name = paste("qmean",vname,sep=".")
         qsdev.name = paste("qsdev",vname,sep=".")
         thismean = mult0 * (aws[[qmean.name]]+add0)
         thissdev = mult0 *  aws[[qsdev.name]]
         #---------------------------------------------------------------------------------#



         #------ Organise months and time of day. -----------------------------------------#
         mmms  = rownames(thismean)
         mons  = mmm2mon(mmms)
         nmons = length(mons)
         hday  = as.numeric(colnames(thismean))
         #---------------------------------------------------------------------------------#


         #----- Fix limits for y axis so all months use the same scale. -------------------#
         if (plotsd){
            ylimit = pretty.xylim(u=c(thismean-thissdev,thismean+thissdev),fracexp=fracexp)
         }else{
            ylimit = pretty.xylim(u=thismean,fracexp=fracexp)
         }#end if
         yat     = pretty(ylimit)
         ylabels = sprintf("%g",yat)
         #---------------------------------------------------------------------------------#


         #---------------------------------------------------------------------------------#
         #     Define nice hours for axis marks.                                           #
         #---------------------------------------------------------------------------------#
         xlimit  = c(0,24)
         xat     = pretty.elapsed(x=xlimit,base=12)
         xlabels = sprintf("%g",xat)
         #---------------------------------------------------------------------------------#


         #---------------------------------------------------------------------------------#
         #     Loop through months.                                                        #
         #---------------------------------------------------------------------------------#
         for (mm in sequence(nmons)){
            #----- Handy aliases. ---------------------------------------------------------#
            cmon = sprintf("%2.2i",mons[mm])
            mmm  = month.name[mm]
            #------------------------------------------------------------------------------#


            #----- Load mean and standard deviation. --------------------------------------#
            monmean = thismean[mm,]
            monsdev = thissdev[mm,]
            err.x = c(hday,rev(hday))
            err.y = c(monmean+monsdev,rev(monmean-monsdev))
            #------------------------------------------------------------------------------#

            #----- Make title. ------------------------------------------------------------#
            letitre = paste0(placename,"(",yeara,"-",yearz,")","\n",desc," - ",mmm)
            lex     = "Time [UTC]"
            ley     = desc.unit(desc=desc,unit=unit)
            #------------------------------------------------------------------------------#


            #------------------------------------------------------------------------------#
            #     Loop through all formats.                                                #
            #------------------------------------------------------------------------------#
            for (o in sequence(nout)){
               #----- Open device. --------------------------------------------------------#
               fichier = file.path(outvar,paste0(vname,"-",cmon,".",outform[o]))
               dummy   = open.plot( fichier = fichier
                                  , outform = outform[o]
                                  , size    = sq.size
                                  , ptsz    = ptsz
                                  , depth   = depth
                                  )#end open.plot
               #---------------------------------------------------------------------------#



               #---------------------------------------------------------------------------#
               #     Plot main figure.                                                     #
               #---------------------------------------------------------------------------#
               par(par.user)
               par(mar=c(4.1,4.6,3.1,1.1))
               plot.new()
               plot.window(xlim=xlimit,ylim=ylimit)
               axis(side=1,las=1,at=xat,labels=xlabels)
               axis(side=2,las=1,at=yat,labels=ylabels)
               title(main=letitre,xlab=lex,ylab=ley)
               if (plotgrid) abline(h=yat,v=xat,col="gray62",lwd=0.25,lty="dotted")
               if (plotsd){
                  polygon( x     = err.x
                         , y     = err.y
                         , col   = colerr
                         , angle = angle
                         , dens  = dens
                         , lty   = "solid"
                         , lwd   =  shwd
                         )#end polygon
               }#end if (plotsd)
               lines(x=hday,y=monmean,type="o",pch=16,cex=1.0,lwd=llwd,col=colmean)
               if (plotsd){
                  legend( x       = "topleft"
                        , inset   = 0.01
                        , legend  = c("Mean","1 Sdev")
                        , fill    = c("transparent",colerr)
                        , angle   = angle
                        , density = c(0,dens)
                        , lwd     = llwd
                        , col     = c(colmean,colerr)
                        , bg      = "white"
                        , lty     = c("solid",NA)
                        , cex     = 1.0
                        , pch     = c(16,NA)
                        )#end if
               }#end if (plotsd)
               box()
               dummy = close.plot(outform=outform[o])
               #---------------------------------------------------------------------------#
            }#end for (o in sequence(nout))
            #------------------------------------------------------------------------------#
         }#end for (m in sequence(nmons))
         #---------------------------------------------------------------------------------#
      }#end if (plotit)
      #------------------------------------------------------------------------------------#
   }#end for (tt in sequence(ntvar)
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #    Plot the monthly means.                                                            #
   #---------------------------------------------------------------------------------------#
   cat0("   - Plot monthly means for ",placename,".")
   for (tt in sequence(ntvar)){
      #----- Handy aliases. ---------------------------------------------------------------#
      vname   = tvar$vnam   [tt]
      desc    = tvar$desc   [tt]
      unit    = tvar$unit   [tt]
      if (vname %in% "rain"){
         unit = untab$mmomo
      }else{
         unit = untab[[unit]]
      }#end if
      add0    = tvar$add    [tt]
      mult0   = tvar$mult   [tt]
      colmean = tvar$colmean[tt]
      colerr  = tvar$colerr [tt]
      coledge = tvar$coledge[tt]
      edge    = tvar$edge   [tt]
      plotit  = tvar$plt    [tt]
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #     Proceed with plotting.                                                         #
      #------------------------------------------------------------------------------------#
      if (plotit){
         cat0("     > ",desc,".")


         #----- Find average and standard deviation. --------------------------------------#
         mmean.name = paste("mmean",vname,sep=".")
         msdev.name = paste("msdev",vname,sep=".")
         monmean    = mult0 * (aws[[mmean.name]] + add0 )
         monsdev    = mult0 *  aws[[msdev.name]]
         #---------------------------------------------------------------------------------#




         #------ Organise months and time of day. -----------------------------------------#
         xlimit = c(1,12)
         mons   = sort(unique(aws$mon))
         nmons  = length(mons)
         mlab   = mon2mmm(mons)
         #---------------------------------------------------------------------------------#




         #----- Fix limits for y axis so all months use the same scale. -------------------#
         if ((! plotsd) || vname %in% c("rain","monrain")){
            ylimit = pretty.xylim(u=monmean,fracexp=fracexp)
         }else{
            ylimit = pretty.xylim(u=c(monmean-monsdev,monmean+monsdev),fracexp=fracexp)
         }#end if
         yat     = pretty(ylimit)
         ylabels = sprintf("%g",yat)
         #---------------------------------------------------------------------------------#

         #------ Find polygons. -----------------------------------------------------------#
         err.x = c(mons,rev(mons))
         err.y = c(monmean+monsdev,rev(monmean-monsdev))
         #---------------------------------------------------------------------------------#



         #------ Plot annotation. ---------------------------------------------------------#
         letitre = paste0(placename,"(",yeara,"-",yearz,")","\n",desc)
         lex     = "Month"
         ley     = desc.unit(desc=desc,unit=unit)
         #---------------------------------------------------------------------------------#


         #---------------------------------------------------------------------------------#
         #     Loop through all formats.                                                   #
         #---------------------------------------------------------------------------------#
         for (o in sequence(nout)){
            #----- Open device. -----------------------------------------------------------#
            fichier = file.path(outmmean,paste0(vname,".",outform[o]))
            dummy   = open.plot( fichier = fichier
                               , outform = outform[o]
                               , size    = sq.size
                               , ptsz    = ptsz
                               , depth   = depth
                               )#end open.plot
            #------------------------------------------------------------------------------#



            #------------------------------------------------------------------------------#
            #     Plot main figure.                                                        #
            #------------------------------------------------------------------------------#
            par(par.user)
            par(mar=c(4.1,4.6,3.1,1.1))
            plot.new()
            plot.window(xlim=xlimit,ylim=ylimit)
            axis(side=1,las=1,at=mons,labels=mlab)
            axis(side=2,las=1,at=yat,labels=ylabels)
            title(main=letitre,xlab=lex,ylab=ley)
            if (plotgrid) abline(v=mons,h=yat,lwd=0.25,col="gray62",lty="dotted")
            if (plotsd){
               polygon(x=err.x,y=err.y,col=colerr,angle=angle,dens=dens,lty="solid"
                      ,lwd=shwd)
            }#end if
            lines(x=mons,y=monmean,type="o",pch=16,cex=1.0,lwd=llwd,col=colmean)
            if (plotsd){
               legend( x       = "topleft"
                     , inset   = 0.01
                     , legend  = c("Mean","1 Sdev")
                     , fill    = c("transparent",colerr)
                     , angle   = angle
                     , density = c(0,dens)
                     , lwd     = llwd
                     , col     = c(colmean,colerr)
                     , bg      = "white"
                     , lty     = c("solid",NA)
                     , cex     = 1.0
                     , pch     = c(16,NA)
                     )#end legend
            }#end if
            box()
            dummy = close.plot(outform=outform[o])
            #------------------------------------------------------------------------------#
         }#end for (o in sequence(nout))
         #---------------------------------------------------------------------------------#
      }#end if (plotit)
      #------------------------------------------------------------------------------------#
   }#end for (tt in sequence(ntvar))
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #    Plot the daily means.                                                              #
   #---------------------------------------------------------------------------------------#
   cat0("   - Plot daily means for ",placename,".")
   for (tt in sequence(ntvar)){
      #----- Handy aliases. ---------------------------------------------------------------#
      vname   = tvar$vnam   [tt]
      desc    = tvar$desc   [tt]
      unit    = tvar$unit   [tt]
      if (vname %in% "rain"){
         unit = untab$mmoday
      }else{
         unit = untab[[unit]]
      }#end if
      add0    = tvar$add    [tt]
      mult0   = tvar$mult   [tt]
      colmean = tvar$colmean[tt]
      colerr  = tvar$colerr [tt]
      coledge = tvar$coledge[tt]
      edge    = tvar$edge   [tt]
      plotit  = tvar$plt    [tt]
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #     Proceed with plotting.                                                         #
      #------------------------------------------------------------------------------------#
      if (plotit){
         cat0("     > ",desc,".")


         #----- Find average and standard deviation. --------------------------------------#
         dmean.name = paste("dmean",vname,sep=".")
         dsdev.name = paste("dsdev",vname,sep=".")
         today   = aws$today
         daymean = mult0 * ( aws[[dmean.name]] + add0 )
         daysdev = mult0 *   aws[[dsdev.name]]
         #---------------------------------------------------------------------------------#


         #---------------------------------------------------------------------------------#
         #      Make the time axis.                                                        #
         #---------------------------------------------------------------------------------#
         days     = aws$today
         ndays    = length(days)
         whenplot = pretty.time(days,n=8)
         xlimit   = range(days)
         #---------------------------------------------------------------------------------#




         #----- Fix limits for y axis so all months use the same scale. -------------------#
         ylimit  = pretty.xylim(u=daymean,fracexp=fracexp)
         yat     = pretty(ylimit)
         ylabels = sprintf("%g",yat)
         #---------------------------------------------------------------------------------#



         #------ Plot annotation. ---------------------------------------------------------#
         letitre = paste0(placename,"(",yeara,"-",yearz,")","\n",desc)
         lex     = "Dates"
         ley     = desc.unit(desc=desc,unit=unit)
         #---------------------------------------------------------------------------------#


         #---------------------------------------------------------------------------------#
         #     Loop through all formats.                                                   #
         #---------------------------------------------------------------------------------#
         for (o in sequence(nout)){
            #----- Open device. -----------------------------------------------------------#
            fichier = file.path(outdmean,paste0(vname,".",outform[o]))
            dummy   = open.plot( fichier = fichier
                               , outform = outform[o]
                               , size    = sq.size
                               , ptsz    = ptsz
                               , depth   = depth
                               )#end open.plot
            #------------------------------------------------------------------------------#



            #------------------------------------------------------------------------------#
            #     Plot main figure.                                                        #
            #------------------------------------------------------------------------------#
            par(par.user)
            par(mar=c(4.1,4.6,3.1,1.1))
            plot.new()
            plot.window(xlim=xlimit,ylim=ylimit)
            axis(side=1,at=whenplot$levels,labels=whenplot$labels,padj=whenplot$padj)
            axis(side=2,las=1,at=yat,labels=ylabels)
            title(main=letitre,xlab=lex,ylab=ley)
            if (plotgrid) abline(v=whenplot$levels,h=yat,lwd=0.25,col="gray62",lty="dotted")
            lines  (x=days,y=daymean,type="o",pch=16,cex=1.0,lwd=llwd,col=colmean)
            box()
            dummy = close.plot(outform=outform[o])
            #------------------------------------------------------------------------------#
         }#end for (o in sequence(nout))
         #---------------------------------------------------------------------------------#
      }#end if (plotit)
      #------------------------------------------------------------------------------------#
   }#end for (tt in sequence(ntvar))
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #    Plot the time series of monthly means.                                             #
   #---------------------------------------------------------------------------------------#
   cat0("   - Plot time series of monthly means for ",placename,".")
   for (tt in sequence(ntvar)){
      #----- Handy aliases. ---------------------------------------------------------------#
      vname   = tvar$vnam   [tt]
      desc    = tvar$desc   [tt]
      unit    = tvar$unit   [tt]
      if (vname %in% "rain"){
         unit = untab$mmomo
      }else{
         unit = untab[[unit]]
      }#end if
      add0    = tvar$add    [tt]
      mult0   = tvar$mult   [tt]
      colmean = tvar$colmean[tt]
      colerr  = tvar$colerr [tt]
      coledge = tvar$coledge[tt]
      edge    = tvar$edge   [tt]
      plotit  = tvar$plt    [tt]
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #     Proceed with plotting.                                                         #
      #------------------------------------------------------------------------------------#
      if (plotit){
         cat0("     > ",desc,".")


         #----- Find average and standard deviation. --------------------------------------#
         emean.name = paste("emean",vname,sep=".")
         esdev.name = paste("esdev",vname,sep=".")
         mtsmean    = mult0 * ( aws[[emean.name]] + add0 )
         mtssdev    = mult0 *   aws[[esdev.name]]
         tomonth    = aws$tomonth
         #---------------------------------------------------------------------------------#


         #---------------------------------------------------------------------------------#
         #      Make the time axis.                                                        #
         #---------------------------------------------------------------------------------#
         nmon     = length(tomonth)
         whenplot = pretty.time(tomonth,n=8)
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #     Find the limits for the y-axis.                                             #
         #---------------------------------------------------------------------------------#
         if ((! plotsd) || vname == "rain" || vname == "monrain"){
            ylimit = pretty.xylim(u=mtsmean,fracexp=fracexp)
         }else{
            ylimit = pretty.xylim(u=c(mtsmean-mtssdev,mtsmean+mtssdev),fracexp=fracexp)
         }#end if ((! plotsd) || vname == "rain" || vname == "monrain")
         #---------------------------------------------------------------------------------#




         #------ Find polygons. -----------------------------------------------------------#
         err.x = c(tomonth,rev(tomonth))
         err.y = c(mtsmean+mtssdev,rev(mtsmean-mtssdev))
         #---------------------------------------------------------------------------------#



         #------ Plot annotation. ---------------------------------------------------------#
         letitre = paste0(placename,"(",yeara,"-",yearz,")","\n",desc)
         lex     = "Time"
         ley     = desc.unit(desc=desc,unit=unit)
         #---------------------------------------------------------------------------------#


         #---------------------------------------------------------------------------------#
         #     Loop through all formats.                                                   #
         #---------------------------------------------------------------------------------#
         for (o in sequence(nout)){
            #----- Open device. -----------------------------------------------------------#
            fichier = file.path(outemean,paste0(vname,".",outform[o]))
            dummy   = open.plot( fichier = fichier
                               , outform = outform[o]
                               , size    = sq.size
                               , ptsz    = ptsz
                               , depth   = depth
                               )#end open.plot
            #------------------------------------------------------------------------------#



            #------------------------------------------------------------------------------#
            #     Plot main figure.                                                        #
            #------------------------------------------------------------------------------#
            par(par.user)
            par(mar=c(4.1,4.6,3.1,1.1))
            plot.new()
            plot.window(xlim=xlimit,ylim=ylimit)
            axis(side=1,at=whenplot$levels,labels=whenplot$labels,padj=whenplot$padj)
            axis(side=2,las=1,at=yat,labels=ylabels)
            title(main=letitre,xlab=lex,ylab=ley)
            if (plotgrid) abline(v=whenplot$levels,h=yat,lwd=0.25,col="gray62",lty="dotted")
            if (plotsd){
               polygon(x=err.x,y=err.y,col=colerr,angle=angle,dens=dens,lty="solid"
                      ,lwd=shwd)
            }#end if
            lines(x=tomonth,y=mtsmean,type="l",cex=1.0,lwd=llwd,col=colmean)
            if (plotsd){
               legend( x       = "topleft"
                     , inset   = 0.01
                     , legend  = c("Mean","1 Sdev")
                     , fill    = c("transparent",colerr)
                     , angle   = angle
                     , density = c(0,dens)
                     , lwd     = llwd
                     , col     = c(colmean,colerr)
                     , bg      = "white"
                     , lty     = c("solid",NA)
                     , cex     = 1.0
                     )#end legend
            }#end if
            box()
            dummy = close.plot(outform=outform[o])
            #------------------------------------------------------------------------------#
         }#end for (o in sequence(nout))
         #---------------------------------------------------------------------------------#
      }#end if (plotit)
      #------------------------------------------------------------------------------------#
   }#end for (tt in sequence(ntvar))
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #    Plot entire time series.                                                           #
   #---------------------------------------------------------------------------------------#
   cat0("   - Plot hourly data for ",placename,".")
   for (tt in sequence(ntvar)){
      #----- Handy aliases. ---------------------------------------------------------------#
      vname   = tvar$vnam   [tt]
      desc    = tvar$desc   [tt]
      unit    = tvar$unit   [tt]
      if (vname %in% "rain"){
         unit = untab$mmoday
      }else{
         unit = untab[[unit]]
      }#end if
      add0    = tvar$add    [tt]
      mult0   = tvar$mult   [tt]
      colmean = tvar$colmean[tt]
      colerr  = tvar$colerr [tt]
      coledge = tvar$coledge[tt]
      edge    = tvar$edge   [tt]
      plotit  = tvar$plt    [tt]
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #     Proceed with plotting.                                                         #
      #------------------------------------------------------------------------------------#
      if (plotit){
         cat0("     > ",desc,".")
         if (vname %in% "rain"){
            varnow  = mult0 * ( aws[[vname]] + add0 ) / 24.
         }else{
            varnow  = mult0 * ( aws[[vname]] + add0 )
         }#end if
         when    = aws$when
         #---------------------------------------------------------------------------------#


         #---------------------------------------------------------------------------------#
         #      Make the time axis.                                                        #
         #---------------------------------------------------------------------------------#
         nwhen    = length(when)
         whenplot = pretty.time(when,n=8)
         xlimit   = range(when)
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #     Find the limits for the y-axis.                                             #
         #---------------------------------------------------------------------------------#
         ylimit  = pretty.xylim(u=varnow,fracexp=fracexp)
         yat     = pretty(ylimit)
         ylabels = sprintf("%g",yat)
         #---------------------------------------------------------------------------------#



         #------ Plot annotation. ---------------------------------------------------------#
         letitre = paste0(placename,"(",yeara,"-",yearz,")","\n",desc)
         lex     = "Dates"
         ley     = desc.unit(desc=desc,unit=unit)
         #---------------------------------------------------------------------------------#


         #---------------------------------------------------------------------------------#
         #     Loop through all formats.                                                   #
         #---------------------------------------------------------------------------------#
         for (o in sequence(nout)){
            #----- Open device. -----------------------------------------------------------#
            fichier = file.path(outhours,paste0(vname,".",outform[o]))
            dummy   = open.plot( fichier = fichier
                               , outform = outform[o]
                               , size    = sq.size
                               , ptsz    = ptsz
                               , depth   = depth
                               )#end open.plot
            #------------------------------------------------------------------------------#



            #------------------------------------------------------------------------------#
            #     Plot main figure.                                                        #
            #------------------------------------------------------------------------------#
            par(par.user)
            par(mar=c(4.1,4.6,3.1,1.1))
            plot.new()
            plot.window(xlim=xlimit,ylim=ylimit)
            axis(side=1,at=whenplot$levels,labels=whenplot$labels,padj=whenplot$padj)
            axis(side=2,las=1,at=yat,labels=ylabels)
            title(main=letitre,xlab=lex,ylab=ley)
            if (plotgrid) abline(v=whenplot$levels,h=yat,lwd=0.25,col="gray62",lty="dotted")
            lines  (x=when,y=varnow,lwd=1,col=colmean)
            box()
            dummy = close.plot(outform=outform[o])
            #------------------------------------------------------------------------------#
         }#end for (o in sequence(nout))
         #---------------------------------------------------------------------------------#
      }#end if (plotit)
      #------------------------------------------------------------------------------------#
   }#end for (tt in sequence(ntvar))
   #---------------------------------------------------------------------------------------#
}#end for (s in sequence(nsites))
#==========================================================================================#
#==========================================================================================#
