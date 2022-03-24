#------------------------------------------------------------------------------------------#
#     This script reads in the version 2.0 data from Natalia Coupe and creates the         #
# meteorological drivers.                                                                  #
#------------------------------------------------------------------------------------------#
rm(list=ls())
options(warn=0)
gc()
graphics.off()
#------------------------------------------------------------------------------------------#



#----- Default settings. ------------------------------------------------------------------#
here    = getwd()                  # Current directory.
outpath = file.path(here,"sites")  # Output path
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#    List of places, with the name, first and last full years, and the output variable .   #
#                                                                                          #
# name     - Name of the site for input file name (which should be filled_<name>.csv).     #
# longname - Longer site name description (for labels).                                    #
# lon      - Site longitude                                                                #
# lat      - Site latitude                                                                 #
# h5pref   - Prefix for meteorological driver file names (output files)                    #
# height   - Reference height (make sure to be above the height of the tallest             #
#            possible tree)                                                                #
# dtdat    - Time interval for output                                                      #
# imetavg  - What is the meaning of the time stamps in the input file. The examples below  #
#            are for hourly data sets and what the time stamp "2022-03-23 16:00 UTC"       #
#            corresponds to:                                                               #
#            1 - End at time stamp       (2022-03-23 15:00 UTC to 2022-03-23 16:00 UTC)    #
#            2 - Beginning at time stamp (2022-03-23 16:00 UTC to 2022-03-23 17:00 UTC)    #
#            3 - Middle at time stamp    (2022-03-23 15:30 UTC to 2022-03-23 16:30 UTC)    #
#                                                                                          #
# *** IMPORTANT NOTES ***                                                                  #
#                                                                                          #
# 1 - Longitude and latitude are important for the radiation components model, set both    #
#     coordinates as accurately as possible,                                               #
# 2 - This script assumes that all times are in UTC. If they aren't, you must preprocess   #
#     the code to make them in UTC                                                         #
# 3 - This script will not fill gaps in the meteorological drivers. This must be done      #
#     externally, before using this script.                                                #
#------------------------------------------------------------------------------------------#
n           = 0
place       = list()
n           = n + 1
place[[n]]  = list( name     = "santarem_km83"
                  , longname = "Santarem - Km 83"
                  , lon      = -54.971
                  , lat      =  -3.018
                  , h5pref   = "Santarem_Km83"
                  , height   = 64
                  , dtdat    = 3600.
                  , imetavg  = 1
                  )#end if
#------------------------------------------------------------------------------------------#



#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#     No need to change anything beyond this point unless you are developing the code.     #
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#



#----- Find the number of places to make. -------------------------------------------------#
nplaces      = length(place)
#------------------------------------------------------------------------------------------#



#----- Avoid unecessary and extremely annoying beeps. -------------------------------------#
options(locatorBell=FALSE)
#------------------------------------------------------------------------------------------#



#----- Load packages. ---------------------------------------------------------------------#
isok.chron = require(chron)
isok.rhdf5 = require(rhdf5)
if (! all(c(isok.chron,isok.rhdf5))){
   cat(" Package \"chron\" was found: ",isok.chron,".","\n",sep="")
   cat(" Package \"rhdf5\" was found: ",isok.rhdf5,".","\n",sep="")
   stop("Install missing packages before running this script.")
}#end if (! all(c(isok.chron,isok.rhdf5)))
#------------------------------------------------------------------------------------------#



#----- Load useful functions. -------------------------------------------------------------#
source("timeutils.r")
source("radutils.r")
#------------------------------------------------------------------------------------------#


#----- Check that output path exists. -----------------------------------------------------#
if (! file.exists(outpath)) dir.create(outpath)
#------------------------------------------------------------------------------------------#



for (p in sequence(nplaces)){
   #----- Copy structure to local variables. ----------------------------------------------#
   info       = place[[p]]
   file.in    = file.path(here,paste0("filled_",info$name,".txt"))
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #     Read the data.                                                                    #
   #---------------------------------------------------------------------------------------#
   datum = read.table(file=file.in,sep=",",na="NaN",header=TRUE,skip=47)
   cat(" + Process data from ",info$longname,".","\n",sep="")
   #---------------------------------------------------------------------------------------#
   #=======================================================================================#
   #=======================================================================================#


   #---------------------------------------------------------------------------------------#
   #     Find the dates.                                                                   #
   #---------------------------------------------------------------------------------------#
   datum$when    = chron( paste(datum$month,datum$day,datum$year,sep="/")
                        , paste(datum$hour ,datum$min,datum$sec ,sep=":") )
   datum$today   = chron( paste(datum$month,datum$day,datum$year,sep="/") )
   datum$tomonth = chron( paste(datum$month,        1,datum$year,sep="/") )
   #---------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   #     Find the zenith angle, to split radiation into components.                        #
   #---------------------------------------------------------------------------------------#
   zen           = ed.zen(lon=info$lon,lat=info$lat,when=datum$when,ed21=TRUE
                         ,zeronight=FALSE,meanval=TRUE,imetavg=info$imetavg,nmean=60)
   datum$cosz    = zen$cosz
   #---------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   #     Split the incoming radiation into components, or estimate NIR from PAR if PAR     #
   # measurements are available.                                                           #
   #---------------------------------------------------------------------------------------#
   prss = mean(datum$atm.prss)
   rad  = rshort.bdown(rad.in=datum$rshort.in,atm.prss=datum$atm.prss,cosz=datum$cosz)
   datum$par.beam = rad$par.beam
   datum$par.diff = rad$par.diff
   datum$nir.beam = rad$nir.beam
   datum$nir.diff = rad$nir.diff
   #---------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   #      Decompose wind.                                                                  #
   #---------------------------------------------------------------------------------------#
   trigo          = (270. - datum$atm.vdir) * pio180
   datum$atm.uspd = datum$atm.vels * cos(trigo)
   datum$atm.vspd = datum$atm.vels * sin(trigo)
   #---------------------------------------------------------------------------------------#






   #=======================================================================================#
   #=======================================================================================#
   #=======================================================================================#
   #=======================================================================================#
   #=======================================================================================#
   #=======================================================================================#
   #=======================================================================================#
   #=======================================================================================#
   #     Make ED output.                                                                   #
   #---------------------------------------------------------------------------------------#
   cat(" + Make ED-friendly output files.\n")

   #----- Make sure that the output directory exists, and if not, create it. --------------#
   siteroot = file.path(outpath,info$h5pref)
   if (! file.exists(siteroot)) dir.create(siteroot)
   #---------------------------------------------------------------------------------------#

   #----- List all possible unique month/year combinations. -------------------------------#
   unique.tomonth   = unique(datum$tomonth)
   n.unique.tomonth = length(unique.tomonth)
   #---------------------------------------------------------------------------------------#

   for (um in sequence(n.unique.tomonth)){
      monyear.now = unique.tomonth[um]

      #----- Get current month and year, 3-letter month, and number of days in the month. -#
      month.now   = nummonths(monyear.now)
      year.now    = numyears (monyear.now)
      daymax.now  = daymax   (month=month.now,year=year.now)
      month.label = toupper(month.abb[month.now])
      year.label  = sprintf("%4.4i",year.now)
      #------------------------------------------------------------------------------------#


      #----- Print banner to entertain the user. ------------------------------------------#
      cat("   - Check data from ",month.name[month.now]," ",year.now,".","\n",sep="")
      #------------------------------------------------------------------------------------#


      #----- Find the indices of data that belong to this month and year. -----------------#
      sel    = datum$month == month.now & datum$year == year.now
      #------------------------------------------------------------------------------------#


      #----- Check that all data are there. -----------------------------------------------#
      nsel      = sum(sel)
      nexpected = daymax.now * day.sec / info$dtdat
      #------------------------------------------------------------------------------------#

      #----- If the data are complete, make the output file. ------------------------------#
      if (nsel == nexpected){
         cat("     * Data time series is complete, making the arrays.","\n")
         #---------------------------------------------------------------------------------#
         #      Create the matrices that will have the data.  These will use ED/NCEP name  #
         # convention, and will have a fake 2x2 matrix with the same data because the      #
         # tower may be needed in a site that is nearby but not with the same longitude    #
         # and latitude.                                                                   #
         #---------------------------------------------------------------------------------#
         dlwrf = array(datum$rlong.in[sel],dim=c(nsel,2,2))
         nbdsf = array(datum$nir.beam[sel],dim=c(nsel,2,2))
         nddsf = array(datum$nir.diff[sel],dim=c(nsel,2,2))
         prate = array(datum$rain    [sel],dim=c(nsel,2,2))
         pres  = array(datum$atm.prss[sel],dim=c(nsel,2,2))
         sh    = array(datum$atm.shv [sel],dim=c(nsel,2,2))
         tmp   = array(datum$atm.tmp [sel],dim=c(nsel,2,2))
         ugrd  = array(datum$atm.uspd[sel],dim=c(nsel,2,2))
         vbdsf = array(datum$par.beam[sel],dim=c(nsel,2,2))
         vddsf = array(datum$par.diff[sel],dim=c(nsel,2,2))
         vgrd  = array(datum$atm.vspd[sel],dim=c(nsel,2,2))
         #---------------------------------------------------------------------------------#

         #---------------------------------------------------------------------------------#
         #     Make sure that all data are there.  If not, crash!                          #
         #---------------------------------------------------------------------------------#
         if (any(is.na(dlwrf))) stop("dlwrf has missing values!")
         if (any(is.na(nbdsf))) stop("nbdsf has missing values!")
         if (any(is.na(nddsf))) stop("nddsf has missing values!")
         if (any(is.na(prate))) stop("prate has missing values!")
         if (any(is.na(pres ))) stop("pres  has missing values!")
         if (any(is.na(sh   ))) stop("sh    has missing values!")
         if (any(is.na(tmp  ))) stop("tmp   has missing values!")
         if (any(is.na(ugrd ))) stop("ugrd  has missing values!")
         if (any(is.na(vbdsf))) stop("vbdsf has missing values!")
         if (any(is.na(vddsf))) stop("vddsf has missing values!")
         if (any(is.na(vgrd ))) stop("vgrd  has missing values!")
         #---------------------------------------------------------------------------------#


         #---------------------------------------------------------------------------------#
         #     Make the dataset.                                                           #
         #---------------------------------------------------------------------------------#
         h5.metdriv = paste0(info$h5pref,"_",year.label,month.label,".h5")
         h5.metdriv = file.path(siteroot,h5.metdriv)
         cat("     * Save data to ",basename(h5.metdriv),".","\n",sep="")
         if (file.exists(h5.metdriv)) dummy   = file.remove(h5.metdriv)
         dummy   = h5save( dlwrf,nbdsf,nddsf,prate,pres,sh,tmp,ugrd,vbdsf,vddsf,vgrd
                         , file=h5.metdriv
                         )#end h5save
         dummy   = H5close()
         #---------------------------------------------------------------------------------#
      }else{
         stop (" ---> Could not find all the data for this month!")
      }#end if
   }#end for
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #     Make the header for this experiment.                                              #
   #---------------------------------------------------------------------------------------#
   header  = file.path(siteroot,"ED_MET_DRIVER_HEADER")
   info    = c("# See README at the bottom of this file."
              ,"1"
              ,file.path(siteroot,paste0(info$h5pref,"_"))
              ,paste (2,2,1.0,1.0,info$lon,info$lat)
              ,"12"
              ,paste0("   'hgt'   'tmp'  'pres'    'sh'  'ugrd'  'vgrd' 'prate'"
                     ," 'dlwrf' 'nbdsf' 'nddsf' 'vbdsf' 'vddsf'")
              ,paste (c(info$height,rep(info$dtdat,times=11)),collapse=" ")
              ,paste (" 4 1 1 1 1 1 0 1 1 1 1 1")
              ,paste0(" ")
              ,paste0("!===========================================================!")
              ,paste0("! README                                                    !")
              ,paste0("!===========================================================!")
              ,paste0("!     The header of the meteorological driver must contain  !")
              ,paste0("! the following lines:                                      !")
              ,paste0("!                                                           !")
              ,paste0("! Line  1 : Banner, it will not be read;                    !")
              ,paste0("! Line  2 : Number of file formats, hereafter N;            !")
              ,paste0("! Lines 3+: For each of the N formats, add the following    !")
              ,paste0("!           lines, going through a-f for the first format,  !")
              ,paste0("!           then through a-f for the second format and so   !")
              ,paste0("!            on:                                            !")
              ,paste0("!    a. Prefixes of the file format;                        !")
              ,paste0("!    b. nlon, nlat, deltalon, deltalat, lon0, lat0.  If     !")
              ,paste0("!       lon and lat are also variables, only nlon and nlat  !")
              ,paste0("!       will be used;                                       !")
              ,paste0("!    c. Number of variables contained in this format;       !")
              ,paste0("!    d. List of variables for each format (see Table 1);    !")
              ,paste0("!    e. Frequency at which vares are updated, or the        !")
              ,paste0("!       constant value if the variable type is 4;           !")
              ,paste0("!    f. Variable type (see Table 2);                        !")
              ,paste0("!                                                           !")
              ,paste0("!===========================================================!")
              ,paste0("! Table 1. Variable names recognized by ED.                 !")
              ,paste0("!===========================================================!")
              ,paste0("! -> lon    -  Longitude                        [      deg] !")
              ,paste0("! -> lat    -  Latitude                         [      deg] !")
              ,paste0("! -> land*  -  Land fraction                    [      0-1] !")
              ,paste0("! -> hgt    -  Reference height                 [    m AGL] !")
              ,paste0("! -> tmp    -  Air temperature                  [        K] !")
              ,paste0("! -> pres   -  Pressure                         [       Pa] !")
              ,paste0("! -> sh     -  Specific humidity                [    kg/kg] !")
              ,paste0("! -> ugrd   -  Zonal wind                       [      m/s] !")
              ,paste0("! -> vgrd   -  Zonal wind                       [      m/s] !")
              ,paste0("! -> prate  -  Precipitation rate               [  kg/m2/s] !")
              ,paste0("! -> dlwrf  -  Downward long wave radiation     [     W/m2] !")
              ,paste0("! -> nbdsf  -  Near-IR beam radiation           [     W/m2] !")
              ,paste0("! -> nddsf  -  Near-IR diffuse radiation        [     W/m2] !")
              ,paste0("! -> vbdsf  -  Visible beam radiation           [     W/m2] !")
              ,paste0("! -> vddsf  -  Visible beam radiation           [     W/m2] !")
              ,paste0("! -> co2**  -  CO2 mixing ratio                 [ umol/mol] !")
              ,paste0("!-----------------------------------------------------------!")
              ,paste0("!  *  Land may be absent, in which case all grid cells are  !")
              ,paste0("!     considered land cells for the met driver and may be   !")
              ,paste0("!     selected to drive the simulations.  Also check vari-  !")
              ,paste0("!     able met_land_min in ed_params.f90.                   !")
              ,paste0("!  ** CO2 may be absent, in which case it must be specified !")
              ,paste0("!     in ED2IN (NL%INITIAL_CO2).                            !")
              ,paste0("!===========================================================!")
              ,paste0("!                                                           !")
              ,paste0("!===========================================================!")
              ,paste0("! Table 2. Variable types recognized by ED.                 !")
              ,paste0("!===========================================================!")
              ,paste0("!                                                           !")
              ,paste0("! 0. Read gridded data - no time interpolation;             !")
              ,paste0("! 1. Read gridded data - with time interpolatation;         !")
              ,paste0("! 2. Read gridded data that is constant in time.            !")
              ,paste0("!    If any of this is lon or lat, then deltalon, deltalat  !")
              ,paste0("!    lon0, and lat0 will be ignored;                        !")
              ,paste0("! 3. Read one value representing the whole grid, no time    !")
              ,paste0("!   interpolation;                                          !")
              ,paste0("! 4. Specify a constant for all polygons, constant in time. !")
              ,paste0("!    In this case, give the constant value at line 'e'      !")
              ,paste0("!    instead of the frequency.                              !")
              ,paste0("! 5. Specify a constant for all polygons, with time         !")
              ,paste0("!    interpolation (normally used for CO2).                 !")
              ,paste0("!===========================================================!")
              )#end c
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Write the header.                                                                 #
   #---------------------------------------------------------------------------------------#
   write (x=info,file=header,ncolumns=1,append=FALSE,sep=" ")
   #---------------------------------------------------------------------------------------#
   #=======================================================================================#
   #=======================================================================================#

   rm(datum  )
   rm(dummy  )
}#end for
#------------------------------------------------------------------------------------------#
