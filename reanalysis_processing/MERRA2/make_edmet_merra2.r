#----- Leave this as the first command, this will reset your R session. -------------------#
rm(list=ls())
graphics.off()
gc()
options(warn=0)
#------------------------------------------------------------------------------------------#



#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#        THE FOLLOWING BLOCK ALLOWS YOU TO CONTROL MOST SETTINGS FOR GAP FILLING.          #
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#

#----- Define some paths. -----------------------------------------------------------------#
here    = getwd()                                       # Current path 
srcdir  = file.path(here,"Rsc")                         #
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#       Plot options.                                                                      #
#------------------------------------------------------------------------------------------#
outform         = c("png")        # Formats for output file.  Supported formats are:
                                  #   - "X11"    - for printing on Linux screen
                                  #   - "quartz" - for printing on Mac screen
                                  #   - "eps"    - for postscript printing
                                  #   - "tif"    - for TIFF printing
                                  #   - "png"    - for PNG printing
                                  #   - "pdf"    - for PDF printing
depth           = 200             # PNG/TIFF resolution, in pixels per inch
paper           = "square"        # Paper size, to define the plot shape
ptsz            = 16              # Font size.
n.colourbar     = 100             # Number of colours for colour palette
f.key           = 2/15             # Fraction of plot size for key.
miss.colour     = "grey36"        # Colour for missing data 
#------------------------------------------------------------------------------------------#




#------------------------------------------------------------------------------------------#
#     ED-2.1 control.                                                                      #
#------------------------------------------------------------------------------------------#
edprefix = "MERRA2_ATLFORRJ"                    # ED output file prefix
edroot   = file.path(here,edprefix)           # ED-2.1 output path
outroot  = file.path(here,"figures",edprefix) # Output directory
#------------------------------------------------------------------------------------------#




#------------------------------------------------------------------------------------------#
#     Input paths.  These correspond to the names of the MERRA2 subsets of variables.      #
# The variables needed for ED2 are scattered across four subsets.  Check varread below to  #
# identify the variables needed and cross them with the following directories to identify  #
# which variables to download from each subset.                                            #
#------------------------------------------------------------------------------------------#
m2paths   = c("tavg1_2d_lnd_Nx","tavg1_2d_flx_Nx","tavg1_2d_rad_Nx","tavg1_2d_slv_Nx")
m2paths   = file.path(here,m2paths)
#------------------------------------------------------------------------------------------#


#----- General settings. ------------------------------------------------------------------#
montha   =    1           # First month to use (set it to 1 except if debugging)
yeara    = 1980           # First year to use
monthz   =   12           # Last month to use (set it to 12 except if debugging)
yearz    = 2017           # Last year to use
outhours = sequence(24)-1 # Input/Output hours
height   = 60.            # Reference height
wlon     =  -45.          # Longitude of the westernmost point
elon     =  -40.          # Longitude of the easternmost point
slat     =  -24.          # Latitude of the southernmost point
nlat     =  -20.          # Latitude of the northernmost point
land.min = 0.25           # Minimum fraction to be considered land (dummy for MERRA-2)
#------------------------------------------------------------------------------------------#


#----- Magnitude above which data can be discarded. ---------------------------------------#
mag.max   = 1.e9
undef.out = -9999.99
#------------------------------------------------------------------------------------------#


#----- Flags to overwrite files. ----------------------------------------------------------#
overwrite.hdf5   = c(FALSE,TRUE)[1]
overwrite.header = c(FALSE,TRUE)[1]
#------------------------------------------------------------------------------------------#


#------- Make time table. -----------------------------------------------------------------#
varread  = list( list(vnam="NIRDF"      ,ved2="nddsf" ,desc="Incoming NIR (diffuse)"     )
               , list(vnam="NIRDR"      ,ved2="nbdsf" ,desc="Incoming NIR (direct)"      )
               , list(vnam="PARDFLAND"  ,ved2="vddsf" ,desc="Incoming PAR (diffuse)"     )
               , list(vnam="PARDRLAND"  ,ved2="vbdsf" ,desc="Incoming PAR (direct)"      )
               , list(vnam="PRECTOTCORR",ved2="prate" ,desc="Precipitation rate"         )
               , list(vnam="EMIS"       ,ved2="emis"  ,desc="Emissivity"                 )
               , list(vnam="LWGAB"      ,ved2="dlwra" ,desc="LW absorbed by surface"     )
               , list(vnam="LWGEM"      ,ved2="dlwre" ,desc="LW emitted by surface"      )
               , list(vnam="LWGNT"      ,ved2="dlwrn" ,desc="Net Incoming LW"            )
               , list(vnam="PS"         ,ved2="pres"  ,desc="Surface pressure"           )
               , list(vnam="QV10M"      ,ved2="sh"    ,desc="Specific humidity"          )
               , list(vnam="T10M"       ,ved2="tmp"   ,desc="Temperature"                )
               , list(vnam="U50M"       ,ved2="ugrd"  ,desc="Zonal wind"                 )
               , list(vnam="V50M"       ,ved2="vgrd"  ,desc="Meridional wind"            )
               )#end list
#------------------------------------------------------------------------------------------#



#----- Variables for which we plot monthly averages. --------------------------------------#
n            = 0
outvars      = list()
n            = n + 1
outvars[[n]] = list( vnam = "pres"
                   , desc = "Surface pressure"
                   , unit = "hpa"
                   , dset = "merra2"
                   , add  = "0."
                   , mult = "0.01"
                   , aggr = "mean"
                   , vlwr = 550
                   , vupr = 1050
                   , csch = "iylgnbu"
                   )#end list
n            = n + 1
outvars[[n]] = list( vnam = "sh"
                   , desc = "Specific humidity"
                   , unit = "gwokg"
                   , dset = "merra2"
                   , add  = "0."
                   , mult = "1000."
                   , aggr = "mean"
                   , vlwr = 0
                   , vupr = 25.
                   , csch = "pubugn"
                   )#end list
n            = n + 1
outvars[[n]] = list( vnam = "tmp"
                   , desc = "Air Temperature"
                   , unit = "degC"
                   , dset = "merra2"
                   , add  = "-t00"
                   , mult = "1."
                   , aggr = "mean"
                   , vlwr = -10.
                   , vupr =  35.
                   , csch = "magma"
                   )#end list
n            = n + 1
outvars[[n]] = list( vnam = "vels"
                   , desc = "Wind speed"
                   , unit = "mos"
                   , dset = "merra2"
                   , add  = "0."
                   , mult = "1."
                   , aggr = "mean"
                   , vlwr = 0.
                   , vupr = 10.
                   , csch = "iylgnbu"
                   )#end list
n            = n + 1
outvars[[n]] = list( vnam = "rh"
                   , desc = "Relative humidity"
                   , unit = "pc"
                   , dset = "merra2"
                   , add  = "0."
                   , mult = "100."
                   , aggr = "mean"
                   , vlwr = 0.
                   , vupr = 100.
                   , csch = "ylgnbu"
                   )#end list
n            = n + 1
outvars[[n]] = list( vnam = "dlwrf"
                   , desc = "Incident LW radiation"
                   , unit = "wom2"
                   , dset = "merra2"
                   , add  = "0."
                   , mult = "1."
                   , aggr = "mean"
                   , vlwr = 150
                   , vupr = 550
                   , csch = "ylorrd"
                   )#end list
n            = n + 1
outvars[[n]] = list( vnam = "nbdsf"
                   , desc = "Incoming NIR (Direct)"
                   , unit = "wom2"
                   , dset = "merra2"
                   , add  = "0."
                   , mult = "1."
                   , aggr = "mean"
                   , vlwr = 0.
                   , vupr = 200.
                   , csch = "ylorrd"
                   )#end list
n            = n + 1
outvars[[n]] = list( vnam = "nddsf"
                   , desc = "Incoming NIR (Diffuse)"
                   , unit = "wom2"
                   , dset = "merra2"
                   , add  = "0."
                   , mult = "1."
                   , aggr = "mean"
                   , vlwr = 0.
                   , vupr = 200.
                   , csch = "ylorrd"
                   )#end list
n            = n + 1
outvars[[n]] = list( vnam = "vbdsf"
                   , desc = "Incoming PAR (Direct)"
                   , unit = "wom2"
                   , dset = "merra2"
                   , add  = "0."
                   , mult = "1."
                   , aggr = "mean"
                   , vlwr = 0.
                   , vupr = 200.
                   , csch = "ylorrd"
                   )#end list
n            = n + 1
outvars[[n]] = list( vnam = "vddsf"
                   , desc = "Incoming PAR (Diffuse)"
                   , unit = "wom2"
                   , dset = "merra2"
                   , add  = "0."
                   , mult = "1."
                   , aggr = "mean"
                   , vlwr = 0.
                   , vupr = 200.
                   , csch = "ylorrd"
                   )#end list
n            = n + 1
outvars[[n]] = list( vnam = "prate"
                   , desc = "Precipitation"
                   , unit = "mmomo"
                   , dset = "merra2"
                   , add  = "0."
                   , mult = "hr.sec"
                   , aggr = "sum"
                   , vlwr = 0
                   , vupr = 500
                   , csch = "bupu"
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
#    Load packages and functions.                                                          #
#------------------------------------------------------------------------------------------#
srcdir = (srcdir[file.exists(srcdir)])[1]
source(file.path(srcdir,"load.everything.r"))
#------------------------------------------------------------------------------------------#



#------- Make time table. -----------------------------------------------------------------#
years  = seq(from=yeara,to=yearz,by=1)
nyears = length(years)
when   = data.frame( year             = rep(years,each=12)
                   , month            = rep(sequence(12),times=nyears)
                   , mlab             = rep(toupper(month.abb),times=nyears)
                   , stringsAsFactors = FALSE
                   )#end data.frame
del    = ( (when$year == yeara & when$month < montha )
         | (when$year == yearz & when$month > monthz )
         )#end del
when   = when[! del,,drop=FALSE]
nwhen  = nrow(when)
#------------------------------------------------------------------------------------------#


#------- Convert variable list to data frame. ---------------------------------------------#
varread  = list.2.data.frame(varread)
nvarread = nrow(varread)
#------------------------------------------------------------------------------------------#


#------- Number of paths to search for data sets. -----------------------------------------#
nm2paths = length(m2paths)
#------------------------------------------------------------------------------------------#



#----- Make data frame. -------------------------------------------------------------------#
outvars  = list.2.data.frame(outvars)
noutvars = nrow(outvars)
#------------------------------------------------------------------------------------------#


#----- Make data frame with scenario settings. --------------------------------------------#
noutform  = length(outform)
#------------------------------------------------------------------------------------------#


#----- Number of hours per day, and default time step. ------------------------------------#
nouthours = length(outhours)
dtout     = day.sec / nouthours
#------------------------------------------------------------------------------------------#


#------ Create output paths in case they doen't exist. ------------------------------------#
dummy = dir.create(outroot,showWarnings=FALSE,recursive=TRUE)
dummy = dir.create(edroot ,showWarnings=FALSE,recursive=TRUE)
#------------------------------------------------------------------------------------------#




#==========================================================================================#
#==========================================================================================#
#     Make ED output.                                                                      #
#------------------------------------------------------------------------------------------#
cat0(" + Loop through all times, and make ED-friendly output files.")


#------------------------------------------------------------------------------------------#
#     Read in the data.
#------------------------------------------------------------------------------------------#
any.written = FALSE
for (ww in sequence(nwhen)){
   #---- Aliases for useful variables. ----------------------------------------------------#
   year.now   = when$year [ww]
   month.now  = when$month[ww]
   mlab.now   = when$mlab [ww]
   show.now   = paste (capwords(tolower(mlab.now)),sprintf("%4.4i",year.now))
   suffix.in  = paste0(sprintf("%4.4i",year.now),sprintf("%2.2i",month.now))
   suffix.out = paste0(sprintf("%4.4i",year.now),mlab.now)
   suffix.fig = paste0(sprintf("%4.4i",year.now),"-",sprintf("%2.2i",month.now))
   h5.merra2  = file.path(edroot,paste0(edprefix,"_",suffix.out,".h5"))
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #     Check whether to process time.  In case you want to reprocess files, delete them  #
   # before running R or set overwrite.hdf5 to TRUE.                                       #
   #---------------------------------------------------------------------------------------#
   if (file.exists(h5.merra2) && (! overwrite.hdf5)){
      cat0(" + HDF5 for ",show.now," has already been created.  Skip this time.")
   }else{
      cat0(" + Process data for ",show.now,".")
      #------------------------------------------------------------------------------------#


      #----- List all possible times. -----------------------------------------------------#
      ndays     = daymax(month=month.now,year=year.now)
      ntimes    = nouthours * ndays
      dates.now = paste(month.now,rep(sequence(ndays),each=nouthours),year.now,sep="/")
      times.now = paste(rep(outhours,times=ndays),0,0,sep=":")
      when.now  = chron(dates=dates.now,times=times.now)
      nwhen.now = length(when.now)
      dates.io  = paste(month.now,sequence(ndays),year.now,sep="/")
      when.io   = chron(dates=dates.io,times="00:00:00")
      nwhen.io  = length(when.io)
      dtdat     = mean(diff(outhours)) * hr.sec
      #------------------------------------------------------------------------------------#


      first = TRUE
      #------------------------------------------------------------------------------------#
      #     Loop through times.                                                            #
      #------------------------------------------------------------------------------------#
      for (wio in sequence(nwhen.io)){
         yyyy = sprintf("%4.4i",numyears (when.io[wio]))
         mm   = sprintf("%2.2i",nummonths(when.io[wio]))
         dd   = sprintf("%2.2i",numdays  (when.io[wio]))
         hh   = sprintf("%2.2i",hours    (when.io[wio]))
         wlabel.now = paste0(yyyy,mm,dd)
         wna        = match(when.io[wio],when.now)
         wn         = seq(from=wna,length.out=nouthours,by=1)

         #---------------------------------------------------------------------------------#
         #     Check every path, as the script must read multiple files.                   #
         #---------------------------------------------------------------------------------#
         for (mp in sequence(nm2paths)){
            m2path      = m2paths[mp]
            nc.file     = list.files( path       = m2path
                                    , pattern    = paste0(wlabel.now,"\\.SUB\\.nc$")
                                    , full.names = TRUE
                                    )#end list.files
            nc.file.bz2 = list.files( path       = m2path
                                    , pattern    = paste0(wlabel.now,"\\.SUB\\.nc\\.bz2$")
                                    , full.names = TRUE
                                    )#end list.files
            nc.file.gz  = list.files( path       = m2path
                                    , pattern    = paste0(wlabel.now,"\\.SUB\\.nc\\.gz$")
                                    , full.names = TRUE
                                    )#end list.files
            #------------------------------------------------------------------------------#


            #------------------------------------------------------------------------------#
            #     Use temporary file to avoid losses and to avoid uncompressing all files. #
            #------------------------------------------------------------------------------#
            if (length(nc.file) > 0){
               nc.file     = nc.file[1]
               nc.base     = basename(nc.file)
               temp.file   = file.path(tempdir(),nc.base)
               dummy       = file.copy(from=nc.file,to=temp.file,overwrite=TRUE)
            }else if (length(nc.file.bz2) > 0){
               nc.file.bz2 = nc.file.bz2[1]
               nc.base     = basename(nc.file.bz2)
               nc.base     = gsub(pattern="\\.bz2$",replacement="",x=nc.base)
               dummy       = bunzip2( filename  = nc.file.bz2
                                    , destname  = temp.file
                                    , overwrite = TRUE
                                    , remove    = FALSE
                                    )#end bunzip2
            }else if (length(nc.file.gz) > 0){
               nc.file.gz  = nc.file.gz[1]
               nc.base     = basename(nc.file.gz)
               nc.base     = gsub(pattern="\\.gz$",replacement="",x=nc.base)
               dummy       = gunzip( filename  = nc.file.gz
                                   , destname  = temp.file
                                   , overwrite = TRUE
                                   , remove    = FALSE
                                   )#end gunzip
            }else{
               stop(" Neither the expanded nor the compressed files were found!")
            }#end if (length(nc.file) > 0)
            #------------------------------------------------------------------------------#


            #------------------------------------------------------------------------------#
            #       Retrieve information.                                                  #
            #------------------------------------------------------------------------------#
            if (mp %in% 1) cat0("   - Read file ",nc.base,".")
            nc.conn  = nc_open(temp.file)
            nc.nvars = nc.conn$nvars
            nc.vlist = rep(NA_character_,times=nc.nvars)
            for (v in sequence(nc.nvars)) nc.vlist[v] = nc.conn$var[[v]]$name
            #------------------------------------------------------------------------------#


            #------------------------------------------------------------------------------#
            #     Retrieve variables.                                                      #
            #------------------------------------------------------------------------------#
            for (v in which(varread$vnam %in% nc.vlist)){
               #----- Aliases for current variable. ---------------------------------------#
               v.vnam = varread$vnam[v]
               v.ved2 = varread$ved2[v]
               v.desc = varread$desc[v]
               v.curr = nc.vlist[which(grepl(pattern=v.vnam,x=nc.vlist,ignore.case=TRUE))]
               #---------------------------------------------------------------------------#


               #----- Extract data. -------------------------------------------------------#
               nc.dat  = ncvar_get(nc.conn,v.curr)
               #---------------------------------------------------------------------------#

               #---------------------------------------------------------------------------#
               #     In case this is the first step, make the grid array.                  #
               #---------------------------------------------------------------------------#
               if (first){
                  #------ Find out global grid size, and determine xy resolution. ---------#
                  nx.glob  = dim(nc.dat)[1]
                  ny.glob  = dim(nc.dat)[2]
                  nt.glob  = dim(nc.dat)[3]
                  #------------------------------------------------------------------------#


                  #----- Correct longitude and latitude grid (180W-180E;90N;90S). ---------#
                  lon1d  = ifelse( test = nc.conn$dim$lon$vals >= 180.
                                 , yes  = nc.conn$dim$lon$vals - 360.
                                 , no   = nc.conn$dim$lon$vals
                                 )#end ifelse
                  lat1d  = nc.conn$dim$lat$vals
                  olon   = order(lon1d,decreasing=FALSE)
                  olat   = order(lat1d,decreasing=TRUE )
                  lon1d  = lon1d[olon]
                  lat1d  = lat1d[olat]
                  #------------------------------------------------------------------------#





                  #----- Build the 2-D global grid meshs. ---------------------------------#
                  lon2d  = matrix(rep(lon1d,times=ny.glob),nrow=nx.glob,ncol=ny.glob)
                  lat2d  = matrix(rep(lat1d,each =nx.glob),nrow=nx.glob,ncol=ny.glob)
                  #------------------------------------------------------------------------#


                  #----- Indices for sub-domain (latitude goes from north to south). ------#
                  xykeep = arrayInd( ind  = which( lon2d %wr% c(wlon,elon) 
                                                 & lat2d %wr% c(slat,nlat)
                                                 )#end which
                                   , .dim = dim(lon2d)
                                   )#end arrayInd
                  xkeep  = seq( from = min(xykeep[,1],na.rm=TRUE)
                              , to   = max(xykeep[,1],na.rm=TRUE)
                              , by   = 1
                              )#end seq
                  ykeep  = seq( from = min(xykeep[,2],na.rm=TRUE)
                              , to   = max(xykeep[,2],na.rm=TRUE)
                              , by   = 1
                              )#end seq
                  #------------------------------------------------------------------------#


                  #----- Crop longitude and latitude. -------------------------------------#
                  lon1d  = lon1d[xkeep]
                  lat1d  = lat1d[ykeep]
                  lon2d  = lon2d[xkeep,ykeep]
                  lat2d  = lat2d[xkeep,ykeep]
                  #------------------------------------------------------------------------#


                  #----- Find sub-domain SW corner and size. ------------------------------#
                  nx     = length(lon1d)
                  ny     = length(lat1d)
                  xsw    = min(lon1d)
                  ysw    = min(lat1d)
                  dx     = mean(diff(sort(lon1d)))
                  dy     = mean(diff(sort(lat1d)))
                  #------------------------------------------------------------------------#


                  #------------------------------------------------------------------------#
                  #     Initialise land mask (assume that all points can be used, then     #
                  # discard those points that have undefined numbers.                      #
                  #------------------------------------------------------------------------#
                  land2d = 1 + 0 * lon2d
                  #------------------------------------------------------------------------#


                  #----- Initialise array structure with all data. ------------------------#
                  template      = array(data=NA_real_,dim=c(nx,ny,nwhen.now))
                  merra2        = replicate(nvarread,list(template))
                  names(merra2) = varread$ved2
                  #------------------------------------------------------------------------#
               }#end if (first)
               #---------------------------------------------------------------------------#


               #----- Correct longitude and latitude grid (180W-180E;90N;90S). ------------#
               nc.dat = nc.dat[olon,olat,,drop=FALSE]
               #---------------------------------------------------------------------------#


               #----- Crop domain. --------------------------------------------------------#
               nc.dat                 = nc.dat[xkeep,ykeep,,drop=FALSE]
               discard                = ! (abs(nc.dat) %<% mag.max)
               nc.dat[discard]        = NA_real_
               merra2[[v.ved2]][,,wn] = nc.dat
               #---------------------------------------------------------------------------#


               #----- Find points that have invalid data. ---------------------------------#
               land.now = apply(X=nc.dat,MARGIN=c(1,2),FUN=function(x) all(is.finite(x)))
               if (first){
                  #----- First time, make the land mask. ----------------------------------#
                  first  = FALSE
                  land2d = ifelse(test=land.now,yes=land2d,no=0.)
                  #------------------------------------------------------------------------#
               }else{
                  #----- Other times, make sure that the land mask is consistent. ---------#
                  flood  = any((land2d > 0.) & (! land.now))
                  if (any(flood)) stop ("Undefined points over land.")
                  #------------------------------------------------------------------------#
               }#end if (first)
               #---------------------------------------------------------------------------#
               rm(nc.dat)
               gc()
               #---------------------------------------------------------------------------#
            }#end for (v in sequence(nvarread))
            #------------------------------------------------------------------------------#


            #----- Close and purge the file. ----------------------------------------------#
            dummy   = nc_close(nc.conn)
            dummy   = file.remove(temp.file)
            #------------------------------------------------------------------------------#
         }#end for (mp in sequence(nm2paths)){
         #---------------------------------------------------------------------------------#
      }#end for (wn in sequence(nwhen.now)){
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #     Discard data over water.                                                       #
      #------------------------------------------------------------------------------------#
      cat0("   - Discard data over water.")
      mask3d = array( data = ifelse(test=land2d > land.min,yes=1,no=NA)
                    , dim  = dim(merra2[[varread$ved2[1]]])
                    )#end array
      for (v in sequence(nvarread)){
         #----- Aliases for current variable. ---------------------------------------------#
         v.vnam = varread$vnam[v]
         v.ved2 = varread$ved2[v]
         v.desc = varread$desc[v]
         merra2[[v.ved2]] = merra2[[v.ved2]] * mask3d
         #---------------------------------------------------------------------------------#
      }#end for (v in sequence(nvarread))

      #------------------------------------------------------------------------------------#
      #     Make sure air is never too dry or too wet.                                     #
      #------------------------------------------------------------------------------------#
      cat0("   - Make sure air is not bone-dry or supersaturated.")
      merra2$ssh   = with(merra2,qslif (pres,tmp   )) + 0. * merra2$tmp
      merra2$shmin = 0.01  * merra2$ssh
      merra2$shmax = 0.995 * merra2$ssh
      merra2$sh    = pmax(merra2$shmin,pmin(merra2$shmax,merra2$sh)) + 0. * merra2$sh
      merra2$rh    = with(merra2,rehuil(pres,tmp,sh)) + 0. * merra2$sh
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #     Find incident longwave radiation.                                              #
      #------------------------------------------------------------------------------------#
      cat0("   - Find incident longwave radiation.")
      merra2$dlwrf = merra2$dlwra / merra2$emis
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #    Find wind magnitude (plotting only).                                            #
      #------------------------------------------------------------------------------------#
      cat0("   - Make wind magnitude.")
      merra2$vels = with(merra2,sqrt(ugrd * ugrd + vgrd * vgrd))
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #     Make the met driver file.                                                      #
      #------------------------------------------------------------------------------------#
      #----- 1st file, all drivers except for rainfall. -----------------------------------#
      cat0("   - Replace NA with missing flags.")
      lon   = aperm(a=lon2d ,perm=c(  2,1))
      lat   = aperm(a=lat2d ,perm=c(  2,1))
      land  = aperm(a=land2d,perm=c(  2,1))
      dlwrf = aperm(a=ifelse(is.finite(merra2$dlwrf),merra2$dlwrf,undef.out),perm=c(3,2,1))
      nbdsf = aperm(a=ifelse(is.finite(merra2$nbdsf),merra2$nbdsf,undef.out),perm=c(3,2,1))
      nddsf = aperm(a=ifelse(is.finite(merra2$nddsf),merra2$nddsf,undef.out),perm=c(3,2,1))
      prate = aperm(a=ifelse(is.finite(merra2$prate),merra2$prate,undef.out),perm=c(3,2,1))
      pres  = aperm(a=ifelse(is.finite(merra2$pres ),merra2$pres ,undef.out),perm=c(3,2,1))
      sh    = aperm(a=ifelse(is.finite(merra2$sh   ),merra2$sh   ,undef.out),perm=c(3,2,1))
      tmp   = aperm(a=ifelse(is.finite(merra2$tmp  ),merra2$tmp  ,undef.out),perm=c(3,2,1))
      ugrd  = aperm(a=ifelse(is.finite(merra2$ugrd ),merra2$ugrd ,undef.out),perm=c(3,2,1))
      vbdsf = aperm(a=ifelse(is.finite(merra2$vbdsf),merra2$vbdsf,undef.out),perm=c(3,2,1))
      vddsf = aperm(a=ifelse(is.finite(merra2$vddsf),merra2$vddsf,undef.out),perm=c(3,2,1))
      vgrd  = aperm(a=ifelse(is.finite(merra2$vgrd ),merra2$vgrd ,undef.out),perm=c(3,2,1))
      cat0("   - Save 6-hourly met drivers to ",basename(h5.merra2),".")
      if (file.exists(h5.merra2)) dummy   = file.remove(h5.merra2)
      dummy = h5save(lon,lat,land,dlwrf,nbdsf,nddsf,prate,pres,sh,tmp,ugrd,vbdsf,vddsf,vgrd
                    ,file=h5.merra2)
      dummy = H5close()
      #------------------------------------------------------------------------------------#


      #----- Update any.written so we know if it's possible to write a header file. -------#
      any.written = TRUE
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #     Find dimensions for plot.                                                      #
      #------------------------------------------------------------------------------------#
      limlon  = range(lon1d,finite=TRUE)
      limlat  = range(lat1d,finite=TRUE)
      f.ext   = f.key  / (1. - f.key )
      ex.size = plotsize( proje     = TRUE
                        , limlon    = limlon
                        , limlat    = limlat
                        , extendfc  = "lon"
                        , extfactor = f.ext
                        , paper     = paper
                        )#end plotsize
      #------------------------------------------------------------------------------------#



      #------ Labels for the longitude and latitude. --------------------------------------#
      lonplot = pretty.lonlat(x=limlon,n=6,type="lon")
      latplot = pretty.lonlat(x=limlat,n=6,type="lat")
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #      Go through variables, and make maps of the monthly mean.                      #
      #------------------------------------------------------------------------------------#
      for (v in sequence(noutvars)){
         #---- Handy aliases. -------------------------------------------------------------#
         v.vnam = outvars$vnam[v]
         v.desc = outvars$desc[v]
         v.unit = outvars$unit[v]
         v.dset = outvars$dset[v]
         v.add  = outvars$add [v]
         v.mult = outvars$mult[v]
         v.aggr = outvars$aggr[v]
         v.vlwr = outvars$vlwr[v]
         v.vupr = outvars$vupr[v]
         v.csch = outvars$csch[v]
         v.show = ! (v.vnam %in% "pr.gpcc")
         v.off  = sqrt(.Machine$double.eps) * (v.vupr-v.vlwr)
         #---------------------------------------------------------------------------------#

         #---------------------------------------------------------------------------------#
         #    Make sure variable is supposed to be shown.                                  #
         #---------------------------------------------------------------------------------#
         if (v.show){
            #----- Path for this variable. ------------------------------------------------#
            cat0("   - Plot monthly mean of ",v.desc,".")
            outvar = file.path(outroot,v.vnam)
            if (! file.exists(outvar)) dir.create(outvar)
            #------------------------------------------------------------------------------#


            #----- Retrieve current variable. ---------------------------------------------#
            arrlab = paste0(v.add," + ",v.mult," * ",v.dset,"$",v.vnam)
            arrnow = eval(parse(text=arrlab))
            mmean  = apply(X=arrnow,MARGIN=c(1,2),FUN=match.fun(v.aggr),na.rm=TRUE)
            mmean  = 0*mmean + pmax(v.vlwr+v.off,pmin(v.vupr-v.off,mmean))
            mmean  = ifelse(land2d %>=% 0.5,mmean,NA)
            #------------------------------------------------------------------------------#

            #----- Title. -----------------------------------------------------------------#
            v.main = desc.unit(desc=paste0(show.now," - ",v.desc),unit=untab[[v.unit]])
            #------------------------------------------------------------------------------#

            #----- Colour ramp. -----------------------------------------------------------#
            v.lim     = c(v.vlwr,v.vupr)
            v.at      = pretty(v.lim)
            v.labels  = sprintf("%g",v.at)
            v.levels  = pretty(v.lim,n=n.colourbar)
            v.levels  = c( min(v.levels) - mean(diff(v.levels))
                            , v.levels
                            , max(v.levels) + mean(diff(v.levels))
                            )#end c
            v.nlevels = length(v.levels)
            #------------------------------------------------------------------------------#



            #----- Loop over formats. -----------------------------------------------------#
            for (o in sequence(noutform)){
               #------ Open file. ---------------------------------------------------------#
               fichier = file.path( outvar
                                  , paste0("mmean_",v.vnam,"_",suffix.fig,".",outform[o])
                                  )#end file.path
               dummy   = open.plot( fichier = fichier
                                  , outform = outform[o]
                                  , size    = ex.size
                                  , depth   = depth
                                  , ptsz    = ptsz
                                  )#end open.plot
               #---------------------------------------------------------------------------#


               #----- Map. ----------------------------------------------------------------#
               image.map( x                = lon2d
                        , y                = lat2d
                        , z                = mmean
                        , colour.palette   = match.fun(v.csch)
                        , levels           = v.levels
                        , nlevels          = v.nlevels
                        , na.col           = miss.colour
                        , xlim             = limlon
                        , ylim             = limlat
                        , x.axis.options   = with( lonplot
                                                 , list(side=1,las=1,at=at,labels=labels)
                                                 )#end with
                        , y.axis.options   = with( latplot
                                                 , list(side=2,las=1,at=at,labels=labels)
                                                 )#end with
                        , main.title       = list(text=v.main,font=2,cex.main=0.8*cex.main)
                        , key.log          = FALSE
                        , key.axis.options = list(side=4,las=1,at=v.at,labels=v.labels)
                        , plot.after       = list(map=list(database="worldHires",add=TRUE))
                        , f.key            = f.key
                        , mar.main         = c(3.1,4.1,1.6,0.3)
                        , mar.key          = c(3.1,0.3,1.6,3.1)
                        )#end image.map
               #---------------------------------------------------------------------------#



               #----- Close the device. ---------------------------------------------------#
               dummy = close.plot(outform=outform[o])
               #---------------------------------------------------------------------------#
            }#end for (o in sequence(noutform))     
            #------------------------------------------------------------------------------#


            #----- Free memory. -----------------------------------------------------------#
            rm(arrlab,mmean)
            #------------------------------------------------------------------------------#
         }#end if (v.show)
         #---------------------------------------------------------------------------------#
      }#end for (v in sequence(noutvars))
      #------------------------------------------------------------------------------------#


      #---- Free memory. ------------------------------------------------------------------#
      rm(merra2,lon,lat,land,dlwrf,nbdsf,nddsf,prate,pres,sh,tmp,ugrd,vbdsf,vddsf,vgrd)
      gc()
      #------------------------------------------------------------------------------------#
   }#end if (file.exists(h5.merra2) && (! overwrite.hdf5))
   #---------------------------------------------------------------------------------------#
}#end for (w in sequence(nwhen))
#------------------------------------------------------------------------------------------#





#------------------------------------------------------------------------------------------#
#     Make the header for this experiment.                                                 #
#------------------------------------------------------------------------------------------#
header.merra2 = file.path(edroot,paste0(edprefix,"_HEADER"))
if ( any.written && ( (! file.exists(header.merra2)) || overwrite.header ) ){
   cat0(" + Write header for hourly variables (",basename(header.merra2),").")
   info    = c("# See README at the bottom of this file."
              ,"1"
              ,file.path(edroot,paste0(edprefix,"_"))
              ,paste(nx,ny,dx,dy,xsw,ysw)
              ,"15"
              ,paste0("   'hgt'   'lon'   'lat'  'land' 'dlwrf' 'nbdsf' 'nddsf' 'prate'"
                     ,"  'pres'    'sh'   'tmp'  'ugrd' 'vbdsf' 'vddsf'  'vgrd'")
              ,paste (c(height,rep(dtdat,times=14)),collapse=" ")
              ,paste (c(4,rep(2,times= 3),rep(1,times=3),2,rep(1,times=7),collapse=" "))
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
   write (x=info,file=header.merra2,ncolumns=1,append=FALSE,sep=" ")
}else if(! any.written){
   #----- Do not write header.  We can't because we don't have all the information. -------#
   cat0(" File ",basename(header.merra2)," not written.  The script requires")
   cat0("    at least one HDF5 file to be written to generate the header file.")
   #---------------------------------------------------------------------------------------#
}else if(! any.written){
   #----- Do not write header.  We can't because we don't have all the information. -------#
   cat0(" File ",basename(header.merra2)," already exists and was not re-created.")
   #---------------------------------------------------------------------------------------#
}#end if ( any.written && ( (! file.exists(header.merra2)) || overwrite.header ) )
#------------------------------------------------------------------------------------------#
