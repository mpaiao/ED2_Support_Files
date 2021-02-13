#==========================================================================================#
#==========================================================================================#
#     Script to generate HDF5 files of soil depth compatible with the ED2 model.  This     #
# assumes that input files are in GeoTIFF, and contain soil depth to bedrock in metres.    #
# In this example, the data were derived from the Pelletier et al. (2016,                  #
# doi:10.1002/2015MS000526) data set (by adding "soil/sedimentary deposit thickness" and   #
# "upland hill slope regolith thickness").  But this script should work for any soil       #
# depth data once the input file is properly read and it comes with the expected units.    #
#                                                                                          #
#  Developed by Marcos Longo.                                                              #
#------------------------------------------------------------------------------------------#



#==========================================================================================#
#==========================================================================================#
#     Reset the session.                                                                   #
#------------------------------------------------------------------------------------------#
rm(list=ls())
graphics.off()
options(warn=0)
#==========================================================================================#
#==========================================================================================#




#==========================================================================================#
#==========================================================================================#
#    Main settings.                                                                        #
#------------------------------------------------------------------------------------------#



#----- File and path information. ---------------------------------------------------------#
here          = getwd()                                   # Working directory
depth.path    = file.path(here,"Input")                   # Path with input data
edprefix      = "Pelletier"                               # Prefix for output data set
edroot        = file.path(here,edprefix)                  # Path with R objects
#------------------------------------------------------------------------------------------#

#----- Prefix of the original data. -------------------------------------------------------#
orig.base = "Pelletier_soildepth_010deg.tif"
#------------------------------------------------------------------------------------------#


#------ Define block length. --------------------------------------------------------------#
lona    = -95
lonz    = -30
lata    = -60
latz    =  30
dblock  =   5
#------------------------------------------------------------------------------------------#

#------ Global soil information. ----------------------------------------------------------#
zundef  = 255
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#     In case this is a debug test, set the block for testing that the file put the data   #
# in the right order.  Hint: when possible, select a block along the coast, with well      #
# known geographic features, such as bays, sounds or capes.  The code will plot the        #
# results for this block and they should show the expected shapes assuming north at the    #
# top.                                                                                     #
#------------------------------------------------------------------------------------------#
run.dbg = FALSE  # Run the debugging test (TRUE|FALSE)
dbg.lon = -45    # Longitude of the SW corner of the block used for debugging
dbg.lat = -25    # Latitude of the SW corner of the block used for debugging
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



#----- Load some useful functions. --------------------------------------------------------#
source(file.path(here,"stext_utils.r"))
#------------------------------------------------------------------------------------------#


#------ Get rid of the extremely annoying and unnecessary bell. ---------------------------#
options(locatorBell=FALSE)
#------------------------------------------------------------------------------------------#



#----- Load some useful packages. ---------------------------------------------------------#
loaded.package            = list()
loaded.package[["chron" ]] = require(chron )
loaded.package[["raster"]] = require(raster)
loaded.package[["rhdf5" ]] = require(rhdf5 )
#---- Make sure all packages are loaded fine. ---------------------------------------------#
loaded.package = unlist(loaded.package)
if (! all(loaded.package)){
   miss = which(! loaded.package)
   cat0(" You must install the following packages before using the scripts:")
   for (m in miss) cat0(" -> ",names(loaded.package)[m],".")
   stop("Missing packages!!!")
}#end if
#------------------------------------------------------------------------------------------#


#------ In case the output directory doesn't exist, create it. ----------------------------#
dummy = dir.create(edroot,recursive=TRUE,showWarnings=FALSE)
#------------------------------------------------------------------------------------------#


#------ Read bedrock depth file. =---------------------------------------------------------#
cat0(" + Read bedrock depth content from file ",orig.base,".")
orig.file = file.path(depth.path,orig.base)
datum     = raster(orig.file)
#------------------------------------------------------------------------------------------#



#----- Find coordinates. ------------------------------------------------------------------#
lonbnds = seq(from=lona,to=lonz,by=dblock)
latbnds = seq(from=lata,to=latz,by=dblock)
lonbnds = lonbnds[lonbnds %wr% c(datum@extent@xmin,datum@extent@xmax)]
latbnds = latbnds[latbnds %wr% c(datum@extent@ymin,datum@extent@ymax)]
west    = lonbnds[-length(lonbnds)]
east    = lonbnds[              -1]
south   = latbnds[-length(latbnds)]
north   = latbnds[              -1]
nlon    = length(west)
nlat    = length(south)
#------------------------------------------------------------------------------------------#


#----- Set x and y loop (typically everything, unless this is a debugging run). -----------#
if (run.dbg){
  xloop = which(west  %in% dbg.lon)
  yloop = which(south %in% dbg.lat)
}else{
  xloop = sequence(nlon)
  yloop = sequence(nlat)
}#end if (run.dbg)
#------------------------------------------------------------------------------------------#



#----- Extract data. ----------------------------------------------------------------------#
premier = TRUE
for (y in yloop){
   #----- Current latitude bounds. --------------------------------------------------------#
   snow   = south[y]
   nnow   = north[y]
   latlab = paste0(sprintf("%2.2i",abs(snow)),ifelse(test=snow < 0,yes="S",no="N"))
   #---------------------------------------------------------------------------------------#


   for (x in xloop){
      #----- Current longitude bounds. ----------------------------------------------------#
      wnow   = west[x]
      enow   = east[x]
      lonlab = paste0(sprintf("%3.3i",abs(wnow)),ifelse(test=wnow < 0,yes="W",no="E"))
      cat0(" + Process block ",lonlab,"; ",latlab,".")
      #------------------------------------------------------------------------------------#

      #----- Crop sand and clay. ----------------------------------------------------------#
      cat0("   - Crop depth.")
      block.now = extent(wnow,enow,snow,nnow)
      depth.now = crop(datum,block.now)
      #------------------------------------------------------------------------------------#


      #------ Combine texture classes. ----------------------------------------------------#
      cat0("   - Reshape soil depth.")
      nxnow = depth.now@ncols
      nynow = depth.now@nrows
      yrev  = rev(sequence(nynow))
      depth = matrix(data=getValues(depth.now),nrow=nxnow,ncol=nynow)[,yrev,drop=FALSE]
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #     Do not bother to write HDF5 in case everything is missing.                     #
      #------------------------------------------------------------------------------------#
      if (any(is.finite(depth))){
         #----- Ensure files don't have any NA. -------------------------------------------#
         if (! run.dbg) depth = ifelse(test=is.finite(depth),yes=depth*100L,no=0L)
         #---------------------------------------------------------------------------------#


         #----- Create file name. ---------------------------------------------------------#
         edbase   = paste0(edprefix,"_",latlab,lonlab,".h5")
         cat0("   - Write output file (",edbase,").")
         h5.block = file.path(edroot,edbase)
         if (file.exists(h5.block)) dummy   = file.remove(h5.block)
         dummy   = h5save(depth,file=h5.block)
         dummy   = H5close()
         #---------------------------------------------------------------------------------#



         #----- Read and plot HDF5. -------------------------------------------------------#
         if (run.dbg){
            fine = h5read(file=h5.block,name="depth")
            fine = 0L*fine + ifelse(test=fine == zundef,yes=NA_integer_,no=fine)
            gridded.plot(fine,colour.palette=magma)
            locator(n=1)
            file.remove(h5.block)
         }#end if (run.dbg)
         #---------------------------------------------------------------------------------#
      }else{
         #----- Create file name. ---------------------------------------------------------#
         edbase   = paste0(edprefix,"_",latlab,lonlab,".h5")
         cat0("   - Skip file (",edbase,"); this block is entirely over water.")
         #---------------------------------------------------------------------------------#
      }#end if (any(is.finite(fao)))
      #------------------------------------------------------------------------------------#


      #----- Find the header information. -------------------------------------------------#
      if (premier){
         cat0("     ~ Write header file.")
         premier  = FALSE
         lonwalls = seq( from       = depth.now@extent@xmin
                       , to         = depth.now@extent@xmax
                       , length.out = depth.now@ncols + 1
                       )#end seq
         latwalls = seq( from       = depth.now@extent@ymin
                       , to         = depth.now@extent@ymax
                       , length.out = depth.now@nrows + 1
                       )#end seq
         dlon      = mean(diff(lonwalls))
         dlat      = mean(diff(latwalls))
         nperdegx  = round(1./dlon)
         nperdegy  = round(1./dlat)
         if (nperdegx != nperdegy){
            cat0(" - Longitude resolution: ",dlon,"; # points per degree: ",nperdegx)
            cat0(" - Latitude  resolution: ",dlat,"; # points per degree: ",nperdegy)
            stop(" Grid count must be the same for x and y.")
         }else{
            first = paste0(sprintf("%6i",nxnow   ),",",sprintf("%6i",nynow),","
                          ,sprintf("%6i",nperdegx),","
                          ,sprintf("%6i",lona    ),",",sprintf("%6i",lata    )
                          )#end paste0
            second = "------------------------------------------------------------------"
            third  = "   nio,   njo,   nperdeg,     iwbego,     isbego,"
            header = c(first,second,third)
            head.file = file.path(edroot,paste0(edprefix,"_HEADER"))
            dummy     = write(x=header,file=head.file,ncolumns=1,append=FALSE)
         }#end if (nperdegx != nperdegy)
         #---------------------------------------------------------------------------------#
      }#end if ((x == 1) && (y == 1))
      #------------------------------------------------------------------------------------#
   }#end for (x in xloop)
   #---------------------------------------------------------------------------------------#
}#end for (y in yloop)
#------------------------------------------------------------------------------------------#
