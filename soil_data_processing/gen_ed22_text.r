#==========================================================================================#
#==========================================================================================#
#     Script to generate HDF5 files of soil texture compatible with the ED2 model, based   #
# on SoilGrids data.                                                                       #
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
stext.path    = file.path(here,"Input")                   # Path with input data
edprefix      = "SoilGrids20"                             # Prefix for output data set
edroot        = file.path(here,edprefix)                  # Path with R objects
#------------------------------------------------------------------------------------------#


#----- Global soil information. -----------------------------------------------------------#
zlayer        = "005-015cm"
zstat         = "mean"
zundef        = -1
#------------------------------------------------------------------------------------------#


#------ Texture files. --------------------------------------------------------------------#
n             = 0
textures      = list( list(vnam="sand",desc="Percent Sand",ved2="sand")
                    , list(vnam="silt",desc="Percent Silt",ved2="silt")
                    , list(vnam="clay",desc="Percent Clay",ved2="clay")
                    )#end list
#------------------------------------------------------------------------------------------#


#------ Define block length. --------------------------------------------------------------#
lona    = -94
lonz    = -30
lata    = -60
latz    =  30
dblock  =   2
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#     In case this is a debug test, set the block for testing that the file put the data   #
# in the right order.  Hint: when possible, select a block along the coast, with well      #
# known geographic features, such as bays, sounds or capes.  The code will plot the        #
# results for this block and they should show the expected shapes assuming north at the    #
# top.                                                                                     #
#------------------------------------------------------------------------------------------#
run.dbg = FALSE  # Run the debugging test (TRUE|FALSE)
dbg.lon = -44    # Longitude of the SW corner of the block used for debugging
dbg.lat = -24    # Latitude of the SW corner of the block used for debugging
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



#------ Site list.  Check whether to exclude thin sites. ----------------------------------#
textures  = list.2.data.frame(textures)
ntextures = nrow(textures)
#------------------------------------------------------------------------------------------#


#------ In case the output directory doesn't exist, create it. ----------------------------#
dummy = dir.create(edroot,recursive=TRUE,showWarnings=FALSE)
#------------------------------------------------------------------------------------------#


#------ Loop through texture files. -------------------------------------------------------#
datum = replicate(ntextures,list())
names(datum) = textures$text
for (tx in sequence(ntextures)){
   #----- Handy aliases. ------------------------------------------------------------------#
   tx.vnam = textures$vnam[tx]
   tx.desc = textures$desc[tx]
   tx.ved2 = textures$ved2[tx]
   tx.base = paste0(tx.vnam,"_",zlayer,"_",zstat,".tif")
   tx.file = file.path(stext.path,tx.base)
   cat0(" + Read ",tx.desc," content from file ",tx.base,".")
   #---------------------------------------------------------------------------------------#


   #----- Load the raster. ----------------------------------------------------------------#
   datum[[tx.ved2]] = raster(tx.file)
   #---------------------------------------------------------------------------------------#
}#end for (tx in sequence(ntextures))
#------------------------------------------------------------------------------------------#




#----- Make sure sand, clay and silt add up to 100%. --------------------------------------#
cat0(" + Make sure sand, silt, and clay contents add up to 100%.")
total      = datum$sand + datum$clay + datum$silt
datum$sand = as.integer(round(1000 * datum$sand / total))
datum$clay = as.integer(round(1000 * datum$clay / total))
datum$silt = 1000 - datum$clay - datum$sand
#------------------------------------------------------------------------------------------#




#----- Find coordinates. ------------------------------------------------------------------#
lonbnds = seq(from=lona,to=lonz,by=dblock)
latbnds = seq(from=lata,to=latz,by=dblock)
lonbnds = lonbnds[lonbnds %wr% c(datum$sand@extent@xmin,datum$sand@extent@xmax)]
latbnds = latbnds[latbnds %wr% c(datum$sand@extent@ymin,datum$sand@extent@ymax)]
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





#------------------------------------------------------------------------------------------#
#      Loop through all variables, we will create multiple data sets.                      #
#------------------------------------------------------------------------------------------#
premier = TRUE
#----- Handy aliases. ---------------------------------------------------------------------#
tx.vnam = "text"
tx.desc = "Texture"
tx.ved2 = "text"
tx.sand = file.path(stext.path,paste0("sand_",zlayer,"_",zstat,".tif"))
tx.clay = file.path(stext.path,paste0("clay_",zlayer,"_",zstat,".tif"))
cat0(" + Generate ED2 files for ",tx.desc,".")
#------------------------------------------------------------------------------------------#


#----- Loop through longitudes. -----------------------------------------------------------#
for (x in xloop){
   #----- Current longitude bounds. -------------------------------------------------------#
   wnow   = west[x]
   enow   = east[x]
   lonlab = paste0(sprintf("%3.3i",abs(wnow)),ifelse(test=wnow < 0,yes="W",no="E"))
   #---------------------------------------------------------------------------------------#

   for (y in yloop){
      #----- Current latitude bounds. -----------------------------------------------------#
      snow   = south[y]
      nnow   = north[y]
      latlab = paste0(sprintf("%2.2i",abs(snow)),ifelse(test=snow < 0,yes="S",no="N"))
      cat0("   - Process block ",lonlab,"; ",latlab,".")
      #------------------------------------------------------------------------------------#




      #----- Crop data. -------------------------------------------------------------------#
      cat0("     ~ Crop textures.")
      block.now = extent(wnow,enow,snow,nnow)
      sand.now  = crop(datum$sand,block.now)
      clay.now  = crop(datum$clay,block.now)
      ved2.now  = sand.now
      #------------------------------------------------------------------------------------#


      #------ Combine texture classes. ----------------------------------------------------#
      cat0("     ~ Reshape texture classes.")
      nxnow = sand.now@ncols
      nynow = sand.now@nrows
      yrev  = rev(sequence(nynow))
      sand  = matrix(data=getValues(sand.now),nrow=nxnow,ncol=nynow)
      sand  = sand[,yrev,drop=FALSE]
      clay  = matrix(data=getValues(clay.now),nrow=nxnow,ncol=nynow)
      clay  = clay[,yrev,drop=FALSE]
      ved2  = sand
      #------------------------------------------------------------------------------------#


      #------ Reshape results to a matrix. ------------------------------------------------#
      cat0("     ~ Find texture classes.")
      fao = mapply( FUN      = function(sandfrac,clayfrac){
                                  if (is.finite(sandfrac) & is.finite(clayfrac)){
                                     ans = sclass(sandfrac,clayfrac)
                                  }else{
                                     ans = NA
                                  }#end if
                                  return(ans)
                               }#end function
                  , sandfrac = c(sand/1000.)
                  , clayfrac = c(clay/1000.)
                  )#end mapply
      fao = matrix(data=c(fao),nrow=nxnow,ncol=nynow)
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #     Do not bother to write HDF5 in case everything is missing.                     #
      #------------------------------------------------------------------------------------#
      edbase = paste0(edprefix,"_",latlab,lonlab,".h5")
      if (any(is.finite(ved2))){
         #----- Fill in the fao data set, and convert it to integer. ----------------------#
         cat0("     ~ Fill missing values using nearest valid point.")
         if (! run.dbg) fao = populate(fao)
         fao = matrix(data=as.integer(round(fao)),nrow=nrow(fao),ncol=ncol(fao))
         #---------------------------------------------------------------------------------#



         #----- Create file name. ---------------------------------------------------------#
         cat0("     ~ Write output file (",edbase,").")
         h5.block = file.path(edroot,edbase)
         if (file.exists(h5.block)) dummy   = file.remove(h5.block)
         dummy   = h5save(fao,file=h5.block)
         dummy   = H5close()
         #---------------------------------------------------------------------------------#



         #----- Read and plot HDF5. -------------------------------------------------------#
         if (run.dbg){
            fine = h5read(file=h5.block,name="fao")
            fine = 0L*fine + ifelse(test=fine == zundef,yes=NA_integer_,no=fine)
            gridded.plot(fine,colour.palette=magma)
            locator(n=1)
            file.remove(h5.block)
         }#end if (run.dbg)
         #---------------------------------------------------------------------------------#
      }else{
         #----- Create file name. ---------------------------------------------------------#
         cat0("     ~ Skip file (",edbase,"); this block is entirely over water.")
         #---------------------------------------------------------------------------------#
      }#end if (any(is.finite(fao)))
      #------------------------------------------------------------------------------------#



      #----- Find the header information. -------------------------------------------------#
      if (premier){
         cat0("     ~ Write header file.")
         premier  = FALSE
         lonwalls = seq( from       = ved2.now@extent@xmin
                       , to         = ved2.now@extent@xmax
                       , length.out = ved2.now@ncols + 1
                       )#end seq
         latwalls = seq( from       = ved2.now@extent@ymin
                       , to         = ved2.now@extent@ymax
                       , length.out = ved2.now@nrows + 1
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
            second = "-----------------------------------------------------------------"
            third  = "   nio,   njo,   nperdeg,     iwbego,     isbego,"
            header = c(first,second,third)
            head.file = file.path(edroot,paste0(edprefix,"_HEADER"))
            dummy     = write(x=header,file=head.file,ncolumns=1,append=FALSE)
         }#end if (nperdegx != nperdegy)
         #---------------------------------------------------------------------------------#
      }#end if (premier)
      #------------------------------------------------------------------------------------#

   }#end for (y in sequence(nlat))
   #---------------------------------------------------------------------------------------#
}#end for (x in sequence(nlon))
#------------------------------------------------------------------------------------------#
