#----- Refresh everything. ----------------------------------------------------------------#
rm(list=ls())
options(warn=0)
gc()
graphics.off()
#------------------------------------------------------------------------------------------#


#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#
#     Define the parameters.                                                               #
#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#

#----- 1. Paths. --------------------------------------------------------------------------#
here    = getwd()                                         # Main path (current path)
inpdir  = file.path(here,"inputs_luh")                    # Path with inputs
histdir = file.path(inpdir,"historical")                  # Input path
#----- Future land use scenarios. ---------------------------------------------------------#
rcp     = c("rcp26_image","rcp45_minicam","rcp60_aim","rcp85_message")[3]
#------------------------------------------------------------------------------------------#

#------------------------------------------------------------------------------------------#
# 2. Define the initial and final coordinates for the output blocks.                       #
#------------------------------------------------------------------------------------------#
lona    = -56.0 # Westernmost edge (first block longitude will be lona + dblock/2)
lonz    = -52.0 # Easternmost edge (the closest that creates regular blocks will be used)
lata    =  -4.0 # Southernmost edge (first block latitude will be lata + dblock/2)
latz    =   0.0 # Northernmost edge (the closest that creates regular blocks will be used)
dblock  =  0.50 # Size of block (both longitude and latitude)
#------------------------------------------------------------------------------------------#

#------------------------------------------------------------------------------------------#
#     The following variables are PFT-dependent, so they can all be vectors.  Current ED   #
# PFTs:                                                                                    #
# 1  - C4 Grass                                                                            #
# 2  - Early successional tropical broadleaf tree (thin bark)                              #
# 3  - Mid successional tropical broadleaf tree (thin bark)                                #
# 4  - Late successional tropical broadleaf tree (thin bark)                               #
# 5  - C3 grass                                                                            #
# 6  - Northern North American pines                                                       #
# 7  - Southern North American pines                                                       #
# 8  - Late successional conifers                                                          #
# 9  - Early successional cold deciduous broadleaf tree                                    #
# 10 - Mid successional cold deciduous broadleaf tree                                      #
# 11 - Late successional cold deciduous broadleaf tree                                     #
# 12 - Early successional tropical broadleaf tree (thick bark)                             #
# 13 - Mid successional tropical broadleaf tree (thick bark)                               #
# 14 - Late successional tropical broadleaf tree (thick bark)                              #
# 15 - Araucaria                                                                           #
# 16 - Tropical C3 grass                                                                   #
# 17 - Liana                                                                               #
#                                                                                          #
# npft.harvest      - Number of PFTs that can be harvested.  If this is set to zero, this  #
#                     means that the information on the bottom shall not be used.          #
# harvest.pft       - Which PFTs are harvested?                                            #
# minimum.dbh.slog  - Minimum DBH to be (selectively) logged in primary forest.            #
#                     This is used only when landuse[12] is set to -1.0                    #
# harvest.prob.slog - Probability of harvesting a tree whose DBH exceeds the minimum DBH.  #
#                     This is used only when landuse[12] is set to -1.0                    #
# minimum.dbh.fplt  - Minimum DBH to be logged in forest plantation.                       #
#                     This is used only when landuse[14] is set to -1.0                    #
# harvest.prob.fplt - Probability of harvesting a tree whose DBH exceeds the minimum DBH.  #
#                     This is used only when landuse[14] is set to -1.0                    #
#------------------------------------------------------------------------------------------#
npft.harvest  = 13
harvest.pft   = c(2,3,4,6,7,8,9,10,11,12,13,14,15)
tropical      = harvest.pft %in% c(2,3,4,12,13,14,15)
mindbh.slog   = ifelse(tropical,50.0, 0.0)
harvprob.slog = ifelse(tropical,0.25,0.25)
mindbh.fplt   = ifelse(tropical,50.0, 0.0)
harvprob.fplt = ifelse(tropical,0.25,0.25)
#------------------------------------------------------------------------------------------#




#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#
#     You should not need to change anything beyond this point unless you are developing   #
# the script...                                                                            #
#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#


#----- Load some scripts. -----------------------------------------------------------------#
source(file.path(here,"solid.angle.r"))
#------------------------------------------------------------------------------------------#


#----- Load some packages. ----------------------------------------------------------------#
#------------------------------------------------------------------------------------------#


#----- Set some constants. ----------------------------------------------------------------#
pio180    <<- pi/ 180.        # Pi/180 (deg -> rad)                             [      ---]
erad      <<- 6370997.        # Earth radius                                    [        m]
#------------------------------------------------------------------------------------------#


#----- Set future land use input, and the output path. ------------------------------------#
rcpdir  = file.path(inpdir,rcp)
extdir  = file.path(inpdir,paste(rcp,"ext",sep="_"))
outpref = paste("luh-1.1",rcp,sep="+")
outmain = file.path(here,outpref)
if (dblock == 1.0){
   outdir = file.path(outmain,"one")
}else if (dblock == 0.5){
   outdir = file.path(outmain,"half")
}else if (dblock == 0.2){
   outdir = file.path(outmain,"fifth")
}else{
   stop ("Invalid dblock")
}#end if
dummy = dir.create(outdir,recursive=TRUE,showWarnings=FALSE)
#------------------------------------------------------------------------------------------#



#----- Find the longitude and latitude of the edges. --------------------------------------#
wlon = seq(from=lona,to=lonz,by=dblock)
elon = wlon[-1]
wlon = wlon[-length(wlon)]
clon = 0.5 * (wlon+elon)

slat = seq(from=lata,to=latz,by=dblock)
nlat = slat[-1]
slat = slat[-length(slat)]
clat = 0.5 * (slat+nlat)
#------------------------------------------------------------------------------------------#


#----- Find the number of longitude and latitude blocks. ----------------------------------#
nblon = length(clon)
nblat = length(clat)
#------------------------------------------------------------------------------------------#


#----- Create list with blocks. -----------------------------------------------------------#
wlon  = rep(wlon,times = nblat)
elon  = rep(elon,times = nblat)
slat  = rep(slat,each  = nblon)
nlat  = rep(nlat,each  = nblon)
clon  = rep(clon,times = nblat)
clat  = rep(clat,each  = nblon)
block = data.frame( wlon   = wlon
                  , elon   = elon
                  , slat   = slat
                  , nlat   = nlat
                  , clon   = clon
                  , clat   = clat
                  , area   = solid.angle( sw      = cbind(wlon,slat)
                                        , ne      = cbind(elon,nlat)
                                        , degrees = TRUE
                                        , radius  = erad
                                        )#end solid.angle
                  , edarea = (erad * pio180)^2 * cos(pio180 * clat) 
                  )#end data.frame
#------------------------------------------------------------------------------------------#


#----- Free memory. -----------------------------------------------------------------------#
rm(wlon,elon,clon,slat,nlat,clat)
#------------------------------------------------------------------------------------------#



#----- Define the number of blocks to be read. --------------------------------------------#
nblocks   = nrow(block)
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#    Big loop over all blocks.                                                             #
#------------------------------------------------------------------------------------------#
for (b in sequence(nblocks)){

   #----- Format longitudes and latitudes depending on the grid resolution. ---------------#
   if (dblock %in% c(0.5)){
      charclon  = sprintf("%.2f",block$clon[b])
      charclat  = sprintf("%.2f",block$clat[b])
   }else{
      charclon  = sprintf("%.1f",block$clon[b])
      charclat  = sprintf("%.1f",block$clat[b])
   }#end if (dblock %in% c(0.5))
   #---------------------------------------------------------------------------------------#


   #----- Create file names. --------------------------------------------------------------#
   histfile  = file.path(histdir,paste0("lat"         ,charclat,"lon",charclon,".lu"))
   rcpfile   = file.path(rcpdir ,paste0("lat"         ,charclat,"lon",charclon,".lu"))
   extfile   = file.path(extdir ,paste0("lat"         ,charclat,"lon",charclon,".lu"))
   outfile   = file.path(outdir ,paste0(outpref,"-lat",charclat,"lon",charclon,".lu"))
   #---------------------------------------------------------------------------------------#

   #---------------------------------------------------------------------------------------#
   #    Not all files exist, we must check it before to avoid crashes..                    #
   #---------------------------------------------------------------------------------------#
   if (rcp %in% "rcp85_message"){
      fine = file.exists(histfile) && file.exists(rcpfile) && file.exists(extfile)
   }else{
      fine = file.exists(histfile) && file.exists(rcpfile)
   }#end if
   #---------------------------------------------------------------------------------------#

   #---------------------------------------------------------------------------------------#
   #    Generate the ED-friendly file in case all input files exist.                       #
   #---------------------------------------------------------------------------------------#
   if (fine){
      cat(" + Process block centred around ","(",block$clon[b],",",block$clat[b],")."
             ,"  File: ",basename(histfile),".","\n",sep="")

      #----- Read historical and scenario, then bind them. --------------------------------#
      hist.disturb   = read.table(file=histfile,header=TRUE,row.names=1)
      rcp.disturb    = read.table(file=rcpfile ,header=TRUE,row.names=1)
      if (rcp %in% "rcp85_message"){
         ext.disturb = read.table(file=extfile ,header=TRUE,row.names=1)
         disturb     = rbind(hist.disturb,rcp.disturb,ext.disturb)
      }else{
         disturb     = rbind(hist.disturb,rcp.disturb)
      }#end if
      #------------------------------------------------------------------------------------#

      #----- Grab years from row names. ---------------------------------------------------#
      years     = as.integer(dimnames(disturb)[[1]])
      #------------------------------------------------------------------------------------#

      #------------------------------------------------------------------------------------#
      #     The new product has biomass in MgC, and has a new target biomass from          #
      # secondary non-forest patches.  To keep it consistent with the old definition, we   #
      # merge sbh2 and sbh3.  There is no point to keep them separated since they will be  #
      # merged in ED2 anyway.                                                              #
      #------------------------------------------------------------------------------------#
      disturb$sbh    = disturb$sbh  * 1000.
      disturb$vbh    = disturb$vbh  * 1000.
      disturb$sbh2   = ( disturb$sbh2 + disturb$sbh3 ) * 1000.
      disturb$vbh2   = disturb$vbh2 * 1000.
      disturb$f_sbh2 = disturb$f_sbh2 + disturb$f_sbh3
      bye_three      = names(disturb) %in% c("sbh3","f_sbh3")
      disturb        = disturb[,! bye_three]
      #------------------------------------------------------------------------------------#


      #----- Make the table and header look nice. -----------------------------------------#
      yeara       = years[1]
      yearz       = years[length(years)]
      charwlon    = paste(  "WEST.LONGITUDE =",sprintf("%.3f",block$wlon[b]))
      charelon    = paste(  "EAST.LONGITUDE =",sprintf("%.3f",block$elon[b]))
      charslat    = paste(  "SOUTH.LATITUDE =",sprintf("%.3f",block$slat[b]))
      charnlat    = paste(  "NORTH.LATITUDE =",sprintf("%.3f",block$nlat[b]))
      chararea    = paste(  "BLOCK.AREA     =",block$area[b])
      charyeara   = paste(  "FIRST.LUYEAR   =",yeara)
      charyearz   = paste(  "LAST.LUYEAR    =",yearz)
      charnpft    = paste(  "N.PFT.HARVEST  =",npft.harvest)
      charpfth    = paste(c("HARVEST.PFT    =",harvest.pft  ),collapse=" ")
      chardbhslog = paste(c("MINDBH.SLOG    =",mindbh.slog  ),collapse=" ")
      charprbslog = paste(c("HARVPROB.SLOG  =",harvprob.slog),collapse=" ")
      chardbhfplt = paste(c("MINDBH.FPLT    =",mindbh.fplt  ),collapse=" ")
      charprbfplt = paste(c("HARVPROB.FPLT  =",harvprob.fplt),collapse=" ")
      charlabel   = paste(c("year",dimnames(disturb)[[2]]),collapse=" ")
      disturb     = format(disturb)
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #     Add a header to include the area that is associated with this block, and the   #
      # edges.  This makes the data set more generalised.                                  #
      #------------------------------------------------------------------------------------#
      dum     = write(x=charwlon   , file=outfile, append = FALSE)
      dum     = write(x=charelon   , file=outfile, append = TRUE )
      dum     = write(x=charslat   , file=outfile, append = TRUE )
      dum     = write(x=charnlat   , file=outfile, append = TRUE )
      dum     = write(x=chararea   , file=outfile, append = TRUE )
      dum     = write(x=charyeara  , file=outfile, append = TRUE )
      dum     = write(x=charyearz  , file=outfile, append = TRUE )
      dum     = write(x=charnpft   , file=outfile, append = TRUE )
      dum     = write(x=charpfth   , file=outfile, append = TRUE )
      dum     = write(x=chardbhslog, file=outfile, append = TRUE )
      dum     = write(x=charprbslog, file=outfile, append = TRUE )
      dum     = write(x=chardbhfplt, file=outfile, append = TRUE )
      dum     = write(x=charprbfplt, file=outfile, append = TRUE )
      dum     = write(x=charlabel,file=outfile, append = TRUE )
      dum     = write.table( x         = disturb
                           , file      = outfile
                           , append    = TRUE
                           , quote     = FALSE
                           , sep       = " "
                           , row.names = TRUE
                           , col.names = FALSE
                           )#end write.table
      #------------------------------------------------------------------------------------#

   }else{

      #----- Skip block. ------------------------------------------------------------------#
      cat(" + Skip block centred around ","(",block$clon[b],",",block$clat[b],")."
             ,"  File: ",basename(histfile)," not found.","\n",sep="")
      #------------------------------------------------------------------------------------#
   }#end if (file.exists(infile))
   #---------------------------------------------------------------------------------------#
}#end for (b in sequence(nblocks))
#------------------------------------------------------------------------------------------#
