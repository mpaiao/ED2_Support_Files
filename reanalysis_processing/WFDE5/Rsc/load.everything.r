#==========================================================================================#
#==========================================================================================#
#     This script loads all other scripts in this path, and also loads all the necessary   #
# packages.                                                                                #
#------------------------------------------------------------------------------------------#
if ("srcdir" %in% ls()){
   srcdir <<- srcdir
}else{
   srcdir <<- getwd()
}#end if
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#      Find which major version of R is calling this script.                               #
#------------------------------------------------------------------------------------------#
R.major <<- as.numeric(R.version$major)
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#      Make the screen output as wide as the screen permits.                               #
#------------------------------------------------------------------------------------------#
ncstring = as.integer(Sys.getenv("COLUMNS"))
if (! is.na(ncstring)){
   if (ncstring > 80 & ncstring < 500) options(width=ncstring)
}#end if (! is.na(ncstring))
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#     Fix the colours according to the current background.                                 #
#------------------------------------------------------------------------------------------#
if (! "ibackground" %in% ls()) ibackground = 0
if (ibackground == 0){
  foreground    = "black"
  background    = "white"
}else if (ibackground == 1){
   foreground    <<- "white"
   background    <<- "black"
}else if (ibackground == 2){
   foreground    <<- "white"
   background    <<- "#282828"
}else{
   stop(paste0(" Invalid ibackground value (",ibackground,")"))
}#end if
#------------------------------------------------------------------------------------------#




#------------------------------------------------------------------------------------------#
#      Define the size of the titles and axes.                                             #
#------------------------------------------------------------------------------------------#
if (! "ptsz" %in% ls()){
   ptsz <<- 16
}else{
   ptsz <<- ptsz
}#end if
cex.ptsz <<- 1.0 * min(1.0,ptsz / 15)
cex.main <<- 1.1 * min(1.0,ptsz / 14)
cex.lab  <<- 1.0 * min(1.0,ptsz / 14)
#------------------------------------------------------------------------------------------#




#------------------------------------------------------------------------------------------#
#      Define the integration period for photosynthesis-related variables.                 #
#------------------------------------------------------------------------------------------#
if (! "iint.photo" %in% ls()){
   iint.photo <<- 0
}else{
   iint.photo <<- iint.photo
}#end if
#------------------------------------------------------------------------------------------#






#----- Create the default plotting settings for R. ----------------------------------------#
par.user <<- list( bg       = "transparent"
                 , col      = foreground
                 , col.axis = foreground
                 , col.lab  = foreground
                 , col.main = foreground
                 , col.sub  = foreground
                 , fg       = foreground
                 , cex.main = cex.main
                 , cex.lab  = cex.lab
                 , family   = "Helvetica"
                 , mar      = c(5.1,4.4,4.1,2.1)
                 , mgp      = c(2.25,0.25,0)
                 , tcl      = +0.25
                 )#end list
#------------------------------------------------------------------------------------------#


#----- Wrapper for loading packages without pop ups. --------------------------------------#
discreet.require <<- function(...){
   dummy = suppressPackageStartupMessages(suppressWarnings(require(...)))
   return(dummy)
}#end discreet.require
#------------------------------------------------------------------------------------------#

#-----------------------------------------------------------------------------------------#
#     Load all packages needed.                                                            #
#------------------------------------------------------------------------------------------#
loaded.package = list()
loaded.package[["chron"       ]] = discreet.require(chron       )
loaded.package[["data.table"  ]] = discreet.require(data.table  )
loaded.package[["grDevices"   ]] = discreet.require(grDevices   )
loaded.package[["maps"        ]] = discreet.require(maps        )
loaded.package[["mapdata"     ]] = discreet.require(mapdata     )
loaded.package[["maptools"    ]] = discreet.require(maptools    )
loaded.package[["ncdf4"       ]] = discreet.require(ncdf4       )
loaded.package[["rworldmap"   ]] = discreet.require(rworldmap   )
loaded.package[["rhdf5"       ]] = discreet.require(rhdf5       )
loaded.package[["R.utils"     ]] = discreet.require(R.utils     )
#---- Packages that must be loaded at the end. --------------------------------------------#
#loaded.package[["forecast"    ]] = discreet.require(forecast    )
#------------------------------------------------------------------------------------------#


#---- Make sure all packages are loaded fine. ---------------------------------------------#
loaded.package = unlist(loaded.package)
if (! all(loaded.package)){
   miss = which(! loaded.package)
   cat(" You must install the following packages before using the scripts:","\n")
   for (m in miss) cat(" -> ",names(loaded.package)[m],"\n",sep="")
   risky = readline(" Are you sure you want to proceed [y|N]? ")
   risky = tolower(risky)
   if (! risky %in% c("y","yes")) stop("Missing packages!!!")
}#end if
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#  SHADY BUSINESS...  We must unlock grav from package boot and replace by our good        #
#                     old value from rconstants.r.                                         #
#------------------------------------------------------------------------------------------#
#envir = as.environment("package:boot")
#try(unlockBinding("grav",envir),silent=TRUE)
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#  SHADY BUSINESS...  We must unlock trim from package R.oo and raster and replace by our  #
#                     function that has more options than the package one.                 #
#------------------------------------------------------------------------------------------#
envir = as.environment("package:R.oo")
try(unlockBinding("trim",envir),silent=TRUE)
#envir = as.environment("package:raster")
#try(unlockBinding("trim",envir),silent=TRUE)
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#  SHADY BUSINESS...  We must unlock RGB from package raster and replace by our function.  #
#------------------------------------------------------------------------------------------#
#envir = as.environment("package:raster")
#try(unlockBinding("RGB",envir),silent=TRUE)
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#  SHADY BUSINESS...  We must unlock tobin from package survival and replace by our        #
# function.                                                                                #
#------------------------------------------------------------------------------------------#
#envir = as.environment("package:survival")
#try(unlockBinding("tobin",envir),silent=TRUE)
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#  SHADY BUSINESS...  We must unlock pspline from package survival and replace by our      #
# function.                                                                                #
#------------------------------------------------------------------------------------------#
#envir = as.environment("package:survival")
#try(unlockBinding("pspline",envir),silent=TRUE)
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#  SHADY BUSINESS...  We must unlock %>% from package forecast and replace by our          #
# function.                                                                                #
#------------------------------------------------------------------------------------------#
#envir = as.environment("package:forecast")
#try(unlockBinding("%>%",envir),silent=TRUE)
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#     Organise the files so we load them in the right order.                               #
#------------------------------------------------------------------------------------------#
at.first      = c("colour.utils.r","rconstants.r","globdims.r","unitlist.r")
at.end        = c("iso3166.r")
myself        = c("load.everything.r")
all.scripts   = sort(list.files(path=srcdir,pattern="\\.[RrSsQq]$"))
back.up       = sort(list.files(path=srcdir,pattern="^[~]"))
keep          = ! ( all.scripts %in% at.first
                  | all.scripts %in% at.end
                  | all.scripts %in% myself
                  | all.scripts %in% back.up
                  )#end
middle        = all.scripts[keep]
order.scripts = c(at.first,middle,at.end)
nscripts      = length(order.scripts)
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#     Load all files, in order.  Here we replace the warnings by errors, just to make sure #
# that all the functions are clean.                                                        #
#------------------------------------------------------------------------------------------#
warn.orig = getOption("warn")
options(warn=2)
cat(" + Load scripts from ",srcdir,".","\n",sep="")
for (iscript in sequence(nscripts)){
   script.now  = order.scripts[iscript]
   full        = file.path(srcdir,script.now)
   isok        = try(source(full),silent=TRUE)
   if ("try-error" %in% is(isok)){
      options(warn=warn.orig)
      cat("   - Script ",script.now," has bugs!  Check the errors/warnings: ","\n",sep="")
      source(full)
      stop("Source code problem")
   }#end if
}#end for
options(warn=warn.orig)
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#       Get rid of the extremely annoying and unnecessary bell.  Also, force the system to #
# use Helvetica as the default font family.                                                #
#------------------------------------------------------------------------------------------#
options(locatorBell=FALSE,family="Helvetica")
#------------------------------------------------------------------------------------------#

#==========================================================================================#
#==========================================================================================#
