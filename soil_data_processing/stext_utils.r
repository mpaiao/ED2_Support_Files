#==========================================================================================#
#==========================================================================================#
#     Auxilliary functions used by the scripts that generate soil data for ED2.  No need   #
# to load them by hand, scripts gen_ed22_depth.r and gen_ed22_text.r will load this file.  #
#                                                                                          #
#  Developed by Marcos Longo.                                                              #
#------------------------------------------------------------------------------------------#




#==========================================================================================#
#==========================================================================================#
#     This function converts lists into data frames, preserving original types (as long as #
# the types are defined here, you may need to add other types.                             #
#------------------------------------------------------------------------------------------#
list.2.data.frame <<- function(x){
   if (! is.list(x)){
      stop("x must be a list!")
   }else if(length(x) == 0){
      out = data.frame()
   }else{
      xlen          = sapply(X = x, FUN = length)

      #------ Check that the lists all match. ---------------------------------------------#
      if (length(unique(xlen)) != 1){
         stop("All list elements must have the same length")
      }else{
         xlen = unique(xlen)
      }#end if
      #------------------------------------------------------------------------------------#


      #------ If lists are empty, bind but make it empty. ---------------------------------#
      if (xlen == 0){
         out           = data.frame(bye=numeric(length(x)))
         rownames(out) = names(x)
         out           = out[,-1]
      }else{
         #----- Make sure that the names all match. ---------------------------------------#
         names.ok  = apply(X=sapply(X=x,FUN=names),MARGIN=1,FUN=length.unique)
         if (any(names.ok) != 1){
            stop("All list must be the same (same names, and same order)")
         }#end if
         #---------------------------------------------------------------------------------#




         #----- Keep the type. ------------------------------------------------------------#
         xwhat = sapply(X = x[[1]], FUN = typeof, simplify = FALSE)
         xchron        = ( sapply(X= x[[1]], FUN = is.chron) 
                         | sapply(X= x[[1]], FUN = is.dates)
                         )#end 
         xwhat[xchron] = "chron"
         #---------------------------------------------------------------------------------#



         #----- Bind the lists and create data drame. -------------------------------------#
         out        = rbind( apply(X=rbind(sapply(X=x,FUN=c)),MARGIN=1,FUN=unlist)
                           , deparse.level = 1
                           )#end rbind
         names.out  = colnames(out)
         out        = split(x=out,f=col(out))
         names(out) = names.out
         #---------------------------------------------------------------------------------#


         #----- Make sure that variable types are preserved. ------------------------------#
         out = as.data.frame( mapply(FUN=as,object=out,Class=xwhat,SIMPLIFY=FALSE)
                            , stringsAsFactors = FALSE
                            )#end as.data.frame
         rownames(out) = names(x)
         #---------------------------------------------------------------------------------#
      }#end if
      #------------------------------------------------------------------------------------#
   }#end if
   #---------------------------------------------------------------------------------------#

   return(out)
}#end function list.2.data.frame
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#     This function counts how many unique values exist in a vector.                       #
#------------------------------------------------------------------------------------------#
length.unique <<- function(x){
   ans = length(unique(x))
   return(ans)
}#end function lit.sample
#==========================================================================================#
#==========================================================================================#




#==========================================================================================#
#==========================================================================================#
#      Find the nearest neighbour.                                                         #
#------------------------------------------------------------------------------------------#
idx.nearest <<- function(here,there){
   xdist = (here[1] - there[,1])
   ydist = (here[2] - there[,2])
   dist  = xdist*xdist + ydist*ydist
   near  = which.min(dist)
   ans   = there[near,]
   return(ans)
}#end idx.nearest
#==========================================================================================#
#==========================================================================================#



#==========================================================================================#
#==========================================================================================#
#      Function to eliminate missing values from the raster file.                          #
#------------------------------------------------------------------------------------------#
populate <<- function(mat,gmin=1,gby=1){
   #----- Do nothing in case the matrix is empty. -----------------------------------------#
   if (all(is.na(mat)) || (! any(is.na(mat)))) return(mat)
   #---------------------------------------------------------------------------------------#


   #----- Find the gap sizes to try filling. ----------------------------------------------#
   xmat  = row(mat)
   ymat  = col(mat)
   miss  = which(is.na(c(mat)))
   fine  = which(is.finite(c(mat)))
   #---------------------------------------------------------------------------------------#


   #----- Create arrays with the coordinates of missing and fine data. --------------------#
   xy.miss   = cbind(xmat[miss],ymat[miss])
   xy.miss   = split(x=xy.miss,f=row(xy.miss))
   xy.fine   = cbind(xmat[fine],ymat[fine])
   xy.near   = t(sapply(X=xy.miss,FUN=idx.nearest,there=xy.fine))
   ans       = mat
   ans[miss] = mat[xy.near]
   #---------------------------------------------------------------------------------------#

   return(ans)
}#end populate
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#       Mix of cat and past0 (hence cat0).                                                 #
#------------------------------------------------------------------------------------------#
cat0 <<- function(...,file="",fill=FALSE,labels=NULL,append=FALSE){
   cat(...,"\n",file=file,fill=fill,sep="",labels=labels,append=append)
}#end cat0
#==========================================================================================#
#==========================================================================================#





#==========================================================================================#
#==========================================================================================#
#     This function determines the soil class based on the fraction of sand, clay, and     #
# silt separates.                                                                          #
#------------------------------------------------------------------------------------------#
sclass <<- function(sandfrac,clayfrac){
    
   #----- Define the percentage of sand, clay, and silt. ----------------------------------#
   sand = 100. * sandfrac
   clay = 100. * clayfrac
   silt = 100. - sand - clay
   #---------------------------------------------------------------------------------------#


   #----- If the fractions are vectors, use mapply to get the answer. ---------------------#
   if (length(sand) > 1){
      l.sand = as.list(0.01*sand)
      l.clay = as.list(0.01*clay)
      ans    = mapply( FUN      = sclass
                     , sandfrac = l.sand
                     , clayfrac = l.clay
                     , SIMPLIFY = TRUE
                     )#end mapply
      return(ans)
   }#end if
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #     Here there is not much we can do other than explore where in the triangle space   #
   # we are.                                                                               #
   #---------------------------------------------------------------------------------------#
   if ( is.na(silt) || is.na(sand) || is.na(clay) ){
      mysoil = NA_integer_
   }else if (silt > 100. | silt < 0. | sand > 100. | sand < 0. | clay > 100. | clay < 0. ) {
      print("---------------------------------------------------")
      print(" At least one of your percentages is screwy...")
      print(paste0("SAND = ",sprintf("%.2f",sand),"%"))
      print(paste0("CLAY = ",sprintf("%.2f",clay),"%"))
      print(paste0("SILT = ",sprintf("%.2f",silt),"%"))
      print("---------------------------------------------------")
      stop ("This soil doesn''t fit into any category...")
      
   }else if(sand > (85.0 + 0.5 * clay)) {
      mysoil =  1 #----- Sand. ------------------------------------------------------------#
   }else if(sand > (70.0 + clay)) {
      mysoil =  2 #----- Loamy sand. ------------------------------------------------------#
   }else if((clay <= 20.0 & sand >= 52.5) | (clay <= 7.5 & silt <= 50.0)) {
      mysoil =  3 #----- Sandy loam. ------------------------------------------------------#
   }else if((clay <= 27.5 & silt > 50.0 & silt <= 80.0) | (silt >  80.0 & clay > 12.5)) {
      mysoil =  4 #----- Silt loam. -------------------------------------------------------#
   }else if(clay > 7.5 & clay <= 27.5 & silt > 27.5 & silt <= 50.0 & sand <= 52.5) {
      mysoil =  5 #----- Loam. ------------------------------------------------------------#
   }else if(clay > 20.0 & clay <= 35.0 & silt <= 27.5 & sand > 45.0) {
      mysoil =  6 #----- Sandy clay loam. -------------------------------------------------#
   }else if(clay > 27.5 & clay <= 40.0 & sand <= 20.0) {
      mysoil =  7 #----- Silty clay loam. -------------------------------------------------#
   }else if(clay > 27.5 & clay <= 40.0 & sand > 20.0 & sand <= 45.0) {
      mysoil =  8 #----- Clayey loam. -----------------------------------------------------#
   }else if(clay > 35.0 & sand > 45.0) {
      mysoil =  9 #----- Sandy clay. ------------------------------------------------------#
   }else if(clay > 40.0 & silt > 40.0) {
      mysoil = 10 #----- Silty clay. ------------------------------------------------------#
   }else if(clay <= 70.0 & sand <= 30.0 & silt <= 30.0) {
      mysoil = 11 #----- Clay. ------------------------------------------------------------#
   }else if( silt > 80.0 & clay <= 12.5) {
      mysoil = 14 #----- Silt. ------------------------------------------------------------#
   }else if( clay > 70.0) {
      mysoil = 15 #----- Heavy clay. ------------------------------------------------------#
   }else if( clay > 40.0 & sand > 30.0 & sand <= 45.0) {
      mysoil = 16 #----- Clayey sand. -----------------------------------------------------#
   }else if( clay > 40.0 & silt > 30.0 & silt <= 40.0) {
      mysoil = 17 #----- Clayey silt. -----------------------------------------------------#
  }else{
      print("---------------------------------------------------")
      print(paste0("SAND = ",sprintf("%.2f",sand),"%"))
      print(paste0("CLAY = ",sprintf("%.2f",clay),"%"))
      print(paste0("SILT = ",sprintf("%.2f",silt),"%"))
      print("---------------------------------------------------")
      stop ("This soil doesn''t fit into any category...")
  }#end if

  return(mysoil)
}#end function
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#     Function that checks whether an object is chron.                                     #
#------------------------------------------------------------------------------------------#
is.dates <<- function(x) inherits(x,"dates")
is.times <<- function(x) inherits(x,"times")
is.time  <<- function(x){ is.chron(x) || is.dates(x) || is.times(x)}
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#     Function gridded.plot                                                                #
#                                                                                          #
#     This function plots a gridded plot, with the colour scheme given.  This is similar   #
# to sombreado, but structured more like image.map.  Currently only one panel is allowed,  #
# though this may change soon.                                                             #
#------------------------------------------------------------------------------------------#
gridded.plot <<- function( x                = seq(from=0,to=1,len=nrow(z))
                         , y                = seq(from=0,to=1,len=ncol(z))
                         , z
                         , xlim             = range(unlist(x),finite=TRUE)
                         , ylim             = range(unlist(y),finite=TRUE)
                         , zlim             = range(unlist(z),finite=TRUE)
                         , xlog             = FALSE
                         , ylog             = FALSE
                         , levels           = if (key.log){
                                                 sort(unique(pretty.log(x=zlim,n=nlevels
                                                                       ,forcelog=TRUE)))
                                              }else{
                                                 sort(unique(pretty    (x=zlim,n=nlevels)))
                                              }#end if
                         , nlevels          = 20
                         , colour.palette   = cm.colors
                         , col              = colour.palette(length(levels)-1)
                         , na.col           = "grey94"
                         , key.log          = FALSE
                         , key.vertical     = TRUE
                         , x.axis.options   = NULL
                         , y.axis.options   = NULL
                         , key.axis.options = NULL
                         , key.options      = NULL
                         , main.title       = NULL
                         , key.title        = NULL
                         , plot.after       = NULL
                         , matrix.plot      = FALSE
                         , legend.options   = NULL
                         , edge.axes        = FALSE
                         , mar              = NULL
                         , mar.main         = mar
                         , mar.key          = NULL
                         , oma              = NULL
                         , omd              = NULL
                         , f.key            = 1/6
                         , f.leg            = 1/6
                         , xaxs             = "i"
                         , yaxs             = "i"
                         , smidgen          = 0
                         , useRaster        = TRUE
                         , ...
                         ){

   #----- Check which kind of input was given. --------------------------------------------#
   if (missing(z)) {
      #----- No z was given x must be a list or the user didn't provide any axis... -------# 
      if (!missing(x)) {
         if (is.list(x)) {
            #----- X is a list, copy the elements to variables. ---------------------------#
            z = x$z
            y = x$y
            x = x$x
         }else{
            #----- x is an array, make up some x axis. ------------------------------------#
            z = x
            x = seq(from = 0, to = 1, length.out = nrow(z))
            y = seq(from = 0, to = 1, length.out = ncol(z))
         }#end if
         #---------------------------------------------------------------------------------#
      }else{
         #----- Bad setting. -------------------------------------------------------------#
         stop("no `z' matrix specified")
         #---------------------------------------------------------------------------------#
      }#end if
      #------------------------------------------------------------------------------------#
   }else if (is.list(x)) {
      #----- Z is there, just need to check whether x and y were given as a list... -------#
      y = x$y
      x = x$x
      #------------------------------------------------------------------------------------#
   }#end if
   #---------------------------------------------------------------------------------------#

   #----- Check whether the z matrix makes sense or not. ----------------------------------#
   if (! (is.matrix(z) || is.data.frame(z))){
      stop("no proper `z' matrix specified")
   }#end if
   #---------------------------------------------------------------------------------------#


   #----- No messed-up axes are allowed, they must increase. ------------------------------#
   if (any(diff(x) %le% 0) || any(diff(y) %le% 0)){
       stop("increasing x and y values expected")
   }#end if
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #     If legend is to be plotted, key.vertical has to be TRUE.  In case the user said   #
   # otherwise, return a warning.  Also, define offsets for X and Y according to the       #
   # legends and keys.                                                                     #
   #---------------------------------------------------------------------------------------#
   plot.legend = ! is.null(legend.options)
   if ( plot.legend && (! key.vertical)){
      warning(" key.vertical=FALSE ignored due to the legend.")
      key.vertical = TRUE
   }#end if
   #---------------------------------------------------------------------------------------#


   #----- Save the margins to avoid losing the data. --------------------------------------#
   par.orig = par(no.readonly=TRUE)
   mar.orig = par.orig$mar
   on.exit(par(par.orig))
   foreground = par("fg")
   par( bg       = "transparent"
      , col      = foreground
      , col.axis = foreground
      , col.lab  = foreground
      , col.main = foreground
      , col.sub  = foreground
      , fg       = foreground
      , cex.main = 1.0
      , cex.lab  = 1.0
      , family   = "Helvetica"
      , mar      = c(5.1,4.4,4.1,2.1)
      , mgp      = c(2.25,0.25,0)
      , tcl      = +0.25
      )#end par
   #---------------------------------------------------------------------------------------#





   #----- Check for outer margins. --------------------------------------------------------#
   if ( (! is.null(oma)) && (! is.null(omd))){
      stop ("You cannot provide both oma and omd!")
   }else if (is.null(oma) && is.null(omd)){
      par(oma=c(0,0,0,0))
   }else if (is.null(omd)){
      par(oma=oma)
   }else{
      par(omd=omd)
   }#end if
   #---------------------------------------------------------------------------------------#



   #----- Split the screen into multiple pieces (legend, key, plots...) -------------------#
   fh.panel = 1. - f.key
   fv.panel = 1. - f.leg
   if (plot.legend){
      layout( mat     = rbind(c(3, 2),c(1,0))
            , heights = c(fv.panel,f.leg)
            , widths  = c(fh.panel,f.key)
            )#end layout
   }else if (key.vertical){
      layout(mat=cbind(2, 1), widths = c(fh.panel,f.key))
   }else{
      layout( mat     = rbind(2,1)
            , heights = c(fh.panel,f.key)
            )#end layout
   }#end if
   #---------------------------------------------------------------------------------------#
   #=======================================================================================#
   #=======================================================================================#






   #=======================================================================================#
   #=======================================================================================#
   #      First plot: the legend.                                                          #
   #---------------------------------------------------------------------------------------#
   if (plot.legend){
      par(mar = c(0.1,0.1,0.1,0.1))
      plot.new()
      plot.window(xlim=c(0,1),ylim=c(0,1))
      do.call(what="legend",args=legend.options)
   }#end if
   #=======================================================================================#
   #=======================================================================================#






   #=======================================================================================#
   #=======================================================================================#
   #      Second plot: the key scale.                                                      #
   #---------------------------------------------------------------------------------------#
      if (! is.null(mar.key)){
         par(mar = mar.key)
      }else if (key.vertical){
         par(mar = pretty.box(n=1)$mar.key)
      }else{
         par(mar = c(2.1,4.6,1.6,2.1))
      }#end if
      plot.new()
      #------------------------------------------------------------------------------------#
      #     Plot in the horizontal or vertical depending on where the scale is going to    #
      # be plotted.                                                                        #
      #------------------------------------------------------------------------------------#
      if (key.vertical){
         #----- Decide whether the scale is logarithmic or not. ---------------------------#
         if (key.log){
            plot.window(xlim=c(0,1),ylim=range(levels),xaxs="i",yaxs="i",log="y")
         }else{
            plot.window(xlim=c(0,1),ylim=range(levels),xaxs="i",yaxs="i")
         }#end if
         #---------------------------------------------------------------------------------#

         #----- Draw the colour bar. ------------------------------------------------------#
         rect(xleft=0,ybottom=levels[-length(levels)],xright=1,ytop=levels[-1],col=col
             ,border=col)
         #---------------------------------------------------------------------------------#

         #----- Check whether there are specific instructions for plotting the key axis. --#
         if (missing(key.axis.options)) {
            key.now = list(side=4,las=1,...)
         }else{
            key.now = modifyList(x=key.axis.options,val=list(side=4,las=1))
         }#end if
         do.call (what="axis",args=key.now)
         #---------------------------------------------------------------------------------#
      }else{
         #----- Decide whether the scale is logarithmic or not. ---------------------------#
         if (key.log){
            plot.window(xlim=range(levels),ylim=c(0,1),xaxs="i",yaxs="i",las=1,log="x")
         }else{
            plot.window(xlim=range(levels),ylim=c(0,1),xaxs="i",yaxs="i",las=1)
         }#end if
         #---------------------------------------------------------------------------------#


         #----- Draw the colour bar. ------------------------------------------------------#
         rect(xleft=levels[-length(levels)],ybottom=0,xright=levels[-1],ytop=1
             ,col=col,border=col)
         #---------------------------------------------------------------------------------#


         #----- Check whether there are specific instructions for plotting the key axis. --#
         if (missing(key.axis.options)) {
            key.now = list(side=1,las=1,...)
         }else{
            key.now = modifyList(x=key.axis.options,val=list(side=1,las=1))
         }#end if
         do.call (what="axis",args=key.now)
         #---------------------------------------------------------------------------------#
      }#end if
      #------------------------------------------------------------------------------------#


      #----- Draw box. --------------------------------------------------------------------#
      box()
      #------------------------------------------------------------------------------------#


      #----- Plot the title. --------------------------------------------------------------#
      if (! is.null(key.title)) do.call(what="title",args=key.title)
      #------------------------------------------------------------------------------------#
   #=======================================================================================#
   #=======================================================================================#


   #----- Set the window. -----------------------------------------------------------------#
   if (! is.null(mar.main)){
      mar.now = mar.main
   }else if (key.vertical){
      mar.now = c(5.1,4.1,4.1,2.1)
   }else{
      mar.now = c(4.1,4.1,4.1,2.1)
   }#end if
   plog = ""
   if (xlog) plog = paste(plog,"x",sep="")
   if (ylog) plog = paste(plog,"y",sep="")
   par(mar = mar.now)
   plot.new()
   plot.window(xlim=xlim,ylim=ylim,log=plog,xaxs=xaxs,yaxs=yaxs,...)
   #---------------------------------------------------------------------------------------#


   #----- Plot the field. -----------------------------------------------------------------#
   xyz = list(x=x,y=y,z=z)
   image(xyz,breaks=levels,col=col,add=TRUE,useRaster=useRaster)
   #---------------------------------------------------------------------------------------#




   #---- Plot the X axis. -----------------------------------------------------------------#
   if (! is.null(x.axis.options)){
      x.axis.now = modifyList(x=x.axis.options,val=list(side=1))
   }else{
      x.axis.now = list(side=1,las=1)
   }#end if
   do.call(what="axis",args=x.axis.now)
   #---------------------------------------------------------------------------------------#




   #---- Plot the Y axis. -----------------------------------------------------------------#
   if (! is.null(y.axis.options)){
      y.axis.now = modifyList(x=y.axis.options,val=list(side=2))
   }else{
      y.axis.now = list(side=2,las=1)
   }#end if
   do.call(what="axis",args=y.axis.now)
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Plot other options.                                                               #
   #---------------------------------------------------------------------------------------#
   n.after = length(plot.after)
   for (a in sequence(n.after)){
       do.call(what=names(plot.after)[a],args=plot.after[[a]])
   }#end for
   #---------------------------------------------------------------------------------------#


   #----- Lastly, add the box (so it stays on top). ---------------------------------------#
   box()
   #---------------------------------------------------------------------------------------#



   #----- Make sure we get the main text. -------------------------------------------------#
   if (! is.list(main.title)){
      main.title=list(main=main.title)
   }else if (! "main" %in% names(main.title)){
      names(main.title)[[1]] = "main"
   }#end if
   do.call(what="title",args=main.title)
   #---------------------------------------------------------------------------------------#


   #----- Don't bother the user with messages. --------------------------------------------#
   invisible()
   #---------------------------------------------------------------------------------------#
}#end function gridded.plot
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#    This function splits the number into number of rows and columns that will make the    #
# prettiest (prettiest here means the closest to a golden ratio rectangle).                #
#------------------------------------------------------------------------------------------#
pretty.box = function(n,horizontal=TRUE,byrow=TRUE,f.ncol=NULL,f.nrow=NULL
                     ,angle.crit=atan2((1.+sqrt(5))/2.,1)*180./pi){


   #---------------------------------------------------------------------------------------#
   #      Check the size of n.  If it has a single value, then we must find the number of  #
   # rows and columns.  If it has lenght 2, then the first value becomes the number of     #
   # columns and the second becomes the number of rows.                                    #
   #---------------------------------------------------------------------------------------#
   n.n = length(n)
   if (n.n == 0){
      #------------------------------------------------------------------------------------#
      #      Empty! Make a single box.                                                     #
      #------------------------------------------------------------------------------------#
      nbrow  = 1
      nbcol  = 1
      nbox   = 1
      mat    = matrix(1,ncol=1,nrow=1)
      xangle = 45
      yangle = 45
      #------------------------------------------------------------------------------------#

   }else if (n.n == 1){
      #------------------------------------------------------------------------------------#
      #      Single site.  We must test whether we can split in nows and columns that      #
      # looks nice by the closest to square root.  If the number is prime, we must add one #
      # so it can be split.   We keep adding empty boxes until we get a number that has a  #
      # good ratio.                                                                        #
      #------------------------------------------------------------------------------------#
      iterate = TRUE
      n       = max(n,1)
      nbox    = n
      while (iterate){

         nbrow.pot = seq(from=1,to=floor(nbox),by = 1)
         nbcol.pot = nbox / nbrow.pot

         #----- Select only the cases that make a 2-D grid with the right orientation. ----#
         if (horizontal){
            ok  = nbcol.pot == as.integer(nbcol.pot) & nbcol.pot >= nbrow.pot
         }else{
            ok  = nbcol.pot == as.integer(nbcol.pot) & nbcol.pot <= nbrow.pot
         }#end if
         #---------------------------------------------------------------------------------#


         #----- Discard bad options. ------------------------------------------------------#
         nbrow.pot = nbrow.pot[ok]
         nbcol.pot = as.integer(nbcol.pot[ok])
         #---------------------------------------------------------------------------------#


         #---------------------------------------------------------------------------------#
         #      Choose the one that is the closest to the square.                          #
         #---------------------------------------------------------------------------------#
         nuse  = which.min(abs(nbrow.pot-nbcol.pot))
         nbrow = nbrow.pot[nuse]
         nbcol = nbcol.pot[nuse]
         #---------------------------------------------------------------------------------#


         #---------------------------------------------------------------------------------#
         #      Make a matrix for the layout.                                              #
         #---------------------------------------------------------------------------------#
         mat          = matrix(sequence(nbrow*nbcol),ncol=nbcol,nrow=nbrow,byrow=byrow)
         mat[mat > n] = 0
         #---------------------------------------------------------------------------------#


         #---------------------------------------------------------------------------------#
         #      Check whether the angle is acceptable...                                   #
         #---------------------------------------------------------------------------------#
         xangle = atan2(nbrow,nbcol) * 180. / pi
         yangle = atan2(nbcol,nbrow) * 180. / pi
         iterate = nbox > 2 && (xangle > angle.crit | yangle > angle.crit)
         if (iterate) nbox = nbox + 1
         #---------------------------------------------------------------------------------#
      }#end while
      #------------------------------------------------------------------------------------#

   }else if (n.n == 2){
      #------------------------------------------------------------------------------------#
      #      Two values were given.  No iteration needed.                                  #
      #------------------------------------------------------------------------------------#
      nbcol  = max(n[1],1)
      nbrow  = max(n[2],1)
      nbox   = nbcol*nbrow
      mat    = matrix(sequence(nbrow*nbcol),ncol=nbcol,nrow=nbrow,byrow=byrow)
      xangle = atan2(nbrow,nbcol) * 180. / pi
      yangle = atan2(nbcol,nbrow) * 180. / pi
      #------------------------------------------------------------------------------------#

   }else{
      #------------------------------------------------------------------------------------#
      #      3-D boxes are not available...                                                #
      #------------------------------------------------------------------------------------#
      cat (" n           = (",paste(n,sep="; "),")","\n")
      cat (" length of n =  ",length(n)            ,"\n")
      stop(" Invalid n! It must have length between 0 and 2")
      #------------------------------------------------------------------------------------#
   }#end if
   #---------------------------------------------------------------------------------------#


   #----- Make a matrix with offset of one in case the user appends a legend. -------------#
   sel           = mat == 0
   mat.off       = mat + 1
   mat.off[sel]  = 0
   mat.off2      = mat + 2
   mat.off2[sel] = 0
   #---------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   #     Create the margins for each box assuming that axes should be plotted only on the  #
   # bottom and left panels.                                                               #
   #---------------------------------------------------------------------------------------#
   n.seq = sequence(nbox)
   if (byrow){
      left   = ( n.seq %% nbcol ) == 1 | nbcol == 1
      right  = ( n.seq %% nbcol ) == 0
      top    = n.seq <=   nbcol
      bottom = n.seq >  ( nbrow - 1 ) * nbcol
   }else{
      left   = n.seq <= nbrow
      right  = n.seq > ( nbcol - 1 ) * nbrow
      top    = ( n.seq %% nbrow ) == 1 | nbrow == 1
      bottom = ( n.seq %% nbrow ) == 0
   }#end if
   centre = ! ( left | right )
   middle = ! ( bottom | top )
   if (nbcol == 1){
      mar.left  = rep(5.1,times=nbox)
      mar.right = rep(1.1,times=nbox)
   }else if (nbcol == 2){
      mar.left  = 2.1 + 2 * left
      mar.right = 0.1 + 2 * right
   }else{
      mar.left  = 3.1 + 0 * left
      mar.right = 0.6 + 0. * right
   }#end if
   if (nbrow == 1){
      mar.bottom = rep(5.1,times=nbox)
      mar.top    = rep(4.1,times=nbox)
   }else if (nbrow == 2){
      mar.bottom = 2.1 + 2. * bottom
      mar.top    = 2.1 + 2. * top
   }else{
      mar.bottom = 3.1 + 0. * bottom
      mar.top    = 2.1 + 0. * bottom
   }#end if
   mar                = cbind(mar.bottom,mar.left,mar.top,mar.right)
   dimnames(mar)[[2]] = c("bottom","left","top","right")

   #----- Margins for when all boxes get axis plots. --------------------------------------#
   if (nbcol == 1){
      mar0.left  = 5.1
      mar0.right = 1.1
   }else{
      mar0.left  = 4.1
      mar0.right = 0.1
   }#end if
   if (nbrow == 1){
      mar0.bottom = 5.1
      mar0.top    = 4.1
   }else{
      mar0.bottom = 4.1
      mar0.top    = 2.1
   }#end if
   mar0        = c(mar0.bottom,mar0.left,mar0.top,mar0.right)
   names(mar0) = c("bottom","left","top","right")
   #---------------------------------------------------------------------------------------#




   #----- Key margins for when all boxes get axis plots. ----------------------------------#
   key.left  = 0.6
   key.right = 4.1
   if (nbrow == 1){
      key.bottom = 5.1
      key.top    = 4.1
   }else{
      key.bottom = 4.1
      key.top    = 2.1
   }#end if
   mar.key        = c(key.bottom,key.left,key.top,key.right)
   names(mar.key) = c("bottom","left","top","right")
   #---------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   #     Return a list with the information.                                               #
   #---------------------------------------------------------------------------------------#
   ans = list( nrow       = nbrow
             , ncol       = nbcol
             , nbox       = nbox
             , angle      = horizontal * yangle + (1 - horizontal) * xangle
             , mat        = mat
             , mat.off    = mat.off
             , mat.off2   = mat.off2
             , bottom     = bottom
             , left       = left
             , top        = top
             , right      = right
             , mar        = mar
             , mar0       = mar0
             , mar.key    = mar.key
             )#end list
   return(ans)
   #---------------------------------------------------------------------------------------#
}#end function
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#     Function that creates a pretty scale in log co-ordinates.                            #
#------------------------------------------------------------------------------------------#
pretty.log = function(x,base=10,n=5,forcelog=FALSE){
   log.neat  = pretty(x=log(x,base=base),n=n)
   dlog.neat = median(diff(log.neat))
   neat      = base^log.neat
   nact      = length(neat)

   #---------------------------------------------------------------------------------------#
   #     In case it is a base 10, make it even prettier.  Also, in case the log scale is   #
   # very close to the actual scale, forget about the log scale and use regular pretty, it #
   # gives a more legible scale.                                                           #
   #---------------------------------------------------------------------------------------#
   if (base == 10 && dlog.neat %wr% c(0.1,0.5) && (! forcelog)){
      sel  = abs(log.neat - as.integer(log.neat)) %lt% (0.5 * dlog.neat)
      tens = sort(c(neat[sel],base^(c(floor(min(log.neat)),ceiling(max(log.neat))))))

      #------------------------------------------------------------------------------------#
      #      Fix scale so it looks nicer if dtens is 2 or 3.                               #
      #------------------------------------------------------------------------------------#
      if (dlog.neat %lt% 0.15){
         mult    = c(1,2,3,5,7)
      }else if (dlog.neat %lt% 0.35){
         mult    = c(1,2,5)
      }else{
         mult    = c(1,3)
      }#end if
      vlevels   = sort(unique(c(mult %o% tens)))
      aa        = min(which(vlevels > min(neat)))-1
      zz        = max(which(vlevels < max(neat)))+1
      vlevels   = vlevels[aa:zz]
      still.log = TRUE
      #------------------------------------------------------------------------------------#
   }else if(dlog.neat %lt% 0.1 && (! forcelog)){
      #-----  The plot is hardly log, use normal units instead. ---------------------------#
      vlevels = pretty(x=x,n=n)
      still.log = FALSE
      #------------------------------------------------------------------------------------#
   }else{
      vlevels   = neat
      still.log = TRUE
   }#end if
   #---------------------------------------------------------------------------------------#



   if (still.log){
      power.base  = base^(floor(log(x=vlevels,base=base)))
      npow        = length(unique(power.base))
      ndigits     = ceiling(log(nact/npow,base=base))
      vlevels     = round(vlevels / power.base,digits=ndigits) * power.base
   }#end if

   return(vlevels)
}#end function
#==========================================================================================#
#==========================================================================================#







#==========================================================================================#
#==========================================================================================#
#     Operators to test whether the values are within range or outside range.              #
#                                                                                          #
# x %wr% y -- TRUE when x is within range of y     (x exactly at bounds returns TRUE )     #
#------------------------------------------------------------------------------------------#
'%wr%' <<- function(x,y){
   ans = (! (is.na(x) | is.nan(x))) & ( x >= min(y,na.rm=TRUE) & x <= max(y,na.rm=TRUE) )
   return(ans)
}#end function
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#      Safe logical operators.  These will always return FALSE if x or y are not finite.   #
#------------------------------------------------------------------------------------------#
'%eq%' <<- function(x,y){
   if (any(c(FALSE,is.numeric(x) & is.numeric(y)),na.rm=TRUE)){
      ans = is.finite(unlist(x)) & is.finite(unlist(y)) & x == y
   }else{
      ans = ! is.na(x) & ! is.na(y) & ! is.nan(x) & ! is.nan(y) & x == y
   }#end if
   return(ans)
}#end function
'%ne%' <<- function(x,y){
   if (any(c(FALSE,is.numeric(x) & is.numeric(y)),na.rm=TRUE)){
      ans = is.finite(unlist(x)) & is.finite(unlist(y)) & x != y
   }else{
      ans = ! is.na(x) & ! is.na(y) & ! is.nan(x) & ! is.nan(y) & x != y
   }#end if
   return(ans)
}#end function
'%gt%' <<- function(x,y) is.finite(unlist(x)) & is.finite(unlist(y)) & x  > y
'%lt%' <<- function(x,y) is.finite(unlist(x)) & is.finite(unlist(y)) & x  < y
'%ge%' <<- function(x,y) is.finite(unlist(x)) & is.finite(unlist(y)) & x >= y
'%le%' <<- function(x,y) is.finite(unlist(x)) & is.finite(unlist(y)) & x <= y
#==========================================================================================#
#==========================================================================================#







#==========================================================================================#
#==========================================================================================#
#      Magma.                                                                              #
#------------------------------------------------------------------------------------------#
magma <<- function(n){
   #----- Color entries. ------------------------------------------------------------------#
   nodes     = c( "#000004", "#1A1042", "#4A1079", "#792282", "#AA337D", "#D9466B"
                , "#F7725C", "#FEAA74", "#FDE2A3", "#FCFDBF")
   #---------------------------------------------------------------------------------------#

   #----- Call the interpolator. ----------------------------------------------------------#
   ans   = colour.interpol(nodes=nodes,n=n)
   return(ans)
   #---------------------------------------------------------------------------------------#
}#end function magma
#----- Inverse colour palette. ------------------------------------------------------------#
imagma <- function(n,alpha=1.0){ rev(magma(n=n))}
#------------------------------------------------------------------------------------------#
#==========================================================================================#
#==========================================================================================#





#==========================================================================================#
#==========================================================================================#
#      General function that interpolates nodes by using HSV.                              #
#------------------------------------------------------------------------------------------#
colour.interpol <<- function(nodes,n){
   #----- Decompose colours. --------------------------------------------------------------#
   rgbmat = data.frame(t(col2rgb(nodes)))
   #---------------------------------------------------------------------------------------#


   #----- Find the reference points for the interpolation functions. ----------------------#
   rgbmat$x = seq(from=0,to=1,length.out=nrow(rgbmat))
   #---------------------------------------------------------------------------------------#


   #----- Create functions to interpolate. ------------------------------------------------#
   f.red   = splinefun(x=rgbmat$x,y=rgbmat$red  ,method="monoH.FC")
   f.green = splinefun(x=rgbmat$x,y=rgbmat$green,method="monoH.FC")
   f.blue  = splinefun(x=rgbmat$x,y=rgbmat$blue ,method="monoH.FC")
   #---------------------------------------------------------------------------------------#


   #----- Create output table. ------------------------------------------------------------#
   rgbout       = data.frame( x = seq(from=0,to=1,length.out=n) )
   rgbout$red   = round(pmax(0,pmin(255,f.red  (rgbout$x))))
   rgbout$green = round(pmax(0,pmin(255,f.green(rgbout$x))))
   rgbout$blue  = round(pmax(0,pmin(255,f.blue (rgbout$x))))
   #---------------------------------------------------------------------------------------#
   

   #------ Colour rainbow. ----------------------------------------------------------------#
   ans = with(rgbout,rgb(red=red,green=green,blue=blue,maxColorValue=255))
   return(ans)
   #---------------------------------------------------------------------------------------#
}#end colour.interpol
#==========================================================================================#
#==========================================================================================#
