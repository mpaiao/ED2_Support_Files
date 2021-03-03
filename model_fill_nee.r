#==========================================================================================#
#     Close all devices, and delete all variables.                                         #
#==========================================================================================#
rm(list=ls())
graphics.off()
options(warn=0)
#==========================================================================================#


#==========================================================================================#
#           Set up                                                                         #
#==========================================================================================#
#------------- Plotting on/off ------------------------------------------------------------#
plot.par.curves = FALSE
plot.peaks      = FALSE
plot.others     = TRUE
#------------- Load a previous session if the RData exists? -------------------------------#
reload         = FALSE
#------------- Paths ----------------------------------------------------------------------#
here           = getwd()
srcdir         = file.path(here,"Rsc")
outroot        = file.path(here,"figures")
#------------- Site-specific information. -------------------------------------------------#
iata           = "tb0"                   # 3-letter identification for the site
asciipref      = "tb0"                   # Prefix for the output ascii file
place          = "Tanguro (Control), MT" # Full name (plot title)
yeara          = 2008                    # First year with data
yearz          = 2016                    # Last year with data
dry.beg        = "08/01"                 # Average onset of the dry season (format "mm/dd")
dry.end        = "12/01"                 # Average end of the dry season   (format "mm/dd")
yeara.ascii    = 2008                    # First year for the output
yearz.ascii    = 2016                    # Last year for the output
#------- Miscellaneous set up. ------------------------------------------------------------#
nrunave        = 50             # Number of days for the running average
#------- Minimum number of consecutive days with no data to be considered blackout. -------#
min.blackout   = 70000
#------------------------------------------------------------------------------------------#
#    The following two variables control the u* filter:                                    #
# - ust.filter - all u* thresholds you would like to try.                                  #
# - ust.output - which u* to use in the output.                                            #
#------------------------------------------------------------------------------------------#
ust.filter   = seq(from=0.00,to=0.55,by=0.01)
ust.output   = 0.15
#------------------------------------------------------------------------------------------#

#----- Variables for nocturnal fluxes. ----------------------------------------------------#
ust.xlim  = c(0,0.8)  # X Limits.  NULL means that R decides it.
flux.ylim = c(-3,12)  # Y Limits.  NULL means that R decides it.
n.boot    = 1000      # Number of bootstrap samples
#------------------------------------------------------------------------------------------#


#------ Determine wet-dry seasons dynamically (not used anywhere in the code). ------------#
find.dyn.season = FALSE
#------------------------------------------------------------------------------------------#


#------------------------------------------------------------------------------------------#
#    The following variables control the photosynthetic capacity:                          #
# - par.phcap.min - Minimum PAR.                                                           #
# - par.phcap.max - Maximum PAR.                                                           #
#------------------------------------------------------------------------------------------#
par.phcap.min = 500
par.phcap.max = 700
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#     Plot options.                                                                        #
#------------------------------------------------------------------------------------------#
outform         = c("pdf")  # Formats for output file (case insensitive):
                            #   - "X11" - for printing on screen
                            #   - "eps" - for postscript printing
                            #   - "png" - for PNG printing
                            #   - "pdf" - for PDF printing
depth           = 96        # PNG resolution, in pixels per inch
ncolours.im     = 200       # Number of colours for image.map plots
paper           = "square"  # Paper size, to define the plot shape
ptsz            = 21        # Font size.
fracexp         = 0.4       # Expansion to fit the legend
f.leg           = 1/6       # Fraction for xy 
#------------------------------------------------------------------------------------------#
#==========================================================================================#
#==========================================================================================#




#==========================================================================================#
#==========================================================================================#
#    List of variables for the output.                                                     #
#------------------------------------------------------------------------------------------#
cflxvar       = list()
cflxvar[[ 1]] = list( name = "NEE"
                    , desc = "Net ecosystem exchange"
                    , unit = "kgcom2oyr"
                    , mult = "umols.2.kgCyr"
                    , col  = "chartreuse"
                    )
cflxvar[[ 2]] = list( name = "GEP"
                    , desc = "Gross ecosystem exchange"
                    , unit = "kgcom2oyr"
                    , mult = "umols.2.kgCyr"
                    , col  = "forestgreen"
                    )
cflxvar[[ 3]] = list( name = "RECO"
                    , desc = "Ecosystem respiration"
                    , unit = "kgcom2oyr"
                    , mult = "umols.2.kgCyr"
                    , col  = "orangered"
                    )
cflxvar[[ 4]] = list( name = "CCNT"
                    , desc = "Gap-filled fraction"
                    , unit = "pc"
                    , mult = NA_character_
                    , col  = "grey16"
                    )
#------------------------------------------------------------------------------------------#




#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#        Changes beyond this point are for developing the script only.                     #
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#





#==========================================================================================#
#     Load libraries and scripts.                                                          #
#==========================================================================================#
srcdir = (srcdir[file.exists(srcdir)])[1]
source(file.path(srcdir,"load.everything.r"))
#==========================================================================================#



#----- Determine how many formats we must output. -----------------------------------------#
outform = tolower(outform)
nout    = length(outform)
#------------------------------------------------------------------------------------------#


#----- Create the path for figures in case it doesn't exist. ------------------------------#
lieu = poilist$short[match(iata,poilist$iata)]
outplace = file.path(outroot,lieu)
if (! file.exists(outroot )) dir.create(outroot )
if (! file.exists(outplace)) dir.create(outplace)
#------------------------------------------------------------------------------------------#



#==========================================================================================#
#    Length of the flux variables (for averaging and plotting).                            #
#==========================================================================================#
cflxvar   = list.2.data.frame(cflxvar)
ncflxvar  = nrow(cflxvar)
#==========================================================================================#



#==========================================================================================#
#    Beginning and end of the time series (including black out periods).                   #
#==========================================================================================#
whena = chron(paste( 1, 1,yeara,sep="/"))
whenz = chron(paste(12,31,yearz,sep="/"))
#==========================================================================================#


#==========================================================================================#
#     Load dataset.                                                                        #
#==========================================================================================#
metfill.file = file.path(here,paste0("metfill_",iata,"_",yeara,"-",yearz,".RData"))
filled.file  = file.path(here,paste0("filled_" ,iata,"_",yeara,"-",yearz,".RData"))
#==========================================================================================#



#----- Define plot window size ------------------------------------------------------------#
f.ext = f.leg / (1. - f.leg)
sq.size = plotsize(proje=FALSE,paper=paper)
ex.size = plotsize(proje=FALSE,paper=paper,extendfc="lon",extfactor=f.ext)
ey.size = plotsize(proje=FALSE,paper=paper,extendfc="lat",extfactor=f.ext)
#------------------------------------------------------------------------------------------#



#==========================================================================================#
#     Check whether to reload the previous calculation or run it again.                    #
#==========================================================================================#
if ( reload && file.exists(filled.file)){
   #=======================================================================================#
   #=======================================================================================#
   #      Reload the previous calculation.                                                 #
   #=======================================================================================#
   #=======================================================================================#
   cat0(" + Retrieve data from ",basename(filled.file),".")
   load(filled.file)
   eft    = get(iata)
   ustall = get(paste("ust",iata,sep="."))
   nwhen  = length(eft$when)
   #---------------------------------------------------------------------------------------#



   #=======================================================================================#
   #    Retrieve the tables and filters from the previous session.                         #
   #=======================================================================================#
   ust.filter = ustall$ust.filter
   nustar     = length(ust.filter)
   ubest      = ustall$ubest
   ualt       = ustall$ualt
   ust.alt    = ust.filter[ualt]
   NEE        = ustall$nee
   GEP        = ustall$gpp
   GEE        = ustall$gee
   RECO       = ustall$reco
   CCNT       = ustall$ccnt
   flag.NEE   = ustall$gfflg.nee
   flag.GEP   = ustall$gfflg.gpp
   flag.RECO  = ustall$gfflg.reco
   #=======================================================================================#

}else{
   #=======================================================================================#
   #=======================================================================================#
   #      Find the productivity and respiration again.                                     #
   #=======================================================================================#
   #=======================================================================================#


   #=======================================================================================#
   #     Load the metfill file.                                                            #
   #=======================================================================================#
   cat0(" + Load data from ",basename(metfill.file),".")
   load(metfill.file)
   eft = get(iata)



   #=======================================================================================#
   #    Length of the dataset.                                                             #
   #=======================================================================================#
   nwhen     = length(eft$when)
   #=======================================================================================#



   #=======================================================================================#
   #      empty is just a vector of length nwhen full of NA.                               #
   #=======================================================================================#
   empty = rep(NA,times=nwhen)
   #=======================================================================================#




   #=======================================================================================#
   #    Create some dummy variables if they aren't already part of the dataset.            #
   #=======================================================================================#
   cat0(" + Find the dew point temperature and setting other T variables to NA.")
   if (! "atm.tdew" %in% names(eft)) eft$atm.tdew = tslif(eft$atm.prss,eft$atm.pvap,"pvap")
   if (! "atm.tamb" %in% names(eft)) eft$atm.tamb = empty
   if (! "atm.tson" %in% names(eft)) eft$atm.tson = empty
   #----- Same goes for the gap filling flags. --------------------------------------------#
   if (! "gfflg.tdew" %in% names(eft)) eft$gfflg.atm.tdew = empty
   if (! "gfflg.tamb" %in% names(eft)) eft$gfflg.atm.tamb = empty
   if (! "gfflg.tson" %in% names(eft)) eft$gfflg.atm.tson = empty
   #=======================================================================================#



   #-------------------------- Convert W/m^2 to umol/m^2/s --------------------------------#
   cat0(" + Convert PAR to umol/m2/s.")
   fit.par.in   = eft$par.in   * Watts.2.Ein * 1e6
   fit.par.out  = eft$par.out  * Watts.2.Ein * 1e6
   fit.par.beam = eft$par.beam * Watts.2.Ein * 1e6
   fit.par.diff = eft$par.diff * Watts.2.Ein * 1e6
   #---------------------------------------------------------------------------------------#


   #=======================================================================================#
   # Useful time vectors 
   #=======================================================================================#
   day       = eft$daytime
   night     = eft$nighttime
   notday    = ! eft$daytime
   months    = match(eft$tomonth,unique(eft$tomonth))
   hhmm      = 100 * eft$hour + eft$minu

   when.drya = chron(paste(dry.beg,numyears(eft$when),sep="/"))
   when.dryz = chron(paste(dry.end,numyears(eft$when),sep="/"))
   syear     = numyears(eft$when) + as.numeric(eft$today >= when.dryz)
   dry       = eft$when >= when.drya & eft$when < when.dryz
   #=======================================================================================#



   #=======================================================================================#
   #      nee.raw is the NEE that is not filtered.                                         #
   #=======================================================================================#
   nee.raw = eft$fco2 + eft$storco2
   #=======================================================================================#



   #=======================================================================================#
   #     Find the u* using the t- and F-tests.                                             #
   #---------------------------------------------------------------------------------------#
   ust.alt = with(eft,Ft.ustar(ustar=ustar,cflxca=fco2,cflxst=storco2,nighttime=nighttime))
   ust.alt = round(ust.alt,2)
   #---------------------------------------------------------------------------------------#


   #=======================================================================================#
   #    Make sure the output u* is one of the u* that we are solving.  Also, find the      #
   # number of tries and which one is the favourite for output.                            #
   #=======================================================================================#
   ust.filter = sort(unique(c(ust.output,ust.alt,ust.filter)))
   ubest      = which(ust.filter == ust.output)
   ualt       = which(ust.filter == ust.alt   )
   nustar     = length(ust.filter)
   #=======================================================================================#




   #=======================================================================================#
   #      Find the periods that are considered bad, and flag them.                         #
   #=======================================================================================#
   cat0(" + Find all the blackout periods.")
   fine = is.finite(nee.raw)
   today.fine     = eft$today[fine]
   nfine          = length(today.fine)
   whena.blackout = c(min(eft$today),chron(today.fine+1))
   whenz.blackout = c(chron(today.fine-1),max(eft$today))
   long.blackout  = (whenz.blackout - whena.blackout) > min.blackout
   whena.blackout = whena.blackout[long.blackout]
   whenz.blackout = whenz.blackout[long.blackout]
   n.blackout     = length(whena.blackout)
   sel.bo         = rep(FALSE,times=nwhen)
   if (n.blackout > 0){
      whena.blackout = chron(whena.blackout)
      whenz.blackout = chron(whenz.blackout)
      for (b in sequence(n.blackout)){
         sel         = ( eft$when >= chron(whena.blackout[b])
                       & eft$when <  chron(whenz.blackout[b]) )
         sel.bo[sel] = TRUE
      }#end for
   }#end if
   #=======================================================================================#



   #=======================================================================================#
   #     Initialise the matrices with GPP and NEE as functions of u* filter.               #
   #=======================================================================================#
   NEE       = matrix(  ncol=nustar,nrow=nwhen,dimnames=list(eft$when,ust.filter))
   NEE.turb  = matrix(  ncol=nustar,nrow=nwhen,dimnames=list(eft$when,ust.filter))
   GEE       = matrix(  ncol=nustar,nrow=nwhen,dimnames=list(eft$when,ust.filter))
   GEP       = matrix(  ncol=nustar,nrow=nwhen,dimnames=list(eft$when,ust.filter))
   RECO      = matrix(  ncol=nustar,nrow=nwhen,dimnames=list(eft$when,ust.filter))
   CCNT      = matrix(0,ncol=nustar,nrow=nwhen,dimnames=list(eft$when,ust.filter))
   flag.NEE  = matrix(  ncol=nustar,nrow=nwhen,dimnames=list(eft$when,ust.filter))
   flag.RECO = matrix(  ncol=nustar,nrow=nwhen,dimnames=list(eft$when,ust.filter))
   flag.GEE  = matrix(  ncol=nustar,nrow=nwhen,dimnames=list(eft$when,ust.filter))
   #=======================================================================================#



   #=======================================================================================#
   #     Loop over all u* filters.                                                         #
   #=======================================================================================#
   cat0(" + Apply different u* filters:")
   for (u in sequence(nustar)){
      cat0("   - U* = ",ust.filter[u],"m/s.")

      #====================================================================================#
      #      Select what to discard.                                                       #
      #====================================================================================#
      keep    = eft$ustar %>=% ust.filter[u]
      #====================================================================================#



      #====================================================================================#
      #     NEE: use raw fco2 + filled canopy storage                                      #
      #====================================================================================#
      NEE     [  keep,u] = nee.raw[keep]
      NEE.turb[  keep,u] = nee.raw[keep]
      CCNT    [! keep,u] = 1
      #====================================================================================#



      #====================================================================================#
      #     RECO: use nighttime NEE                                                        #
      #====================================================================================#
      RECO[notday,u]          = NEE[notday,u]
      reco.index              = is.finite(RECO[,u])
      flag.RECO[reco.index,u] = 0
      #====================================================================================#




      #------------------------------------------------------------------------------------#
      #         Create index "ngood" describing how many good hourly values of R           #
      #           have been measured before the current hour                               #
      #------------------------------------------------------------------------------------#
      cat0("     * Find REco running average over ",nrunave," good nighttime hours.")
      ngood         = pmax(0,cumsum(reco.index)-1)


      #----- Ensure averaging doesn't happen over blank period. ---------------------------#
      flux.times = eft$when >= whena & eft$when < whenz & ! sel.bo
      for (b in sequence(n.blackout)){
         sel        = eft$when >= whenz.blackout[b]
         ngood[sel] = ngood[sel] + 2 * nrunave
      }#end for
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #      Using for loop.    
      #------------------------------------------------------------------------------------#
      use          = is.finite(RECO[,u])
      reco.sparse  = RECO [use,u]        # cut down on computing time by 
      ngood.sparse = ngood[use]          #   leaving NA's out of calculation

      reco.runavg  = empty

      for (j in which(flux.times)){
         index          = ( ngood.sparse >= (ngood[j] - nrunave/2) 
                          & ngood.sparse <  (ngood[j] + nrunave/2))
         reco.runavg[j] = mean(reco.sparse[index],na.rm=TRUE)
      }#end for
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #     Print the 5-number statistics of the gap filling.                              #
      #------------------------------------------------------------------------------------#
      fill              = is.na(RECO[,u]) & ! sel.bo
      cat0("     * Five-number summary (RECO): ")
      cat0(c("     > RAW  : ",sprintf("%7.2f",fivenum(RECO       [! fill,u])),"."))
      cat0(c("     > NOT  : ",sprintf("%7.2f",fivenum(reco.runavg[! fill  ])),"."))
      cat0(c("     > ADD  : ",sprintf("%7.2f",fivenum(reco.runavg[  fill  ])),"."))
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #      Gap Fill R with Running Average
      #------------------------------------------------------------------------------------#
      RECO[fill,u]      = reco.runavg[fill]
      flag.RECO[fill,u] = 1
      #------------------------------------------------------------------------------------#
      #====================================================================================#



      #====================================================================================#
      #           GEE                                                                      #
      #====================================================================================#


      #----- Fill in the points in which both NEE and RECO exist. -------------------------#
      GEE[,u] = NEE[,u] - RECO[,u]
      flag.GEE[is.finite(GEE[,u]) | notday,u] = 0
      flag.GEE[is.na(GEE[,u])     &  day  ,u] = 1
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #            vector of seasons to have GEE modeled                                   #
      #------------------------------------------------------------------------------------#
      season      = syear + 0.5 * as.numeric(dry)
      last        = season == yearz+1 & eft$year == yearz
      #------------------------------------------------------------------------------------#


      #----- Temporary fix to allow the end of the last year to be processed. -------------#
      season[last]      = yearz + 0.5          
      fill              = is.na(RECO[,u]) & season == yearz + 0.5
      RECO[fill,u]      = mean (RECO[season== yearz + 0.5,u],na.rm=TRUE)
      flag.RECO[fill,u] = 1
      #------------------------------------------------------------------------------------#


      #----- Select the seasons that have data. -------------------------------------------#
      seasons = unique(season[! sel.bo])
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #------------- Loop over seasons to fit in GEE with light curve ---------------------#
      #------------------------------------------------------------------------------------#
      cat0("     * Fit GEE Light Response Curves.")
      outcurve = file.path(outplace,"par_curves")
      if (! file.exists(outcurve)) dir.create(outcurve)
      for (j in 1:length(seasons)){
         par.d     = empty
         gee.d     = empty
         vpd.d     = empty
         tmp.d     = empty
         gee.model = empty
         index     = day  & season==seasons[j]

         cat0("       > Season: ",seasons[j],".")
         gee.d[index ] = NEE[index,u] - RECO[index,u]
         par.d[index ] = fit.par.in[index ]
         vpd.d[index ] = eft$atm.vpd[index]
         tmp.d[index ] = eft$atm.tmp[index]

         #---------------------------------------------------------------------------------#
         #      Non-linear fitting.                                                        #
         #---------------------------------------------------------------------------------#
         wgts    = 1/((par.d*0.002+3))
         fit.par = nls(gee.d~(a1+(a2*par.d)/(a3+par.d))
                      ,weights = wgts
                      ,start   = list(a1 = 1,a2 = -40, a3 = 500)
                      ,control = list(maxiter=500,minFactor=2^-26) )
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #      Find some of the statistics.                                               #
         #---------------------------------------------------------------------------------#
         summ.par   = summary(fit.par)
         aa         = summ.par$parameters
         coeff.par  = summ.par$parameters[,1]
         pvalue.par = summ.par$parameters[,4]
         #----- Find the adjusted R2. -----------------------------------------------------#
         ss.par    = sum(summ.par$residuals^2)
         df.par    = summ.par$df[2]
         mean.y    = mean(gee.d[index],na.rm=TRUE)
         ss.tot    = sum((gee.d[index]-mean.y)^2,na.rm=TRUE)
         df.tot    = sum(is.finite(gee.d[index])) - 1
         r2.par    = 1.0 - ss.par * df.tot / (ss.tot * df.par)
         pred.par  = predict(object=fit.par
                            ,newdata = data.frame(par.d=fit.par.in,vpd.d=eft$atm.pvap))
         #---------------------------------------------------------------------------------#




         #---------------------------------------------------------------------------------#
         #     Print the statistics of the gap-filling model.                              #
         #---------------------------------------------------------------------------------#
         cat0  ("         # Five-number summary (GEE):")
         cat0(c("           + RAW  : ",sprintf("%7.2f",fivenum(gee.d     [index])),"."))
         cat0(c("           + FIT  : ",sprintf("%7.2f",fivenum(pred.par  [index])),"."))
         cat0  ("         # Model statistics \n")
         cat0(c("           + COEFF   : ",sprintf("%9.3f",coeff.par ),"."))
         cat0(c("           + P-VALUE : ",sprintf("%9.2e",pvalue.par),"."))
         cat0(c("           + R^2     : ",sprintf("%9.3f",r2.par    ),"."))
         #---------------------------------------------------------------------------------#


         #----- Create the text for the model fitting, so we attach it to the plot. -------#
         aa     = sprintf( "%.3f",summ.par$parameters[1,1])
         bb     = sprintf("%.3f" ,summ.par$parameters[2,1])
         cc     = sprintf("%.3f" ,summ.par$parameters[3,1])
         paa    = sprintf("%9.2e",summ.par$parameters[1,4])
         pbb    = sprintf("%9.2e",summ.par$parameters[2,4])
         pcc    = sprintf("%9.2e",summ.par$parameters[3,4])
         rrfit  = sprintf(" %.3f",r2.par                  )
         xyfit  = substitute( expr = y == ( aa + bb * PAR / ( cc + PAR ) ) 
                            , env  = list(aa = aa, bb = bb, cc = cc) )
         r2text = substitute( expr = R ^ 2 == r2   , env = list(r2=rrfit)      )
         px     = substitute( p(a,b,c) == p(paa,pbb,pcc)
                            , env  = list( paa=paa,pbb=pbb,pcc=pcc))
         #---------------------------------------------------------------------------------#


         #---------------------------------------------------------------------------------#
         #------------------- Generate plot within loop for this season -------------------#
         #---------------------------------------------------------------------------------#
         if (plot.par.curves){
            for (o in sequence(nout)){
               #----- Filename. -----------------------------------------------------------#
               lab.ust = sprintf("%.2f",ust.filter[u])
               fichier = file.path(outcurve
                                  ,paste0("par_curve_",lab.ust,"season_",seasons[j],"."
                                         ,outform[o])
                                  )
               dummy   = open.plot( fichier = fichier
                                  , outform = outform[o]
                                  , size    = sq.size
                                  , ptsz    = ptsz
                                  , depth   = depth
                                  )#end open.plot
               #---------------------------------------------------------------------------#


               #---------------------------------------------------------------------------#
               #     Title.                                                                #
               #---------------------------------------------------------------------------#
               if (seasons[j] %% 1 == 0.5){
                  drywet = "Dry"
               }else{
                  drywet = "Wet"
               }#end if
               letitre = paste0(place,"\n U* = ",lab.ust," m/s; "
                              ,"Period: ",drywet," season ",floor(seasons[j]))
               lex     = desc.unit(desc="PAR",unit=untab$umolom2os)
               ley     = desc.unit(desc="GEE",unit=untab$umolom2os)
               #---------------------------------------------------------------------------#


               #---------------------------------------------------------------------------#
               #     Plot the observed and modelled GEE as functions of PAR.               #
               #---------------------------------------------------------------------------#
               xlimit = c(0,2000)
               ylimit = c(-50,20)
               plot(x=xlimit,y=ylimit,main=letitre,xlab=lex,ylab=ley,type="n")
               abline(h=0,v=0,lwd=2,col="black")
               points(x=par.d,y=gee.d,pch=16,col="grey33",cex=0.8)
               points(par.d,pred.par,col="red3",pch=16,cex=0.6)
               #---------------------------------------------------------------------------#


               #---------------------------------------------------------------------------#
               #      Find some good position to plot the fit.                             #
               #---------------------------------------------------------------------------#
               xt = xlimit[1] + 1.0*diff(xlimit); yt = ylimit[1] + 1.00*diff(ylimit)
               text(xt,yt,xyfit,adj=c(1,0.5),cex=0.8)
               xt = xlimit[1] + 1.0*diff(xlimit); yt = ylimit[1] + 0.95*diff(ylimit)
               text(xt,yt,r2text,adj=c(1,0.5),cex=0.8)
               xt = xlimit[1] + 1.0*diff(xlimit); yt = ylimit[1] + 0.90*diff(ylimit)
               text(xt,yt,px,adj=c(1,0.5),cex=0.8)
               #---------------------------------------------------------------------------#


               #---------------------------------------------------------------------------#
               #      Close output device.                                                 #
               #---------------------------------------------------------------------------#
               dummy = close.plot(outform=outform[o])
               #---------------------------------------------------------------------------#
            }#end for
         }#end plotting
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #------------------------ Fill in Missing Daytime GEE ----------------------------#
         #---------------------------------------------------------------------------------#
         fill.index        = day & season==seasons[j] & is.na(GEE[,u])
         GEE[fill.index,u] = pred.par[fill.index]
         fill.index        = notday & season == seasons[j] 
         GEE[fill.index,u] = 0 
         #---------------------------------------------------------------------------------#
      }#end for
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #     Make sure none of the fluxes will be defined during the black out periods.     #
      #  Then find NEE.                                                                    #
      #------------------------------------------------------------------------------------#
      #GEE      [sel.bo,u] = NA
      #RECO     [sel.bo,u] = NA
      #flag.GEE [sel.bo,u] = NA
      #flag.RECO[sel.bo,u] = NA
      NEE [,u] =  GEE[,u] + RECO[,u]
      GEP [,u] = -GEE[,u]
      #------------------------------------------------------------------------------------#


      #====================================================================================#
      #       Create Flags for NEE                                                         #
      #====================================================================================#
      flag.NEE[is.finite(NEE.turb[,u])   , u]  = 0
      flag.NEE[flag.RECO[,u]  ==1 & !day , u]  = 1
      flag.NEE[flag.GEE[,u] == 1         , u]  = 2
      #====================================================================================#


      #------------------------------------------------------------------------------------#
      #            Revert seasons.                                                         #
      #------------------------------------------------------------------------------------#
      season      = syear + 0.5 * as.numeric(dry)
      #------------------------------------------------------------------------------------#
   }#end for
   #=======================================================================================#
   #=======================================================================================#




   #=======================================================================================#
   #    Save variables to the data structure.                                              #
   #=======================================================================================#
   cat0(" + Copy the flux of choice to the structure.")
   eft$nee.raw        = nee.raw
   eft$nee            = NEE      [,ubest]
   eft$nee.turb       = NEE.turb [,ubest]
   eft$nep            = -NEE     [,ubest]
   eft$gpp            = -GEE     [,ubest]
   eft$reco           = RECO     [,ubest]
   eft$gfflg.nee      = flag.NEE [,ubest]
   eft$gfflg.reco     = flag.RECO[,ubest]
   eft$gfflg.gpp      = flag.GEE [,ubest]
   eft$nee.alt        = NEE      [,ualt ]
   eft$nee.turb.alt   = NEE.turb [,ualt ]
   eft$nep.alt        = -NEE     [,ualt ]
   eft$gpp.alt        = -GEE     [,ualt ]
   eft$reco.alt       = RECO     [,ualt ]
   eft$gfflg.nee.alt  = flag.NEE [,ualt ]
   eft$gfflg.reco.alt = flag.RECO[,ualt ]
   eft$gfflg.gpp.alt  = flag.GEE [,ualt ]
   #=======================================================================================#




   #=======================================================================================#
   #    Save variables to the data structure.                                              #
   #=======================================================================================#
   cat0(" + Save all u* tests to a list.")
   ustall            = list()
   ustall$when       = eft$when
   ustall$ust.filter = ust.filter
   ustall$ubest      = ubest
   ustall$ualt       = ualt
   ustall$nee        = NEE
   ustall$nep        = -NEE
   ustall$gpp        = GEP
   ustall$gee        = GEE
   ustall$reco       = RECO
   ustall$ccnt       = CCNT
   ustall$gfflg.nee  = flag.NEE
   ustall$gfflg.gpp  = flag.GEE
   ustall$gfflg.reco = flag.RECO
   #=======================================================================================#





   #=======================================================================================#
   #      Dynamic Season Variable: Occurs at the source/sink transition
   #        i.e. max's and min's in the Cumulative NEE
   #=======================================================================================#
   if (find.dyn.season){
      cat0(" + Find the dynamic seasons.")
      years        = unique(eft$year)
      nyears       = length(years)
      un.seasons   = unique(season)
      n.un.seasons = length(years)

      #------------------- Peaks in Cumulative NEE ----------------------------------------#
      peak.wet = chron(NULL)
      peak.dry = chron(NULL)
      outpeak  = file.path(outplace,"peaks")
      if (! file.exists(outpeak)) dir.create(outpeak)
      for (j in 1:nyears){
         t.start = chron(paste(5,1,years[j]  ,sep="/"))
         t.end   = chron(paste(3,1,years[j]+1,sep="/"))

         index = eft$when > t.start & eft$when < t.end & !is.na(NEE[,ubest])
         if (any(index)){
            sum.nee = cumsum(NEE[index,ubest])
            days    = eft$today[index]

            sum.day = tapply(sum.nee,days,mean)
            tt.day  = unique(days)
            mon     = nummonths(tt.day)

            max.j = tt.day[sum.day==max(sum.day[mon    > 4    ])]
            min.j = tt.day[sum.day==min(sum.day[tt.day > max.j])]
            peak.dry = c(peak.dry,max.j)
            peak.wet = c(peak.wet,min.j)

            if (plot.peaks){
               for (o in sequence(nout)){

                   #----- Filename. -------------------------------------------------------#
                   lab.ust = sprintf("%.2f",ust.filter[u])
                   fichier = file.path(outpeak
                                      ,paste0("peaks_",years[j],"-",iata,".",outform[o])
                   dummy   = open.plot( fichier = fichier
                                      , outform = outform[o]
                                      , size    = sq.size
                                      , ptsz    = ptsz
                                      , depth   = depth
                                      )#end open.plot
                   #-----------------------------------------------------------------------#


                   #-----------------------------------------------------------------------#
                   #     Title.                                                            #
                   #-----------------------------------------------------------------------#
                   letitre = paste0("Peak - ",place," - Year: ",floor(years[j]))
                   lex     = "Day"
                   ley     = "NEE sum"
                   #-----------------------------------------------------------------------#


                   plot(tt.day,sum.day,main = letitre,xlab=lex,ylab=ley)
                   points(min.j, sum.day[tt.day==min.j],col="red3",pch=16,cex=2)
                   points(max.j, sum.day[tt.day==max.j],col="red3",pch=16,cex=2)


                   #-----------------------------------------------------------------------#
                   #      Close output device.                                             #
                   #-----------------------------------------------------------------------#
                   dummy = close.plot(outform=outform[o])
                   #-----------------------------------------------------------------------#
               }#end for
               #---------------------------------------------------------------------------#
            }#end plotting
            #------------------------------------------------------------------------------#
         }else{
            cat0("   - Not enough points to define the season, use default season.")
            peak.dry = c(peak.dry,chron(paste(dry.beg,years[j],sep="/")))
            peak.wet = c(peak.wet,chron(paste(dry.end,years[j],sep="/")))
         }#end if(any(index))
         #---------------------------------------------------------------------------------#
      }#end for
      #------------------------------------------------------------------------------------#



      #------------------------- Build "dyn.seas" -----------------------------------------#
      cat0(" + Build the vector with dynamic seasons.")
      peaks             = sort(c(peak.wet,peak.dry))
      peak.brks         = c(-Inf,peaks)
      cut.when          = cut(eft$today,peak.brks)
      lev.when          = levels(cut.when)
      idx.when          = match(cut.when,lev.when)
      idx.sel           = is.finite(idx.when)
      dyn.seas          = rep(0,nwhen)
      dyn.seas[idx.sel] = un.seasons[idx.when[idx.sel]]
      for (ss in un.seasons){
         cat0("   - Length of season ",ss," = ",sum(dyn.seas == ss),".")
      }#end for
      #====================================================================================#
   }else{
      dyn.seas          = rep(0,nwhen)
   }#end if
   #=======================================================================================#



   #=======================================================================================#
   # Additional Variables
   #=======================================================================================#
   cat0(" + Find some additional variables.")
   # Evapotranspiration
   et = ( eft$fh2o                 # kg/m2/day
        / 24.      )               # output/day
                                   # =  kg/m2/hr = mm/hr
   # Latent Heat
   lambda    = alvli(eft$atm.tmp)      # J/kg
   LHdry     = et * (1/3600) * lambda  # J/m2/s = W/m2

   # fh2o in mmol/m2/s
   fh2o.mmol = et * 1e3 / mmh2o
   #=======================================================================================#


   #=======================================================================================#
   #         Reformat Hourly Flux Table
   #=======================================================================================#
   #---------------------- Reformat the table ---------------------------------------------#
   flux = data.frame(
                when       = eft$when                             # TIMES
               ,month      = eft$month                            #
               ,year       = eft$year                             #
               ,day        = eft$day                              #
               ,hour       = eft$hour                             #
               ,minu       = eft$minu                             #
               ,today      = eft$today                            #
               ,tomonth    = eft$tomonth                          #
               ,doy        = eft$doy                              #
               ,cosz       = eft$cosz                             # SUN MODEL
               ,sunhgt     = eft$sunhgt                           #
               ,par.pot    = eft$par.pot     * Watts.2.Ein * 1.e6 #
               ,nighttime  = eft$nighttime                        #
               ,daytime    = eft$daytime                          #
               ,twilight   = eft$twilight                         #
               ,season                                            #
               ,dyn.seas                                          #
               ,temp       = eft$atm.tmp - t00                    # ATMOSPHERE
               ,prss       = eft$atm.prss                         #
               ,windspd    = eft$atm.vels                         #
               ,winddir    = eft$atm.vdir                         #
               ,tdew       = eft$atm.tdew                         #
               ,pvap       = eft$atm.pvap                         #
               ,vpd        = eft$atm.vpd                          #
               ,tson       = eft$atm.tson                         #
               ,tamb       = eft$atm.tamb                         #
               ,sh         = eft$atm.shv                          #
               ,rh         = eft$atm.rhv                          #
               ,ustar      = eft$ustar                            #
               ,rain       = eft$rain     * 3600.                 # FLUXES
               ,fh2o       = fh2o.mmol                            # 
               ,fheat      = eft$fheat                            #
               ,et                                                #
               ,H          = eft$fsens                            #
               ,LHdry                                             #
               ,fco2       = eft$fco2                             #
               ,fco2.op    = eft$fco2.op                          #
               ,fco2.cp    = eft$fco2.cp                          #
               ,storage    = eft$storco2                          #
               ,co2        = eft$co2                              #
               ,NEE.dbraw  = eft$nee.dbraw                        #
               ,NEE.db     = eft$nee.db                           #
               ,R.db       = eft$reco.db                          #
               ,GEE.db     = -eft$gpp.db                          #
               ,NEE.raw    = nee.raw                              #
               ,NEE.turb   = NEE.turb [,ubest]                    #
               ,NEE        = NEE      [,ubest]                    #
               ,R          = RECO     [,ubest]                    #
               ,GEE        = GEE      [,ubest]                    #
               ,par.in     = eft$par.in      * Watts.2.Ein * 1.e6 # RADIATION
               ,par.out    = eft$par.out     * Watts.2.Ein * 1.e6 #
               ,par.beam   = eft$par.beam    * Watts.2.Ein * 1.e6 #
               ,par.diff   = eft$par.diff    * Watts.2.Ein * 1.e6 #
               ,rshort.in  = eft$rshort.in                        #
               ,rshort.out = eft$rshort.out)                      # 
   #----------------- Save Filling Flags as a Seperate Table ------------------------------#
   flags =  data.frame(
                when       = eft$when                             # TIMES
               ,month      = eft$month                            #
               ,year       = eft$year                             #
               ,day        = eft$day                              #
               ,hour       = eft$hour                             #
               ,minu       = eft$minu                             #
               ,today      = eft$today                            #
               ,tomonth    = eft$tomonth                          #
               ,doy        = eft$doy                              #
               ,temp       = eft$gfflg.atm.tmp                    # ATMOSPHERE
               ,prss       = eft$gfflg.atm.prss                   #
               ,windspd    = eft$gfflg.atm.vels                   #
               ,pvap       = eft$gfflg.atm.pvap                   #
               ,tson       = eft$gfflg.atm.tson                   #
               ,tamb       = eft$gfflg.atm.tamb                   #
               ,ustar      = eft$gfflg.ustar                      #
               ,rain       = eft$gfflg.rain                       # FLUXES
               ,fh2o       = eft$gfflg.fh2o                       # 
               ,fheat      = eft$gfflg.fheat                      #
               ,storage    = eft$gfflg.storco2                    #
               ,NEE        = flag.NEE [,ubest]                    #
               ,R          = flag.RECO[,ubest]                    #
               ,GEE        = flag.GEE [,ubest]                    #
               ,par.in     = eft$gfflg.par.in                     # RADIATION
               ,par.out    = eft$gfflg.par.out)                   #
   #---------------------------------------------------------------------------------------#
   #   Simplify Flags:      1. Objective Analysis                                          #
   #                        2. Harmonic  Analysis                                          #
   #                        3. Filled from relationship w/ another measurement             #
   #                        4. Misc. (variable-specific)                                   #
   #---------------------------------------------------------------------------------------#
   flags[flags==-1] = 3
   #---------------------------------------------------------------------------------------#
   cat0(" + Made Hourly Table.")





   #=======================================================================================#
   #          Create a Monthly data frame 
   #=======================================================================================#
   sel.phcap = flux$par.in > par.phcap.min & flux$par.in < par.phcap.max
   monthly = data.frame(
         year         = tapply(flux$year            ,months        ,mean ,na.rm=T)
        ,month        = tapply(flux$month           ,months        ,mean ,na.rm=T)
        ,tomonth      = tapply(flux$tomonth         ,months        ,mean ,na.rm=T)
        ,temp         = tapply(flux$temp            ,months        ,mean ,na.rm=T)
        ,temp.day     = tapply(flux$temp [day]      ,months[day]   ,mean ,na.rm=T)
        ,temp.night   = tapply(flux$temp [night]    ,months[night] ,mean ,na.rm=T)
        ,prss         = tapply(flux$prss            ,months        ,mean ,na.rm=T)
        ,sh           = tapply(flux$sh              ,months        ,mean ,na.rm=T)
        ,rh           = tapply(flux$rh              ,months        ,mean ,na.rm=T)
        ,vpd          = tapply(flux$vpd             ,months        ,mean ,na.rm=T)
        ,ustar        = tapply(flux$ustar           ,months        ,mean ,na.rm=T)
        ,rain         = tapply(flux$rain            ,months        ,sum  ,na.rm=T)
        ,fh2o         = tapply(flux$fh2o            ,months        ,mean ,na.rm=T)
        ,et           = tapply(flux$et              ,months        ,sum  ,na.rm=T)
        ,fheat        = tapply(flux$fheat           ,months        ,mean ,na.rm=T)
        ,H            = tapply(flux$H               ,months        ,mean ,na.rm=T)
        ,fco2         = tapply(flux$fco2            ,months        ,mean ,na.rm=T)
        ,fco2.op      = tapply(flux$fco2.op         ,months        ,mean ,na.rm=T)
        ,fco2.cp      = tapply(flux$fco2.cp         ,months        ,mean ,na.rm=T)
        ,tson         = tapply(flux$tson            ,months        ,mean ,na.rm=T)
        ,par.day      = tapply(flux$par.in  [day]   ,months[day]   ,mean ,na.rm=T)
        ,par.pot.day  = tapply(flux$par.pot [day]   ,months[day]   ,mean ,na.rm=T)
        ,par.beam.day = tapply(flux$par.beam[day]   ,months[day]   ,mean ,na.rm=T)
        ,par.diff.day = tapply(flux$par.diff[day]   ,months[day]   ,mean ,na.rm=T)   
        ,NEE.turb     = rowMeans(tapply(flux$NEE.turb,list(months,hhmm),mean ,na.rm=T)
                                ,na.rm=T)
        ,NEE          = tapply(flux$NEE             ,months        ,mean ,na.rm=T)
        ,R            = tapply(flux$R       [night] ,months[night] ,mean ,na.rm=T)
        ,GEE          = tapply(flux$GEE     [day]   ,months[day]   ,mean ,na.rm=T)
        ,GEE.turb     = rowMeans(tapply(flux$NEE.turb[day]-flux$R[day]
                                       ,list(months[day],hhmm[day]),mean ,na.rm=T)
                                ,na.rm=T)
        ,NEE.dbraw    = tapply(flux$NEE.dbraw       ,months        ,mean ,na.rm=T)
        ,NEE.db       = tapply(flux$NEE.db          ,months        ,mean ,na.rm=T)
        ,GEE.db       = tapply(flux$GEE.db          ,months        ,mean ,na.rm=T)
        ,R.db         = tapply(flux$R.db            ,months        ,mean ,na.rm=T)
        ,Pc.raw       = tapply(-flux$NEE.turb[sel.phcap],months[sel.phcap],mean ,na.rm=T)
        )#end data.frame

   bo.tomonth = rep(TRUE,times=length(monthly$tomonth))
   if (n.blackout > 0){
      for (b in sequence(n.blackout)){
         sel = monthly$tomonth %wr% c(whena.blackout[b],whenz.blackout[b])
         bo.tomonth[sel] = FALSE
      }#end for
   }#end for
   #monthly$GEE[bo.tomonth] = NA
   #monthly$Pc [bo.tomonth] = NA
   #monthly$et [bo.tomonth] = NA
   cat0(" + Made Monthly Table.")
   #=======================================================================================#






   #=======================================================================================#
   #    Make table for public use                                                          #
   #=======================================================================================#
   public = data.frame(
          GMT           = as.character(sprintf("%.5f",as.numeric(flux$when)))
         ,year          = flux$year
         ,month         = as.character(sprintf("%02.0f"  ,flux$month    ))
         ,day           = as.character(sprintf("%02.0f"  ,flux$day      ))
         ,hour          = as.character(sprintf("%02.0f"  ,flux$hour     ))
         ,minu          = as.character(sprintf("%02.0f"  ,flux$minu     ))
         ,DOY           = as.character(sprintf("%009.5f" ,flux$doy      ))
         ,sunhgt        = as.character(sprintf("%06.3f"  ,flux$sunhgt   ))
         ,daytime       = as.numeric(flux$daytime  )
         ,nighttime     = as.numeric(flux$nighttime)
         ,twilight      = as.numeric(flux$twilight )
         ,season        = as.character(sprintf("%6.1f"   ,flux$season   ))
         ,temp          = as.character(sprintf("%.3f"    ,flux$temp     ))
         ,prss          = as.character(sprintf("%010.3f" ,flux$prss     ))
         ,windspd       = as.character(sprintf("%06.3f"  ,flux$windspd  ))
         ,winddir       = as.character(sprintf("%.3f"    ,flux$winddir  ))
         ,SH            = as.character(sprintf("%.5f"    ,flux$sh       ))
         ,RH            = as.character(sprintf("%.3f"    ,flux$rh       ))
         ,VPD           = as.character(sprintf("%0008.3f",flux$vpd      ))
         ,PAR           = as.character(sprintf("%08.3f"  ,flux$par.in   ))
         ,ustar         = as.character(sprintf("%.3f"    ,flux$ustar    ))
         ,rain          = as.character(sprintf("%.3f"    ,flux$rain     ))
         ,fh2o          = as.character(sprintf("%.3f"    ,flux$fh2o     ))
         ,fheat         = as.character(sprintf("%.3f"    ,flux$fheat    ))
         ,ET            = as.character(sprintf("%.6f"    ,flux$et       ))
         ,H             = as.character(sprintf("%.3f"    ,flux$H        ))
         ,LE            = as.character(sprintf("%.3f"    ,flux$LHdry    ))
         ,fco2          = as.character(sprintf("%.3f"    ,flux$fco2     ))
         ,fco2.op       = as.character(sprintf("%.3f"    ,flux$fco2.op  ))
         ,fco2.cp       = as.character(sprintf("%.3f"    ,flux$fco2.cp  ))
         ,storage       = as.character(sprintf("%.3f"    ,flux$storage  ))
         ,NEE.raw       = as.character(sprintf("%.3f"    ,flux$NEE.raw  ))
         ,NEE.turb      = as.character(sprintf("%.3f"    ,flux$NEE.turb ))
         ,NEE.fill      = as.character(sprintf("%.3f"    ,flux$NEE      ))
         ,R             = as.character(sprintf("%.3f"    ,flux$R        ))
         ,GEE           = as.character(sprintf("%.3f"    ,flux$GEE      ))
         ,NEE.dbraw     = as.character(sprintf("%.3f"    ,flux$NEE.dbraw))
         ,NEE.db        = as.character(sprintf("%.3f"    ,flux$NEE.db   ))
         ,R.db          = as.character(sprintf("%.3f"    ,flux$R.db     ))
         ,GEE.db        = as.character(sprintf("%.3f"    ,flux$GEE.db   ))
         ,flag.temp     = eft$gfflg.atm.tmp
         ,flag.prss     = eft$gfflg.atm.prss
         ,flag.windspd  = eft$gfflg.atm.vels
         ,flag.winddir  = eft$gfflg.atm.vels
         ,flag.SH       = eft$gfflg.atm.pvap
         ,flag.PAR      = eft$gfflg.par.in
         ,flag.ustar    = eft$gfflg.ustar
         ,flag.rain     = eft$gfflg.rain
         ,flag.fh2o     = eft$gfflg.fh2o
         ,flag.fheat    = eft$gfflg.fheat
         ,flag.storage  = eft$gfflg.storco2
         ,flag.NEE.fill = flag.NEE [,ubest]
         ,flag.R        = flag.RECO[,ubest]
         ,flag.GEE      = flag.GEE [,ubest]
         )#end data.frame
   #=======================================================================================#



   #=======================================================================================#
   #    Make public output.                                                                #
   #=======================================================================================#
   #----- Eliminate data that doesn't belong to us nor are public. ------------------------#
   #----- Simplify PAR flag. --------------------------------------------------------------#
   public$flag.PAR [public$flag.PAR==4] =  1
   public = public [flux$year >= yeara.ascii & flux$year <= yearz.ascii,]
   #----- Write the ASCII file. -----------------------------------------------------------#
   asciifile = file.path(here ,paste0(asciipref,".eddyflux.v2." ,yeara.ascii,"."
                                     ,yearz.ascii,".txt")
   write.table(format(public,drop0trailing=F),asciifile,sep=" ",row.names=F,quote=F)
   cat0(" + Saved Hourly ASCII Table.")
   #=======================================================================================#




   #=======================================================================================#
   #    Save the flux variables (MNH style).                                               #
   #=======================================================================================#
   save(list=c("flux","flags","monthly"), file = file.path(here,"flux_filled.RData"))
   cat0(" + Saved R Workspace (MNH).")
   #=======================================================================================#



   #=======================================================================================#
   #    Save the flux variables (MLO style).                                               #
   #=======================================================================================#
   assign(x=iata,value=eft)
   ust.iata = paste("ust",iata,sep=".")
   assign(x=ust.iata,value=ustall)
   save(list=c(iata,ust.iata), file = filled.file)
   cat0(" + Saved R Workspace (MLO).")
   #=======================================================================================#
}#end if (reload && file.exists(filled.file)
#==========================================================================================#
#==========================================================================================#





#==========================================================================================#
#==========================================================================================#
#     Quit the script if the user doesn't want plots.                                      #
#------------------------------------------------------------------------------------------#
if (plot.others){
   #=======================================================================================#
   #=======================================================================================#
   #     Aggregate NEE, GEE, and GEP by month, and compare the differences.                #
   #---------------------------------------------------------------------------------------#
   cat(" + Compute the means.")
   means = list( today   = unique(eft$today)
               , tomonth = unique(eft$tomonth)
               , tohour  = unique(eft$hour + eft$minu / 60)
               )#end list
   hhmm  = 100 * eft$hour + eft$minu
   for (cv in 1:ncflxvar){
      this.name  = cflxvar[[cv]]$name
      this.desc  = cflxvar[[cv]]$desc
      this.unit  = untab[[cflxvar[[cv]]$unit]]
      this.mult  = cflxvar[[cv]]$mult
      this.col   = cflxvar[[cv]]$col

      cat0("   - Variable:",this.name,".")
      this.var                 = get(this.name)

      if (this.name %in% "CCNT"){
         means[[this.name]]$dmean = 100. * ( qapply( X     = this.var
                                                   , DIM   = 1
                                                   , INDEX = eft$today
                                                   , FUN   = sum
                                                   , na.rm = TRUE
                                                   )#end qapply
                                           / qapply( X     = is.finite(this.var)
                                                   , DIM   = 1
                                                   , INDEX = eft$today
                                                   , FUN   = sum
                                                   , na.rm = TRUE
                                                   )#end qapply
                                           ) #end dmean
         means[[this.name]]$emean = 100. * ( qapply( X     = this.var
                                                   , DIM   = 1
                                                   , INDEX = eft$tomonth
                                                   , FUN   = sum
                                                   , na.rm = TRUE
                                                   )#end qapply
                                           / qapply( X     = is.finite(this.var)
                                                   , DIM   = 1
                                                   , INDEX = eft$tomonth
                                                   , FUN   = sum
                                                   , na.rm = TRUE
                                                   )#end qapply
                                           ) #end dmean
         means[[this.name]]$mmean = 100. * ( qapply( X     = this.var
                                                   , DIM   = 1
                                                   , INDEX = eft$month
                                                   , FUN   = sum
                                                   , na.rm = TRUE
                                                   )#end qapply
                                           / qapply( X     = is.finite(this.var)
                                                   , DIM   = 1
                                                   , INDEX = eft$month
                                                   , FUN   = sum
                                                   , na.rm = TRUE
                                                   )#end qapply
                                           ) #end dmean
         means[[this.name]]$qmean = 100. * ( qapply( X     = this.var
                                                   , DIM   = 1
                                                   , INDEX = list(hhmm,eft$month )
                                                   , FUN   = sum
                                                   , na.rm = TRUE
                                                   )#end qapply
                                           / qapply( X     = is.finite(this.var)
                                                   , DIM   = 1
                                                   , INDEX = list(hhmm,eft$month )
                                                   , FUN   = sum
                                                   , na.rm = TRUE
                                                   )#end qapply
                                           ) #end dmean
         means[[this.name]]$gmean = 100. * ( colSums (this.var,na.rm=TRUE)
                                           / colSums (is.finite(this.var),na.rm=TRUE)
                                           )#end gmean
      }else{

         means[[this.name]]$dmean = qapply( X     = this.var
                                          , DIM   = 1
                                          , INDEX = eft$today
                                          , FUN   = mean
                                          , na.rm = TRUE
                                          )#end qapply
         means[[this.name]]$emean = qapply( X     = this.var
                                          , DIM   = 1
                                          , INDEX = eft$tomonth
                                          , FUN   = mean
                                          , na.rm = TRUE
                                          )#end qapply
         means[[this.name]]$mmean = qapply( X     = this.var
                                          , DIM   = 1
                                          , INDEX = eft$month
                                          , FUN   = mean
                                          , na.rm = TRUE
                                          )#end qapply
         means[[this.name]]$qmean = qapply( X     = this.var
                                          , DIM   = 1
                                          , INDEX = list(hhmm,eft$month )
                                          , FUN   = mean
                                          , na.rm = TRUE
                                          )#end qapply
         means[[this.name]]$gmean = colMeans(this.var,na.rm=TRUE)
      }#end if

   }#end for
   #=======================================================================================#
   #=======================================================================================#





   #=======================================================================================#
   #=======================================================================================#
   #       Now we make a plot similar to figure 3 of Bonal et al. (2008).                  #
   #---------------------------------------------------------------------------------------#
   cat0(" + Find the fluxes as functions of u* bins.")

   outustar = file.path(outplace,"ustar_filter")
   if (! file.exists(outustar)) dir.create(outustar)

   #---------------------------------------------------------------------------------------#
   #     Make sure that we don't use any gap filled data.                                  #
   #---------------------------------------------------------------------------------------#
   sel      = ( eft$gfflg.fco2  %==% 0 & eft$gfflg.storco2 %==% 0
              & eft$gfflg.ustar %==% 0 & eft$nighttime            )
   ustar    = eft$ustar   [sel]
   fco2     = eft$fco2    [sel]
   storco2  = eft$storco2 [sel]
   nee      = fco2 + storco2
   #---------------------------------------------------------------------------------------#



   #----- Split the data into u* bins. ----------------------------------------------------#
   cat0("   - Split u* into bins.")
   quant.breaks  = seq(from=0.0,to=1.0,by=0.05)
   ustar.breaks  = quantile(ustar,probs=quant.breaks,na.rm=TRUE)
   sel           = ! duplicated(ustar.breaks)
   ustar.breaks  = ustar.breaks[sel]
   quant.breaks  = quant.breaks[sel]
   quant.midst   = mid.points(quant.breaks)
   ustar.midst   = quantile(ustar,probs=quant.midst,na.rm=TRUE)
   n.ustar      = length(ustar.breaks) - 1 
   ustar.cut    = cut(ustar, breaks = ustar.breaks )
   #---------------------------------------------------------------------------------------#




   #----- Find the mean fluxes for each u* bin. -------------------------------------------#
   cat0("   - Find the mean fluxes for each u* bin.")
   nee.ust     = tapply(X = nee    , INDEX = ustar.cut, FUN = boot.six.summary, R = n.boot)
   cflxca.ust  = tapply(X = fco2   , INDEX = ustar.cut, FUN = boot.six.summary, R = n.boot)
   cflxst.ust  = tapply(X = storco2, INDEX = ustar.cut, FUN = boot.six.summary, R = n.boot)
   nee.ust     = sapply(X = nee.ust   , FUN = c)
   cflxca.ust  = sapply(X = cflxca.ust, FUN = c)
   cflxst.ust  = sapply(X = cflxst.ust, FUN = c)
   k           = match(c("expected","ci.lower","ci.upper"),six.summary.names)
   nee.ust     = nee.ust   [k,]
   cflxca.ust  = cflxca.ust[k,]
   cflxst.ust  = cflxst.ust[k,]

   ustar.step  = rep(ustar.breaks,each=2)[-c(1,2*(n.ustar+1))]
   nee.step    = t(apply(X=nee.ust   ,MARGIN=1,FUN=rep,each=2))
   cflxca.step = t(apply(X=cflxca.ust,MARGIN=1,FUN=rep,each=2))
   cflxst.step = t(apply(X=cflxst.ust,MARGIN=1,FUN=rep,each=2))
   #---------------------------------------------------------------------------------------#


   #----- Find the plot limits. -----------------------------------------------------------#
   cat0("   - Find the plot ranges.")
   if (is.null(ust.xlim)){
      xlim = range(ustar.breaks)
   }else{
      xlim = ust.xlim
   }#end if
   ylim = range(c(nee.ust,cflxca.ust,cflxst.ust),na.rm=TRUE)
   if (is.null(flux.ylim)){
      ylim = pretty.xylim(u=ylim,fracexp=0.0,is.log=FALSE)
   }else{
      ylim = flux.ylim
   }#end if
   #---------------------------------------------------------------------------------------#


   #----- Plot annotation. ----------------------------------------------------------------#
   letitre = paste0("Average night time fluxes by u* quantiles - ",place)
   lex     = desc.unit(desc="u* class",unit=untab$mos)
   ley     = desc.unit(desc="Flux",unit=untab$umolcom2os)
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #    Plot the average night time fluxes as a function of u*.                            #
   #---------------------------------------------------------------------------------------#
   cat0(" + Plot the averages.")
   for (o in 1:nout){
      #----- File name and format. --------------------------------------------------------#
      fichier = file.path(outustar,paste0(iata,"-nightflux_by_ustar.",outform[o]))
      dummy   = open.plot( fichier = fichier
                         , outform = outform[o]
                         , size    = ey.size
                         , ptsz    = ptsz
                         , depth   = depth
                         )#end open.plot
      #------------------------------------------------------------------------------------#



      #----- Split area into two. ---------------------------------------------------------#
      par(par.user)
      layout(mat=rbind(2,1),heights=c(1.-f.leg,f.leg))
      #------------------------------------------------------------------------------------#




      #------------------------------------------------------------------------------------#
      #     Plot legend.                                                                   #
      #------------------------------------------------------------------------------------#
      par(mar=c(0.1,0.1,0.1,0.1))
      plot.new()
      plot.window(xlim=c(0,1),ylim=c(0,1))
      legend( x       = "bottom"
            , inset   = 0.0
            , legend  = parse( text = c( "N*E*E"
                                       , "C*O[2]*phantom(1)*f*l*u*x"
                                       , "C*O[2]*phantom(1)*s*t*o*r*a*g*e"
                                       , "L*i*t*e*r*a*t*u*r*e*phantom(1)*u[m*i*n]^symbol(\"\\052\")"
                                       , "A*l*t*e*r*n*a*t*i*v*e*phantom(1)*u[m*i*n]^symbol(\"\\052\")"
                                       )#end c
                             )#end parse
            , fill    = c("#A3CC52","#3B24B3","#990F0F","transparent","transparent")
            , border  = c("#A3CC52","#3B24B3","#990F0F","transparent","transparent")
            , angle   = c(-45,45,90,0, 0)
            , density = c(20,20,20,20,20)
            , col     = c("#4B6614","#160959","#4D0404","violetred","dodgerblue")
            , pch     = c(rep(16,3),NA,NA)
            , lty     = c(rep("solid",3),"dotdash","longdash")
            , lwd     = 2
            , ncol    = 2
            , cex     = 0.7 * cex.ptsz
            , xpd     = TRUE
            , bty     = "n"
            )#end if
      #------------------------------------------------------------------------------------#





      #------------------------------------------------------------------------------------#
      #     Plot the averaged fluxes.                                                      #
      #------------------------------------------------------------------------------------#
      par(mar=par.user$mar)
      plot.new()
      plot.window(xlim=xlim,ylim=ylim)
      #----- Open the plot window. --------------------------------------------------------#
      axis(side=1)
      axis(side=2,las=1)
      title(xlab=lex,ylab=ley,main=letitre)
      #----- Background grid. -------------------------------------------------------------#
      grid(col=grid.colour,lty="solid")
      #----- Plot the u* threshold. -------------------------------------------------------#
      abline(v=ust.output,col="violetred" ,lty="dotdash" ,lwd=2)
      abline(v=ust.alt   ,col="dodgerblue",lty="longdash",lwd=2)
      #----- Plot the binned averages. ----------------------------------------------------#
      epolygon( x       = c(ustar.step,rev(ustar.step))
              , y       = c(cflxca.step[2,],rev(cflxca.step[3,]))
              , col     = "#3B24B3"
              , angle   = 45
              , density = 20
              )#end epolygon
      epolygon( x       = c(ustar.step,rev(ustar.step))
              , y       = c(cflxst.step[2,],rev(cflxst.step[3,]))
              , col     = "#990F0F"
              , angle   = 90
              , density = 20
              )#end epolygon
      epolygon( x       = c(ustar.step,rev(ustar.step))
              , y       = c(nee.step[2,],rev(nee.step[3,]))
              , col     = "#A3CC52"
              , angle   = -45
              , density = 20
              )#end epolygon
      lines (x = ustar.step,y = cflxca.step[1,], type="l",col="#160959",lwd=2.0)
      lines (x = ustar.step,y = cflxst.step[1,], type="l",col="#4D0404",lwd=2.0)
      lines (x = ustar.step,y = nee.step   [1,], type="l",col="#4B6614",lwd=2.0)
      points(x = ustar.midst,y = cflxca.ust[1,], type="p",col="#160959",pch=16 ,cex=0.8)
      points(x = ustar.midst,y = cflxst.ust[1,], type="p",col="#4D0404",pch=16 ,cex=0.8)
      points(x = ustar.midst,y = nee.ust   [1,], type="p",col="#4B6614",pch=16 ,cex=0.8)
      box()
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #      Close the device.                                                             #
      #------------------------------------------------------------------------------------#
      dummy = close.plot(outform=outform[o])
      #------------------------------------------------------------------------------------#
   }#end for
   #---------------------------------------------------------------------------------------#
   #=======================================================================================#
   #=======================================================================================#






   #=======================================================================================#
   #=======================================================================================#
   #       Similar to the figure above, but plotting the entire distribution.              #
   #---------------------------------------------------------------------------------------#
   cat0(" + Find the box plots as functions of u* bins.")

   outustar = file.path(outplace,"ustar_filter")
   if (! file.exists(outustar)) dir.create(outustar)

   #---------------------------------------------------------------------------------------#
   #     Make sure that we don't use any gap filled data.                                  #
   #---------------------------------------------------------------------------------------#
   sel            = ( ( is.finite(eft$gfflg.fco2)    & eft$gfflg.fco2    == 0 )
                    & ( is.finite(eft$gfflg.storco2) & eft$gfflg.storco2 == 0 )
                    & ( is.finite(eft$gfflg.ustar  ) & eft$gfflg.ustar   == 0 )
                    & eft$nighttime          )
   ustar    = eft$ustar   [sel]
   fco2     = eft$fco2    [sel]
   storco2  = eft$storco2 [sel]
   nee      = fco2 + storco2
   #---------------------------------------------------------------------------------------#



   #----- Split the data into u* bins. ----------------------------------------------------#
   cat0("   - Split u* into bins.")
   ustar.breaks = quantile(ustar,probs=seq(from=0.0,to=1.0,by=0.10))
   sel           = ! duplicated(ustar.breaks)
   ustar.breaks  = ustar.breaks[sel]
   n.ustar      = length(ustar.breaks) - 1 
   ustar.mid    = ustar.breaks[1:n.ustar] + 0.5 * diff(ustar.breaks)
   ustar.cut    = cut(ustar, breaks = ustar.breaks )
   ustar.lev    = levels(ustar.cut)
   ustar.idx    = match(ustar.cut,ustar.lev)
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #      Create the list from which we will draw the box plots.                           #
   #---------------------------------------------------------------------------------------#
   bplot = list()
   for (u in sequence(n.ustar)){
      sel            = ustar.idx == u 
      off            = (u-1) * 3
      bplot[[1+off]] = nee    [sel]
      bplot[[2+off]] = fco2   [sel]
      bplot[[3+off]] = storco2[sel]
   }#end for
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #      Assign the colours.                                                              #
   #---------------------------------------------------------------------------------------#
   bpcols = rep(c("chartreuse","steelblue","gold"),times=n.ustar)
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #      Define the labels for the x axis.                                                #
   #---------------------------------------------------------------------------------------#
   xat    = seq(from=0,to=3*n.ustar,by=3) + 0.5
   xlabel = sprintf("%.2f",ustar.breaks)
   xlimit = range(xat)
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #      Define the range for the y axis, and expand it to fit the legend.                #
   #---------------------------------------------------------------------------------------#
   ylimit = range(sapply(X=bplot,FUN=range,simplify=TRUE,na.rm=TRUE),na.rm=TRUE)
   ylimit = pretty.xylim(u=ylimit,fracexp=0.0,is.log=FALSE)
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #    Plot the average night time fluxes as a function of u*.                            #
   #---------------------------------------------------------------------------------------#
   for (o in sequence(nout)){
      #----- File name and format. --------------------------------------------------------#
      fichier = paste(outustar,paste0("boxplot-night.nee-",iata,".",outform[o]))
      dummy   = open.plot( fichier = fichier
                         , outform = outform[o]
                         , size    = ey.size
                         , ptsz    = ptsz
                         , depth   = depth
                         )#end open.plot
      #------------------------------------------------------------------------------------#



      #----- Split device into two. -------------------------------------------------------#
      par(par.user)
      layout(mat=rbind(2,1),heights=c(1.-f.leg,f.leg))
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #     Plot legend.                                                                   #
      #------------------------------------------------------------------------------------#
      par(mar=c(0.1,0.1,0.1,0.1))
      plot.new()
      plot.window(xlim=c(0,1),ylim=c(0,1))
      legend( x      = "bottom"
            , inset  = 0.0
            , legend = c("NEE","CO2 flux","CO2 storage")
            , fill   = c("chartreuse","steelblue","gold")
            , xpd    = TRUE
            , bty    = "n"
            )#end legend
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #     Plot the box plots.                                                            #
      #------------------------------------------------------------------------------------#
      par(mar=par.user$mar)
      #----- Open the plot window. --------------------------------------------------------#
      plot.new()
      plot.window(xlim=xlimit,ylim=ylimit)
      title(xlab=desc.unit(desc="u* intervals",unit=untab$mos)
           ,ylab=desc.unit(desc="Fluxes",unit=untab$umolcom2os)
           ,main=paste("Boxplot of night time fluxes by u* quantiles",place,sep=" - ")
           )#end title
      #----- X axis. ----------------------------------------------------------------------#
      axis(side=1,at=xat,labels=xlabel)
      axis(side=2,las=1)
      #----- Background grid. -------------------------------------------------------------#
      abline(h=axTicks(2),v=xat,col=grid.colour,lty="solid")
      #----- Plot the binned averages. ----------------------------------------------------#
      boxplot(x=bplot,col=bpcols,add=TRUE,show.names=FALSE)
      box()
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #      Close the device.                                                             #
      #------------------------------------------------------------------------------------#
      dummy = close.plot(outform=outform[o])
      #------------------------------------------------------------------------------------#
   }#end for
   #---------------------------------------------------------------------------------------#
   #=======================================================================================#
   #=======================================================================================#







   #=======================================================================================#
   #=======================================================================================#
   #     Plot the variables as functions of the ustar filter and the month.                #
   #---------------------------------------------------------------------------------------#
   cat0(" + Plot global means.")
   outustar = file.path(outplace,"ustar_filter")
   if (! file.exists(outustar)) dir.create(outustar)
   ylimit = 0    # To make sure zero is included
   leglab = NULL
   legcol = NULL
   cvuse  = NULL
   for (cv in 1:ncflxvar){
      this.name  = cflxvar[[cv]]$name
      if (this.name != "CCNT"){
         ylimit     = range(c(ylimit,means[[this.name]]$gmean),na.rm=TRUE)
         leglab     = c(leglab,cflxvar[[cv]]$name)
         legcol     = c(legcol,cflxvar[[cv]]$col)
         cvuse      = c(cvuse,cv)
      }#end if
   }#end for
   #---------------------------------------------------------------------------------------#


   #----- Make y axis longer so it fits a legend. -----------------------------------------#
   ylimit = pretty.xylim(u=ylimit,fracexp=0.0,is.log=FALSE)
   #---------------------------------------------------------------------------------------#


   #----- Title and axes. -----------------------------------------------------------------#
   letitre = paste("Mean carbon balance - ",place,"\n Period : ",yeara,"-",yearz,sep="")
   lex     = desc.unit(desc="U* thresold",unit=untab$mos)
   ley     = desc.unit(desc="Mean Carbon fluxes",unit=untab$umolcom2os)
   #---------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   #     Plot the graphs.                                                                  #
   #---------------------------------------------------------------------------------------#
   for (o in sequence(nout)){
      #----- Filename. --------------------------------------------------------------------#
      fichier = file.path(outustar,paste0("gmean_cflx_",iata,".",outform[o]))
      dummy   = open.plot( fichier = fichier
                         , outform = outform[o]
                         , size    = ey.size
                         , ptsz    = ptsz
                         , depth   = depth
                         )#end open.plot
      #------------------------------------------------------------------------------------#



      #----- Split device. ----------------------------------------------------------------#
      par(par.user)
      layout(mat=rbind(2,1),heights=c(1.-f.leg,f.leg))
      #------------------------------------------------------------------------------------#



      #------ Plot the legend. ------------------------------------------------------------#
      par(mar=c(0.1,0.1,0.1,0.1))
      plot.new()
      plot.window(xlim=c(0,1),ylim=c(0,1))
      legend(x="bottom",inset=0.0,legend=leglab,col=legcol,lwd=2,pch=16,bty="n",xpd=TRUE)
      #------------------------------------------------------------------------------------#




      #----- Plot the main frame. ---------------------------------------------------------#
      par(mar=par.user$mar)
      plot.new()
      plot.window(xlim=range(ust.filter),ylim=ylimit)
      axis(side=1)
      axis(side=2,las=1)
      title(main=letitre,xlab=lex,ylab=ley,ylim=ylimit)
      grid(col=grid.colour,lty="solid")
      abline(h=0,v=ust.output,col="black",lwd=2,lty="dotdash")
      for (cv in cvuse){
         this.name = cflxvar[[cv]]$name
         lines(x=ust.filter,y=means[[this.name]]$gmean,type="o",col=cflxvar[[cv]]$col
              ,pch=16,lwd=2)
      }#end for
      box()
      #------------------------------------------------------------------------------------#




      #------------------------------------------------------------------------------------#
      #      Close output device.                                                          #
      #------------------------------------------------------------------------------------#
      dummy = close.plot(outform=outform[o])
      #------------------------------------------------------------------------------------#
   }#end for
   #---------------------------------------------------------------------------------------#

   #=======================================================================================#
   #=======================================================================================#






   #=======================================================================================#
   #=======================================================================================#
   #     Plot the monthly season means as a function of the u* filter.                     #
   #---------------------------------------------------------------------------------------#
   cat0(" + Plot monthly means.")
   outustar = file.path(outplace,"ustar_filter")
   if (! file.exists(outustar)) dir.create(outustar)

   #---------------------------------------------------------------------------------------#
   #    Find the time to plot.                                                             #
   #---------------------------------------------------------------------------------------#
   mon    = c(0.5,sequence(12),12.5)
   monat  = sequence(12)
   monlab = substring(month.abb,1,1)
   nmon   = length(mon)
   mon2d  = matrix(data=mon       ,nrow=nmon,ncol=nustar,byrow=FALSE)
   ust2d  = matrix(data=ust.filter,nrow=nmon,ncol=nustar,byrow=TRUE )
   ust.at = pretty(ust.filter,n=10)
   #---------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   #     Plot defaults.                                                                    #
   #---------------------------------------------------------------------------------------#
   x.axis.options = list(side=1,at=monat,labels=monlab)
   y.axis.options = list(side=2,at=ust.at)
   plot.after     = list( abline = list(v= monat,h=ust.at,col=grid.colour,lty="dotted")
                        , abline = list(h=ust.output,col="black",lwd=2,lty="dotdash")
                        )#end list
   #---------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   #     Plot the graphs.                                                                  #
   #---------------------------------------------------------------------------------------#
   for (cv in sequence(ncflxvar)){
      #------------------------------------------------------------------------------------#
      this.name  = cflxvar$name[cv]
      this.desc  = cflxvar$desc[cv]
      this.unit  = untab[[cflxvar$unit[cv]]]
      this.mult  = cflxvar$mult[cv]
      this.col   = cflxvar$col [cv]
      if (is.na(this.mult)) this.mult = 1. else this.mult = eval(parse(text=this.mult))
      cat0("   - ",this.name,".")
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #    Get the variable, and repeat Jan and Dec so we can see these months.            #
      #------------------------------------------------------------------------------------#
      this.var = means[[this.name]]$mmean * this.mult
      this.var = rbind(this.var[12,],this.var,this.var[1,])
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #     Plot the graphs.                                                               #
      #------------------------------------------------------------------------------------#
      for (o in sequence(nout)){
         #----- Filename. -----------------------------------------------------------------#
         fichier = file.path(outustar
                            ,paste0("mmean_",this.name,"_ustar-",iata,".",outform[o])
                            )#end file.path
         dummy   = open.plot( fichier = fichier
                            , outform = outform[o]
                            , size    = ex.size
                            , ptsz    = ptsz
                            , depth   = depth
                            )#end open.plot
         #---------------------------------------------------------------------------------#


         #----- Title and axes. -----------------------------------------------------------#
         letitre = paste0(this.desc," - ",place,"\n Period : ",yeara,"-",yearz)
         lacle   = desc.unit(desc=NULL,unit=this.unit)
         lex     = desc="Month"
         ley     = desc.unit(desc="u* threshold",unit=untab$mos)
         #---------------------------------------------------------------------------------#


         #----- Plot the main frame. ------------------------------------------------------#
         par(par.user)
         image.map( x              = mon2d
                  , y              = ust2d
                  , z              = this.var
                  , colour.palette = panoply
                  , nlevels        = ncolours.im
                  , na.col         = miss.colour
                  , x.axis.options = x.axis.options
                  , y.axis.options = y.axis.options
                  , main.title     = list(main=letitre,cex.main=cex.main)
                  , key.title      = list(main=lacle,cex.main=0.8*cex.main,line=1)
                  , key.log        = FALSE
                  , plot.after     = plot.after
                  , f.key          = f.leg
                  )#end image.map
         #---------------------------------------------------------------------------------#


         #---------------------------------------------------------------------------------#
         #      Close output device.                                                       #
         #---------------------------------------------------------------------------------#
         dummy = close.plot(outform=outform[o])
         #---------------------------------------------------------------------------------#
      }#end for
      #------------------------------------------------------------------------------------#
   }#end for
   #---------------------------------------------------------------------------------------#

   #=======================================================================================#
   #=======================================================================================#






   #=======================================================================================#
   #=======================================================================================#
   #     Plot the monthly season means as a function of the u* filter.                     #
   #---------------------------------------------------------------------------------------#
   cat0(" + Plotting time series of the monthly means.")
   outustar = file.path(outplace,"ustar_filter")
   if (! file.exists(outustar)) dir.create(outustar)

   #---------------------------------------------------------------------------------------#
   #    Find the time to plot.                                                             #
   #---------------------------------------------------------------------------------------#
   whenplot = pretty.time(means$tomonth,n=6)
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #    Find the time to plot.                                                             #
   #---------------------------------------------------------------------------------------#
   nmon   = length(means$tomonth)
   mon2d  = matrix(data=means$tomonth,nrow=nmon,ncol=nustar,byrow=FALSE)
   ust2d  = matrix(data=ust.filter   ,nrow=nmon,ncol=nustar,byrow=TRUE )
   ust.at = pretty(ust.filter,n=10)
   #---------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   #     Plot defaults.                                                                    #
   #---------------------------------------------------------------------------------------#
   x.axis.options = list(side=1,at=whenplot$levels,labels=whenplot$labels
                        ,padj=whenplot$padj)
   y.axis.options = list(side=2,at=ust.at)
   plot.after     = list( abline = list(v=whenplot$levels,h=ust.at,col=grid.colour
                                       ,lty="dotted")
                        , abline = list(h=ust.output,col="black",lwd=2,lty="dotdash")
                        )#end list
   #---------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   #     Plot the graphs.                                                                  #
   #---------------------------------------------------------------------------------------#
   for (cv in sequence(ncflxvar)){
      #------------------------------------------------------------------------------------#
      this.name  = cflxvar$name[cv]
      this.desc  = cflxvar$desc[cv]
      this.unit  = untab[[cflxvar$unit[cv]]]
      this.col   = cflxvar$col [cv]
      this.mult  = cflxvar$mult[cv]
      cat0("   - ",this.name,".")
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #    Get the variable.                                                               #
      #------------------------------------------------------------------------------------#
      this.var = means[[this.name]]$emean * this.mult
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #     Plot the graphs.                                                               #
      #------------------------------------------------------------------------------------#
      for (o in sequence(nout)){
         #----- Filename. -----------------------------------------------------------------#
         fichier = file.path(outustar
                            ,paste0("emean_",this.name,"_ustar-",iata,".",outform[o])
                            )#end file.path
         dummy   = open.plot( fichier = fichier
                            , outform = outform[o]
                            , size    = ex.size
                            , ptsz    = ptsz
                            , depth   = depth
                            )#end open.plot
         #---------------------------------------------------------------------------------#


         #----- Title and axes. -----------------------------------------------------------#
         letitre = paste0(this.desc," - ",place)
         lacle   = desc.unit(desc=NULL,unit=this.unit)
         lex     = desc="Month"
         ley     = desc.unit(desc="u* threshold",unit=untab$mos)
         #---------------------------------------------------------------------------------#





         #----- Plot the main frame. ------------------------------------------------------#
         par(par.user)
         image.map( x              = mon2d
                  , y              = ust2d
                  , z              = this.var
                  , colour.palette = panoply
                  , nlevels        = ncolours.im
                  , na.col         = miss.colour
                  , x.axis.options = x.axis.options
                  , y.axis.options = y.axis.options
                  , main.title     = list(main=letitre,cex.main=cex.main)
                  , key.title      = list(main=lacle,cex.main=0.8*cex.main,line=1)
                  , key.log        = FALSE
                  , plot.after     = plot.after
                  , f.key          = f.leg
                  )#end image.map
         #---------------------------------------------------------------------------------#


         #---------------------------------------------------------------------------------#
         #      Close output device.                                                       #
         #---------------------------------------------------------------------------------#
         dummy = close.plot(outform=outform[o])
         #---------------------------------------------------------------------------------#
      }#end for
      #------------------------------------------------------------------------------------#
   }#end for
   #---------------------------------------------------------------------------------------#

   #=======================================================================================#
   #=======================================================================================#






   #=======================================================================================#
   #=======================================================================================#
   #     Plot the monthly season means as a function of the u* filter.                     #
   #---------------------------------------------------------------------------------------#
   cat0(" + Plot monthly mean of the diurnal cycle.")
   outqmean = file.path(outplace,"qmean")
   if (! file.exists(outqmean)) dir.create(outqmean)

   #---------------------------------------------------------------------------------------#
   #    Find the time to plot.                                                             #
   #---------------------------------------------------------------------------------------#
   hhmm   = means$tohour
   nhhmm  = length(hhmm)
   dhhmm  = median(diff(hhmm))
   hhmm   = c(hhmm[1]-dhhmm,hhmm,hhmm[nhhmm]+dhhmm)
   nhhmm  = length(hhmm)
   hhmmat = seq(from=0,to=24,by=3)
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #    Find the time to plot.                                                             #
   #---------------------------------------------------------------------------------------#
   hhmm2d = matrix(data=hhmm      ,nrow=nhhmm,ncol=nustar,byrow=FALSE)
   ust2d  = matrix(data=ust.filter,nrow=nhhmm,ncol=nustar,byrow=TRUE )
   ust.at = pretty(ust.filter,n=10)
   #---------------------------------------------------------------------------------------#




   #---------------------------------------------------------------------------------------#
   #     Plot defaults.                                                                    #
   #---------------------------------------------------------------------------------------#
   x.axis.options = list(side=1,at=hhmmat,labels=hhmmat)
   y.axis.options = list(side=2,at=ust.at)
   plot.after     = list( abline = list(v=hhmmat,h=ust.at,col=grid.colour,lty="dotted")
                        , abline = list(h=ust.output,col="black",lwd=2,lty="dotdash")
                        )#end list
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #     Plot the graphs.                                                                  #
   #---------------------------------------------------------------------------------------#
   for (cv in sequence(ncflxvar)){
      #------------------------------------------------------------------------------------#
      this.name  = cflxvar$name[cv]
      this.desc  = cflxvar$desc[cv]
      this.unit  = untab[[cflxvar$unit[cv]]]
      this.col   = cflxvar$col [cv]
      this.mult  = cflxvar$mult[cv]
      cat0("   - ",this.name,".")
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #     Output for this variable.                                                      #
      #------------------------------------------------------------------------------------#
      outvar = file.path(outqmean,this.name)
      if (! file.exists(outvar)) dir.create(outvar)
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #    Get the variable.                                                               #
      #------------------------------------------------------------------------------------#
      this.var = means[[this.name]]$qmean * this.mult
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #    Find the range (keep the same scale for all months).                            #
      #------------------------------------------------------------------------------------#
      zlim     = range(this.var,na.rm=TRUE)
      zlevels  = pretty(zlim,n=200)
      nzlevels = length(zlevels)
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #    Loop over all months.                                                           #
      #------------------------------------------------------------------------------------#
      for (m in sequence(12)){
         cmon = sprintf("%2.2i",m)
         cat0("     * ",mlist[m],".")


         #---------------------------------------------------------------------------------#
         #     Select the current month, then repeat the first and the last hour so we can #
         # see them.                                                                       #
         #---------------------------------------------------------------------------------#
         mon.var = this.var[,m,]
         mon.var = rbind(mon.var[dim(mon.var)[1],],mon.var,mon.var[1,])
         #---------------------------------------------------------------------------------#


         #---------------------------------------------------------------------------------#
         #     Plot the graphs.                                                            #
         #---------------------------------------------------------------------------------#
         for (o in sequence(nout)){
            #----- Filename. --------------------------------------------------------------#
            fichier = file.path( outvar
                               , paste0("qm_",this.name,"_ustar-month_",cmon,"-",iata,"."
                                       ,outform[o])
                               )#end file.path
            dummy   = open.plot( fichier = fichier
                               , outform = outform[o]
                               , size    = ex.size
                               , ptsz    = ptsz
                               , depth   = depth
                               )#end open.plot
            #------------------------------------------------------------------------------#


            #----- Title and axes. --------------------------------------------------------#
            letitre = paste0(this.desc," - ",place," - ",mlist[m]
                            ,"\n Period : ",yeara,"-",yearz)
            lacle   = desc.unit(desc=NULL,unit=this.unit)
            lex     = desc.unit(desc="Time",unit=untab$utc)
            ley     = desc.unit(desc="u* threshold",unit=untab$mos)
            #------------------------------------------------------------------------------#



            #----- Plot the main frame. ---------------------------------------------------#
            par(par.user)
            image.map( x              = hhmm2d
                     , y              = ust2d
                     , z              = mon.var
                     , colour.palette = panoply
                     , nlevels        = ncolours.im
                     , na.col         = miss.colour
                     , x.axis.options = x.axis.options
                     , y.axis.options = y.axis.options
                     , main.title     = list(main=letitre,cex.main=cex.main)
                     , key.title      = list(main=lacle,cex.main=0.8*cex.main,line=1)
                     , key.log        = FALSE
                     , plot.after     = plot.after
                     , f.key          = key.frac
                     )#end image.map
            #------------------------------------------------------------------------------#


            #------------------------------------------------------------------------------#
            #      Close output device.                                                    #
            #------------------------------------------------------------------------------#
            dummy = close.plot(outform=outform[o])
            #------------------------------------------------------------------------------#
         }#end for
         #---------------------------------------------------------------------------------#
      }#end for
      #------------------------------------------------------------------------------------#
   }#end for
   #---------------------------------------------------------------------------------------#
}#end if
#==========================================================================================#
#==========================================================================================#
