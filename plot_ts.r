#----- Leave this as the first command, this will reset your R session. -------------------#
rm(list=ls())
graphics.off()
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
here        = getwd()                   # Current path 
srcdir      = file.path(here,"Rsc"    ) # Path with additional scripts
ibackground = 0 
outroot     = file.path(here,"figures")
#------------------------------------------------------------------------------------------#




#------------------------------------------------------------------------------------------#
#     Long name, for plotting titles.                                                      #
#------------------------------------------------------------------------------------------#
longname = "Tanguro (Control), MT"  # Nice name for the site
iata     = "tb0"                    # 3-letter/number code used to identify this site.
yeara    = 2008                     # First year to use
yearz    = 2016                     # Last year to use
lat      = -13.082                  # Latitude (degrees north)
lon      = -52.376                  # Longitude (degrees east)
alt      = 377.                     # Altitude  (metres above sea level)
#------------------------------------------------------------------------------------------#




#------------------------------------------------------------------------------------------#
#     Session information, for saving/loading data already read/analysed.                  #
#------------------------------------------------------------------------------------------#
gflevel        = "metfill"
eft.rdata      = paste0(gflevel,"_",iata,"_",yeara,"-",yearz,".RData")
#------------------------------------------------------------------------------------------#




#------------------------------------------------------------------------------------------#
#     Plot options.                                                                        #
#------------------------------------------------------------------------------------------#
outform         = c("pdf")        # Formats for output file (case insensitive):
                                  #   - "X11" - for printing on screen
                                  #   - "eps" - for postscript printing
                                  #   - "png" - for PNG printing
                                  #   - "pdf" - for PDF printing
byeold          = TRUE            # Remove old files of the given format?
depth           = 96              # PNG resolution, in pixels per inch
paper           = "square"        # Paper size, to define the plot shape
ptsz            = 18              # Font size.
plotgrid        = TRUE            # Should I plot the grid in the background?
legwhere        = "topleft"       # Where should I place the legend?
inset           = 0.01            # inset distance between legend and edge of plot region.
legbg           = "white"         # Legend background colour.
plotsd          = FALSE           # Plot the standard deviation (T/F)
angle           = 45              # Angle for shaded 1-SD
dens            = 40              # Density of lines
llwd            = 2.5             # Line width
fracexp         = 0.2             # Expansion to fit the legend
shwd            = 1.0             # Width of shaded lines
plot.blackout   = FALSE           # Should I plot a gray box during the "blackout" period
col.blackout    = "bisque"        # Colour for the blackout period
whena.blackout  = c("01/01/2001") # Beginning of the black-out periods
whenz.blackout  = c("01/01/2002") # End of the black-out periods
f.leg          = 1/6              # Fraction of the plotting area for legend.
#------------------------------------------------------------------------------------------#




#------------------------------------------------------------------------------------------#
#    Plot the gaps.  The variables on the list are defined as:                             #
#                                                                                          #
# vari:     The variable name                                                              #
# desc:     Variable description (for plotting)                                            #
# unit:     Variable units (for plotting)                                                  #
# harm:     Use harmonic analysis for filling (TRUE/FALSE)                                 #
# form:     In case harm = FALSE, the linear model to be used.                             #
#           If form=NA and harm=FALSE, then the gaps will not be filled unless you hard-   #
#           code a calculation.                                                            #
#------------------------------------------------------------------------------------------#
n         = 0
tvar      = list()
n         = n + 1
tvar[[n]] = list( vnam      = "rshort.in"
                , desc      = "Incident SW"
                , unit      = "wom2"
                , add       = NA_character_
                , mult      = NA_character_
                , colmean   = "darkorange3"
                , colerr    = "gold"
                , coledge   = "gray31"
                , plt       = TRUE
                , edge      = TRUE
                , out.hour  = TRUE
                , out.daily = TRUE
                )#end list
n         = n + 1
tvar[[n]] = list( vnam      = "rshort.out"
                , desc      = "Outgoing SW"
                , unit      = "wom2"
                , add       = NA_character_
                , mult      = NA_character_
                , colmean   = "darkorange3"
                , colerr    = "gold"
                , coledge   = "gray31"
                , plt       = TRUE
                , edge      = TRUE
                , out.hour  = TRUE
                , out.daily = TRUE
                )#end list
n         = n + 1
tvar[[n]] = list( vnam      = "par.in"
                , desc      = "Incident PAR"
                , unit      = "umolom2os"
                , add       = NA_character_
                , mult      = "Watts.2.Ein * 1.e6"
                , colmean   = "darkorange3"
                , colerr    = "gold"
                , coledge   = "gray31"
                , plt       = TRUE
                , edge      = TRUE
                , out.hour  = TRUE
                , out.daily = TRUE
                )#end list
n         = n + 1
tvar[[n]] = list( vnam      = "par.out"
                , desc      = "Outgoing PAR"
                , unit      = "umolom2os"
                , add       = NA_character_
                , mult      = "Watts.2.Ein * 1.e6"
                , colmean   = "darkorange3"
                , colerr    = "gold"
                , coledge   = "gray31"
                , plt       = TRUE
                , edge      = TRUE
                , out.hour  = TRUE
                , out.daily = TRUE
                )#end list
n         = n + 1
tvar[[n]] = list( vnam      = "rlong.in"
                , desc      = "Incident LW"
                , unit      = "wom2"
                , add       = NA_character_
                , mult      = NA_character_
                , colmean   = "midnightblue"
                , colerr    = "steelblue"
                , coledge   = "gray31"
                , plt       = TRUE
                , edge      = TRUE
                , out.hour  = TRUE
                , out.daily = TRUE
                )#end list
n         = n + 1
tvar[[n]] = list( vnam      = "rlong.out"
                , desc      = "Outgoing LW"
                , unit      = "wom2"
                , add       = NA_character_
                , mult      = NA_character_
                , colmean   = "midnightblue"
                , colerr    = "steelblue"
                , coledge   = "gray31"
                , plt       = TRUE
                , edge      = TRUE
                , out.hour  = TRUE
                , out.daily = TRUE
                )#end list
n         = n + 1
tvar[[n]] = list( vnam      = "rnet"
                , desc      = "Net radiation"
                , unit      = "wom2"
                , add       = NA_character_
                , mult      = NA_character_
                , colmean   = "midnightblue"
                , colerr    = "steelblue"
                , coledge   = "gray31"
                , plt       = TRUE
                , edge      = TRUE
                , out.hour  = TRUE
                , out.daily = TRUE
                )#end list
n         = n + 1
tvar[[n]] = list( vnam      = "atm.tmp"
                , desc      = "Air temperature"
                , unit      = "degC"
                , add       = "-t00"
                , mult      = NA_character_
                , colmean   = "darkorange3"
                , colerr    = "gold"
                , coledge   = "gray31"
                , plt       = TRUE
                , edge      = TRUE
                , out.hour  = TRUE
                , out.daily = TRUE
                )#end list
n         = n + 1
tvar[[n]] = list( vnam      = "atm.rhv"
                , desc      = "Relative humidity"
                , unit      = "pc"
                , add       = NA_character_
                , mult      = "100."
                , colmean   = "midnightblue"
                , colerr    = "steelblue"
                , coledge   = "gray31"
                , plt       = TRUE
                , edge      = TRUE
                , out.hour  = TRUE
                , out.daily = TRUE
                )#end list
n         = n + 1
tvar[[n]] = list( vnam      = "atm.pvap"
                , desc      = "Vapour pressure"
                , unit      = "hpa"
                , add       = NA_character_
                , mult      = "0.01"
                , colmean   = "midnightblue"
                , colerr    = "steelblue"
                , coledge   = "gray31"
                , plt       = TRUE
                , edge      = TRUE
                , out.hour  = TRUE
                , out.daily = TRUE
                )#end list
n         = n + 1
tvar[[n]] = list( vnam      = "rain"
                , desc      = "Precipitation rate"
                , unit      = "kgwom2ohr"
                , add       = NA_character_
                , mult      = "hr.sec"
                , colmean   = "midnightblue"
                , colerr    = "steelblue"
                , coledge   = "gray31"
                , plt       = TRUE
                , edge      = FALSE
                , out.hour  = FALSE
                , out.daily = FALSE
                )#end list
n         = n + 1
tvar[[n]] = list( vnam      = "atm.vels"
                , desc      = "Wind speed"
                , unit      = "mos"
                , add       = NA_character_
                , mult      = NA_character_
                , colmean   = "purple4"
                , colerr    = "mediumpurple1"
                , coledge   = "gray31"
                , plt       = TRUE
                , edge      = TRUE
                , out.hour  = TRUE
                , out.daily = TRUE
                )#end list
n         = n + 1
tvar[[n]] = list( vnam      = "atm.prss"
                , desc      = "Atmospheric Pressure"
                , unit      = "hpa"
                , add       = NA_character_
                , mult      = "0.01"
                , colmean   = "purple4"
                , colerr    = "mediumpurple1"
                , coledge   = "gray31"
                , plt       = TRUE
                , edge      = TRUE
                , out.hour  = TRUE
                , out.daily = TRUE
                )#end list
n         = n + 1
tvar[[n]] = list( vnam      = "ustar"
                , desc      = "Friction velocity"
                , unit      = "mos"
                , add       = NA_character_
                , mult      = NA_character_
                , colmean   = "purple4"
                , colerr    = "mediumpurple1"
                , coledge   = "gray31"
                , plt       = TRUE
                , edge      = TRUE
                , out.hour  = TRUE
                , out.daily = TRUE
                )#end list
n         = n + 1
tvar[[n]] = list( vnam      = "fco2"
                , desc      = "CO2 flux"
                , unit      = "umolom2os"
                , add       = NA_character_
                , mult      = NA_character_
                , colmean   = "forestgreen"
                , colerr    = "chartreuse"
                , coledge   = "gray31"
                , plt       = TRUE
                , edge      = TRUE
                , out.hour  = TRUE
                , out.daily = TRUE
                )#end list
n         = n + 1
tvar[[n]] = list( vnam      = "storco2"
                , desc      = "CO2 storage"
                , unit      = "umolom2os"
                , add       = NA_character_
                , mult      = NA_character_
                , colmean   = "forestgreen"
                , colerr    = "chartreuse"
                , coledge   = "gray31"
                , plt       = TRUE
                , edge      = TRUE
                , out.hour  = TRUE
                , out.daily = TRUE
                )#end list
n         = n + 1
tvar[[n]] = list( vnam      = "nee"
                , desc      = "Net ecosystem exchange"
                , unit      = "umolom2os"
                , add       = NA_character_
                , mult      = NA_character_
                , colmean   = "forestgreen"
                , colerr    = "chartreuse"
                , coledge   = "gray31"
                , plt       = TRUE
                , edge      = TRUE
                , out.hour  = TRUE
                , out.daily = TRUE
                )#end list
n         = n + 1
tvar[[n]] = list( vnam      = "fsens"
                , desc      = "Sensible heat flux"
                , unit      = "wom2"
                , add       = NA_character_
                , mult      = NA_character_
                , colmean   = "darkorange3"
                , colerr    = "gold"
                , coledge   = "gray31"
                , plt       = TRUE
                , edge      = TRUE
                , out.hour  = TRUE
                , out.daily = TRUE
                )#end list
n         = n + 1
tvar[[n]] = list( vnam      = "fh2o"
                , desc      = "Water vapour flux"
                , unit      = "kgwom2oday"
                , add       = NA_character_
                , mult      = "day.sec"
                , colmean   = "midnightblue"
                , colerr    = "steelblue"
                , coledge   = "gray31"
                , plt       = TRUE
                , edge      = TRUE
                , out.hour  = TRUE
                , out.daily = TRUE
                )#end list
n         = n + 1
tvar[[n]] = list( vnam      = "atm.shv"
                , desc      = "Specific humidity"
                , unit      = "gokg"
                , add       = NA_character_
                , mult      = "1000."
                , colmean   = "midnightblue"
                , colerr    = "steelblue"
                , coledge   = "gray31"
                , plt       = TRUE
                , edge      = TRUE
                , out.hour  = TRUE
                , out.daily = TRUE
                )#end list
n         = n + 1
tvar[[n]] = list( vnam      = "atm.vpd"
                , desc      = "Vapour pressure deficit"
                , unit      = "hpa"
                , add       = NA_character_
                , mult      = "0.01"
                , colmean   = "midnightblue"
                , colerr    = "steelblue"
                , coledge   = "gray31"
                , plt       = TRUE
                , edge      = TRUE
                , out.hour  = TRUE
                , out.daily = TRUE
                )#end list
n         = n + 1
tvar[[n]] = list( vnam      = "atm.rhos"
                , desc      = "Air density"
                , unit      = "kgom3"
                , add       = NA_character_
                , mult      = NA_character_
                , colmean   = "purple4"
                , colerr    = "mediumpurple1"
                , coledge   = "gray31"
                , plt       = TRUE
                , edge      = TRUE
                , out.hour  = TRUE
                , out.daily = TRUE
                )#end list
n         = n + 1
tvar[[n]] = list( vnam      = "par.beam"
                , desc      = "Direct PAR"
                , unit      = "umolom2os"
                , add       = NA_character_
                , mult      = "Watts.2.Ein * 1.e6"
                , colmean   = "darkorange3"
                , colerr    = "gold"
                , coledge   = "gray31"
                , plt       = TRUE
                , edge      = TRUE
                , out.hour  = TRUE
                , out.daily = TRUE
                )#end list
n         = n + 1
tvar[[n]] = list( vnam      = "par.diff"
                , desc      = "Diffuse PAR"
                , unit      = "umolom2os"
                , add       = NA_character_
                , mult      = "Watts.2.Ein * 1.e6"
                , colmean   = "midnightblue"
                , colerr    = "steelblue"
                , coledge   = "gray31"
                , plt       = TRUE
                , edge      = TRUE
                , out.hour  = TRUE
                , out.daily = TRUE
                )#end list
n         = n + 1
tvar[[n]] = list( vnam      = "nep"
                , desc      = "Net ecosystem productivity"
                , unit      = "umolom2os"
                , add       = NA_character_
                , mult      = NA_character_
                , colmean   = "forestgreen"
                , colerr    = "chartreuse"
                , coledge   = "gray31"
                , plt       = TRUE
                , edge      = TRUE
                , out.hour  = TRUE
                , out.daily = TRUE
                )#end list
n         = n + 1
tvar[[n]] = list( vnam      = "gpp"
                , desc      = "Gross primary productivity"
                , unit      = "umolom2os"
                , add       = NA_character_
                , mult      = NA_character_
                , colmean   = "forestgreen"
                , colerr    = "chartreuse"
                , coledge   = "gray31"
                , plt       = TRUE
                , edge      = TRUE
                , out.hour  = TRUE
                , out.daily = TRUE
                )#end list
n         = n + 1
tvar[[n]] = list( vnam      = "reco"
                , desc      = "Ecosystem respiration"
                , unit      = "umolom2os"
                , add       = NA_character_
                , mult      = NA_character_
                , colmean   = "forestgreen"
                , colerr    = "chartreuse"
                , coledge   = "gray31"
                , plt       = TRUE
                , edge      = TRUE
                , out.hour  = TRUE
                , out.daily = TRUE
                )#end list
n         = n + 1
tvar[[n]] = list( vnam      = "nee.pmb"
                , desc      = "Net ecosystem exchange - Brando"
                , unit      = "umolom2os"
                , add       = NA_character_
                , mult      = NA_character_
                , colmean   = "forestgreen"
                , colerr    = "chartreuse"
                , coledge   = "gray31"
                , plt       = TRUE
                , edge      = TRUE
                , out.hour  = TRUE
                , out.daily = TRUE
                )#end list
n         = n + 1
tvar[[n]] = list( vnam      = "gpp.pmb"
                , desc      = "Gross primary productivity - Brando"
                , unit      = "umolom2os"
                , add       = NA_character_
                , mult      = NA_character_
                , colmean   = "forestgreen"
                , colerr    = "chartreuse"
                , coledge   = "gray31"
                , plt       = TRUE
                , edge      = TRUE
                , out.hour  = TRUE
                , out.daily = TRUE
                )#end list
n         = n + 1
tvar[[n]] = list( vnam      = "reco.pmb"
                , desc      = "Ecosystem respiration - Brando"
                , unit      = "umolom2os"
                , add       = NA_character_
                , mult      = NA_character_
                , colmean   = "forestgreen"
                , colerr    = "chartreuse"
                , coledge   = "gray31"
                , plt       = TRUE
                , edge      = TRUE
                , out.hour  = TRUE
                , out.daily = TRUE
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



#------------------------------------------------------------------------------------------#
#    Load data and transform so we can tack on monthly means.                              #
#------------------------------------------------------------------------------------------#
cat0(" + Load data from ",basename(eft.rdata),".")
load(file=eft.rdata)
eft = as.list(eft)
#------------------------------------------------------------------------------------------#



#----- Create output directory in case it doesn't exist. ----------------------------------#
if (! file.exists(outroot)) dir.create(outroot)
#------------------------------------------------------------------------------------------#



#----- Define plot window size ------------------------------------------------------------#
f.ext   = f.leg / (1. - f.leg)
sq.size = plotsize(proje=FALSE,paper=paper)
ex.size = plotsize(proje=FALSE,paper=paper,extendfc="lat",extfactor=f.ext)
#------------------------------------------------------------------------------------------#



#----- Determine how many variables we must plot. -----------------------------------------#
tvar    = list.2.data.frame(tvar)
ntvar   = nrow(tvar)
#------------------------------------------------------------------------------------------#




#----- Determine how many formats we must output. -----------------------------------------#
outform = tolower(outform)
nout    = length(outform)
#------------------------------------------------------------------------------------------#




#----- Determine how many black-out periods there were, and convert their times to chron. -#
nblackout      = length(whena.blackout)
whena.blackout = chron(whena.blackout)
whenz.blackout = chron(whenz.blackout)
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#     Adjust some data according to the gap-filling level.                                 #
#------------------------------------------------------------------------------------------#
if (! gflevel %in% c("metfill","filled")){
   cat0(" + Adjust and find some variables that may not be available.")
   eft$atm.rhos = idealdenssh(eft$atm.prss,eft$atm.tmp,eft$atm.shv)
   eft$atm.vpd  = eslif(eft$atm.tmp) - eft$atm.pvap
   eft$par.beam = eft$rshort.in * NA
   eft$par.diff = eft$rshort.in * NA
   eft$nir.beam = eft$rshort.in * NA
   eft$nir.diff = eft$rshort.in * NA
   eft$atm.shv  = ep * eft$atm.pvap / (eft$atm.prss - (1. - ep) * eft$atm.pvap)
}#end if
#------------------------------------------------------------------------------------------#


#------------------------------------------------------------------------------------------#
#      Do not show ecosystem fluxes unless this is the final stage.                        #
#------------------------------------------------------------------------------------------#
if (! gflevel %in% c("filled")){
   cat0(" + Censor ecosystem fluxes at the ",gflevel," stage (they are not ready).")
   eft$nep  = NA + eft$cosz
   eft$nee  = NA + eft$cosz
   eft$gpp  = NA + eft$cosz
   eft$reco = NA + eft$cosz
}#end if
#------------------------------------------------------------------------------------------#


#------------------------------------------------------------------------------------------#
#      Find the aggregated variables.                                                      #
#------------------------------------------------------------------------------------------#
cat0(" + Find daily and monthly time vectors.")
eft$today   = chron(unique(eft$today))
eft$tomonth = chron(unique(dates(eft$when[eft$day == 1])))
#------------------------------------------------------------------------------------------#


#------------------------------------------------------------------------------------------#
#     Find several auxiliary variables for daily and monthly statistics...                 #
#------------------------------------------------------------------------------------------#
cat0(" + Find time-related variables for the daily and monthly statistics.")
mm    = eft$month
hh    = eft$hour
yy    = eft$year
hrs   = sort(unique(hh))
mons  = tolower(month.abb[sort(unique(mm))])
yrs   = unique(eft$year)
nhrs  = length(hrs)
nmons = length(mons)
nyrs  = length(yrs)
#----- Find the maximum number of days for normal and leap years. -------------------------#
emean.dmax = daymax(nummonths(eft$tomonth),numyears(eft$tomonth))
mmean.dmax = tapply(X=emean.dmax,INDEX=nummonths(eft$tomonth),FUN=mean,na.rm=TRUE)
#------------------------------------------------------------------------------------------#


#------------------------------------------------------------------------------------------#
#     Obtain the mean diurnal cycle.                                                       #
#------------------------------------------------------------------------------------------#
cat0(" + Find all sorts of averages:")
for (v in sequence(ntvar)){
    varname   = tvar$vnam[v]
    vardesc   = tvar$desc[v]
    varval   = eft[[varname]]
    varval2  = varval * varval
    cat0("    * ",vardesc,".")

    okvar    = is.finite(varval)
    #----- Mean diurnal cycle. ------------------------------------------------------------#
    qmean    = aggregate(varval ,by=list(mm,hh),FUN=mean,na.rm=TRUE,simplify=TRUE)
    qmean    = matrix   (qmean[,3],nrow=nmons,ncol=nhrs,byrow=FALSE)
    dimnames(qmean) = list(mons,hrs)
    #----- Mean sum of squares (diurnal cycle). -------------------------------------------#
    qmsqu    = aggregate(varval2,by=list(mm,hh),FUN=mean,na.rm=TRUE,simplify=TRUE)
    qmsqu    = matrix   (qmsqu[,3],nrow=nmons,ncol=nhrs,byrow=FALSE)
    dimnames(qmsqu) = list(mons,hrs)
    #----- Standard deviation. ------------------------------------------------------------#
    qsdev    = aggregate(varval  ,by=list(mm,hh),FUN=sd  ,na.rm=TRUE,simplify=TRUE)
    qsdev    = matrix   (qsdev[,3],nrow=nmons,ncol=nhrs,byrow=FALSE)
    dimnames(qsdev) = list(mons,hrs)
    #----- Valid data. --------------------------------------------------------------------#
    qccnt    = aggregate(okvar    ,by=list(mm,hh),FUN=sum  ,na.rm=TRUE,simplify=TRUE)
    qccnt    = matrix   (qccnt[,3],nrow=nmons,ncol=nhrs,byrow=FALSE)
    dimnames(qccnt) = list(mons,hrs)
    #--------------------------------------------------------------------------------------#


    #----- Time series of monthly means. --------------------------------------------------#
    a.mean = aggregate(x=varval ,by=list(hh,mm,yy),FUN=mean,na.rm=TRUE,simplify=TRUE)
    a.msqu = aggregate(x=varval2,by=list(hh,mm,yy),FUN=mean,na.rm=TRUE,simplify=TRUE)
    a.ccnt = aggregate(x=okvar  ,by=list(hh,mm,yy),FUN=sum ,na.rm=TRUE,simplify=TRUE)
    nagg   = length(a.mean[,4])
    ind    = paste(substring(10000+a.mean[,3],2,5),substring(100+a.mean[,2],2,3)
                  ,sep="",concatenate="")
    emean   = tapply(X=a.mean[,4],INDEX=ind,FUN=mean,na.rm=FALSE)
    emsqu   = tapply(X=a.msqu[,4],INDEX=ind,FUN=mean,na.rm=FALSE)
    eccnt   = tapply(X=a.ccnt[,4],INDEX=ind,FUN=sum ,na.rm=FALSE)
    srnorm1 = sqrt(1./(1. - 1./eccnt))
    esdev   = sqrt((emsqu - emean*emean)*srnorm1)
    #--------------------------------------------------------------------------------------#



    #----- Daily mean and standard deviation. ---------------------------------------------#
    today = dates(eft$when)
    dmean = tapply(X=varval ,INDEX=today,FUN=mean,na.rm=FALSE)
    dmsqu = tapply(X=varval2,INDEX=today,FUN=mean,na.rm=FALSE)
    dsdev = tapply(X=varval ,INDEX=today,FUN=sd  ,na.rm=FALSE)
    dccnt = tapply(X=okvar  ,INDEX=today,FUN=sum ,na.rm=FALSE)
    #--------------------------------------------------------------------------------------#



    #----- Monthly mean and standard deviation. -------------------------------------------#
    mmean   = rowMeans(qmean,na.rm=TRUE)
    mmsqu   = rowMeans(qmsqu,na.rm=TRUE)
    mccnt   = rowMeans(qccnt,na.rm=TRUE)
    srnorm1 = sqrt(1./(1. - 1./qccnt))
    msdev   = rowMeans(sqrt(qmsqu-qmean*qmean)*srnorm1,na.rm=TRUE)
    #--------------------------------------------------------------------------------------#



    #----- Change units for precipitation (mm/month) for monthly means. -------------------#
    if (varname %in% "rain"){
       dmean = dmean              * day.sec
       dmsqu = dmsqu              * day.sec
       dsdev = dsdev              * day.sec
       mmean = mmean * mmean.dmax * day.sec
       mmsqu = mmsqu * mmean.dmax * day.sec
       msdev = msdev * mmean.dmax * day.sec
       emean = emean * emean.dmax * day.sec
       emsqu = emsqu * emean.dmax * day.sec
       esdev = esdev * emean.dmax * day.sec
    }#end if
    #--------------------------------------------------------------------------------------#



    #--------------------------------------------------------------------------------------#
    #     Make unique names for the statistics.                                            #
    #--------------------------------------------------------------------------------------#
    qmean.name = paste0("qmean.",varname)
    qmsqu.name = paste0("qmsqu.",varname)
    qsdev.name = paste0("qsdev.",varname)
    qccnt.name = paste0("qccnt.",varname)
    dmean.name = paste0("dmean.",varname)
    dmsqu.name = paste0("dmsqu.",varname)
    dsdev.name = paste0("dsdev.",varname)
    dccnt.name = paste0("dccnt.",varname)
    mmean.name = paste0("mmean.",varname)
    mmsqu.name = paste0("mmsqu.",varname)
    msdev.name = paste0("msdev.",varname)
    mccnt.name = paste0("mccnt.",varname)
    emean.name = paste0("emean.",varname)
    emsqu.name = paste0("emsqu.",varname)
    esdev.name = paste0("esdev.",varname)
    eccnt.name = paste0("eccnt.",varname)
    #--------------------------------------------------------------------------------------#



    #--------------------------------------------------------------------------------------#
    #     Append the statistics to the list.                                               #
    #--------------------------------------------------------------------------------------#
    eft[[qmean.name]] = qmean
    eft[[qmsqu.name]] = qmsqu
    eft[[qsdev.name]] = qsdev
    eft[[qccnt.name]] = qccnt
    eft[[dmean.name]] = dmean
    eft[[dmsqu.name]] = dmsqu
    eft[[dsdev.name]] = dsdev
    eft[[dccnt.name]] = dccnt
    eft[[mmean.name]] = mmean
    eft[[mmsqu.name]] = mmsqu
    eft[[msdev.name]] = msdev
    eft[[mccnt.name]] = mccnt
    eft[[emean.name]] = emean
    eft[[emsqu.name]] = emsqu
    eft[[esdev.name]] = esdev
    eft[[eccnt.name]] = eccnt
    #--------------------------------------------------------------------------------------#
}#end for (v in sequence(ntvar)
#===========================================================================================#
#===========================================================================================#






#===========================================================================================#
#===========================================================================================#
#   Plotting output.                                                                        #
#-------------------------------------------------------------------------------------------#
cat0(" + Plot averages for ",longname,".")

#----- Copy structure to local variables. -------------------------------------------------#
outlevel  = file.path(outroot,gflevel)


outmdcyc   = file.path(outlevel,"qmean")
outmmean   = file.path(outlevel,"mmean")
outdmean   = file.path(outlevel,"dmean")
outemean   = file.path(outlevel,"emean")
outhours   = file.path(outlevel,"hours")
outycomp   = file.path(outlevel,"ycomp")

if (! file.exists(outlevel)) dir.create(outlevel)
if (! file.exists(outmdcyc)) dir.create(outmdcyc)
if (! file.exists(outmmean)) dir.create(outmmean)
if (! file.exists(outdmean)) dir.create(outdmean)
if (! file.exists(outemean)) dir.create(outemean)
if (! file.exists(outhours)) dir.create(outhours)

#------------------------------------------------------------------------------------------#
#    Plot the diurnal cycle by month.                                                      #
#------------------------------------------------------------------------------------------#
for (tt in sequence(ntvar)){
   #------ Handy aliases. -----------------------------------------------------------------#
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
   #---------------------------------------------------------------------------------------#


   #------ Convert unit conversion into actual numbers. -----------------------------------#
   if (is.na(add0 )) add0  = 0. else add0  = eval(parse(text=add0 ))
   if (is.na(mult0)) mult0 = 1. else mult0 = eval(parse(text=mult0))
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #      Plot averages.                                                                   #
   #---------------------------------------------------------------------------------------#
   if (plotit){
      cat0("   - Variable ",desc,".")
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
      #    Mean diurnal cycle.                                                             #
      #------------------------------------------------------------------------------------#
      cat0("     > Mean diurnal cycle.")

      #----- Create path in case it doesn't exist. ----------------------------------------#
      outvar = file.path(outmdcyc,vname)
      if (! file.exists(outvar)) dir.create(outvar)
      #------------------------------------------------------------------------------------#

      #----- Retrieve mean and standard deviation. ----------------------------------------#
      qmean.name = paste0("qmean.",vname)
      qsdev.name = paste0("qsdev.",vname)
      thismean = mult0 * (eft[[qmean.name]]+add0)
      thissdev = mult0 *  eft[[qsdev.name]]
      #------------------------------------------------------------------------------------#


      #----- Retrieve hour and month information. -----------------------------------------#
      mmms  = rownames(thismean)
      mons  = match(mmms,tolower(month.abb))
      nmons = length(mons)
      hday  = as.numeric(colnames(thismean))
      #------------------------------------------------------------------------------------#


      #----- Set nice x and y limits. -----------------------------------------------------#
      xlimit  = pretty.xylim(u=hday,fracexp=0.0)
      xat     = pretty.elapsed(xlimit,base=24)
      xlabels = sprintf("%g",xat)
      if (plotsd){
         ylimit = pretty.xylim(c(thismean-thissdev,thismean+thissdev),fracexp=fracexp)
      }else{
         ylimit = pretty.xylim(thismean)
      }#end if (plotsd)
      yat = pretty(ylimit)
      ylabels = sprintf("%g",yat)
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #     Loop through months.                                                           #
      #------------------------------------------------------------------------------------#
      for (mm in sequence(nmons)){
         #----- Handy aliases. ------------------------------------------------------------#
         cmon = sprintf("%2.2i",mons[mm])
         mmm  = month.name[mm]
         #---------------------------------------------------------------------------------#

         #----- Build plot data and polygons. ---------------------------------------------#
         monmean = thismean[mm,]
         monsdev = thissdev[mm,]
         err.x = c(hday,rev(hday))
         err.y = c(monmean+monsdev,rev(monmean-monsdev))
         #---------------------------------------------------------------------------------#



         #----- Plot annotation. ----------------------------------------------------------#
         letitre = paste0(longname,"(",yeara,"-",yearz,")","\n",desc," - ",mmm)
         lex     = desc.unit(desc="Time",unit=untab$utc)
         ley     = desc.unit(desc=desc  ,unit=unit)
         #---------------------------------------------------------------------------------#


         #---------------------------------------------------------------------------------#
         #     Loop through formats.                                                       #
         #---------------------------------------------------------------------------------#
         for (o in sequence(nout)){
            #----- Open plotting device. --------------------------------------------------#
            fichier = file.path(outvar,paste0(vname,"-",cmon,".",outform[o]))
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
            #------------------------------------------------------------------------------#
         }#end for (o in sequence(nout))
         #---------------------------------------------------------------------------------#
      }#end for (m in sequence(nmonths))
      #------------------------------------------------------------------------------------#
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#





      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
      #    Monthly means.                                                                  #
      #------------------------------------------------------------------------------------#
      cat0("     > ",desc,".")


      #----- Find average and standard deviation. -----------------------------------------#
      mmean.name = paste("mmean",vname,sep=".")
      msdev.name = paste("msdev",vname,sep=".")
      if (vname %in% c("rain","monrain")){
         monmean    = eft[[mmean.name]]
         monsdev    = eft[[msdev.name]]
      }else{
         monmean    = mult0 * (eft[[mmean.name]] + add0 )
         monsdev    = mult0 *  eft[[msdev.name]]
      }#end if (vname %in% c("rain","monrain"))
      #------------------------------------------------------------------------------------#




      #------ Organise months and time of day. --------------------------------------------#
      xlimit = c(1,12)
      mons   = sort(unique(eft$mon))
      nmons  = length(mons)
      mlab   = mon2mmm(mons)
      #------------------------------------------------------------------------------------#




      #----- Fix limits for y axis so all months use the same scale. ----------------------#
      if ((! plotsd) || vname %in% c("rain","monrain")){
         ylimit = pretty.xylim(u=monmean,fracexp=fracexp)
      }else{
         ylimit = pretty.xylim(u=c(monmean-monsdev,monmean+monsdev),fracexp=fracexp)
      }#end if
      yat     = pretty(ylimit)
      ylabels = sprintf("%g",yat)
      #------------------------------------------------------------------------------------#

      #------ Find polygons. --------------------------------------------------------------#
      err.x = c(mons,rev(mons))
      err.y = c(monmean+monsdev,rev(monmean-monsdev))
      #------------------------------------------------------------------------------------#



      #------ Plot annotation. ------------------------------------------------------------#
      letitre = paste0(longname,"(",yeara,"-",yearz,")","\n",desc)
      lex     = "Month"
      if (vname %in% c("rain","monrain")){
         ley     = desc.unit(desc=desc,unit=untab$mmomo)
      }else{
         ley     = desc.unit(desc=desc,unit=unit)
      }#end if (vname %in% c("rain","monrain"))
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #     Loop through all formats.                                                      #
      #------------------------------------------------------------------------------------#
      for (o in sequence(nout)){
         #----- Open device. --------------------------------------------------------------#
         fichier = file.path(outmmean,paste0(vname,".",outform[o]))
         dummy   = open.plot( fichier = fichier
                            , outform = outform[o]
                            , size    = sq.size
                            , ptsz    = ptsz
                            , depth   = depth
                            )#end open.plot
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #     Plot main figure.                                                           #
         #---------------------------------------------------------------------------------#
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
         #---------------------------------------------------------------------------------#
      }#end for (o in sequence(nout))
      #------------------------------------------------------------------------------------#
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#





      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
      #     Daily means.                                                                   #
      #------------------------------------------------------------------------------------#
      cat0("     > Daily averages.")


      #----- Find average and standard deviation. -----------------------------------------#
      dmean.name = paste("dmean",vname,sep=".")
      dsdev.name = paste("dsdev",vname,sep=".")
      today   = eft$today
      if (vname %in% c("rain","monrain")){
         daymean    = eft[[dmean.name]]
         daysdev    = eft[[dsdev.name]]
      }else{
         daymean    = mult0 * (eft[[dmean.name]] + add0 )
         daysdev    = mult0 *  eft[[dsdev.name]]
      }#end if (vname %in% c("rain","monrain"))
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #      Make the time axis.                                                           #
      #------------------------------------------------------------------------------------#
      days     = eft$today
      ndays    = length(days)
      whenplot = pretty.time(days,n=8)
      xlimit   = range(days)
      #------------------------------------------------------------------------------------#




      #----- Fix limits for y axis so all months use the same scale. ----------------------#
      ylimit  = pretty.xylim(u=daymean,fracexp=fracexp)
      yat     = pretty(ylimit)
      ylabels = sprintf("%g",yat)
      #------------------------------------------------------------------------------------#



      #------ Plot annotation. ------------------------------------------------------------#
      letitre = paste0(longname,"(",yeara,"-",yearz,")","\n",desc)
      lex     = "Dates"
      if (vname %in% c("rain","monrain")){
         ley     = desc.unit(desc=desc,unit=untab$mmoday)
      }else{
         ley     = desc.unit(desc=desc,unit=unit)
      }#end if (vname %in% c("rain","monrain"))
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #     Loop through all formats.                                                      #
      #------------------------------------------------------------------------------------#
      for (o in sequence(nout)){
         #----- Open device. --------------------------------------------------------------#
         fichier = file.path(outdmean,paste0(vname,".",outform[o]))
         dummy   = open.plot( fichier = fichier
                            , outform = outform[o]
                            , size    = sq.size
                            , ptsz    = ptsz
                            , depth   = depth
                            )#end open.plot
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #     Plot main figure.                                                           #
         #---------------------------------------------------------------------------------#
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
         #---------------------------------------------------------------------------------#
      }#end for (o in sequence(nout))
      #------------------------------------------------------------------------------------#
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#





      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
      #     Time series of monthly means.                                                  #
      #------------------------------------------------------------------------------------#
      cat0("     > Time series of monthly averages.")


      #----- Find average and standard deviation. -----------------------------------------#
      emean.name = paste("emean",vname,sep=".")
      esdev.name = paste("esdev",vname,sep=".")
      if (vname %in% c("rain","monrain")){
         mtsmean    = eft[[emean.name]]
         mtssdev    = eft[[esdev.name]]
      }else{
         mtsmean    = mult0 * ( eft[[emean.name]] + add0 )
         mtssdev    = mult0 *   eft[[esdev.name]]
      }#end if (vname %in% c("rain","monrain"))
      tomonth    = eft$tomonth
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #      Make the time axis.                                                           #
      #------------------------------------------------------------------------------------#
      nmon     = length(tomonth)
      whenplot = pretty.time(tomonth,n=8)
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #     Find the limits for the y-axis.                                                #
      #------------------------------------------------------------------------------------#
      if ((! plotsd) || vname == "rain" || vname == "monrain"){
         ylimit = pretty.xylim(u=mtsmean,fracexp=fracexp)
      }else{
         ylimit = pretty.xylim(u=c(mtsmean-mtssdev,mtsmean+mtssdev),fracexp=fracexp)
      }#end if ((! plotsd) || vname == "rain" || vname == "monrain")
      yat = pretty(ylimit)
      ylabels = sprintf("%g",yat)
      #------------------------------------------------------------------------------------#




      #------ Find polygons. --------------------------------------------------------------#
      err.x = c(tomonth,rev(tomonth))
      err.y = c(mtsmean+mtssdev,rev(mtsmean-mtssdev))
      #------------------------------------------------------------------------------------#



      #------ Plot annotation. ------------------------------------------------------------#
      letitre = paste0(longname,"(",yeara,"-",yearz,")","\n",desc)
      lex     = "Time"
      if (vname %in% c("rain","monrain")){
         ley     = desc.unit(desc=desc,unit=untab$mmomo)
      }else{
         ley     = desc.unit(desc=desc,unit=unit)
      }#end if (vname %in% c("rain","monrain"))
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #     Loop through all formats.                                                      #
      #------------------------------------------------------------------------------------#
      for (o in sequence(nout)){
         #----- Open device. --------------------------------------------------------------#
         fichier = file.path(outemean,paste0(vname,".",outform[o]))
         dummy   = open.plot( fichier = fichier
                            , outform = outform[o]
                            , size    = sq.size
                            , ptsz    = ptsz
                            , depth   = depth
                            )#end open.plot
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #     Plot main figure.                                                           #
         #---------------------------------------------------------------------------------#
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
         #---------------------------------------------------------------------------------#
      }#end for (o in sequence(nout))
      #------------------------------------------------------------------------------------#
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#





      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
      #     Hourly time series.                                                            #
      #------------------------------------------------------------------------------------#
      cat0("     > Time series of hourly values.")
      if (vname %in% "rain"){
         varnow  = mult0 * ( eft[[vname]] + add0 ) / 24.
      }else{
         varnow  = mult0 * ( eft[[vname]] + add0 )
      }#end if
      when    = eft$when
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #      Make the time axis.                                                           #
      #------------------------------------------------------------------------------------#
      nwhen    = length(when)
      whenplot = pretty.time(when,n=8)
      xlimit   = range(when)
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #     Find the limits for the y-axis.                                                #
      #------------------------------------------------------------------------------------#
      ylimit  = pretty.xylim(u=varnow,fracexp=fracexp)
      yat     = pretty(ylimit)
      ylabels = sprintf("%g",yat)
      #------------------------------------------------------------------------------------#



      #------ Plot annotation. ------------------------------------------------------------#
      letitre = paste0(longname,"(",yeara,"-",yearz,")","\n",desc)
      lex     = "Dates"
      ley     = desc.unit(desc=desc,unit=unit)
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #     Loop through all formats.                                                      #
      #------------------------------------------------------------------------------------#
      for (o in sequence(nout)){
         #----- Open device. --------------------------------------------------------------#
         fichier = file.path(outhours,paste0(vname,".",outform[o]))
         dummy   = open.plot( fichier = fichier
                            , outform = outform[o]
                            , size    = sq.size
                            , ptsz    = ptsz
                            , depth   = depth
                            )#end open.plot
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #     Plot main figure.                                                           #
         #---------------------------------------------------------------------------------#
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
         #---------------------------------------------------------------------------------#
      }#end for (o in sequence(nout))
      #------------------------------------------------------------------------------------#
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
   }#end if (plotit)
   #---------------------------------------------------------------------------------------#
}#end for (tt in sequence(ntvar))
#------------------------------------------------------------------------------------------#
