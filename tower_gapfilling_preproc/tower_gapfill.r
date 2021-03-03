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
here       = getwd()                         # Current path 
srcdir     = file.path(here,"Rsc"          ) # Path with additional scripts
outroot    = file.path(here,"gaps"         ) # Output directory
other.path = file.path(here,"other_sites"  ) # Path with RData for additional sites.
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
dtdat    = 3600                     # Tower period between two records in seconds
                                    #    This is used to correct rainfall so it is in 
                                    #    kg/m2/s.  Everything else uses hourly.
ref.prss = 95000.                   # Reference pressure, for potential
                                    #    radiation.  Use a value that
                                    #    is typical for the site, 
                                    #    usually a value in the low end
                                    #    (but not outliers, typical low
                                    #    pressure values).  Values in
                                    #    Pa, please.
off.utc   = -4                      # Offset for UTC.  If data are in local time, the time
                                    #    stamp corresponds to the end of the averaging
                                    #    period, and the time zone is UTC-4, set
                                    #    off.utc=-4.  If the data are already in UTC, set 
                                    #    off.utc=0.  If the time is UTC but the time stamp
                                    #    is  the beginning of the hourly average period, 
                                    #    set off.utc=-1 
#------------------------------------------------------------------------------------------#




#------------------------------------------------------------------------------------------#
#     File information.  Here we save several levels of QA/QC/GF, so we can check the time #
# series in their various stages (the stages are sequential, so each level includes the    #
# previous levels).                                                                        #
# * raw     -- Raw data.  The original dataset, shifted so the time stamp is at the end of #
#              the hour and standardised units, but otherwise no QA/QC/GF applied.         #
# * clean   -- Clean data.  We remove known bad periods (hardcoded...), outliers, correct  #
#              drift and bias, and gap fill some variables with redundant measurements     #
# * selfgf  -- Partially gap filled data.  We use similar measurements from the site       #
#              itself to fill in data.                                                     #
# * objana  -- Partially gap filled data.  We use objective analysis from neighbour AWS.   #
# * harmana -- Partially gap filled data.  We use harmonic analysis for some variables     #
# * fspline -- Partially gap filled data.  We use spline to fill in any minor remaining    #
#              gap (recommended only for those variables filled using objective analysis   #
#              and harmonic analysis).                                                     #
# * metfill -- Parametric gap filling for variables that couldn't be filled with harmonic  #
#              analysis.  All variables but NEE, GPP, and respiration are completely gap   #
#              filled.                                                                     #
# * filled  -- Gap filling of NEE, GPP, and respiration (done in a separate script).       #
#------------------------------------------------------------------------------------------#
rdata.suffix   = paste0(iata,"_",yeara,"-",yearz,".RData")
reload.raw     = TRUE   ;  raw.rdata      = paste0("raw_"    ,rdata.suffix)
reload.clean   = TRUE   ;  clean.rdata    = paste0("clean_"  ,rdata.suffix)
reload.selfgf  = TRUE   ;  selfgf.rdata   = paste0("selfgf_" ,rdata.suffix)
reload.objana  = TRUE   ;  objana.rdata   = paste0("objana_" ,rdata.suffix)
reload.harmana = TRUE   ;  harmana.rdata  = paste0("harmana_",rdata.suffix)
reload.fspline = TRUE   ;  fspline.rdata  = paste0("fspline_",rdata.suffix)
reload.metfill = TRUE   ;  metfill.rdata  = paste0("metfill_",rdata.suffix)
#------------------------------------------------------------------------------------------#




#------------------------------------------------------------------------------------------#
#     Define the site-specific information for reading the dataset.                        #
#------------------------------------------------------------------------------------------#
eddyfile    = "Control_Gap_Filled_MDS_masterFEB2017.csv"  # Eddy flux
raindayfile = "inmet.83270.csv"                           # Daily precipitation file
other.rbase = "raw_othersites-2008-2016.RData"            # R object with neighbour sites
                                                          #    (path will be added later)
whichsep     = ","              # String that defines a separation of fiels
                                #    comma is ",", space is " ", and tab is "\t".
na.strings   = "NA"             # String that identifies missing data (all strings must
                                #    be listed if there is more than one way to flag
                                #    missing data
comment.char = ""               # Flag to denote comment in the input file
                                #    If there is no comment, make this "" and the reading
                                #    will be significantly faster.
frac.cosz    = 2./3.            # Fraction of winter time noon 
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#     List of variables to be read in.  Each file has a different format, so it's unlikely #
# you will read a new dataset without modifying the code below.  Anyway, it's worth to     #
# make the code the most standardised we can...  This is the "default" name for each vari- #
# able, and the units are the recommended output units (you may edit the code to convert   #
# them to the right units).                                                                #
#                                                                                          #
# year           - calendar year                                              [       UTC] #
# month          - calendar month                                             [       UTC] #
# day            - calendar day                                               [       UTC] #
# hour           - hour                                                       [       UTC] #
# minu           - minute                                                     [       UTC] #
# dfrac          - day fraction (32.75 means Feb 1, YYYY - 18:00:00)          [       UTC] #
# doy            - day of year  (32 means Feb 1, YYYY)                        [       UTC] #
# hhmm           - hour and minute                                            [       UTC] #
# fco2           - CO2 flux, not storage-corrected                            [ umol/m2/s] #
# storco2        - CO2 storage                                                [ umol/m2/s] #
# ustar          - Friction velocity                                          [       m/s] #
# zeta           - z/Obukhov length                                           [       m/m] #
# fsens          - Sensible heat flux                                         [      W/m2] #
# storsens       - Sensible heat storage                                      [      W/m2] #
# flatent        - Latent heat flux                                           [   kg/m2/s] #
# storlatent     - Latent heat storage                                        [   kg/m2/s] #
# fh2o           - Water vapour flux                                          [   kg/m2/s] #
# storh2o        - Water vapour storage                                       [   kg/m2/s] #
# par.in         - Incoming PAR                                               [ umol/m2/s] #
# par.diff       - Incoming diffuse PAR                                       [ umol/m2/s] #
# par.out        - Outgoing PAR (2004-2006)                                   [ umol/m2/s] #
# rshort.in      - Incoming solar radiation                                   [      W/m2] #
# rshort.out     - Outgoing solar radiation                                   [      W/m2] #
# rlong.in       - Incoming longwave radiation                                [      W/m2] #
# rlong.out      - Outgoing longwave radiation                                [      W/m2] #
# rnet           - Net radiation                                              [      W/m2] #
# atm.tmp        - Air temperature                                            [         K] #
# atm.rhv        - Air relative humidity                                      [        --] #
# atm.shv        - Air specific humidity                                      [     kg/kg] #
# atm.vpd        - Air vapour pressure deficit                                [        Pa] #
# rain           - Rainfall                                                   [   kg/m2/s] #
# soil.water.XXX - Soil water content (XXX is the depth in cm)                [     m3/m3] #
# soil.temp.XXX  - Soil temperature (XXX is the depth in cm)                  [         K] #
# wdir           - Wind direction                                             [        --] #
# atm.vels       - Wind speed                                                 [       m/s] #
# atm.prss       - Air pressure                                               [        Pa] #
# atm.uspd       - Zonal wind                                                 [       m/s] #
# atm.vspd       - Meridional wind                                            [       m/s] #
# albedo         - Solar albedo                                               [        --] #
#------------------------------------------------------------------------------------------#
n         = 0
epro      = list()
n         = n + 1
epro[[n]] = list( vorig = "Year"
                , vname = "year"
                , type  = "integer"
                , mult  = NA_character_
                , add   = NA_character_
                )#end list
n         = n + 1
epro[[n]] = list( vorig = "Day"
                , vname = "doy"
                , type  = "integer"
                , mult  = NA_character_
                , add   = NA_character_
                )#end list
n         = n + 1
epro[[n]] = list( vorig = "Hour"
                , vname = "hfrac"
                , type  = "numeric"
                , mult  = "1./day.hr"
                , add   = NA_character_
                )#end list
n         = n + 1
epro[[n]] = list( vorig = "NEE_Thres_NEE_f"
                , vname = "nee.pmb"
                , type  = "numeric"
                , mult  = NA_character_
                , add   = NA_character_
                )#end list
n         = n + 1
epro[[n]] = list( vorig = "GPP_Thres_NEE_f"
                , vname = "gpp.pmb"
                , type  = "numeric"
                , mult  = NA_character_
                , add   = NA_character_
                )#end list
n         = n + 1
epro[[n]] = list( vorig = "Reco_Thres_NEE"
                , vname = "reco.pmb"
                , type  = "numeric"
                , mult  = NA_character_
                , add   = NA_character_
                )#end list
n         = n + 1
epro[[n]] = list( vorig = "LE"
                , vname = "flatent"
                , type  = "numeric"
                , mult  = NA_character_
                , add   = NA_character_
                )#end list
n         = n + 1
epro[[n]] = list( vorig = "H"
                , vname = "fsens"
                , type  = "numeric"
                , mult  = NA_character_
                , add   = NA_character_
                )#end list
n         = n + 1
epro[[n]] = list( vorig = "h2o_flux"
                , vname = "fh2o"
                , type  = "numeric"
                , mult  = "0.001*mmh2o"
                , add   = NA_character_
                )#end list
n         = n + 1
epro[[n]] = list( vorig = "Rg"
                , vname = "fgnd"
                , type  = "numeric"
                , mult  = NA_character_
                , add   = NA_character_
                )#end list
n         = n + 1
epro[[n]] = list( vorig = "Ustar"
                , vname = "ustar"
                , type  = "numeric"
                , mult  = NA_character_
                , add   = NA_character_
                )#end list
n         = n + 1
epro[[n]] = list( vorig = "Tau"
                , vname = "fmom"
                , type  = "numeric"
                , mult  = NA_character_
                , add   = NA_character_
                )#end list
n         = n + 1
epro[[n]] = list( vorig = "co2_flux"
                , vname = "fco2"
                , type  = "numeric"
                , mult  = NA_character_
                , add   = NA_character_
                )#end list
n         = n + 1
epro[[n]] = list( vorig = "H_strg"
                , vname = "storsens"
                , type  = "numeric"
                , mult  = NA_character_
                , add   = NA_character_
                )#end list
n         = n + 1
epro[[n]] = list( vorig = "h2o_strg"
                , vname = "storh2o"
                , type  = "numeric"
                , mult  = "0.001*mmh2o"
                , add   = NA_character_
                )#end list
n         = n + 1
epro[[n]] = list( vorig = "co2_strg"
                , vname = "storco2"
                , type  = "numeric"
                , mult  = NA_character_
                , add   = NA_character_
                )#end list
n         = n + 1
epro[[n]] = list( vorig = "sonic_temperature"
                , vname = "atm.tson"
                , type  = "numeric"
                , mult  = NA_character_
                , add   = NA_character_
                )#end list
n         = n + 1
epro[[n]] = list( vorig = "air_temperature"
                , vname = "atm.tmp"
                , type  = "numeric"
                , mult  = NA_character_
                , add   = NA_character_
                )#end list
n         = n + 1
epro[[n]] = list( vorig = "air_pressure"
                , vname = "atm.prss"
                , type  = "numeric"
                , mult  = NA_character_
                , add   = NA_character_
                )#end list
n         = n + 1
epro[[n]] = list( vorig = "e"
                , vname = "atm.pvap"
                , type  = "numeric"
                , mult  = NA_character_
                , add   = NA_character_
                )#end list
n         = n + 1
epro[[n]] = list( vorig = "co2_mixing_ratio"
                , vname = "co2"
                , type  = "numeric"
                , mult  = NA_character_
                , add   = NA_character_
                )#end list
n         = n + 1
epro[[n]] = list( vorig = "specific_humidity"
                , vname = "atm.shv"
                , type  = "numeric"
                , mult  = NA_character_
                , add   = NA_character_
                )#end list
n         = n + 1
epro[[n]] = list( vorig = "Tdew"
                , vname = "atm.tdew"
                , type  = "numeric"
                , mult  = NA_character_
                , add   = NA_character_
                )#end list
n         = n + 1
epro[[n]] = list( vorig = "wind_speed"
                , vname = "atm.vels"
                , type  = "numeric"
                , mult  = NA_character_
                , add   = NA_character_
                )#end list
n         = n + 1
epro[[n]] = list( vorig = "wind_dir"
                , vname = "atm.vdir"
                , type  = "numeric"
                , mult  = NA_character_
                , add   = NA_character_
                )#end list
n         = n + 1
epro[[n]] = list( vorig = "LWin_1_1_1"
                , vname = "rlong.in"
                , type  = "numeric"
                , mult  = NA_character_
                , add   = NA_character_
                )#end list
n         = n + 1
epro[[n]] = list( vorig = "LWout_1_1_1"
                , vname = "rlong.out"
                , type  = "numeric"
                , mult  = NA_character_
                , add   = NA_character_
                )#end list
n         = n + 1
epro[[n]] = list( vorig = "PPFD_1_1_1"
                , vname = "par.in"
                , type  = "numeric"
                , mult  = "Ein.2.Watts*1.e-6"
                , add   = NA_character_
                )#end list
n         = n + 1
epro[[n]] = list( vorig = "P_1_1_1"
                , vname = "rain"
                , type  = "numeric"
                , mult  = "2000./hr.sec"
                , add   = NA_character_
                )#end list
n         = n + 1
epro[[n]] = list( vorig = "SWin_1_1_1"
                , vname = "rshort.in"
                , type  = "numeric"
                , mult  = NA_character_
                , add   = NA_character_
                )#end list
n         = n + 1
epro[[n]] = list( vorig = "SWout_1_1_1"
                , vname = "rshort.out"
                , type  = "numeric"
                , mult  = NA_character_
                , add   = NA_character_
                )#end list
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#     Gap filling options.  The list of variables will be defined in a few lines below.    #
#------------------------------------------------------------------------------------------#
use.similar    = TRUE              # Should I use similar/duplicated measurements to fill 
                                   #   in missing data?
detrend.method = "mean"            # Which method to use to detrend series 
                                   #     (case insensitive)
                                   # - "loess":  local fitting technique; this is good for 
                                   #             most cases but the first and last point 
                                   #             must be defined. 
                                   # - "linear": simple linear detrending (requires 
                                   #             package RSEIS). this is the only option 
                                   #             in case you have gaps at the beginning or 
                                   #             end of the time series.
                                   # - "mean":  remove the mean, no true detrending.  It 
                                   #            may be the best option in case the trends
                                   #            are not clear or the time series is too
                                   #            short.
signal.retain  = 0.80              # Aimed fraction of signal to be maintained.  The ideal
                                   #    number is to retain as much signal as possible, but
                                   #    not adding everything because the faintest signals
                                   #    are typically noise.  A good first guess is 80%
conv.threshold = 0.001             # Tolerance for iterations at any given power
minmod         = 5                 # Maximum number of modes to consider
maxmod         = 600               # Maximum number of modes to consider
maxin          = 50                # Maximum number of inner iteractions
harmverb       = 1                 # How verbose the harmonic filling should be
                                   # 0 - silent
                                   # 1 - outer loop only
                                   # 2 - outer and inner loops
yeara.objana   = c(2008)           # First and last year to consider when obtaining the
yearz.objana   = c(2016)           #    statistics for the objective analysis.  Ideally 
                                   #    this should be the period that all data sets
                                   #    overlap, to avoid introducing biases due to inter-
                                   #    -annual variability.  The numbers are taken as 
                                   # options (first is the preferred period, then others
                                   # are tried in case the previous doesn't have statistics 
                                   # for all times needed.
yeara.harm     = c(2014)           # List periods to be filled before the entire time series
                                   #   is filled.  Splitting the time series may be a good
                                   #   option in case you have blocks with data and large
                                   #   gaps between these blocks (Ideally really you 
                                   #   shouldn't try to use harmonic analysis to fill 
                                   #   long periods because the amplitudes tend to decrease
                                   #   the longer the period without data).
yearz.harm     = c(2016)           #   - yeara.harm, beginning of each period
                                   #   - yearz.harm, end of each period
err.jackknife  = TRUE              # Estimate the error using jackknife?
del.jackknife  = 1./3.             # Fraction of data to be deleted for each jackknife 
                                   #     realisation.
max.jackknife  = 10                # Maximum number of iterations for jackknife error 
                                   #    estimation.
#------------------------------------------------------------------------------------------#




#------------------------------------------------------------------------------------------#
#     Minimum wind speed [m/s] and minimum u*.                                             #
#------------------------------------------------------------------------------------------#
atm.vels.min = 0.05
ustar.min    = 0.001
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#     Plot options.                                                                        #
#------------------------------------------------------------------------------------------#
plotregression = TRUE           # Plot the regressions for method 1 and NEE. 
plotgapfill    = TRUE           # Plot the gap filled time series.
outform        = c("pdf")       # Formats for output file (case insensitive):
                                #   - "X11" - for printing on screen
                                #   - "eps" - for postscript printing
                                #   - "png" - for PNG printing
                                #   - "pdf" - for PDF printing
byeold         = TRUE           # Remove old files of the given format?
depth          = 96             # PNG resolution, in pixels per inch
paper          = "square"       # Paper size, to define the plot shape
ptsz           = 18             # Font size.
plotgrid       = TRUE           # Should I plot the grid in the background?
legwhere       = "topleft"      # Where should I place the legend?
inset          = 0.01           # inset distance between legend and edge of plot box
legbg          = "white"        # Legend background colour.
maxgap         = 120            # Maximum number of plots
                                #    (if more gaps exist, they will be combined)
f.leg          = 1/6            # Fraction of the plotting area for legend.
#------------------------------------------------------------------------------------------#




#------------------------------------------------------------------------------------------#
#     Number of sub-periods for the time series that will be filled by the harmonic        #
# analysis.  The number is always going to added by one last loop in which the entire      #
# series is filled.                                                                        #
#------------------------------------------------------------------------------------------#
nharm = length(yeara.harm) + as.integer( length(yeara.harm) > 1)
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
n           = 0
vargap      = list()
n           = n + 1
vargap[[n]] = list( vari        = "atm.tmp"   
                  , desc        = "Air temperature"
                  , unit        = "degC"
                  , add0        = "-t00"
                  , mult        = NA_character_
                  , harm        = TRUE
                  , harm.tback  = TRUE
                  , oban        = TRUE
                  , fspline     = TRUE
                  , out.hour    = TRUE
                  , out.all     = FALSE
                  , drift       = FALSE
                  , quadratic   = FALSE
                  , sw.typevar  = FALSE
                  , potential   = NA_character_
                  , whena.drift = c("01/01/2014")
                  , whenz.drift = c("01/01/2017")
                  , form        = NA_character_
                  ) #end list
n           = n + 1
vargap[[n]] = list( vari        = "atm.tson"
                  , desc        = "Sonic temperature"
                  , unit        = "degC"
                  , add0        = "-t00"
                  , mult        = NA_character_
                  , harm        = FALSE
                  , harm.tback  = TRUE
                  , oban        = FALSE
                  , fspline     = FALSE
                  , out.hour    = TRUE
                  , out.all     = FALSE
                  , drift       = FALSE
                  , quadratic   = FALSE
                  , sw.typevar  = FALSE
                  , potential   = NA_character_
                  , whena.drift = c("01/01/2014")
                  , whenz.drift = c("01/01/2017")
                  , form        = NA_character_
                  ) #end list
n           = n + 1
vargap[[n]] = list( vari        = "atm.pvap"
                  , desc        = "Vapour pressure"
                  , unit        = "hpa"
                  , add0        = NA_character_
                  , mult        = "0.01"
                  , harm        = TRUE
                  , harm.tback  = TRUE
                  , oban        = TRUE
                  , fspline     = TRUE
                  , out.hour    = TRUE
                  , out.all     = FALSE
                  , drift       = FALSE
                  , quadratic   = FALSE
                  , sw.typevar  = FALSE
                  , potential   = NA_character_
                  , whena.drift = c("01/01/2014")
                  , whenz.drift = c("01/01/2017")
                  , form        = NA_character_
                  ) #end list
n           = n + 1
vargap[[n]] = list( vari        = "atm.prss"
                  , desc        = "Atmospheric pressure"
                  , unit        = "hpa"
                  , add0        = NA_character_
                  , mult        = "0.01"
                  , harm        = TRUE
                  , harm.tback  = TRUE
                  , oban        = TRUE
                  , fspline     = FALSE
                  , out.hour    = TRUE
                  , out.all     = FALSE
                  , drift       = FALSE
                  , quadratic   = FALSE
                  , sw.typevar  = FALSE
                  , potential   = NA_character_
                  , whena.drift = c("01/01/2014")
                  , whenz.drift = c("01/01/2017")
                  , form        = NA_character_
                  ) #end list
n           = n + 1
vargap[[n]] = list( vari        = "atm.uspd"
                  , desc        = "Zonal wind"
                  , unit        = "mos"
                  , add0        = NA_character_
                  , mult        = NA_character_
                  , harm        = TRUE
                  , harm.tback  = TRUE
                  , oban        = FALSE
                  , fspline     = TRUE
                  , out.hour    = FALSE
                  , out.all     = FALSE
                  , drift       = FALSE
                  , quadratic   = FALSE
                  , sw.typevar  = FALSE
                  , potential   = NA_character_
                  , whena.drift = c("01/01/2014")
                  , whenz.drift = c("01/01/2017")
                  , form        = NA_character_
                  ) #end list
n           = n + 1
vargap[[n]] = list( vari        = "atm.vels"
                  , desc        = "Wind speed"
                  , unit        = "mos"
                  , add0        = NA_character_
                  , mult        = NA_character_
                  , harm        = TRUE
                  , harm.tback  = TRUE
                  , oban        = TRUE
                  , fspline     = TRUE
                  , out.hour    = TRUE
                  , out.all     = FALSE
                  , drift       = FALSE
                  , quadratic   = FALSE
                  , sw.typevar  = FALSE
                  , potential   = NA_character_
                  , whena.drift = c("01/01/2014")
                  , whenz.drift = c("01/01/2017")
                  , form        = NA_character_
                  ) #end list
n           = n + 1
vargap[[n]] = list( vari        = "atm.vspd"
                  , desc        = "Meridional wind"
                  , unit        = "mos"
                  , add0        = NA_character_
                  , mult        = NA_character_
                  , harm        = TRUE
                  , harm.tback  = TRUE
                  , oban        = FALSE
                  , fspline     = TRUE
                  , out.hour    = FALSE
                  , out.all     = FALSE
                  , drift       = FALSE
                  , quadratic   = FALSE
                  , sw.typevar  = FALSE
                  , potential   = NA_character_
                  , whena.drift = c("01/01/2014")
                  , whenz.drift = c("01/01/2017")
                  , form        = NA_character_
                  ) #end list
n           = n + 1
vargap[[n]] = list( vari        = "rnet"
                  , desc        = "Net radiation"
                  , unit        = "wom2"
                  , add0        = NA_character_
                  , mult        = NA_character_
                  , harm        = FALSE
                  , harm.tback  = TRUE
                  , oban        = FALSE
                  , fspline     = FALSE
                  , out.hour    = FALSE
                  , out.all     = FALSE
                  , drift       = FALSE
                  , quadratic   = FALSE
                  , sw.typevar  = FALSE
                  , potential   = NA_character_
                  , whena.drift = c("01/01/2014")
                  , whenz.drift = c("01/01/2017")
                  , form        = NA_character_
                  ) #end list
n           = n + 1
vargap[[n]] = list( vari        = "rshort.in"
                  , desc        = "Incoming shortwave radiation"
                  , unit        = "wom2"
                  , add0        = NA_character_
                  , mult        = NA_character_
                  , harm        = TRUE
                  , harm.tback  = TRUE
                  , oban        = TRUE
                  , fspline     = TRUE
                  , out.hour    = TRUE
                  , out.all     = FALSE
                  , drift       = FALSE
                  , quadratic   = FALSE
                  , sw.typevar  = TRUE
                  , potential   = "rshort.pot"
                  , whena.drift = c("01/01/2014")
                  , whenz.drift = c("01/01/2017")
                  , form        = NA_character_
                  ) #end list
n           = n + 1
vargap[[n]] = list( vari        = "rshort.out"
                  , desc        = "Outgoing shortwave radiation"
                  , unit        = "wom2"
                  , add0        = NA_character_
                  , mult        = NA_character_
                  , harm        = FALSE
                  , harm.tback  = TRUE
                  , oban        = FALSE
                  , fspline     = FALSE
                  , out.hour    = TRUE
                  , out.all     = FALSE
                  , drift       = FALSE
                  , quadratic   = TRUE
                  , sw.typevar  = TRUE
                  , potential   = "rshort.pot"
                  , whena.drift = c("01/01/2014")
                  , whenz.drift = c("01/01/2017")
                  , form        = NA_character_
                  ) #end list
n           = n + 1
vargap[[n]] = list( vari        = "rlong.in"
                  , desc        = "Incoming longwave radiation"
                  , unit        = "wom2"
                  , add0        = NA_character_
                  , mult        = NA_character_
                  , harm        = FALSE
                  , harm.tback  = TRUE
                  , oban        = TRUE
                  , fspline     = FALSE
                  , out.hour    = TRUE
                  , out.all     = FALSE
                  , drift       = FALSE
                  , quadratic   = FALSE
                  , sw.typevar  = FALSE
                  , potential   = NA_character_
                  , whena.drift = c("01/01/2014")
                  , whenz.drift = c("01/01/2017")
                  , form        = NA_character_
                  ) #end list
n           = n + 1
vargap[[n]] = list( vari        = "rlong.out"
                  , desc        = "Outgoing longwave radiation"
                  , unit        = "wom2"
                  , add0        = NA_character_
                  , mult        = NA_character_
                  , harm        = FALSE
                  , harm.tback  = TRUE
                  , oban        = FALSE
                  , fspline     = FALSE
                  , out.hour    = TRUE
                  , out.all     = FALSE
                  , drift       = FALSE
                  , quadratic   = FALSE
                  , sw.typevar  = FALSE
                  , potential   = NA_character_
                  , whena.drift = c("01/01/2014")
                  , whenz.drift = c("01/01/2017")
                  , form        = NA_character_
                  ) #end list
n           = n + 1
vargap[[n]] = list( vari        = "rain"
                  , desc        = "Precipitation rate"
                  , unit        = "kgwom2ohr"
                  , add0        = NA_character_
                  , mult        = "hr.sec"
                  , harm        = FALSE
                  , harm.tback  = TRUE
                  , oban        = TRUE
                  , fspline     = FALSE
                  , out.hour    = FALSE
                  , out.all     = FALSE
                  , drift       = FALSE
                  , quadratic   = FALSE
                  , sw.typevar  = FALSE
                  , potential   = NA_character_
                  , whena.drift = c("01/01/2014")
                  , whenz.drift = c("01/01/2017")
                  , form        = NA_character_
                  ) #end list
n           = n + 1
vargap[[n]] = list( vari        = "ustar"
                  , desc        = "Friction velocity"
                  , unit        = "mos"
                  , add0        = NA_character_
                  , mult        = NA_character_
                  , harm        = TRUE
                  , harm.tback  = TRUE
                  , oban        = FALSE
                  , fspline     = FALSE
                  , out.hour    = TRUE
                  , out.all     = FALSE
                  , drift       = FALSE
                  , quadratic   = FALSE
                  , sw.typevar  = FALSE
                  , potential   = NA_character_
                  , whena.drift = c("01/01/2014")
                  , whenz.drift = c("01/01/2017")
                  , form        = NA_character_
                  ) #end list
n           = n + 1
vargap[[n]] = list( vari        = "fco2"
                  , desc        = "CO2 flux"
                  , unit        = "umolom2os"
                  , add0        = NA_character_
                  , mult        = NA_character_
                  , harm        = FALSE
                  , harm.tback  = TRUE
                  , oban        = FALSE
                  , fspline     = FALSE
                  , out.hour    = TRUE
                  , out.all     = FALSE
                  , drift       = FALSE
                  , quadratic   = FALSE
                  , sw.typevar  = FALSE
                  , potential   = NA_character_
                  , whena.drift = c("01/01/2014")
                  , whenz.drift = c("01/01/2017")
                  , form        = NA_character_
                  ) #end list
n           = n + 1
vargap[[n]] = list( vari        = "fsens"
                  , desc        = "Sensible heat flux"
                  , unit        = "wom2"
                  , add0        = NA_character_
                  , mult        = NA_character_
                  , harm        = FALSE
                  , harm.tback  = TRUE
                  , oban        = FALSE
                  , fspline     = FALSE
                  , out.hour    = TRUE
                  , out.all     = FALSE
                  , drift       = FALSE
                  , quadratic   = FALSE
                  , sw.typevar  = FALSE
                  , potential   = NA_character_
                  , whena.drift = c("01/01/2014")
                  , whenz.drift = c("01/01/2017")
                  , form        = NA_character_
                  ) #end list
n           = n + 1
vargap[[n]] = list( vari        = "storsens"
                  , desc        = "Sensible heat storage"
                  , unit        = "wom2"
                  , add0        = NA_character_
                  , mult        = NA_character_
                  , harm        = FALSE
                  , harm.tback  = TRUE
                  , oban        = FALSE
                  , fspline     = FALSE
                  , out.hour    = TRUE
                  , out.all     = FALSE
                  , drift       = FALSE
                  , quadratic   = FALSE
                  , sw.typevar  = FALSE
                  , potential   = NA_character_
                  , whena.drift = c("01/01/2014")
                  , whenz.drift = c("01/01/2017")
                  , form        = NA_character_
                  ) #end list
n           = n + 1
vargap[[n]] = list( vari        = "fh2o"
                  , desc        = "Water flux"
                  , unit        = "kgom2oday"
                  , add0        = NA_character_
                  , mult        = "day.sec"
                  , harm        = FALSE
                  , harm.tback  = TRUE
                  , oban        = FALSE
                  , fspline     = FALSE
                  , out.hour    = TRUE
                  , out.all     = FALSE
                  , drift       = FALSE
                  , quadratic   = FALSE
                  , sw.typevar  = FALSE
                  , potential   = NA_character_
                  , whena.drift = c("01/01/2014")
                  , whenz.drift = c("01/01/2017")
                  , form        = NA_character_
                  ) #end list
n           = n + 1
vargap[[n]] = list( vari        = "storh2o"
                  , desc        = "Water storage"
                  , unit        = "kgom2oday"
                  , add0        = NA_character_
                  , mult        = "day.sec"
                  , harm        = FALSE
                  , harm.tback  = TRUE
                  , oban        = FALSE
                  , fspline     = FALSE
                  , out.hour    = TRUE
                  , out.all     = FALSE
                  , drift       = FALSE
                  , quadratic   = FALSE
                  , sw.typevar  = FALSE
                  , potential   = NA_character_
                  , whena.drift = c("01/01/2014")
                  , whenz.drift = c("01/01/2017")
                  , form        = NA_character_
                  ) #end list
n           = n + 1
vargap[[n]] = list( vari        = "storco2"
                  , desc        = "CO2 Canopy storage"
                  , unit        = "umolom2os"
                  , add0        = NA_character_
                  , mult        = NA_character_
                  , harm        = FALSE
                  , harm.tback  = TRUE
                  , oban        = FALSE
                  , fspline     = FALSE
                  , out.hour    = TRUE
                  , out.all     = FALSE
                  , drift       = FALSE
                  , quadratic   = FALSE
                  , sw.typevar  = FALSE
                  , potential   = NA_character_
                  , whena.drift = c("01/01/2014")
                  , whenz.drift = c("01/01/2017")
                  , form        = NA_character_
                  ) #end list
n           = n + 1
vargap[[n]] = list( vari        = "co2"
                  , desc        = "CO2 concentration"
                  , unit        = "umolomol"
                  , add0        = NA_character_
                  , mult        = NA_character_
                  , harm        = FALSE
                  , harm.tback  = TRUE
                  , oban        = FALSE
                  , fspline     = FALSE
                  , out.hour    = TRUE
                  , out.all     = FALSE
                  , drift       = FALSE
                  , quadratic   = FALSE
                  , sw.typevar  = FALSE
                  , potential   = NA_character_
                  , whena.drift = c("01/01/2014")
                  , whenz.drift = c("01/01/2017")
                  , form        = NA_character_
                  ) #end list
n           = n + 1
vargap[[n]] = list( vari        = "par.in"
                  , desc        = "Incoming PAR"
                  , unit        = "umolom2os"
                  , add0        = NA_character_
                  , mult        = "Watts.2.Ein * 1e+6"
                  , harm        = FALSE
                  , harm.tback  = TRUE
                  , oban        = FALSE
                  , fspline     = FALSE
                  , out.hour    = TRUE
                  , out.all     = FALSE
                  , drift       = FALSE
                  , quadratic   = FALSE
                  , sw.typevar  = TRUE
                  , potential   = "par.pot"
                  , whena.drift = c("01/01/2014")
                  , whenz.drift = c("01/01/2017")
                  , form        = NA_character_
                  ) #end list
n           = n + 1
vargap[[n]] = list( vari        = "par.out"
                  , desc        = "Outgoing PAR"
                  , unit        = "umolom2os"
                  , add0        = NA_character_
                  , mult        = "Watts.2.Ein * 1e+6"
                  , harm        = FALSE
                  , harm.tback  = TRUE
                  , oban        = FALSE
                  , fspline     = FALSE
                  , out.hour    = TRUE
                  , out.all     = FALSE
                  , drift       = FALSE
                  , quadratic   = FALSE
                  , sw.typevar  = TRUE
                  , potential   = "par.pot"
                  , whena.drift = c("01/01/2014")
                  , whenz.drift = c("01/01/2017")
                  , form        = NA_character_
                  ) #end list
n           = n + 1
vargap[[n]] = list( vari        = "nee"
                  , desc        = "Net ecosystem exchange"
                  , unit        = "umolom2os"
                  , add0        = NA_character_
                  , mult        = NA_character_
                  , harm        = FALSE
                  , harm.tback  = TRUE
                  , oban        = FALSE
                  , fspline     = FALSE
                  , out.hour    = FALSE
                  , out.all     = FALSE
                  , drift       = FALSE
                  , quadratic   = FALSE
                  , sw.typevar  = FALSE
                  , potential   = NA_character_
                  , whena.drift = c("01/01/2014")
                  , whenz.drift = c("01/01/2017")
                  , form        = NA_character_
                  ) #end list
n           = n + 1
vargap[[n]] = list( vari        = "gpp.pmb"
                  , desc        = "Gross primary productivity - Brando"
                  , unit        = "umolom2os"
                  , add0        = NA_character_
                  , mult        = NA_character_
                  , harm        = FALSE
                  , harm.tback  = TRUE
                  , oban        = FALSE
                  , fspline     = FALSE
                  , out.hour    = FALSE
                  , out.all     = FALSE
                  , drift       = FALSE
                  , quadratic   = FALSE
                  , sw.typevar  = FALSE
                  , potential   = NA_character_
                  , whena.drift = c("01/01/2014")
                  , whenz.drift = c("01/01/2017")
                  , form        = NA_character_
                  ) #end list
n           = n + 1
vargap[[n]] = list( vari        = "reco.pmb"
                  , desc        = "Ecosystem respiration - Brando"
                  , unit        = "umolom2os"
                  , add0        = NA_character_
                  , mult        = NA_character_
                  , harm        = FALSE
                  , harm.tback  = TRUE
                  , oban        = FALSE
                  , fspline     = FALSE
                  , out.hour    = FALSE
                  , out.all     = FALSE
                  , drift       = FALSE
                  , quadratic   = FALSE
                  , sw.typevar  = FALSE
                  , potential   = NA_character_
                  , whena.drift = c("01/01/2014")
                  , whenz.drift = c("01/01/2017")
                  , form        = NA_character_
                  ) #end list
n           = n + 1
vargap[[n]] = list( vari        = "nee.pmb"
                  , desc        = "Net Ecosystem Exchange - Brando"
                  , unit        = "umolom2os"
                  , add0        = NA_character_
                  , mult        = NA_character_
                  , harm        = FALSE
                  , harm.tback  = TRUE
                  , oban        = FALSE
                  , fspline     = FALSE
                  , out.hour    = FALSE
                  , out.all     = FALSE
                  , drift       = FALSE
                  , quadratic   = FALSE
                  , sw.typevar  = FALSE
                  , potential   = NA_character_
                  , whena.drift = c("01/01/2014")
                  , whenz.drift = c("01/01/2017")
                  , form        = NA_character_
                  ) #end list
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




#----- Count the number of variables to be checked for gap-filling. -----------------------#
vargap   = list.2.data.frame(vargap)
nvargap  = nrow(vargap)
#------------------------------------------------------------------------------------------#



#----- Create output directory in case it doesn't exist. ----------------------------------#
if (! file.exists(outroot)) dir.create(outroot)
#------------------------------------------------------------------------------------------#



#----- Determine how many formats we must output. -----------------------------------------#
outform = tolower(outform)
nout    = length(outform)
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#     List of variables related to shortwave radiation.  These variables must be set to    #
# zero during nighttime.                                                                   #
#------------------------------------------------------------------------------------------#
radvars = c("rshort.in","rshort.out","par.in","par.out")
#------------------------------------------------------------------------------------------#


#------------------------------------------------------------------------------------------#
#    Limits for the combined dataset.                                                      #
#------------------------------------------------------------------------------------------#
whena = chron(dates=paste( 1, 1,yeara,sep="/"),times=paste( 0, 0, 0,sep=":"))
whenz = chron(dates=paste(12,31,yearz,sep="/"),times=paste(23, 0, 0,sep=":"))
#------------------------------------------------------------------------------------------#



#----- Define the R Data file name with raw and gap filled data. --------------------------#
session.raw     = file.path(here,raw.rdata    )
session.clean   = file.path(here,clean.rdata  )
session.selfgf  = file.path(here,selfgf.rdata )
session.objana  = file.path(here,objana.rdata )
session.harmana = file.path(here,harmana.rdata)
session.fspline = file.path(here,fspline.rdata)
session.metfill = file.path(here,metfill.rdata)
#------------------------------------------------------------------------------------------#


#----- Full file name for ancillary data set. ---------------------------------------------#
other.rdata = file.path(other.path,other.base)
#------------------------------------------------------------------------------------------#


#----- Transform lists into data frames. --------------------------------------------------#
epro    = list.2.data.frame(epro )
nepro  = nrow(epro )
#------------------------------------------------------------------------------------------#

#----- Find the number of objective analysis periods to try. ------------------------------#
if (length(yeara.objana) != length(yearz.objana)){
   cat0(" --> Length of yeara.objana: ",length(yeara.objana),".")
   cat0(" --> Length of yearz.objana: ",length(yearz.objana),".")
   stop(" Check settings: 'yeara.objana' and 'yearz.objana' must have the same length!")
}else{
   nobjana = length(yeara.objana)
}#end if (length(yeara.objana) != length(yearz.objana))
#------------------------------------------------------------------------------------------#


#----- Find the cosine threshold for drift correction and gap filling. --------------------#
cosz.thresh   = max(cosz.min, frac.cosz * cos((abs(pi*lat/180)+abs(capri))))
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
#    Load data from nearby towers, so we can run an objective analysis for some of the     #
# variables.                                                                               #
#------------------------------------------------------------------------------------------#
other.full = file.path(here,other.rdata)
if (! file.exists(other.full)){
   stop(paste0(" + File ",other.rdata," was not found!"))
}else{
   cat0(" - Load data from nearby sites.")
   load(file=other.full)
   nother = length(fitz.aws)
}#end if
#---------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#    Read in the data.                                                                     #
#------------------------------------------------------------------------------------------#
if (reload.raw && file.exists(session.raw)){
   cat0(" - Load raw data from file ",basename(session.raw),".")
   load(session.raw)
   raw    = eft
   ntimes = length(raw$when)
}else{
   cat0(" - Read data from file ",eddyfile,".")
   nlines = as.integer( system( command = paste0("cat ",eddyfile," | wc -l")
                              , intern  = TRUE
                              )#end system
                      )#end as.integer
   eftlt  = read.csv( file             = eddyfile
                    , header           = TRUE
                    , sep              = ","
                    , comment.char     = ""
                    , nrows            = nlines
                    , colClasses       = c("character",rep("numeric",64),rep("character",2)
                                          ,rep("numeric",135),rep("character")
                                          ,rep("numeric",7),rep("character",2)
                                          ,rep("numeric",4))
                    , stringsAsFactors = FALSE
                    )#end read.table
   #---------------------------------------------------------------------------------------#

   #----- Define time zero to ensure proper mapping. --------------------------------------#
   zero     = chron(paste(12,31,yeara-1,sep="/"))
   #---------------------------------------------------------------------------------------#


   #----- Throw away columns that we won't use. -------------------------------------------#
   eftlt        = eftlt[,epro$vorig,drop=FALSE]
   names(eftlt) = epro$vname
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #      Standardise scale of variables that must be scaled.                              #
   #---------------------------------------------------------------------------------------#
   cat0(" - Standardise variable scales:")
   for (e in which(epro$type %in% "numeric")){
      e.vnam          = epro$vname[e]
      e.mult          = epro$mult [e]
      e.add           = epro$add  [e]
      if (is.na(e.mult)) e.mult = 1. else e.mult = eval(parse(text=e.mult))
      if (is.na(e.add )) e.add  = 0. else e.add  = eval(parse(text=e.add ))
      eftlt[[e.vnam]] = e.add + e.mult * eftlt[[e.vnam]]
      cat0("     * ",e.vnam,".")
   }#end for (e in which(epro$type %in% "numeric"))
   #---------------------------------------------------------------------------------------#
         
   #---------------------------------------------------------------------------------------#
   #       Create time variable.                                                           #
   #---------------------------------------------------------------------------------------#
   eftlt$when    = chron( chron(paste(12,31,eftlt$year-1,sep="/"))
                        + eftlt$doy + eftlt$hfrac - off.utc / day.hr
                        )#end chron
   dt.scale      = hr.sec / dtdat
   eftlt$elapsed = round(day.hr * dt.scale * (eftlt$when - zero)) / dt.scale
   inp.fac       = ceiling(eftlt$elapsed)
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #     Check that all variables that must exist are assigned in 'raw'.  In case any      #
   # variable is missing, make a dummy vector.                                             #
   #---------------------------------------------------------------------------------------#
   cat0("   - Ensure that all critical variables exist in 'raw'.")
   for (v in which(! vargap$vari %in% names(eftlt))){
      v.must          = vargap$vari[v]
      eftlt[[v.must]] = rep(NA,times=length(eftlt$when))
   }#end for (m in which(! must.exist %in% names(eftlt))
   #---------------------------------------------------------------------------------------#




   #----- Define the time vector that will be used by the full time series. ---------------#
   cat0(" - Define the time vector.")
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
   dt.scale = hr.sec / dtdat
   elapsed  = round(day.hr * dt.scale * (when - zero)) / dt.scale
   #---------------------------------------------------------------------------------------#


   #----- Start populating the data frame which will hold all variables. ------------------#
   cat0(" - Initialise the data frame.")
   eft     = data.frame(when=chron(when),elapsed=elapsed)
   eft.fac = ceiling(eft$elapsed)
   nwhen   = length(eft$when)
   empty   = rep(NA,times=nwhen)
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #      For all the other variables, find the averages.                                  # 
   #---------------------------------------------------------------------------------------#
   cat0(" - Average the hours:")
   avgvar = epro$vname[! epro$vname %in% c("year","dfrac","doy","hhmm")]
   navgs  = length(avgvar)
   for (a in sequence(navgs)){
      v.now = avgvar[a]
      cat0("   * ",v.now,".")
      mean.now               = tapply(X = eftlt[[v.now]], INDEX=inp.fac, FUN=mean)
      fac.now                = as.numeric(names(mean.now))
      idx                    = match(fac.now,eft.fac)
      sel                    = is.finite(idx)
      eft[[v.now]]           = empty
      eft[[v.now]][idx[sel]] = mean.now[sel]
   }#end for (a in sequence(navgs))
   #---------------------------------------------------------------------------------------#



   #----- Recalculate all time stuff. -----------------------------------------------------#
   cat0(" - Make derived time information.")
   eft       = alltimes(datin=eft,lon=lon,lat=lat,ed21=TRUE,zeronight=FALSE,meanval=TRUE
                       ,imetavg=1,nmean=60,na.rm=TRUE)
   #---------------------------------------------------------------------------------------#



   #----- Find auxiliary variables that will help to make the combined dataset. -----------#
   ntimes    = length(eft$when)
   empty     = rep(NA,times=ntimes)
   #---------------------------------------------------------------------------------------#



   #----- Calculate partial pressure of water vapour. -------------------------------------#
   cat0("   - Find vapour pressure.")
   fill         = is.na(eft$atm.pvap)
   if ("atm.tdew" %in% names(eft)){
      eft$atm.pvap[fill] = eslif(eft$atm.tdew[fill])
      fill               = is.na(eft$atm.pvap)
   }#end if("atm.tdew" %in% names(eft))
   if ("atm.rhv" %in% names(eft)){
      eft$atm.pvap[fill] = eft$atm.rhv[fill] * eslif(eft$atm.tmp[fill])
      fill               = is.na(eft$atm.pvap)
   }#end if("atm.rhv" %in% names(eft))
   if ("atm.shv" %in% names(eft)){
      eft$atm.pvap[fill] = ( eft$atm.prss[fill] * eft$atm.shv[fill]
                           / ( ep + (1. - ep)*eft$atm.shv[fill] )
                           )#end atm.pvap
   }#end if ("atm.shv" %in% names(eft))
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #      Find wind speed components (so we can interpolate direction if needed be.        #
   #---------------------------------------------------------------------------------------#
   eft$atm.uspd = - eft$atm.vels * sin(eft$atm.vdir * pio180)
   eft$atm.vspd = - eft$atm.vels * cos(eft$atm.vdir * pio180)
   #---------------------------------------------------------------------------------------#



   #------ Find net radiation. ------------------------------------------------------------#
   eft$rnet = eft$rshort.in + eft$rlong.in - eft$rshort.out - eft$rlong.out
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #    Find the lowest pressure amongst all measurements, and use it as the reference     #
   # pressure, to be on the conservative side.  Then we create a 1-minute time series and  #
   # find the cosine of the zenith angle, then find the maximum possible radiation for     #
   # every minute, then make the hourly average to find the potential maximum.             #
   #---------------------------------------------------------------------------------------#
   cat(" + Find the 1-minute maximum radiation...","\n")
   onemin           = list()
   i.01m            = seq(from=1,to=ntimes*60,by=1)-60
   index            = ceiling((i.01m-0.5)/60)
   onemin$when      = whena + i.01m / day.min
   n.01m            = length(onemin$when)
   onemin           = alltimes(datin=onemin,lon=lon,lat=lat
                              ,ed21=TRUE,zeronight=FALSE,meanval=FALSE,imetavg=1
                              ,nmean=1,na.rm=TRUE)
   onemin$rshort.in = rep(eft$rshort.in,each=60)
   sunny            = rshort.bdown( rad.in   = onemin$rshort.in
                                  , atm.prss = ref.prss + 0. * onemin$cosz
                                  , cosz     = onemin$cosz
                                  , rad.type = "rshort")
   eft$rshort.pot   = tapply(X=sunny$rshort.max,INDEX=index,FUN=mean,na.rm=TRUE)
   eft$par.pot      = tapply(X=sunny$par.max   ,INDEX=index,FUN=mean,na.rm=TRUE)
   rm(onemin,sunny,i.01m,index,n.01m)
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #     Save the compounded dataset so we don't load them again.                          #
   #---------------------------------------------------------------------------------------#
   cat0(" + Save tower data to file ",basename(session.raw),".")
   save(eft,file=session.raw)
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #    Copy eft to raw.  From this point on, the raw data will not be called eft, but     #
   # raw because we always use eft as the structure name that is saved to the RData.       #
   #---------------------------------------------------------------------------------------#
   raw = eft
   #---------------------------------------------------------------------------------------#

}#end if
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
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#     Check whether to run the data QA/QC or load previously computed clean time series.   #
#------------------------------------------------------------------------------------------#
if (reload.clean && file.exists(session.clean)){
   #----- Reload the gap-filling. ---------------------------------------------------------#
   cat0(" + Load clean data from file ",basename(session.clean),".")
   load(session.clean)
   clean   = eft
   ntimes  = length(clean$when)
   #---------------------------------------------------------------------------------------#

}else{

   #---------------------------------------------------------------------------------------#
   #     Clean the raw time series.  In this section we perform the following tasks.       #
   # 1. Eliminate obviously bad data                                                       #
   # 2. Eliminate outliers                                                                 #
   # 3. Correct drift and poor calibration (run outlier elimination again if needed).      #
   # 4. Gap fill using redundant variables (i.e., using variables that are highly          #
   #    autocorrelated and from measurements from the same site... Most gap filling will   #
   #    happen in the next block.                                                          #
   #---------------------------------------------------------------------------------------#


   #----- Initialise the data frame that will contain the clean dataset. ------------------#
   eft     = raw
   ntimes  = length(eft$when)
   #---------------------------------------------------------------------------------------#



   #=======================================================================================#
   #=======================================================================================#
   #     Loop over all the variables that will be gap-filled, and assign an initial gap-   #
   # filling flag and error.  The flags will be always as follow:                          #
   #                                                                                       #
   #  0 -- Same value as in clean                                                          #
   # -1 -- Gap filled using redundant variables.                                           #
   #  1 -- Gap filled using objective analysis.                                            #
   #  2 -- Gap filled using harmonic analysis.                                             #
   #  3 -- Gap filled using parametric analysis.                                           #
   #  4 -- Gap filled using alternative, variable-specific method.                         #
   #  5 -- NEE gap filled using linear, multi-variable model.                              #
   #  6 -- NEE gap filled by nighttime respiration interpolation                           #
   #  7 -- NEE gap filled by harmonic analysis.                                            #
   # NA -- Variable that could not be gap filled.                                          #
   #                                                                                       #
   #    The error will be 0 for observed variable, otherwise it will contain an estimate   #
   # of the error associated with gap filling.  Some of the errors may be skipped so they  #
   # will remain NA.                                                                       #
   #---------------------------------------------------------------------------------------#
   for (vv in sequence(nvargap)){
      #---- Handy aliases. ----------------------------------------------------------------#
      this.vnam       = vargap$vari      [vv]
      this.desc       = vargap$desc      [vv]
      this.unit       = vargap$unit      [vv]
      this.harm       = vargap$harm      [vv]
      this.harm.tback = vargap$harm.tback[vv]
      this.oban       = vargap$oban      [vv]
      this.out.hour   = vargap$out.hour  [vv]
      this.form       = vargap$form      [vv]
      #------------------------------------------------------------------------------------#



      #----- Make a new variable containing both the gap filling flag and the error. ------#
      gfflg = paste0("gfflg.",this.vnam)
      gferr = paste0("gferr.",this.vnam)
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #     Initialise the gap filling flag and estimate by setting all points as NA.      #
      #------------------------------------------------------------------------------------#
      eft[[gfflg]] = rep(NA,times=ntimes)
      eft[[gferr]] = rep(NA,times=ntimes)
      #------------------------------------------------------------------------------------#
   }#end for (vv in sequence(nvargap))
   #=======================================================================================#
   #=======================================================================================#



   #---------------------------------------------------------------------------------------#
   #     Initialise the list that will contain the drift and bias corrections, and the     #
   # simple gap filling using redundant data.                                              #
   #---------------------------------------------------------------------------------------#
   fixlist = list()
   #---------------------------------------------------------------------------------------#



   #=======================================================================================#
   #=======================================================================================#
   #     Discard drifting values.                                                          #
   #---------------------------------------------------------------------------------------#
   cat0(" + Discard obviously bad periods.")
      #------------------------------------------------------------------------------------#
      #   1. Pressure had some weird behaviour in May and June 2014.                       #
      #------------------------------------------------------------------------------------#
      sel = eft$when %wr% chron(c("05/01/2014","07/01/2014"))
      eft$atm.prss [sel] = NA
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #   2. Rain turned out to be zero for most of 2016.  It was a dry year, but zero     #
      #      from March through December is rather unlikely -- INMET sites not too far     #
      #      Tanguro had rain.  Remove all data from 16 Feb 2016 onwards.                  #
      #------------------------------------------------------------------------------------#
      sel = eft$when %>=% chron("02/16/2016")
      eft$rain     [sel] = NA
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #   3. Sensible heat looks weird in August 2014 (too much oscillation).              #
      #------------------------------------------------------------------------------------#
      sel = eft$when %wr% chron(c("07/31/2014","09/10/2014"))
      eft$fsens[sel] = NA
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #   4. CO2 storage looks odd in the first months. CO2 flux is undefined too, so      #
      #      remove storage.  CO                                                              #
      #------------------------------------------------------------------------------------#
      sel = eft$when %<=% chron("06/01/2014") 
      eft$storco2 [sel] = NA
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #   5. CO2 has weird patterns at times.                                              #
      #------------------------------------------------------------------------------------#
      sel = eft$when %<=% chron("06/01/2014")
      eft$co2 [sel] = NA
      sel = eft$when %wr% chron(c("02/01/2015","03/01/2015"))
      eft$co2 [sel] = NA
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #   6. Humidity has a weird spike in one day.                                        #
      #------------------------------------------------------------------------------------#
      sel = eft$when %wr% chron(c("04/05/2016","04/06/2016"))
      eft$atm.pvap[sel] = NA
      eft$atm.rhv [sel] = NA
      eft$atm.shv [sel] = NA
      eft$atm.tdew[sel] = NA
      #------------------------------------------------------------------------------------#
   #=======================================================================================#
   #=======================================================================================#




   #=======================================================================================#
   #=======================================================================================#
   #   Make nighttime radiation temporarily NA so the skewed normal fitting works.  They   #
   # will be replaced by zeroes after the outlier test.                                    #
   #---------------------------------------------------------------------------------------#
   eft$rshort.in      [eft$nighttime] = NA
   eft$rshort.out     [eft$nighttime] = NA
   eft$par.in         [eft$nighttime] = NA
   eft$par.out        [eft$nighttime] = NA
   #=======================================================================================#
   #=======================================================================================#




   #=======================================================================================#
   #=======================================================================================#
   #      Run two filters to reject bad data and make them NA.  Not all variables should   #
   # be checked this way, and we must skip variables whose distribution is too far from    #
   # Gaussian at diurnal cycle or at the daily scale.                                      #
   #---------------------------------------------------------------------------------------#
   cat0(" + Remove outliers by hour and by day.")
   for (vv in sequence(nvargap)){
      #------------------------------------------------------------------------------------#
      #     Copy the settings for this variable to scratch variables.                      #
      #------------------------------------------------------------------------------------#
      this.vnam        = vargap$vari       [vv]
      this.desc        = vargap$desc       [vv]
      this.unit        = vargap$unit       [vv]
      this.harm        = vargap$harm       [vv]
      this.harm.tback  = vargap$harm.tback [vv]
      this.oban        = vargap$oban       [vv]
      this.out.hour    = vargap$out.hour   [vv]
      this.out.all     = vargap$out.all    [vv]
      this.drift       = vargap$drift      [vv]
      this.sw.typevar  = vargap$sw.typevar [vv]
      this.potential   = vargap$potential  [vv]
      this.whena.drift = vargap$whena.drift[vv]
      this.whenz.drift = vargap$whenz.drift[vv]
      this.form        = vargap$form       [vv]
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #     Remove outliers.                                                               #
      #------------------------------------------------------------------------------------#
      cat0("   - Variable:",this.vnam,".")
      thisvar          = eft[[this.vnam]]
      eft[[this.vnam]] = del.outliers( x        = thisvar
                                     , when     = eft$when
                                     , out.hour = this.out.hour
                                     , out.all  = this.out.all
                                     )#end del.outliers
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #    In case vapour pressure turns out to be bad, discard specific humidity,         #
      # relative humidity, and dew point temperature.                                      #
      #------------------------------------------------------------------------------------#
      if (this.vnam %in% c("atm.pvap")){
         discard               = is.finite(thisvar) & is.na(eft[[this.vnam]])
         eft$atm.shv [discard] = NA
         eft$atm.rhv [discard] = NA
         eft$atm.tdew[discard] = NA
      }#end if
      #------------------------------------------------------------------------------------#
   }#end for
   #=======================================================================================#
   #=======================================================================================#






   #=======================================================================================#
   #=======================================================================================#
   #      Delete the wind components in case the wind speed was bad.                       #
   #---------------------------------------------------------------------------------------#
   del               = is.na(eft$atm.vels) & is.finite(eft$atm.vdir)
   eft$atm.uspd[del] = NA
   eft$atm.vspd[del] = NA
   #=======================================================================================#
   #=======================================================================================#






   #=======================================================================================#
   #=======================================================================================#
   #      Delete net radiation in case any of the radiation components have been           #
   # discarded.                                                                            #
   #---------------------------------------------------------------------------------------#
   del           = ( ( is.na(eft$rshort.in ) | is.na(eft$rlong.in )
                     | is.na(eft$rshort.out) | is.na(eft$rlong.out) )
                   & is.finite(eft$rnet) )
   eft$rnet[del] = NA
   #=======================================================================================#
   #=======================================================================================#





   #=======================================================================================#
   #=======================================================================================#
   #=======================================================================================#
   #=======================================================================================#
   #    Some sensors drift away from their calibration, so we check whether to fix the     #
   # drift or not.                                                                         #
   #---------------------------------------------------------------------------------------#
   cat0("  + Correct bogus drifts:")
   for (vv in sequence(nvargap)){
      #------------------------------------------------------------------------------------#
      #     Copy the settings for this variable to scratch variables.                      #
      #------------------------------------------------------------------------------------#
      this.vnam        = vargap$vari      [vv]
      this.desc        = vargap$desc      [vv]
      this.unit        = vargap$unit      [vv]
      this.harm        = vargap$harm      [vv]
      this.harm.tback  = vargap$harm.tback[vv]
      this.oban        = vargap$oban      [vv]
      this.out.hour    = vargap$out.hour  [vv]
      this.out.all     = vargap$out.all   [vv]
      this.drift       = vargap$drift     [vv]
      this.quadratic   = vargap$quadratic [vv]
      this.sw.typevar  = vargap$sw.typevar[vv]
      this.potential   = vargap$potential [vv]
      this.whena       = chron(vargap$whena.drift[vv])
      this.whenz       = chron(vargap$whenz.drift[vv])
      this.form        = vargap$form      [vv]
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #     Check whether to correct the drift ot not.                                     #
      #------------------------------------------------------------------------------------#
      if (this.drift){
         cat0("   - Variable: ",this.vnam,".")

         #----- Initialise a list that will keep information about the drift correction. --#
         fixlist[[this.vnam]] = list()
         #---------------------------------------------------------------------------------#



         #----- Number of independent sub-periods to correct drift. -----------------------#
         nsubp       = length(this.whena)
         this.period = paste0("(",numyears(this.whena),"-",numyears(this.whenz),")")
         #---------------------------------------------------------------------------------#


         #------ Find out which days are complete. ----------------------------------------#
         num.when           = as.numeric(eft$when)
         drift.correction   = rep(1.0,times=ntimes)
         int.tomonth        = as.numeric(eft$tomonth)
         unique.int.tomonth = unique(int.tomonth)
         unique.tomonth     = chron(unique.int.tomonth)
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #     Check whether this is a radiation type of variable.  If it is, then we      #
         # use only the top 5% relative maximum, otherwise we use all data.                #
         #---------------------------------------------------------------------------------#
         if (this.sw.typevar){
            #----- Find the period and which observations to consider. --------------------#
            sel.period    = eft$when %wr% chron(c(this.whena,this.whenz)) 
            sel.fit       = sel.period & eft$cosz > cosz.thresh
            #------------------------------------------------------------------------------#


            #----- Remove data outside the range. -----------------------------------------#
            out           = ! sel.fit
            this.var      = eft[[this.vnam     ]]
            this.pot      = eft[[this.potential]]
            this.var[out] = NA
            this.pot[out] = NA 
            #------------------------------------------------------------------------------#


            #----- Find the fraction of the total. ----------------------------------------#
            this.frac     = this.var / this.pot
            #------------------------------------------------------------------------------#


            #------------------------------------------------------------------------------#
            #     Find the 95% quantile of the variable for each month, but keep only      #
            # the months that have enough information.                                     #
            #------------------------------------------------------------------------------#
            this.frac.month      = tapply( X     = this.frac
                                         , INDEX = int.tomonth
                                         , FUN   = quantile
                                         , prob  = 0.95
                                         , na.rm = TRUE
                                         )#tapply
            ok.month             = tapply( X     = is.finite(this.var)
                                         , INDEX = int.tomonth
                                         , FUN   = sum
                                         , na.rm = TRUE
                                         )#tapply
            tot.month            = tapply( X     = is.finite(this.pot)
                                         , INDEX = int.tomonth
                                         , FUN   = sum
                                         , na.rm = TRUE
                                         )#tapply
            frac.ok              = ok.month / tot.month
            del                  = ! ( frac.ok %>=% 0.50 )
            this.frac.month[del] = NA
            #------------------------------------------------------------------------------#



            #------------------------------------------------------------------------------#
            #     Find the threshold for each month, and remove the lower 95% of the       #
            # data plus data for months with few measurements.                             #
            #------------------------------------------------------------------------------#
            idx.month    = match(int.tomonth,unique.int.tomonth)
            frac.thresh  = this.frac.month[idx.month]
            #------------------------------------------------------------------------------#



            #------------------------------------------------------------------------------#
            #     Remove the bottom 95% of the data.                                       #
            #------------------------------------------------------------------------------#
            this.ts         = this.frac
            del             = ! ( this.frac %>=% frac.thresh)
            this.ts[del]    = NA
            #------------------------------------------------------------------------------#
         }else{
            #------------------------------------------------------------------------------#
            #     Use the entire time series, just scaling by the 95% quantile so the      #
            # method is normalised (uneccesary, just to compare different variables).      #
            #------------------------------------------------------------------------------#
            sel.period    = eft$when %wr% chron(c(this.whena,this.whenz))
            out           = ! sel.period
            this.use      = eft[[this.vnam]]
            this.use[out] = NA
            this.pot      = quantile(this.use,prob=0.95,na.rm=TRUE)
            this.ts       = this.use / this.pot
            #------------------------------------------------------------------------------#
         }#end if
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #     Make the data frame.                                                        #
         #---------------------------------------------------------------------------------#
         data.in = data.frame(when=as.numeric(eft$when),this.ts=this.ts)
         #----- Find the total number of valid entries. -----------------------------------#
         use     = is.finite(this.ts)
         nuse    = length(use)
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #     Fit a curve.  First we try a parabolic curve, then we check the p-values.   #
         # If the quadratic term is significant, keep it, otherwise, simplify to a         #
         # straight line.  If even the straight line is non-significant, then we assume    #
         # that no drift correction is needed.  Regardless the curve, we use a robust      #
         # linear fit, because it is less sensitive to outliers.                           #
         #---------------------------------------------------------------------------------#
         if (this.quadratic){
            fit.in    = rlm( this.ts ~ 1 + when + I(when^2),data=data.in)
            summ.in   = summary(object=fit.in)
            pvalue.in = 2.0 * pt(-abs(summ.in$coefficients[,3]),df=summ.in$df[2])
            polyfit   = 1 + as.integer(pvalue.in[3] <= pval.max)
         }else{
            polyfit   = 1
         }#end if
         #----- The variable shan't be fitted to a quadratic, try linear. -----------------#
         if (polyfit == 1){
            fit.in       = rlm( this.ts ~ 1 + when,data=data.in)
            summ.in      = summary(object=fit.in)
            pvalue.in    = 2*pt(-abs(summ.in$coefficients[,3]),df=summ.in$df[2])
            polyfit      = 0 + as.integer(pvalue.in[2] <= pval.max)
            #----- Linear fit wasn't significant, use a flat line (no drift). -------------#
            if (polyfit == 0){
               fit.in    = rlm( this.ts ~ 1,data=data.in)
               summ.in   = summary(object=fit.in)
               pvalue.in = 2*pt(-abs(summ.in$coefficients[,3]),df=summ.in$df[2])
            }#end if 
         }#end if
         cat0("       > Adjusted polynomial of order ",polyfit,".")
         #---------------------------------------------------------------------------------#



         #----- Find the adjusted R2. -----------------------------------------------------#
         ss.err                    = sum(summ.in$residuals^2)
         df.err                    = summ.in$df[2]
         mean.y                    = mean(this.ts,na.rm=TRUE)
         ss.tot                    = sum((this.ts-mean.y)^2,na.rm=TRUE)
         df.tot                    = sum(is.finite(this.ts)) - 1
         r2.in                     = 1 - ss.err * df.tot / (ss.tot * df.err)
         this.pred                 = predict(object=fit.in,newdata = data.in)
         fixlist[[this.vnam]][[n]] = list( whena     = this.whena[n]
                                         , whenz     = this.whenz[n]
                                         , fit.in    = fit.in
                                         , summ.in   = summ.in
                                         , pvalue.in = pvalue.in
                                         , r2.in     = r2.in
                                         )#end list
         #---------------------------------------------------------------------------------#




         #---------------------------------------------------------------------------------#
         #     Plot the regression line                                                    #
         #---------------------------------------------------------------------------------#
         cat0("     * Plot fitted curve.")
         if (plotregression){

            #----- Make a directory for the output. ---------------------------------------#
            outdrift = file.path(outroot,"drift")
            if (! file.exists(outdrift)) dir.create(outdrift)
            #------------------------------------------------------------------------------#


            #----- Select only the time we will plot. -------------------------------------#
            when.now     = eft$when[use]
            hour.now     = eft$hour[use]
            this.ts.now  = this.ts [use]
            #------------------------------------------------------------------------------#



            #----- Set limits. ------------------------------------------------------------#
            xlimit = pretty.xylim(u=when.now,fracexp=0.2)
            ylimit = pretty.xylim(u=this.ts.now,fracexp=0.2)
            #------------------------------------------------------------------------------#



            #----- Find a nice scale for time. --------------------------------------------#
            whenplot     = pretty.time(chron(xlimit),n=8)
            #------------------------------------------------------------------------------#


            #------------------------------------------------------------------------------#
            #     Set up title,  and axes labels, colour of the data points, and the       #
            # prefix of the file name.                                                     #
            #------------------------------------------------------------------------------#
            if (this.sw.typevar){
               zen.label = sprintf("%2.2i",round(acos(cosz.thresh)*180./pi))
               zen.title = sprintf("%.1f" ,acos(cosz.thresh)*180./pi)
               letitre   = paste0 ("Relative ",this.desc," \n"
                                  ,"Maximum zenith angle = ",zen.title)
               lex       = "Time"
               ley       = paste0 ("Relative ",this.desc)
               #----- Colour-code the scatter-plot with the hour of the maximum. ----------#
               which.hour = sort(unique(hour.now))
               leg.hour   = sprintf("%2.2i",which.hour)
               rbow.hour  = panoply(n=length(which.hour))
               col.now    = rbow.hour[match(hour.now,which.hour)]
               #---------------------------------------------------------------------------#


               #---------------------------------------------------------------------------#
               #    File prefix.                                                           #
               #---------------------------------------------------------------------------#
               this.prefix = paste0(this.vnam,"-zen_",zen.label)
               #---------------------------------------------------------------------------#
            }else{
               letitre   = paste0("Relative ",this.desc)
               lex       = "Time"
               ley       = paste0("Relative ",this.desc)
               
               #----- Use one colour only. ------------------------------------------------# 
               col.now   = "steelblue"
               #---------------------------------------------------------------------------#


               #---------------------------------------------------------------------------#
               #    File prefix.                                                           #
               #---------------------------------------------------------------------------#
               this.prefix = this.vnam
               #---------------------------------------------------------------------------#
            }#end if
            #------------------------------------------------------------------------------#




            #------------------------------------------------------------------------------#
            #     Retrieve R2 and p-values.                                                #
            #------------------------------------------------------------------------------#
            if (polyfit == 2){
               aa     = sprintf( "%.3f",fit.in$coefficients[1])
               bb     = sprintf("%9.2e",fit.in$coefficients[2])
               cc     = sprintf("%9.2e",fit.in$coefficients[3])
               paa    = sprintf("%9.2e",pvalue.in[1]          )
               pbb    = sprintf("%9.2e",pvalue.in[2]          )
               pcc    = sprintf("%9.2e",pvalue.in[3]          )
               rrfit  = sprintf(" %.3f",r2.in                 )
               xyfit  = substitute( expr = y == aa + bb * x + cc * x^2
                                  , env  = list(aa = aa, bb = bb, cc = cc))
               r2text = substitute( expr = R ^ 2 == r2   , env = list(r2=rrfit))
               px     = substitute( p(1,x,x^2) == p(paa,pbb,pcc)
                                  , env  = list( paa=paa,pbb=pbb,pcc=pcc))
            }else if(polyfit == 1){
               aa     = sprintf( "%.3f",fit.in$coefficients[1])
               bb     = sprintf("%9.2e",fit.in$coefficients[2])
               paa    = sprintf("%9.2e",pvalue.in[1]          )
               pbb    = sprintf("%9.2e",pvalue.in[2]          )
               rrfit  = sprintf(" %.3f",r2.in                 )
               xyfit  = substitute( expr = y == aa + bb * x
                                  , env  = list(aa = aa, bb = bb))
               r2text = substitute( expr = R ^ 2 == r2   , env = list(r2=rrfit))
               px     = substitute( p(1,x) == p(paa,pbb)
                                  , env  = list( paa=paa,pbb=pbb))
            }else{
               aa     = sprintf( "%.3f",fit.in$coefficients[1])
               paa    = sprintf("%9.2e",pvalue.in[1]          )
               rrfit  = sprintf(" %.3f",r2.in                 )
               xyfit  = substitute( expr = y == aa
                                  , env  = list(aa = aa))
               r2text = substitute( expr = R ^ 2 == r2   , env = list(r2=rrfit))
               px     = substitute( p(1) == p(paa)
                                  , env  = list( paa=paa))
            }#end if
            #------------------------------------------------------------------------------#



            #------------------------------------------------------------------------------#
            #      Loop through formats.                                                   #
            #------------------------------------------------------------------------------#
            for (o in sequence(nout)){
               #----- File name and format. -----------------------------------------------#
               fichier = file.path(outdrift,paste0(this.prefix,".",outform[o]))
               dummy   = open.file( fichier = fichier
                                  , outform = outform[o]
                                  , size    = sq.size
                                  , ptsz    = ptsz
                                  , depth   = depth
                                  )#end open.file
               #---------------------------------------------------------------------------#


               #----- Open plot window. ---------------------------------------------------#
               par(par.user)
               plot.new()
               plot.window(xlim=xlimit,ylim=ylimit)
               #----- Plot time axis. -----------------------------------------------------#
               axis(side=1,at=whenplot$levels,labels=whenplot$labels,padj=whenplot$padj)
               axis(side=2,las=1)
               #---------------------------------------------------------------------------#
               #     Make a grid that matches the time, and add beginning and end of gap.  #
               #---------------------------------------------------------------------------#
               if (plotgrid){ 
                  abline(v=whenplot$levels,h=axTicks(side=2),col="grey76",lty="dotted")
               }#end if
               #----- Plot the model prediction. ------------------------------------------#
               lines(x=eft$when,y=this.pred,col="grey28",lty="solid",lwd=2)
               #----- Plot the points, colours according to times. ------------------------#
               points(x=when.now,y=this.ts.now,pch=16,cex=1.1,col=col.now)
               #----- Plot the legend. ----------------------------------------------------#
               if (this.sw.typevar){
                  legend(x="topright",inset=0.01,legend=leg.hour,ncol=2
                        ,title="Time [UTC]",pch=16,pt.cex=1.2,col=rbow.hour,bg="white")
               }#end if
               title(main=letitre,xlab=lex,ylab=ley)
               box()
               #---------------------------------------------------------------------------#


               #---------------------------------------------------------------------------#
               #      Find some good position to plot the fit.                             #
               #---------------------------------------------------------------------------#
               xt = xlimit[1] + 0.01*diff(xlimit); yt = ylimit[1] + 0.12*diff(ylimit)
               text(xt,yt,xyfit,adj=c(0,0.5))
               xt = xlimit[1] + 0.01*diff(xlimit); yt = ylimit[1] + 0.06*diff(ylimit)
               text(xt,yt,r2text,adj=c(0,0.5))
               xt = xlimit[1] + 0.01*diff(xlimit); yt = ylimit[1] + 0.00*diff(ylimit)
               text(xt,yt,px,adj=c(0,0.5))
               #---------------------------------------------------------------------------#


               #---------------------------------------------------------------------------#
               #      Close the device.                                                    #
               #---------------------------------------------------------------------------#
               dummy = close.plot(outform=outform[o])
               #---------------------------------------------------------------------------#
            }#end for
            #------------------------------------------------------------------------------#
         }#end if
         #---------------------------------------------------------------------------------#


         #---------------------------------------------------------------------------------#
         #     Find ratio between the actual variable and the "no drift" variable.         #
         #---------------------------------------------------------------------------------#
         cat0("     * Find the drift correction term.")
         new.data                     = data.frame(when=as.numeric(eft$when))
         ref.0                        = match(min(eft$when[sel.period]),eft$when)
         data.0                       = data.frame(when=as.numeric(eft$when[ref.0]))
         month.mult                   = ( predict(object=fit.in,newdata=data.0  )
                                        / predict(object=fit.in,newdata=new.data) )
         drift.correction[sel.period] = month.mult[sel.period]
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #     Plot the drift correction time series                                       #
         #---------------------------------------------------------------------------------#
         cat0("     * Plot the drift correction term.")
         if (plotregression){

            #----- Make a directory for the output. ---------------------------------------#
            outdrift = file.path(outroot,"drift")
            if (! file.exists(outdrift)) dir.create(outdrift)
            #------------------------------------------------------------------------------#




            #----- Set limits. ------------------------------------------------------------#
            xlimit = pretty.xylim(u=eft$when        ,fracexp=0.0)
            ylimit = pretty.xylim(u=drift.correction,fracexp=0.2)
            #------------------------------------------------------------------------------#



            #----- Find a nice scale for time. --------------------------------------------#
            whenplot     = pretty.time(chron(xlimit),n=8)
            #------------------------------------------------------------------------------#


            #------------------------------------------------------------------------------#
            #      Set up title and axes labels, and the file prefix according to the      #
            # type of variable.                                                            #
            #------------------------------------------------------------------------------#
            if (this.sw.typevar){
               zen.label = sprintf("%2.2i",round(acos(cosz.thresh)*180./pi))
               zen.title = sprintf("%.1f" ,acos(cosz.thresh)*180./pi)
               letitre  = paste0( "Drift correction for ",this.desc," \n"
                                , "Maximum zenith angle = ",zen.title)
               lex      = "Time"
               ley      = "Drift multiplication term [---]"
               prefix   = paste0(this.vnam,"-drift_correction-zen_",zen.label)
            }else{
               letitre  = paste0( "Drift correction for ",this.desc)
               lex      = "Time"
               ley      = "Drift multiplication term [---]"
               prefix   = paste0(this.vnam,"-drift_correction")
            }#end if
            #------------------------------------------------------------------------------#

            #------------------------------------------------------------------------------#
            #     Loop through formats.                                                    #
            #------------------------------------------------------------------------------#
            for (o in sequence(nout)){
               #----- File name and format. -----------------------------------------------#
               fichier = file.path(outdrift,paste0(prefix,".",outform[o]))
               dummy   = open.file( fichier = fichier
                                  , outform = outform[o]
                                  , size    = sq.size
                                  , ptsz    = ptsz
                                  , depth   = depth
                                  )#end open.file
               #---------------------------------------------------------------------------#


               #----- Open plot window. ---------------------------------------------------#
               par(par.user)
               plot.new()
               plot.window(xlim=xlimit,ylim=ylimit)
               title(main=letitre,xlab=lex,ylab=ley)
               #----- Plot time axis. -----------------------------------------------------#
               axis(side=1,at=whenplot$levels,labels=whenplot$labels,padj=whenplot$padj)
               axis(side=2,las=1)
               #---------------------------------------------------------------------------#
               #     Make a grid that matches the time, and add beginning and end of gap.  #
               #---------------------------------------------------------------------------#
               if (plotgrid){ 
                  abline(v=whenplot$levels,h=axTicks(side=2),col="grey76",lty="dotted")
               }#end if
               lines(x=eft$when,y=drift.correction,col="forestgreen",lty="solid",lwd=2)
               title(main=letitre,xlab=lex,ylab=ley)
               box()
               #---------------------------------------------------------------------------#


               #---------------------------------------------------------------------------#
               #      Close the device.                                                    #
               #---------------------------------------------------------------------------#
               dummy = close.plot(outform=outform[o])
               #---------------------------------------------------------------------------#
            }#end for (o in sequence(nout))
            #------------------------------------------------------------------------------#
         }#end if (plotregression)
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #     Apply the drift correction to the variable.                                 #
         #---------------------------------------------------------------------------------#
         original   = eft[[this.vnam]]
         detrended  = original * drift.correction
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #     Apply the drift correction to the variable.                                 #
         #---------------------------------------------------------------------------------#
         fixlist[[this.vnam]]$drift.correction = drift.correction
         eft[[this.vnam]]                      = eft[[this.vnam]] * drift.correction
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #     Now that the time series was drift-corrected, we check whether the periods  #
         # must be inter-calibrated.  Again, we decide which statistic to use depending    #
         # on whether it is a radiation type of variable or not.                           #
         #---------------------------------------------------------------------------------#
         cat0("     * Find the calibration correction term:")
         scale.correction     = rep(1.0,times=ntimes)
         #----- Reference period. ---------------------------------------------------------#
         cat("       > Reference period.")
         sel.reference        = eft$when %wr% chron(c(this.whena,this.whenz))
         if (this.sw.typevar){
            #----- Radiation variable.  Pick up the top 5% of the variables. --------------#
            out                   = (! sel.reference ) | eft$cosz <= cosz.thresh
            this.reference        = eft[[this.vnam]] / eft[[this.potential]]
            this.reference[out]   = NA
            top.05.reference      = quantile(this.reference,prob=0.95,na.rm=TRUE)
            top.05                = this.reference >= top.05.reference
            top.05[is.na(top.05)] = FALSE
            mean.reference        = mean(this.reference[top.05],na.rm=TRUE)
            #------------------------------------------------------------------------------#
         }else{
            #----- Other variables.  Take the mean as the reference. ----------------------#
            out                   = (! sel.reference )
            this.reference        = eft[[this.vnam]]
            this.reference[out]   = NA
            mean.reference        = mean(this.reference,na.rm=TRUE)
         }#end if
         scale.correction = ifelse(sel.reference,mean.reference,1.0)
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #     Make a calibration factor.                                                  #
         #---------------------------------------------------------------------------------#
         fixlist[[this.vnam]]$scale.correction = scale.correction
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #     Correct the time series.                                                    #
         #---------------------------------------------------------------------------------#
         eft[[this.vnam]] = eft[[this.vnam]] * scale.correction
         detrended        = eft[[this.vnam]]
         #---------------------------------------------------------------------------------#


         #---------------------------------------------------------------------------------#
         #     Convert units for PAR.                                                      #
         #---------------------------------------------------------------------------------#
         if (this.vnam %in% c("par.in","par.out")){
            dc.orig = original  * Watts.2.Ein * 1E6
            dc.detr = detrended * Watts.2.Ein * 1E6
            dc.unit = "umolom2os"
         }else{
            dc.orig = original
            dc.detr = detrended
            dc.unit = this.unit
         }#end if
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #     Find the range of both variables.                                           #
         #---------------------------------------------------------------------------------#
         span.original   = diff(range(dc.orig,na.rm=TRUE))
         span.detrended  = diff(range(dc.detr,na.rm=TRUE))
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #     Find the plotting limits.                                                   #
         #---------------------------------------------------------------------------------#
         xlimit  = pretty.xylim(u=eft$when          ,fracexp=0.0)
         ylimit  = pretty.xylim(u=c(dc.orig,dc.detr),fracexp=0.0)
         yat     = pretty(ylimit)
         ylabels = sprintf("%g",yat)
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #     Make a nice time axis.                                                      #
         #---------------------------------------------------------------------------------#
         whenplot = pretty.time(eft$when,n=8)
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #     Set title and axes labels.                                                  #
         #---------------------------------------------------------------------------------#
         letitre = paste("Time series of ",this.desc," - ",longname,sep="")
         lexlab  = desc.unit(desc="Time"   ,unit=untab$utc         )
         leylab  = desc.unit(desc=this.desc,unit=untab[[this.unit]])
         #---------------------------------------------------------------------------------#


         #----- Set colours. --------------------------------------------------------------#
         if (span.original >= span.detrended){
            detr.col = "midnightblue"
            orig.col = "deepskyblue"
         }else{
            orig.col = "midnightblue"
            detr.col = "deepskyblue"
         }#end if
         #---------------------------------------------------------------------------------#


         #---------------------------------------------------------------------------------#
         #     Plot the trend versus detrended.  Decide the order for plotting based on    #
         # the one with the widest range (so we can see both plots.                        #
         #---------------------------------------------------------------------------------#
         for (o in sequence(nout)){
            #----- File name and format. --------------------------------------------------#
            fichier = file.path(outdrift,paste0(this.vnam,"_detrend.",outform[o]))
            dummy   = open.file( fichier = fichier
                               , outform = outform[o]
                               , size    = ex.size
                               , ptsz    = ptsz
                               , depth   = depth
                               )#end open.file
            #------------------------------------------------------------------------------#



            #------------------------------------------------------------------------------#
            #     Split device in two.                                                     #
            #------------------------------------------------------------------------------#
            par(par.user)
            layout(mat=rbind(2,1),heights=c(1.-f.leg,f.leg))
            #------------------------------------------------------------------------------#



            #------------------------------------------------------------------------------#
            #     Plot legend.                                                             #
            #------------------------------------------------------------------------------#
            par(mar=c(0.1,4.6,0.1,1.1))
            plot.new()
            plot.window(xlim=c(0,1),ylim=c(0,1))
            legend( x       = "bottom"
                  , inset   = 0
                  , legend  = c("Raw","Drift corrected")
                  , col     = c(orig.col,detr.col)
                  , lwd     = 2
                  , ncol    = 2
                  )#end legend
            #------------------------------------------------------------------------------#





            #------------------------------------------------------------------------------#
            #    Plot time series.                                                         #
            #------------------------------------------------------------------------------#
            par(mar=c(4.1,4.6,4.1,1.1))
            plot.new()
            plot.window(xlim=xlimit,ylim=ylimit)
            axis(side=1,las=1,at=whenplot$levels,labels=whenplot$labels,padj=whenplot$padj)
            axis(side=2,las=1,at=yat,labels=ylabels)
            title(main=letitre,xlab=lexlab,ylab=leylab)
            #------------------------------------------------------------------------------#
            #     Make a grid that matches the time, and add beginning and end of gap.     #
            #------------------------------------------------------------------------------#
            if (plotgrid){ 
               abline(v=whenplot$levels,h=yat,col=grid.colour,lty="dotted")
            }#end if
            #------------------------------------------------------------------------------#



            #------------------------------------------------------------------------------#
            #      Plot time series, with the one with larger range in the background.     #
            #------------------------------------------------------------------------------#
            if (span.original >= span.detrended){
               lines(x=eft$when,y=dc.orig,type="l",col=orig.col)
               lines(x=eft$when,y=dc.detr,type="l",col=detr.col)
            }else{
               lines(x=eft$when,y=dc.detr,type="l",col=detr.col)
               lines(x=eft$when,y=dc.orig,type="l",col=orig.col)
            }#end if
            box()
            #------------------------------------------------------------------------------#


            #------------------------------------------------------------------------------#
            #      Close the device.                                                       #
            #------------------------------------------------------------------------------#
            dummy = close.plot(outform=outform[o])
            #------------------------------------------------------------------------------#
         }#end for (o in sequence(nout))
         #---------------------------------------------------------------------------------#
      }#end if (this.drift)
      #------------------------------------------------------------------------------------#
   }#end for (vv in sequence(nvargap))
   #---------------------------------------------------------------------------------------#
   #=======================================================================================#
   #=======================================================================================#



   #=======================================================================================#
   #=======================================================================================#
   #    Remove data points and/or periods that are very suspicious.                        #
   #---------------------------------------------------------------------------------------#
   cat0(" + Remove suspicious incoming shortwave radiation.")
   eft = del.bad.rshort(dat=eft)
   #=======================================================================================#
   #=======================================================================================#



   #---------------------------------------------------------------------------------------#
   #     Save the compounded dataset so we don't load them again.                          #
   #---------------------------------------------------------------------------------------#
   cat0(" + Save cleaned data to file ",basename(session.clean),".")
   save(eft,fixlist,file=session.clean)
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #    Copy eft to clean.  From this point on, the cleaned data will not be called eft,   #
   # but clean because we always use eft as the structure name that is saved to the RData  #
   # regardless of the level.                                                              #
   #---------------------------------------------------------------------------------------#
   clean = eft
   #---------------------------------------------------------------------------------------#

}#end if
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
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#     Check whether to gap fill using redundant variables.                                 #
#------------------------------------------------------------------------------------------#
if (reload.selfgf && file.exists(session.selfgf)){
   #----- Reload the gap-filling. ---------------------------------------------------------#
   cat0(" + Load self gap-filled data from file ",basename(session.selfgf),".")
   load(session.selfgf)
   selfgf  = eft
   ntimes  = length(selfgf$when)
}else{


   #----- Initialise the data frame that will contain the gap-filled dataset. -------------#
   eft     = clean
   ntimes  = length(eft$when)
   #---------------------------------------------------------------------------------------#




   #=======================================================================================#
   #=======================================================================================#
   #     Gap-filling with redundant variables.  This is the preferred method because it    #
   # uses other data available from the same site.                                         #
   #---------------------------------------------------------------------------------------#
   if (use.similar){

      cat0( " + Preliminary gap-filling using redundant variables from the site.")
      #------------------------------------------------------------------------------------#
      #     Now that PAR was drift-corrected, we find the expected rshort given the        #
      # incoming PAR, and fill in the data with this information.  We also use the ratio   #
      # between the observed shortwave radiation and the obtained by the model to correct  #
      # PAR.                                                                               #
      #------------------------------------------------------------------------------------#
      cat("   - Use PAR to estimate incoming shortwave...","\n")

         #----- Make data frame with variables we judge relevant. -------------------------#
         rsbdown = rshort.bdown( rad.in   = eft$par.in
                               , atm.prss = ref.prss + 0. * eft$cosz
                               , cosz     = eft$cosz
                               , rad.type = "par")

         #----- Find the predicted rshort, and the save the model. ------------------------#
         rshort.in.wn85       = rsbdown$rshort.full
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #     Make the data frames.                                                       #
         #---------------------------------------------------------------------------------#
         both      = is.finite(rshort.in.wn85) & is.finite(eft$rshort.in) & eft$daytime
         nboth     = sum(both,na.rm=TRUE)
         colfit    = cloudy(n=nboth)
         data.in   = data.frame( rshort.in.wn85 = rshort.in.wn85 [both]
                               , cosz           = eft$cosz       [both]
                               , rshort.in.eft  = eft$rshort.in  [both]
                               ) # end data.frame
         idx = order(data.in$rshort.in.wn85,na.last=TRUE,decreasing=FALSE)
         data.pred = data.frame( rshort.in.wn85 = data.in$rshort.in.wn85 [idx]
                               , cosz           = data.in$cosz           [idx]
                               ) # end data.frame
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #     Fit a line using time crossing zero.                                        #
         #---------------------------------------------------------------------------------#
         fit.in    = rlm(rshort.in.eft ~ rshort.in.wn85 * cosz,data=data.in)
         summ.in   = summary(object=fit.in)
         pvalue.in = 2. * pt(-abs(summ.in$coefficients[,3]),df=summ.in$df[2])
         #---------------------------------------------------------------------------------#



         #----- Find the adjusted R2. -----------------------------------------------------#
         ss.err                  = sum(summ.in$residuals^2)
         df.err                  = summ.in$df[2]
         mean.y                  = mean(data.in$rshort.in.eft,na.rm=TRUE)
         ss.tot                  = sum((data.in$rshort.in.eft-mean.y)^2,na.rm=TRUE)
         df.tot                  = sum(is.finite(data.in$rshort.in.eft)) - 1
         r2.in                   = 1.0 - ss.err * df.tot / (ss.tot * df.err)
         data.pred$rshort.in.eft = predict(object=fit.in,newdata = data.pred)
         fixlist$rshort.in       = list()
         fixlist$rshort.in$model = list( fit       = fit.in
                                       , summ      = summ.in
                                       , pvalue.in = pvalue.in
                                       , r2.in     = r2.in
                                       )#end list
         #---------------------------------------------------------------------------------#




         #---------------------------------------------------------------------------------#
         #     Retrieve R2 and p-value of slope.                                           #
         #---------------------------------------------------------------------------------#
         aa     = fit.in$coefficients[1]
         bb     = fit.in$coefficients[2]
         cc     = fit.in$coefficients[3]
         dd     = fit.in$coefficients[4]
         pa     = pvalue.in          [1]
         pb     = pvalue.in          [2]
         pc     = pvalue.in          [3]
         pd     = pvalue.in          [4]
         xyfit  = substitute( expr = y == aa + bb * x + cc * z + dd * x * z
                            , env  = list( aa = sprintf("%9.2e",aa)
                                         , bb = sprintf("%9.2e",bb)
                                         , cc = sprintf("%9.2e",cc)
                                         , dd = sprintf("%9.2e",dd) ) )
         r2text = substitute( expr = R ^ 2 == r2
                            , env  = list(r2=sprintf("%.3f",r2.in)) )
         ppval  = substitute( p(0,x,z,xz) == p(pa,pb,pc,pd)
                            , env  = list( pa = sprintf("%9.2e",pa)
                                         , pb = sprintf("%9.2e",pb)
                                         , pc = sprintf("%9.2e",pc)
                                         , pd = sprintf("%9.2e",pd) ) )
         #---------------------------------------------------------------------------------#


         #---------------------------------------------------------------------------------#
         #       Plot the model fitting.                                                   #
         #---------------------------------------------------------------------------------#
         if (plotregression){


            zen.label = sprintf("%2.2i",round(acos(cosz.thresh)*180./pi))
            zen.title = sprintf("%.1f" ,acos(cosz.thresh)*180./pi)


            #----- Set limits. ------------------------------------------------------------#
            xlimit = pretty.xylim(data.in$rshort.in.wn85,fracexp=0.0)
            ylimit = pretty.xylim(data.in$rshort.in.eft ,fracexp=0.2)
            #------------------------------------------------------------------------------#




            #------------------------------------------------------------------------------#
            #     Set title and axes labels.                                               #
            #------------------------------------------------------------------------------#
            letitre = paste("Incoming shortwave radiation [W/m2]"
                           ,"Time goes from blue to grey"
                           ,paste0("Maximum zenith angle = ",zen.title)
                           ,sep="\n")
            lexlab  = "Drift-corrected, PAR based estimate using Weiss-Norman (1985)"
            leylab  = "Observed data"
            #------------------------------------------------------------------------------#



            #------------------------------------------------------------------------------#
            #    Loop through all formats.                                                 #
            #------------------------------------------------------------------------------#
            for (o in sequence(nout)){
               #----- File name and format. -----------------------------------------------#
               fichier = file.path( outroot
                                  , paste0("rshort-lm-zen_",zen.label,".",outform[o])
                                  )#end file.path
               dummy   = open.file( fichier = fichier
                                  , outform = outform[o]
                                  , size    = sq.size
                                  , ptsz    = ptsz
                                  , depth   = depth
                                  )#end open.file
               #---------------------------------------------------------------------------#


               #----- Plot window and a grid. ---------------------------------------------#
               par(par.user)
               plot.new()
               plot.window(xlim=xlimit,ylim=ylimit)
               axis(side=1,las=1)
               axis(side=2,las=1)
               title(main=letitre,xlab=lexlab,ylab=leylab)
               if (plotgrid) grid(col="grey63",lty="dotted",lwd=0.25)
               #---------------------------------------------------------------------------#



               #---------------------------------------------------------------------------#
               #      Plot points, the 1:1 line, and the linear fit.                       #
               #---------------------------------------------------------------------------#
               points(x=data.in$rshort.in.wn85,y=data.in$rshort.in.eft,type="p",col=colfit
                     ,pch=16)
               abline(a=0,b=1,col="goldenrod",lwd=2,lty="dashed")
               lines (x=data.pred$rshort.in.wn85,y=data.pred$rshort.in.eft
                     ,col="red3",lty="solid",lwd=2)
               #----- Legend for the lines. -----------------------------------------------#
               legend("topleft",inset=0.01,legend=c("1:1 line","Predicted values")
                     ,bg="white",col=c("goldenrod","red3"),lwd=2,lty=c("dashed","solid"))
               box()
               #---------------------------------------------------------------------------#


               #---------------------------------------------------------------------------#
               #      Find some good position to plot the fit.                             #
               #---------------------------------------------------------------------------#
               xt = xlimit[1] + 1.02 * diff(xlimit); yt = ylimit[1] + 0.10 * diff(ylimit)
               text(xt,yt,xyfit,adj=c(1,0.5),cex=0.7)
               xt = xlimit[1] + 1.02 * diff(xlimit); yt = ylimit[1] + 0.05 * diff(ylimit)
               text(xt,yt,r2text,adj=c(1,0.5),cex=0.7)
               xt = xlimit[1] + 1.02 * diff(xlimit); yt = ylimit[1] + 0.00 * diff(ylimit)
               text(xt,yt,ppval,adj=c(1,0.5),cex=0.7)
               #---------------------------------------------------------------------------#


               #---------------------------------------------------------------------------#
               #      Close the device.                                                    #
               #---------------------------------------------------------------------------#
               dummy = close.plot(outform=outform[o])
               #---------------------------------------------------------------------------#
            }#end for (o in sequence(nout))
            #------------------------------------------------------------------------------#
         }#end if (plotregression)
         #---------------------------------------------------------------------------------#




         #------ Predict shortwave radiation from PAR. ------------------------------------#
         data.pred   = data.frame( rshort.in.wn85 = rshort.in.wn85
                                 , cosz           = eft$cosz
                                 ) # end data.frame
         data.pred$rshort.in.eft    = predict(object = fit.in, newdata=data.pred)
         miss                       = ( eft$daytime
                                      & data.pred$rshort.in.eft >= 0.
                                      & data.pred$rshort.in.eft <= eft$rshort.pot
                                      & is.na(eft$rshort.in)
                                      & is.finite(eft$par.in))
         miss[is.na(miss)]          = FALSE
         
         eft$rshort.in      [miss]    = data.pred$rshort.in.eft[miss]
         eft$gfflg.rshort.in[miss]    = -1
         eft$gferr.rshort.in[miss]    = summ.in$sigma
         fixlist$rshort.in            = list()
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #      Delete outliers.                                                           #
         #---------------------------------------------------------------------------------#
         eft$rshort.in[eft$nighttime] = NA
         eft$rshort.in                = del.outliers( x        = eft$rshort.in
                                                    , when     = eft$when
                                                    , out.hour = TRUE
                                                    , out.all  = FALSE
                                                    )#end del.outliers
         eft$rshort.in[eft$nighttime] = 0.
         fixlist$rshort.in$gap.filled = miss & is.finite(eft$rshort.in)
         #---------------------------------------------------------------------------------#
      #------------------------------------------------------------------------------------#







      #------------------------------------------------------------------------------------#
      #     Filling PAR from incoming shortwave radiation.                                 #
      #------------------------------------------------------------------------------------#
      cat0("   - Using incoming shortwave to estimate PAR.")

         #----- Make data frame with variables we judge relevant. -------------------------#
         rsbdown = rshort.bdown( rad.in   = eft$rshort.in
                               , atm.prss = ref.prss + 0. * eft$cosz
                               , cosz     = eft$cosz
                               , rad.type = "rshort")

         #----- Find the predicted PAR, and the save the model. ---------------------------#
         par.in.wn85       = rsbdown$par.full
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #     Make the data frames.                                                       #
         #---------------------------------------------------------------------------------#
         both      = ( is.finite(par.in.wn85)   & is.finite(eft$par.in)      & eft$daytime
                     & is.finite(eft$rshort.in) & is.na(eft$gfflg.rshort.in) )
         nboth     = sum(both,na.rm=TRUE)
         colfit    = cloudy(n=nboth)
         data.in   = data.frame( par.in.wn85    = par.in.wn85    [both]
                               , cosz           = eft$cosz       [both]
                               , par.in.eft     = eft$par.in     [both]
                               ) # end data.frame
         idx = order(data.in$par.in.wn85,na.last=TRUE,decreasing=FALSE)
         data.pred = data.frame( par.in.wn85    = data.in$par.in.wn85    [idx]
                               , cosz           = data.in$cosz           [idx]
                               ) # end data.frame
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #     Fit a line using time crossing zero.                                        #
         #---------------------------------------------------------------------------------#
         fit.in    = rlm(par.in.eft ~ par.in.wn85 * cosz,data=data.in)
         summ.in   = summary(object=fit.in)
         pvalue.in = 2. * pt(-abs(summ.in$coefficients[,3]),df=summ.in$df[2])
         #---------------------------------------------------------------------------------#



         #----- Find the adjusted R2. -----------------------------------------------------#
         ss.err                  = sum(summ.in$residuals^2)
         df.err                  = summ.in$df[2]
         mean.y                  = mean(data.in$par.in.eft,na.rm=TRUE)
         ss.tot                  = sum((data.in$par.in.eft-mean.y)^2,na.rm=TRUE)
         df.tot                  = sum(is.finite(data.in$par.in.eft)) - 1
         r2.in                   = 1.0 - ss.err * df.tot / (ss.tot * df.err)
         data.pred$par.in.eft    = predict(object=fit.in,newdata = data.pred)
         fixlist$par.in          = list()
         fixlist$par.in$model    = list( fit       = fit.in
                                       , summ      = summ.in
                                       , pvalue.in = pvalue.in
                                       , r2.in     = r2.in
                                       )#end list
         #---------------------------------------------------------------------------------#




         #---------------------------------------------------------------------------------#
         #     Retrieve R2 and p-value of slope.                                           #
         #---------------------------------------------------------------------------------#
         aa     = fit.in$coefficients[1]
         bb     = fit.in$coefficients[2]
         cc     = fit.in$coefficients[3]
         dd     = fit.in$coefficients[4]
         pa     = pvalue.in          [1]
         pb     = pvalue.in          [2]
         pc     = pvalue.in          [3]
         pd     = pvalue.in          [4]
         xyfit  = substitute( expr = y == aa + bb * x + cc * z + dd * x * z
                            , env  = list( aa = sprintf("%9.2e",aa)
                                         , bb = sprintf("%9.2e",bb)
                                         , cc = sprintf("%9.2e",cc)
                                         , dd = sprintf("%9.2e",dd) ) )
         r2text = substitute( expr = R ^ 2 == r2
                            , env  = list(r2=sprintf("%.3f",r2.in)) )
         ppval  = substitute( p(0,x,z,xz) == p(pa,pb,pc,pd)
                            , env  = list( pa = sprintf("%9.2e",pa)
                                         , pb = sprintf("%9.2e",pb)
                                         , pc = sprintf("%9.2e",pc)
                                         , pd = sprintf("%9.2e",pd) ) )
         #---------------------------------------------------------------------------------#


         #---------------------------------------------------------------------------------#
         #       Plot the model fitting.                                                   #
         #---------------------------------------------------------------------------------#
         if (plotregression){


            zen.label = sprintf("%2.2i",round(acos(cosz.thresh)*180./pi))
            zen.title = sprintf("%.1f" ,acos(cosz.thresh)*180./pi)


            #----- Set limits. ------------------------------------------------------------#
            xlimit = pretty.xylim(data.in$par.in.wn85,fracexp=0.0)
            ylimit = pretty.xylim(data.in$par.in.eft ,fracexp=0.2)
            #------------------------------------------------------------------------------#


            #------------------------------------------------------------------------------#
            #     Set title and axes labels.                                               #
            #------------------------------------------------------------------------------#
            letitre = paste("Incoming PAR [W/m2]"
                           ,"Time goes from blue to grey"
                           ,paste0("Maximum zenith angle = ",zen.title)
                           ,sep="\n")
            lexlab  = "Drift-corrected, PAR based estimate using Weiss-Norman (1985)"
            leylab  = "Observed data"
            #------------------------------------------------------------------------------#


            #------------------------------------------------------------------------------#
            #    Loop through formats.                                                     #
            #------------------------------------------------------------------------------#
            for (o in sequence(nout)){
               #----- File name and format. -----------------------------------------------#
               fichier = file.path(outroot,paste0("par-lm-zen_",zen.label,".",outform[o]))
               dummy   = open.file( fichier = fichier
                                  , outform = outform[o]
                                  , size    = sq.size
                                  , ptsz    = ptsz
                                  , depth   = depth
                                  )#end open.file
               #---------------------------------------------------------------------------#


               #----- Plot window and a grid. ---------------------------------------------#
               par(par.user)
               plot.new()
               plot.window(xlim=xlimit,ylim=ylimit)
               axis(side=1,las=1)
               axis(side=2,las=1)
               title(main=letitre,xlab=lexlab,ylab=leylab)
               if (plotgrid) grid  (col="grey63",lty="dotted",lwd=0.25)
               #---------------------------------------------------------------------------#



               #---------------------------------------------------------------------------#
               #      Plot points, the 1:1 line, and the linear fit.                       #
               #---------------------------------------------------------------------------#
               points(x=data.in$par.in.wn85,y=data.in$par.in.eft,type="p",col=colfit
                     ,pch=16)
               abline(a=0,b=1,col="goldenrod",lwd=2,lty="dashed")
               lines (x=data.pred$par.in.wn85,y=data.pred$par.in.eft
                     ,col="red3",lty="solid",lwd=2)
               #----- Legend for the lines. -----------------------------------------------#
               legend("topleft",inset=0.01,legend=c("1:1 line","Predicted values")
                     ,bg="white",col=c("goldenrod","red3"),lwd=2,lty=c("dashed","solid"))
               box()
               #---------------------------------------------------------------------------#


               #---------------------------------------------------------------------------#
               #      Find some good position to plot the fit.                             #
               #---------------------------------------------------------------------------#
               xt = xlimit[1] + 1.02 * diff(xlimit); yt = ylimit[1] + 0.10 * diff(ylimit)
               text(xt,yt,xyfit,adj=c(1,0.5),cex=0.7)
               xt = xlimit[1] + 1.02 * diff(xlimit); yt = ylimit[1] + 0.05 * diff(ylimit)
               text(xt,yt,r2text,adj=c(1,0.5),cex=0.7)
               xt = xlimit[1] + 1.02 * diff(xlimit); yt = ylimit[1] + 0.00 * diff(ylimit)
               text(xt,yt,ppval,adj=c(1,0.5),cex=0.7)
               #---------------------------------------------------------------------------#


               #---------------------------------------------------------------------------#
               #      Close the device.                                                    #
               #---------------------------------------------------------------------------#
               dummy = close.plot(outform=outform[o])
               #---------------------------------------------------------------------------#
            }#end for (o in sequence(nout))
            #------------------------------------------------------------------------------#
         }#end if (plotregression)
         #---------------------------------------------------------------------------------#



         #------ Predict shortwave radiation from PAR. ------------------------------------#
         data.pred   = data.frame( par.in.wn85    = par.in.wn85
                                 , cosz           = eft$cosz
                                 ) # end data.frame
         data.pred$rshort.in.eft    = predict(object = fit.in, newdata=data.pred)
         miss                       = ( eft$daytime
                                      & data.pred$par.in.eft >= 0.
                                      & data.pred$par.in.eft <= eft$par.pot
                                      & is.na(eft$par.in)
                                      & is.finite(eft$rshort.in))
         miss[is.na(miss)]          = FALSE
         
         eft$par.in         [miss]    = data.pred$par.in.eft[miss]
         eft$gfflg.par.in   [miss]    = -1
         eft$gferr.par.in   [miss]    = summ.in$sigma
         fixlist$par.in               = list()
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #      Delete outliers.                                                           #
         #---------------------------------------------------------------------------------#
         eft$par.in[eft$nighttime]    = NA
         eft$par.in                   = del.outliers( x        = eft$par.in
                                                    , when     = eft$when
                                                    , out.hour = TRUE
                                                    , out.all  = FALSE
                                                    )#end del.outliers
         eft$par.in[eft$nighttime]    = 0.
         fixlist$par.in$gap.filled    = miss & is.finite(eft$par.in)
         #---------------------------------------------------------------------------------#
      #------------------------------------------------------------------------------------#








      #------------------------------------------------------------------------------------#
      #     Find a line that relates Aspirated temperature to sonic temperature.           #
      #------------------------------------------------------------------------------------#
      cat0("   - Find a linear model relating ambient and sonic temperature.")



         #---------------------------------------------------------------------------------#
         #     Make the data frames.                                                       #
         #---------------------------------------------------------------------------------#
         both      = is.finite(eft$atm.tmp) & is.finite(eft$atm.tson)
         nboth     = sum(both,na.rm=TRUE)
         data.in   = data.frame( atm.tmp  = eft$atm.tmp [both]
                               , atm.tson = eft$atm.tson[both]
                               ) # end data.frame
         data.pred = data.frame( atm.tson = eft$atm.tson)
         data.plot = data.frame( atm.tson = seq( from = min(eft$atm.tson[both])
                                               , to   = max(eft$atm.tson[both])
                                               , length.out = 1000 )
                               ) # end data.frame
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #     Fit a line using time crossing zero.                                        #
         #---------------------------------------------------------------------------------#
         fit.in    = rlm(atm.tmp ~ atm.tson,data=data.in)
         summ.in   = summary(object=fit.in)
         pvalue.in = 2. * pt(-abs(summ.in$coefficients[,3]),df=summ.in$df[2])
         #---------------------------------------------------------------------------------#



         #----- Find the adjusted R2. -----------------------------------------------------#
         ss.err                  = sum(summ.in$residuals^2)
         df.err                  = summ.in$df[2]
         mean.y                  = mean(data.in$atm.tmp,na.rm=TRUE)
         ss.tot                  = sum((data.in$atm.tmp-mean.y)^2,na.rm=TRUE)
         df.tot                  = sum(is.finite(data.in$atm.tmp)) - 1
         r2.in                   = 1.0 - ss.err * df.tot / (ss.tot * df.err)
         data.plot$atm.tmp       = predict(object=fit.in,newdata = data.plot)
         data.pred$atm.tmp       = predict(object=fit.in,newdata = data.pred)
         fixlist$atm.tmp         = list()
         fixlist$atm.tmp$model   = list( fit       = fit.in
                                       , summ      = summ.in
                                       , pvalue.in = pvalue.in
                                       , r2.in     = r2.in
                                       )#end list
         #---------------------------------------------------------------------------------#




         #---------------------------------------------------------------------------------#
         #     Retrieve R2 and p-value of slope.                                           #
         #---------------------------------------------------------------------------------#
         aa     = sprintf("%.3f" ,fit.in$coefficients[1])
         bb     = sprintf("%.3f" ,fit.in$coefficients[2])
         r2     = sprintf("%.3f" ,summ.in$adj.r.squared )
         pa     = sprintf("%9.2e",pvalue.in[1])
         pb     = sprintf("%9.2e",pvalue.in[2])
         xyfit  = substitute( expr = y == aa + bb * x   , env  = list( aa = aa, bb = bb) )
         r2text = substitute( expr = R ^ 2 == r2        , env  = list( r2 = r2.in      ) )
         ppval  = substitute( expr = p(1,x) == p(pa,pb) , env  = list( pa = pa, pb = pb) )
         #---------------------------------------------------------------------------------#


         #---------------------------------------------------------------------------------#
         #       Plot the model fitting.                                                   #
         #---------------------------------------------------------------------------------#
         if (plotregression){

            #----- Set limits. ------------------------------------------------------------#
            xlimit = pretty.xylim(data.in$atm.tson-t00,fracexp=0.0)
            ylimit = pretty.xylim(data.in$atm.tmp-t00 ,fracexp=0.2)
            #------------------------------------------------------------------------------#




            #------------------------------------------------------------------------------#
            #     Set title and axes labels.                                               #
            #------------------------------------------------------------------------------#
            letitre = paste("Air temperature model")
            lexlab  = desc.unit(desc="Ambient temperature",unit=untab$degC)
            leylab  = desc.unit(desc="Sonic temperature"  ,unit=untab$degC)
            #------------------------------------------------------------------------------#




            #------------------------------------------------------------------------------#
            #     Loop through formats.                                                    #
            #------------------------------------------------------------------------------#
            for (o in sequence(nout)){
               #----- File name and format. -----------------------------------------------#
               fichier = file.path(outroot,paste0("tson2tmp-lm.",outform[o]))
               dummy   = open.file( fichier = fichier
                                  , outform = outform[o]
                                  , size    = sq.size
                                  , ptsz    = ptsz
                                  , depth   = depth
                                  )#end open.file
               #---------------------------------------------------------------------------#


               #----- Plot window and a grid. ---------------------------------------------#
               par(par.user)
               par(mar=c(4.1,4.6,2.1,1.6))
               plot.new()
               plot.window(xlim=xlimit,ylim=ylimit)
               axis(side=1,las=1)
               axis(side=2,las=1)
               title(main=letitre,xlab=lexlab,ylab=leylab)
               if (plotgrid) grid(col="grey63",lty="dotted",lwd=0.25)
               #---------------------------------------------------------------------------#



               #---------------------------------------------------------------------------#
               #      Plot points, the 1:1 line, and the linear fit.                       #
               #---------------------------------------------------------------------------#
               points(x=data.in$atm.tson-t00,y=data.in$atm.tmp-t00,type="p",col="sienna"
                     ,pch=16)
               abline(a=0,b=1,col="goldenrod",lwd=2,lty="dashed")
               lines (x=data.plot$atm.tson-t00,y=data.plot$atm.tmp-t00
                     ,col="red3",lty="solid",lwd=2)
               #----- Legend for the lines. -----------------------------------------------#
               legend("topleft",inset=0.01,legend=c("1:1 line","Linear fit"),bg="white"
                     ,col=c("goldenrod","red3"),lwd=2,lty=c("dashed","solid"))
               box()
               #---------------------------------------------------------------------------#


               #---------------------------------------------------------------------------#
               #      Find some good position to plot the fit.                             #
               #---------------------------------------------------------------------------#
               xt = xlimit[1] + 0.99 * diff(xlimit); yt = ylimit[1] + 0.10 * diff(ylimit)
               text(xt,yt,xyfit,adj=c(1,0.5),cex=0.7)
               xt = xlimit[1] + 0.99 * diff(xlimit); yt = ylimit[1] + 0.05 * diff(ylimit)
               text(xt,yt,r2text,adj=c(1,0.5),cex=0.7)
               xt = xlimit[1] + 0.99 * diff(xlimit); yt = ylimit[1] + 0.00 * diff(ylimit)
               text(xt,yt,ppval,adj=c(1,0.5),cex=0.7)
               #---------------------------------------------------------------------------#


               #---------------------------------------------------------------------------#
               #      Close the device.                                                    #
               #---------------------------------------------------------------------------#
               dummy = close.plot(outform=outform[o])
            }#end for (o in sequence(nout))
            #------------------------------------------------------------------------------#
         }#end if (plotregression)
         #---------------------------------------------------------------------------------#



         #------ Apply the linear fit to fill some data, then check for new outliers. -----#
         miss                       = is.na  (eft$atm.tmp) & is.finite(eft$atm.tson)
         eft$atm.tmp      [miss]    = data.pred$atm.tmp[miss]
         eft$gfflg.atm.tmp[miss]    = -1
         eft$gferr.atm.tmp[miss]    = summ.in$sigma
         eft$atm.tmp                = del.outliers( x        = eft$atm.tmp
                                                  , when     = eft$when
                                                  , out.hour = TRUE
                                                  , out.all  = FALSE
                                                  )#end del.outliers
         fixlist$atm.tmp$gap.filled = miss & is.finite(eft$atm.tmp)
         #---------------------------------------------------------------------------------#
      #------------------------------------------------------------------------------------#

   }#end if (use.similar)
   #=======================================================================================#
   #=======================================================================================#



   #=======================================================================================#
   #=======================================================================================#
   #     Remove the flags for variables that have gone missing.                            #
   #---------------------------------------------------------------------------------------#
   cat0("  + Remove flags for variables that have gone missing.")
   for (vv in sequence(nvargap)){
      #------------------------------------------------------------------------------------#
      #     Copy the settings for this variable to scratch variables.                      #
      #------------------------------------------------------------------------------------#
      this.vnam       = vargap$vari      [vv]
      this.desc       = vargap$desc      [vv]
      this.unit       = vargap$unit      [vv]
      this.harm       = vargap$harm      [vv]
      this.harm.tback = vargap$harm.tback[vv]
      this.oban       = vargap$oban      [vv]
      this.out.hour   = vargap$out.hour  [vv]
      this.form       = vargap$form      [vv]
      #------------------------------------------------------------------------------------#



      #----- Make a new variable containing both the gap filling flag and the error. ------#
      gfflg = paste0("gfflg.",this.vnam)
      gferr = paste0("gferr.",this.vnam)
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #     Initialise the gap filling flag and estimate by setting all points as NA.      #
      #------------------------------------------------------------------------------------#
      fine         = is.finite(eft[[this.vnam]]) & is.na(eft[[gfflg]])
      eft[[gfflg]][fine] = 0
      eft[[gferr]][fine] = 0
      #------------------------------------------------------------------------------------#
   }#end for
   #=======================================================================================#
   #=======================================================================================#



   #---------------------------------------------------------------------------------------#
   #     Save the compounded dataset so we don't load them again.                          #
   #---------------------------------------------------------------------------------------#
   cat0(" + Save self gap-filled data to file ",basename(session.selfgf),".")
   save(eft,fixlist,file=session.selfgf)
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #    Copy eft to selfgf.  From this point on, the cleaned data will not be called eft,  #
   # but selfgf because we always use eft as the structure name that is saved to the RData #
   # regardless of the level.                                                              #
   #---------------------------------------------------------------------------------------#
   selfgf = eft
   #---------------------------------------------------------------------------------------#

}#end if
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
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#     Check whether to run gap-filling with objective analysis or load previously computed #
# time series.                                                                             #
#------------------------------------------------------------------------------------------#
if (reload.objana && file.exists(session.objana)){
   #----- Reload the gap-filling. ---------------------------------------------------------#
   cat0(" + Load objective analysis data from file ",basename(session.objana),".")
   load(session.objana)
   objana  = eft
   ntimes  = length(objana$when)
   #---------------------------------------------------------------------------------------#
}else{
   #---------------------------------------------------------------------------------------#
   #     Find the gap-filled time series using objective analysis.                         #
   #---------------------------------------------------------------------------------------#

   #----- Initialise the data frame that will contain the gap-filled dataset. -------------#
   eft     = selfgf
   ntimes  = length(eft$when)
   #---------------------------------------------------------------------------------------#


   #----- Initialise the list that will contain the linear models used here. --------------#
   filllist = list()
   #---------------------------------------------------------------------------------------#






   #=======================================================================================#
   #=======================================================================================#
   #     METHOD 1.  We fill in the gaps using the objective analysis.                      #
   #                                                                                       #
   #     The objective analysis is done on the anomaly so we reduce the biases due to      #
   # altitude.  Alternatively, we could use variables that depend less on height such as   #
   # potential temperature and mean sea level pressure, but that would reduce the number   #
   # of valid points for the objective because these are combined variables and would be   #
   # undefined when one of the measurements were missing.                                  #
   #     The one exception to this is precipitation rate.  Rainfall is way too far from    #
   # a normal distribution, so we must use a combined statistics: we find the number of    #
   # sites that had rainfall, and assume that that was the probability of rainfall at that #
   # hour.  Then we run an objective analysis of the logarithm of the rainfall rate, and   #
   # use a random variable to decide whether to fill in with the objective analysis or set #
   # the rainfall to zero.                                                                 #
   #---------------------------------------------------------------------------------------#
   cat0( " + Find the gap-filled time series using objective analysis.")
   #----- Save the entire time series that comes from the objective analysis. -------------#
   objana.100 = list()
   #---------------------------------------------------------------------------------------#
   #    Now we loop over all variables, and check whether to fill in the data using a      #
   # Barnes' objective analysis.                                                           #
   #---------------------------------------------------------------------------------------#
   for (vv in sequence(nvargap)){
      #----- Copy the structure values to some scratch variables. -------------------------#
      this.vnam       = vargap$vari      [vv]
      this.desc       = vargap$desc      [vv]
      this.unit       = vargap$unit      [vv]
      this.harm       = vargap$harm      [vv]
      this.harm.tback = vargap$harm.tback[vv]
      this.oban       = vargap$oban      [vv]
      this.form       = vargap$form      [vv]
      #------------------------------------------------------------------------------------#




      #----- Make the names of the gap-filling flag and error. ----------------------------#
      this.gfflg = paste0("gfflg.",this.vnam)
      this.gferr = paste0("gferr.",this.vnam)
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #     Fill in with objective analysis if this variable is to be eft.              #
      #------------------------------------------------------------------------------------#
      if (this.oban){
         cat0("  - Fill ",this.desc," with objective analysis.")

         #---------------------------------------------------------------------------------#
         #     Set the hour and minute as an index so we find the mean diurnal cycle and   #
         # its standard deviation.                                                         #
         #---------------------------------------------------------------------------------#
         e.hhmm = paste0(sprintf("%2.2i",eft$month),"-"
                        ,sprintf("%2.2i",eft$hour),sprintf("%2.2i",eft$minu)
                        )#end paste0
         #---------------------------------------------------------------------------------#



         #----- Save data to a scratch variable. ------------------------------------------#
         this.data  = eft[[this.vnam]]
         ok.data    = is.finite(this.data)
         ndata      = length(this.data)
         #---------------------------------------------------------------------------------#


         #---------------------------------------------------------------------------------#
         #    Find the mean diurnal cycle and standard deviation of this variables in case #
         # it is not rainfall.                                                             #
         #---------------------------------------------------------------------------------#
         if (this.vnam %in% "rain"){
            #------------------------------------------------------------------------------#
            #     Rain.  We don't normalise the data but we must find the number of        #
            # finite observations by the hour of the day.                                  #
            #------------------------------------------------------------------------------#
            this.ndcycle = tapply(X=is.finite(this.data),INDEX=e.hhmm,FUN=sum ,na.rm=TRUE)
            #------------------------------------------------------------------------------#



            #----- Match the full time series with the average hour. ----------------------#
            this.idx    = match(e.hhmm,names(this.ndcycle))
            #------------------------------------------------------------------------------#
         }else{
            #------------------------------------------------------------------------------#
            #    Shortwave radiation.  We mustn't use the nighttime.                       #
            #------------------------------------------------------------------------------#
            if (this.vnam %in% radvars){
               this.data[eft$nighttime] = NA
            }#end if
            #------------------------------------------------------------------------------#


            #------------------------------------------------------------------------------#
            #     Find the mean diurnal cycle and the mean variability and number of valid #
            # points of the diel.                                                          #
            #------------------------------------------------------------------------------#
            this.location.dcycle = tapply(X=this.data,INDEX=e.hhmm,FUN=sn.location,na.rm=T)
            this.scale.dcycle    = tapply(X=this.data,INDEX=e.hhmm,FUN=sn.scale   ,na.rm=T)
            this.shape.dcycle    = tapply(X=this.data,INDEX=e.hhmm,FUN=sn.shape   ,na.rm=T)
            this.ndcycle         = tapply(X=ok.data  ,INDEX=e.hhmm,FUN=sum        ,na.rm=T)
            #------------------------------------------------------------------------------#



            #----- Match the full time series with the average hour. ----------------------#
            this.idx    = match(e.hhmm,names(this.ndcycle))
            #------------------------------------------------------------------------------#



            #------------------------------------------------------------------------------#
            #     This.reve will have the points for which the standard deviation is zero, #
            # which will be a singularity and therefore we will not solve using objective  #
            # analysis.  This is likely to happen only for radiation during the night, but #
            # we know it should be always zero.                                            #
            #------------------------------------------------------------------------------#
            this.reve  = ( is.finite(this.scale.dcycle[this.idx]) 
                         & this.scale.dcycle[this.idx] == 0)
            #------------------------------------------------------------------------------#
         }#end if
         #---------------------------------------------------------------------------------#


         #---------------------------------------------------------------------------------#
         #      Bind together the data from other sides.                                   #
         #---------------------------------------------------------------------------------#
         other.dist  = NULL
         other.iata  = NULL
         other.datn  = NULL
         other.reve  = NULL
         other.lon   = NULL
         other.lat   = NULL
         for (o in sequence(nother)){
            #----- Copy the variables to a scratch. ---------------------------------------#
            other.now  = fitz.aws[[o]]
            dist.now   = other.now$dist.ref
            desc.now   = other.now$longname
            iata.now   = other.now$iata
            data.full  = other.now[[this.vnam]]
            when.full  = other.now$when
            year.full  = other.now$year
            month.full = other.now$month
            hour.full  = other.now$hour
            minu.full  = other.now$minu
            lon.now    = other.now$lon
            lat.now    = other.now$lat
            #------------------------------------------------------------------------------#
            cat0("    * Load data from site ",desc.now,".")



            #------------------------------------------------------------------------------#
            #     Find the tag for all times.                                              #
            #------------------------------------------------------------------------------#
            hhmm.full  = paste0(sprintf("%2.2i",month.full),"-"
                               ,sprintf("%2.2i",hour.full),sprintf("%2.2i",minu.full)
                               )#end paste0
            #------------------------------------------------------------------------------#



            #------------------------------------------------------------------------------#
            #    Select the same time period of tower, in case the other sites have longer #
            # time periods.                                                                #
            #------------------------------------------------------------------------------#
            sel.now   = when.full %wr% c(eft$when[1],eft$when[ntimes])
            data.now  = data.full [sel.now]
            when.now  = when.full [sel.now]
            year.now  = when.full [sel.now]
            month.now = month.full[sel.now]
            hour.now  = hour.full [sel.now]
            minu.now  = minu.full [sel.now]
            o.hhmm    = hhmm.full [sel.now]
            #------------------------------------------------------------------------------#


            #------------------------------------------------------------------------------#
            #    Normalise the variable in case it is not rainfall.                        #
            #------------------------------------------------------------------------------#
            if (! ( this.vnam %in% "rain")){
               #----- Find the names of the statistics for this variable. -----------------#
               location.name = paste0("qlocation.",this.vnam)
               scale.name    = paste0("qscale."   ,this.vnam)
               shape.name    = paste0("qshape."   ,this.vnam)
               #---------------------------------------------------------------------------#



               #---------------------------------------------------------------------------#
               #     Loop through all possible base periods for obtaining the statistics.  #
               # In case none of them are suitable, make all statistics NA.                #
               #---------------------------------------------------------------------------#
               stat.fine = FALSE
               oo        = 0
               while ( (! stat.fine) & (oo < nobjana) ){
                  oo       = oo + 1
                  sel.cal  = year.full %wr% c(yeara.objana[oo],yearz.objana[oo])
                  data.cal = data.full[sel.cal]
                  c.hhmm   = hhmm.full[sel.cal]

                  #---- Find the monthly-dependent diel statistics. -----------------------#
                  location.now = tapply(X=data.cal,INDEX=c.hhmm,FUN=sn.location,na.rm=TRUE)
                  scale.now    = tapply(X=data.cal,INDEX=c.hhmm,FUN=sn.scale   ,na.rm=TRUE)
                  shape.now    = tapply(X=data.cal,INDEX=c.hhmm,FUN=sn.shape   ,na.rm=TRUE)
                  names.now    = names(location.now)
                  #------------------------------------------------------------------------#

                  stat.fine = all(is.finite(c(location.now,scale.now,shape.now)))
               }#end while ( (! stat.fine) & (oo < nobjana) )
               #---------------------------------------------------------------------------#


               #----- In case none of the periods worked, discard the statistics. ---------#
               if (stat.fine){
                  cat0( "      > Statistics based on data from "
                      , yeara.objana[oo],"-",yearz.objana[oo],"."
                      )#end cat0
               }else{
                  cat0( "      > Time series doesn't have sufficient data, skip site.")
                  location.now = location.now + NA
                  scale.now    = scale.now    + NA
                  shape.now    = shape.now    + NA
               }#end if
               #---------------------------------------------------------------------------#


               #----- Match the full time series with the average hour. -------------------#
               idx.now     = match(o.hhmm,names.now)
               #---------------------------------------------------------------------------#


               #---------------------------------------------------------------------------#
               #      Make the normalised variable.  Reve now is used only to determine    #
               # whether the standard deviation is absolutely zero.  This would be a       #
               # singularity so we don't trust the objective analysis for these points.    #
               # This is likely to happen only for shortwave radiation during the night,   #
               # which is safe to assume that is always zero (we will ignore the effect    #
               # of fireflies, but I think it is fine...).                                 #
               #---------------------------------------------------------------------------#
               datn.now           = skew2normal( x        = data.now
                                               , location = location.now
                                               , scale    = scale.now
                                               , shape    = shape.now
                                               , idx      = idx.now
                                               )#end skew2normal
               reve.now           = scale.now[idx.now] == 0.
               datn.now[reve.now] = NA
               #---------------------------------------------------------------------------#
            }else{
               #----- Set rain events to 1, no rain events to 0, and missing data to NA. --#
               reve.now = as.integer(data.now > 0.0)
               #---------------------------------------------------------------------------#



               #----- Copy the logarithm of rainfall, removing no-rain events. ------------#
               del             = ! (reve.now %>% 0.0)
               data.now[del]   = NA
               datn.now        = data.now
               #---------------------------------------------------------------------------#
            }#end if
            #------------------------------------------------------------------------------#



            #------------------------------------------------------------------------------#
            #     Concatenate this vector to the working dataset.                          #
            #------------------------------------------------------------------------------#
            other.dist  =     c(other.dist ,dist.now )
            other.iata  =     c(other.iata ,iata.now )
            other.lon   =     c(other.lon  ,lon.now  )
            other.lat   =     c(other.lat  ,lat.now  )
            other.datn  = cbind(other.datn ,datn.now )
            other.reve  = cbind(other.reve ,reve.now )
            #------------------------------------------------------------------------------#
         }#end for
         #---------------------------------------------------------------------------------#


         #----- Name the objects so we have a better way to track troublesome sites. ------#
         names   (other.dist ) = other.iata
         dimnames(other.datn ) = list(NULL,other.iata)
         #---------------------------------------------------------------------------------#


         #------Make a matrix with the coordinates to find the typical site distance . ----#
         other.coord           = cbind(other.lon,other.lat)
         dimnames(other.coord) = list(other.iata,c("lon","lat"))
         #---------------------------------------------------------------------------------#


         #---------------------------------------------------------------------------------#
         #     Run the objective analysis.                                                 #
         #---------------------------------------------------------------------------------#
         if (this.vnam %in% "rain"){
            #------------------------------------------------------------------------------#
            #     Rainfall.  We run two objective analysis, one for the binary flag tell-  #
            # ing whether it rained or not, and another for the log of the precipitation   #
            # rate.  The result of the first objective analysis will be the probability    #
            # that it rained at the reference site, and the second objective analysis will #
            # tell the magnitude of the event in case we decide that it rained given the   #
            # probability.  By doing that we will always use the rainfall estimate if all  #
            # sites available reported rain, and we will always assign zero rainfall if    #
            # none of the sites reported rain.  If some of the sites reported rain whilst  #
            # other sites did not, then we will assign rain only with a certain            #
            # probability.  We use the objective analysis to determine the probability so  #
            # the chances will be biased towards the closest measurements: if sites next   #
            # to the reference reported rain whilst sites further away did not, then it    #
            # should be more likely that rain happened at the reference site than if it    #
            # was the other way round.                                                     #
            #------------------------------------------------------------------------------#
            prob.oban = obj.analysis( ref.lola = t(c(lon,lat))
                                    , rem.lola = other.coord
                                    , rem.dat  = other.reve
                                    , verbose  = TRUE )
            rain.oban = obj.analysis( ref.lola = t(c(lon,lat))
                                    , rem.lola = other.coord
                                    , rem.dat  = other.datn
                                    , verbose  = TRUE ) 
            #------------------------------------------------------------------------------#
            #     Fill in the rainfall events.                                             #
            #------------------------------------------------------------------------------#
            objana.100$rain         = rain.oban$first
            #------------------------------------------------------------------------------#
            #     We "toss the dice" for each event, using a random variable with uniform  #
            # distribution between 0+epsilon and 1-epsilon.  We will assume that it didn't #
            # rain when the values of the dice are greater than the probability.  If all   #
            # sites reported rain then the probability is one so the dice will never be    #
            # greater than the probability.  Conversely, if none of the sites reported     #
            # rain, then the probability will always be less than the dice, so it will     #
            # always assume that no rain fell at the site.                                 #
            #------------------------------------------------------------------------------#
            epsil                   = .Machine$double.neg.eps
            dice                    = runif(n=ndata,min=epsil,max=1.-epsil)
            norain                  = ( is.finite(prob.oban$first) 
                                      & dice > prob.oban$first )
            objana.100$rain[norain] = 0.0
         }else{
            #------------------------------------------------------------------------------#
            #     Other variables.  Run the objective analysis.  The result from the       #
            # objective analysis will be the normalised value, so we convert back to       #
            # actual measurements.                                                         #
            #------------------------------------------------------------------------------#
            this.oban               = obj.analysis( ref.lola = t(c(lon,lat))
                                                  , rem.lola = other.coord
                                                  , rem.dat  = other.datn
                                                  , verbose  = TRUE
                                                  )#end obj.analysis
            objana.100[[this.vnam]] = normal2skew ( xnorm    = this.oban$first
                                                  , location = this.location.dcycle
                                                  , scale    = this.scale.dcycle
                                                  , shape    = this.shape.dcycle
                                                  , idx      = this.idx
                                                  )#end normal2skew
            nused                   = this.oban$nvalid
            #----- Discard data for when standard deviation is zero. ----------------------#
            objana.100[[this.vnam]][this.reve] = NA
            #------------------------------------------------------------------------------#
         }#end if
         #---------------------------------------------------------------------------------#



         #----- Copy the objective analysis to the missing values. ------------------------#
         miss                       = ( is.na    (eft[[this.vnam]]) 
                                      & is.finite(objana.100[[this.vnam]]) )
         eft[[this.vnam ]][miss] = objana.100[[this.vnam]][miss]
         eft[[this.gfflg]][miss] = 1
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #     Estimate the error.  Because we normalised the data by the hour of the day, #
         # we make separate estimates of the error for each hour.                          #
         #---------------------------------------------------------------------------------#
         avail                      = ( is.finite(eft[[this.vnam]]) 
                                      & is.finite(objana.100[[this.vnam]]) )
         gfill.ndcycle              = tapply( X     = is.finite(objana.100[[this.vnam]])
                                            , INDEX = e.hhmm
                                            , FUN   = sum
                                            , na.rm = TRUE
                                            )#tapply
         this.navail                = gfill.ndcycle[this.idx]
         #----- Estimate the error for each hour. -----------------------------------------#
         err.fact                   = ( (objana.100[[this.vnam]] - clean[[this.vnam]])^2 
                                      / (this.navail-3*nused) )
         rmse.cycle                 = sqrt( tapply( X     = err.fact
                                                  , INDEX = e.hhmm
                                                  , FUN   = sum
                                                  , na.rm = TRUE
                                                  )#end tapply
                                          )#end sqrt
         rmse                       = rmse.cycle[this.idx]
         #----- Assign the error for the missing values for the given hour of the day. ----#
         eft[[this.gferr]][miss] = rmse[miss]
         #---------------------------------------------------------------------------------#


         #---------------------------------------------------------------------------------#
         #   If this is pressure, update the vapour pressure, dew point temperature, and   #
         # relative humidity based upon the gap-filled pressure.  Assign a gap-filling     #
         # flag of -2 since it is not "pure".  Later we fix the gap-filling flag for       #
         # specific humidity to 0 if gap filling flag for vapour pressure is -2.           #
         #---------------------------------------------------------------------------------#
         if (this.vnam %in% c("atm.pvap","atm.prss")){
            cat0("    * Fill in derived humidity data.")
            #----- Specific humidity. =----------------------------------------------------#
            miss                     = ( is.finite(eft$atm.pvap) & is.finite(eft$atm.prss)
                                       & is.na(eft$atm.shv) )
            eft$atm.shv       [miss] = ( ep * eft$atm.pvap[miss]
                                       / (eft$atm.prss[miss] - (1.-ep)*eft$atm.pvap[miss])
                                       )#end eft$atm.shv[miss]
            #----- Dew point temperature. -------------------------------------------------#
            miss                     = ( is.finite(eft$atm.pvap) & is.finite(eft$atm.prss)
                                       & is.na(eft$atm.tdew) )
            eft$atm.tdew      [miss] = tslif( pres     = eft$atm.prss[miss]
                                            , hum      = eft$atm.pvap[miss]
                                            , type.hum = "pvap"              )
            eft$gfflg.atm.tdew[miss] = -2
            eft$glerr.atm.tdew[miss] = NA
            #------------------------------------------------------------------------------#


            #------ Relative humidity doesn't directly depend on pressure. ----------------#
            if (this.vnam %in% "atm.pvap"){
               #----- Fill in relative humidity for all points we can do it. --------------#
               miss                   = ( is.finite(eft$atm.pvap) & is.finite(eft$atm.tmp)
                                        & is.na(eft$atm.rhv) )
               eft$atm.rhv    [miss]  = eft$atm.pvap[miss] / eslif(eft$atm.tmp[miss])
               #---------------------------------------------------------------------------#
            }#end if (this.vnam %in% "atm.pvap")
            miss = is.finite(eft$atm.prss) & is.finite(eft$atm.shv)
            #------------------------------------------------------------------------------#
         }#end if
         #---------------------------------------------------------------------------------#
      }#end if
      #------------------------------------------------------------------------------------#
   }#end for
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #      Run radiation sanity check again, to remove bogus intepolated values.            #
   #---------------------------------------------------------------------------------------#
   cat0(" + Remove suspicious incoming shortwave radiation.")
   eft                        = del.bad.rshort(dat=eft)
   #----- Update the gap-filled status in case we delete data. ----------------------------#
   miss                       = is.na(eft$rshort.in)
   eft$gfflg.rshort.in[miss]  = NA
   eft$gferr.rshort.in[miss]  = NA
   miss                       = is.na(eft$rshort.out)
   eft$gfflg.rshort.out[miss] = NA
   eft$gferr.rshort.out[miss] = NA
   #---------------------------------------------------------------------------------------#
   #=======================================================================================#
   #=======================================================================================#



   #----- Find the shortwave radiation breakdown. -----------------------------------------#
   cat0( " + Fill in additional PAR.")
   rsbreak                = rshort.bdown( rad.in   = eft$rshort.in
                                        , atm.prss = ref.prss + 0. * eft$cosz
                                        , cosz     = eft$cosz
                                        , rad.type = "rshort"
                                        )#end rshort.bdown
   miss                   = is.na(eft$par.in) & is.finite(rsbreak$par.full)
   eft$par.in      [miss] = rsbreak$par.full[miss]
   eft$gfflg.par.in[miss] = 1
   eft$gferr.par.in[miss] = ( eft$gferr.rshort.in[miss] 
                            * rsbreak$par.full   [miss] / eft$rshort.in[miss] )
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #     Save the objective analysis dataset so we don't load them again.                  #
   #---------------------------------------------------------------------------------------#
   cat0(" + Save objective analysis to file ",basename(session.objana),".")
   save(eft,objana.100,filllist,file=session.objana)
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #    Copy eft to objana.  From this point on, the objective analysis data will not be   #
   # called eft, but objana because we always use eft as the structure name that is saved  #
   # to the RData regardless of the level.                                                 #
   #---------------------------------------------------------------------------------------#
   objana = eft
   #---------------------------------------------------------------------------------------#
}#end if
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
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#     Check whether to run gap-filling with harmonic analysis or load previously computed  #
# time series.                                                                             #
#------------------------------------------------------------------------------------------#
if (reload.harmana && file.exists(session.harmana)){
   #----- Reload the gap-filling. ---------------------------------------------------------#
   cat0(" + Load harmonic analysis data from file ",basename(session.harmana),".")
   load(session.harmana)
   harmana  = eft
   ntimes   = length(harmana$when)
}else{
   #---------------------------------------------------------------------------------------#
   #     Find the gap-filled time series using harmonic analysis.                          #
   #---------------------------------------------------------------------------------------#

   #----- Initialise the data frame that will contain the gap-filled dataset. -------------#
   eft     = objana
   ntimes  = length(eft$when)
   #---------------------------------------------------------------------------------------#


   #----- Initialise the list that will contain the linear models used here. --------------#
   filllist = list()
   #---------------------------------------------------------------------------------------#




   #=======================================================================================#
   #=======================================================================================#
   #    METHOD 2:  Fill remaining gaps with harmonic analysis (iterative Fourier trans-    #
   #               form).  This method will fill all the remaining gaps using the most     #
   #               significant powers of the Fourier analysis.                             #
   #               This method will produce a continuous and nice looking gap fill when    #
   #               gaps are short, but it tends to dampen the oscillations and make them   #
   #               look overly smoothed.  This happens because we don't use all the powers #
   #               and we don't use all the powers for two important reasons:              #
   #               1.  It takes really long to run with all powers                         #
   #               2.  The weakest powers are more associated with random noise than       #
   #                   actual oscillations, and adding them may actually increase the      #
   #                   errors.  Retaining 75-80% of the signal is usually fine for time    #
   #                   series with strong oscillation patterns.                            #
   #                                                                                       #
   #               Since this method has some shortcomings, we use them only when other    #
   #               methods have been applied and some remaining gaps still exist.          #
   #---------------------------------------------------------------------------------------#
   cat0( " + Run harmonic analyses for some variables.")
   #---- Fourier is the list with the variables that contains the gap-eft time series. -#
   fourier = list()
   for (vv in sequence(nvargap)){
      #----- Copy the structure values to some scratch variables. -------------------------#
      this.vnam        = vargap$vari      [vv]
      this.desc        = vargap$desc      [vv]
      this.unit        = vargap$unit      [vv]
      this.harm        = vargap$harm      [vv]
      this.harm.tback  = vargap$harm.tback[vv]
      this.oban        = vargap$oban      [vv]
      this.form        = vargap$form      [vv]
      this.potential   = vargap$potential [vv]
      this.sw.typevar  = vargap$sw.typevar[vv]
      #------------------------------------------------------------------------------------#



      #----- Make the names of the gap-filling flag and error. ----------------------------#
      this.gfflg = paste0("gfflg.",this.vnam)
      this.gferr = paste0("gferr.",this.vnam)
      #------------------------------------------------------------------------------------#


      #----- Flag the missing values. -----------------------------------------------------#
      miss       = is.na(eft[[this.vnam]])
      nmiss      = sum(miss)
      #------------------------------------------------------------------------------------#

      if ( this.harm & nmiss > 0 ){
         cat0("  - Fill ",this.desc," with harmonic filling.")


         #----- Copy time series to some auxiliary variables. -----------------------------#
         hf.index = seq(from=1,to=ntimes,by=1)
         hf.guess = eft[[this.vnam]]
         hf.error = rep(0,times=ntimes)
         #---------------------------------------------------------------------------------#


         #---------------------------------------------------------------------------------#
         #      Loop over all sub-periods.                                                 #
         #---------------------------------------------------------------------------------#
         for (h in sequence(nharm)){
            #----- Select the period of interest. -----------------------------------------#
            if (h == nharm){
               cat0("    * Entire period","\n")
               sel = rep(TRUE,times=ntimes)
            }else{
               cat("    * Period ",h,": ",yeara.harm[h],"-",yearz.harm[h],"\n")
               sel   = eft$year %wr% c(yeara.harm[h],yearz.harm[h])
            }#end if
            #------------------------------------------------------------------------------#



            #------------------------------------------------------------------------------#
            #      In case this is a sun variable, use the relative value, otherwise use   #
            # the absolute value.   Because sun variables are not available during the     #
            # night, we fill the values with simple interpolation of daytime relative      #
            # values when the gap is not too long.                                         #
            #------------------------------------------------------------------------------#
            if (this.sw.typevar){
               vnnow = ifelse( test = eft$daytime
                             , yes  = eft[[this.vnam]] / eft[[this.potential]]
                             , no   = NA)
               ntnow = na.approx(vnnow,na.rm=FALSE,maxgap=24)
               tsnow = ifelse(test=eft$daytime[sel],yes=vnnow[sel],no=ntnow[sel])
               hfnow = hf.index[sel]
            }else{
               tsnow = eft[[this.vnam]][sel]
               hfnow = hf.index[sel]
            }#end if (this.sw.typevar)
            #------------------------------------------------------------------------------#


            #----- Find the gap-eft time series (and error if sought). -----------------#
            fourier.fit = harmfill( x              = tsnow
                                  , detrend.method = detrend.method
                                  , trend.back     = this.harm.tback
                                  , signal.retain  = signal.retain
                                  , conv.threshold = conv.threshold
                                  , minmod         = minmod
                                  , maxmod         = maxmod
                                  , maxin          = maxin
                                  , verbose        = harmverb
                                  , rmse           = err.jackknife
                                  , del.frac       = del.jackknife
                                  , n.jack         = max.jackknife )
            hf.guess[hfnow] = fourier.fit$xfit
            hf.error[hfnow] = fourier.fit$error
            #------------------------------------------------------------------------------#
         }#end for (h in sequence(nharm))
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #     Replace the remaining NA by the values found by the Fourier analysis.       #
         #---------------------------------------------------------------------------------#
         fourier[[this.vnam]]        = hf.guess
         if (this.sw.typevar){
            #----- Copy the results from the Fourier, scaling by the potential variable. --#
            if (this.harm.tback){
               eft [[this.vnam ]][miss] = hf.guess[miss] * eft[[this.potential]][miss]
               eft [[this.gferr]][miss] = hf.error[miss]
               eft [[this.gfflg]][miss] = 2
            }else{
               eft [[this.vnam ]]       = hf.guess * eft[[this.potential]]
               eft [[this.gferr]][miss] = hf.error[miss]
               eft [[this.gfflg]][miss] = 2
            }#end if
            #------------------------------------------------------------------------------#
         }else{
            #----- Copy the results from the Fourier directly. ----------------------------#
            if (this.harm.tback){
               eft [[this.vnam ]][miss] = hf.guess[miss]
               eft [[this.gferr]][miss] = hf.error[miss]
               eft [[this.gfflg]][miss] = 2
            }else{
               eft [[this.vnam ]]       = hf.guess
               eft [[this.gferr]][miss] = hf.error[miss]
               eft [[this.gfflg]][miss] = 2
            }#end if
            #------------------------------------------------------------------------------#
         }#end if (this.sw.typevar)
         #---------------------------------------------------------------------------------#
      }#end if
      #------------------------------------------------------------------------------------#
   }#end for
   #=======================================================================================#
   #=======================================================================================#



   #---------------------------------------------------------------------------------------#
   #      Run radiation sanity check again, to remove bogus intepolated values.            #
   #---------------------------------------------------------------------------------------#
   cat0(" + Remove suspicious incoming shortwave radiation.")
   eft                        = del.bad.rshort(dat=eft)
   #----- Update the gap-filled status in case we delete data. ----------------------------#
   miss                       = is.na(eft$rshort.in)
   eft$gfflg.rshort.in[miss]  = NA
   eft$gferr.rshort.in[miss]  = NA
   miss                       = is.na(eft$rshort.out)
   eft$gfflg.rshort.out[miss] = NA
   eft$gferr.rshort.out[miss] = NA
   #---------------------------------------------------------------------------------------#
   #=======================================================================================#
   #=======================================================================================#



   #----- Find the shortwave radiation breakdown. -----------------------------------------#
   cat0( " + Fill in additional PAR.")
   rsbreak                = rshort.bdown( rad.in   = eft$rshort.in
                                        , atm.prss = ref.prss + 0. * eft$cosz
                                        , cosz     = eft$cosz
                                        , rad.type = "rshort"
                                        )#end rshort.bdown
   miss                   = is.na(eft$par.in) & is.finite(rsbreak$par.full)
   eft$par.in      [miss] = rsbreak$par.full[miss]
   eft$gfflg.par.in[miss] = 1
   eft$gferr.par.in[miss] = ( eft$gferr.rshort.in[miss] 
                            * rsbreak$par.full   [miss] / eft$rshort.in[miss] )
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Save the harmonic analysis dataset so we don't load them again.                   #
   #---------------------------------------------------------------------------------------#
   cat0(" + Save harmonic analysis to file ",basename(session.harmana),".")
   save(eft,filllist,file=session.harmana)
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #    Copy eft to harmana.  From this point on, the harmonic analysis data will not be   #
   # called eft, but harmana because we always use eft as the structure name that is saved #
   # to the RData regardless of the level.                                                 #
   #---------------------------------------------------------------------------------------#
   harmana = eft
   #---------------------------------------------------------------------------------------#
}#end if
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
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#     Check whether to run gap-filling with spline or load previously computed time        #
# series.                                                                                  #
#------------------------------------------------------------------------------------------#
if (reload.fspline && file.exists(session.fspline)){
   #----- Reload the gap-filling. ---------------------------------------------------------#
   cat0(" + Load spline-filled data from file ",basename(session.fspline),".")
   load(session.fspline)
   fspline  = eft
   ntimes   = length(fspline$when)
}else{
   #---------------------------------------------------------------------------------------#
   #     Find the gap-filled time series using spline.                                     #
   #---------------------------------------------------------------------------------------#

   #----- Initialise the data frame that will contain the gap-filled dataset. -------------#
   eft     = harmana
   ntimes  = length(eft$when)
   #---------------------------------------------------------------------------------------#


   #----- Initialise the list that will contain the linear models used here. --------------#
   filllist = list()
   #---------------------------------------------------------------------------------------#




   #=======================================================================================#
   #=======================================================================================#
   #    METHOD 3.5:  Fill remaining gaps with spline.   We define the function based on    #
   #                 the anomalies of the diurnal cycle in case it is a non-radiation      #
   #                 variable, and the fraction of total incoming radiation in case it is  #
   #                 a radiation variable.  This is not the best method to fill in the     #
   #                 data, and should be used only to fill in some minor gaps in the       #
   #                 absence of any better method.                                         #
   #---------------------------------------------------------------------------------------#
   cat0( " + Run spline filling for some variables.")
   #---- Fourier is the list with the variables that contains the gap-eft time series. -#
   for (vv in sequence(nvargap)){
      #----- Copy the structure values to some scratch variables. -------------------------#
      this.vnam        = vargap$vari      [vv]
      this.desc        = vargap$desc      [vv]
      this.unit        = vargap$unit      [vv]
      this.out.hour    = vargap$out.hour  [vv]
      this.out.all     = vargap$out.all   [vv]
      this.fspline     = vargap$fspline   [vv]
      this.potential   = vargap$potential [vv]
      this.sw.typevar  = vargap$sw.typevar[vv]
      #------------------------------------------------------------------------------------#



      #----- Make the names of the gap-filling flag and error. ----------------------------#
      this.gfflg = paste0("gfflg.",this.vnam)
      this.gferr = paste0("gferr.",this.vnam)
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #      Run an outlier check for this function in case we will run spline.            #
      #------------------------------------------------------------------------------------#
      if (this.fspline){
         cat0("  - Fill ",this.desc," with spline.")


         #---- Check whether objective+harmonic analysis introduced bogus entries. --------#
         cat0("    * Remove outliers by hour and by day.")
         thisvar          = eft[[this.vnam]]
         eft[[this.vnam]] = del.outliers( x        = thisvar
                                        , when     = eft$when
                                        , out.hour = this.out.hour
                                        , out.all  = this.out.all
                                        )#end del.outliers
         #---------------------------------------------------------------------------------#



         #----- Flag the missing values. --------------------------------------------------#
         miss       = is.na(eft[[this.vnam]])
         nmiss      = sum(miss)
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #      Check whether filling is necessary.                                        #
         #---------------------------------------------------------------------------------#
         cat0("    * Standardise time series.")
         if ( nmiss > 0 ){
            #------------------------------------------------------------------------------#
            #      In case this is a sun variable, discard nighttime data.                 #
            #------------------------------------------------------------------------------#
            if (this.sw.typevar){
               #----- Normalise time series, using potential values. ----------------------#
               datf.now = ifelse(test=eft$nighttime,yes=NA,no=eft[[this.vnam]])
            }else{
               datf.now = eft[[this.vnam]]
            }#end if
            #------------------------------------------------------------------------------#


            #------------------------------------------------------------------------------#
            #     Set the hour and minute as an index so we find the mean diurnal cycle    #
            # and its standard deviation.                                                  #
            #------------------------------------------------------------------------------#
            f.hhmm = paste0( sprintf("%2.2i",eft$month),"-"
                           , sprintf("%2.2i",eft$hour),sprintf("%2.2i",eft$minu)
                           )#end paste0
            #------------------------------------------------------------------------------#

            #----- Find the statistics. ---------------------------------------------------#
            location.now = tapply(X=datf.now,INDEX=f.hhmm,FUN=sn.location,na.rm=T)
            scale.now    = tapply(X=datf.now,INDEX=f.hhmm,FUN=sn.scale   ,na.rm=T)
            shape.now    = tapply(X=datf.now,INDEX=f.hhmm,FUN=sn.shape   ,na.rm=T)
            names.now    = names(location.now)
            #------------------------------------------------------------------------------#


            #----- Match the full time series with the average hour. ----------------------#
            idx.now     = match(f.hhmm,names.now)
            #------------------------------------------------------------------------------#


            #------------------------------------------------------------------------------#
            #      Make the normalised variable.                                           #
            #------------------------------------------------------------------------------#
            datn.now           = skew2normal( x        = datf.now
                                            , location = location.now
                                            , scale    = scale.now
                                            , shape    = shape.now
                                            , idx      = idx.now
                                            )#end skew2normal
            #------------------------------------------------------------------------------#



            #------------------------------------------------------------------------------#
            #      Fill in nighttime data with dummy interpolation.                        #
            #------------------------------------------------------------------------------#
            if (this.sw.typevar){
               #----- Normalise time series, using potential values. ----------------------#
               datz.now = na.approx(datn.now,na.rm=FALSE)
               datn.now = ifelse(test=eft$nighttime,yes=datz.now,no=datn.now)
            }#end if
            #------------------------------------------------------------------------------#



            #------------------------------------------------------------------------------#
            #    Select the times with finite entries.                                     #
            #------------------------------------------------------------------------------#
            cat0("    * Derive then apply spline function.")
            sel       = is.finite(datn.now)
            fun.fill  = splinefun( x      = as.numeric(eft$when[sel])
                                 , y      = datn.now[sel]
                                 , method = "monoH.FC"
                                 )#end splinefun
            datn.fill = fun.fill(x=as.numeric(eft$when))
            #------------------------------------------------------------------------------#



            #------------------------------------------------------------------------------#
            #     Back-scale filled data.                                                  #
            #------------------------------------------------------------------------------#
            cat0("    * Re-scale filled data.")
            data.fill = normal2skew( xnorm    = datn.fill
                                   , location = location.now
                                   , scale    = scale.now
                                   , shape    = shape.now
                                   , idx      = idx.now
                                   )#end normal2skew 
            #------------------------------------------------------------------------------#



            #------------------------------------------------------------------------------#
            #      In case this is a shortwave variable, make sure the resulting value     #
            # will not exceed the maximum for the hour or be below zero.                   #
            #------------------------------------------------------------------------------#
            if (this.sw.typevar){
               #----- Normalise time series, using potential values. ----------------------#
               data.fill = pmax(0.,pmin(data.fill,eft[[this.potential]]))
               #---------------------------------------------------------------------------#
            }#end if
            #------------------------------------------------------------------------------#


            #------------------------------------------------------------------------------#
            #     Replace the remaining NA by the values found by the Fourier analysis.    #
            #------------------------------------------------------------------------------#
            eft[[this.vnam ]][miss] = data.fill[miss]
            eft[[this.gferr]][miss] = NA
            eft[[this.gfflg]][miss] = 5
            #------------------------------------------------------------------------------#
         }#end if (nmiss > 0)
         #---------------------------------------------------------------------------------#
      }#end if (this.fspline)
      #------------------------------------------------------------------------------------#
   }#end for (vv in sequence(nvargap))
   #=======================================================================================#
   #=======================================================================================#



   #----- Find the shortwave radiation breakdown. -----------------------------------------#
   cat0( " + Fill in additional PAR.")
   rsbreak                = rshort.bdown( rad.in   = eft$rshort.in
                                        , atm.prss = ref.prss + 0. * eft$cosz
                                        , cosz     = eft$cosz
                                        , rad.type = "rshort"
                                        )#end rshort.bdown
   miss                   = is.na(eft$par.in) & is.finite(rsbreak$par.full)
   eft$par.in      [miss] = rsbreak$par.full[miss]
   eft$gfflg.par.in[miss] = 6
   eft$gferr.par.in[miss] = ( eft$gferr.rshort.in[miss] 
                            * rsbreak$par.full   [miss] / eft$rshort.in[miss] )
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Save the spline filling dataset so we don't load them again.                   #
   #---------------------------------------------------------------------------------------#
   cat0(" + Save spline filling to file ",basename(session.fspline),".")
   save(eft,filllist,file=session.fspline)
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #    Copy eft to fspline.  From this point on, the harmonic analysis data will not be   #
   # called eft, but fspline because we always use eft as the structure name that is saved #
   # to the RData regardless of the level.                                                 #
   #---------------------------------------------------------------------------------------#
   fspline = eft
   #---------------------------------------------------------------------------------------#
}#end if
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
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#     Check whether to run gap-filling with other method or load previously computed       #
# time series.                                                                             #
#------------------------------------------------------------------------------------------#
if (reload.metfill && file.exists(session.metfill)){
   #----- Reload the gap-filling. ---------------------------------------------------------#
   cat0(" + Load gap-filled data (except NEE) from file ",basename(session.metfill),".")
   load(session.metfill)
   metfill  = eft
   ntimes   = length(metfill$when)
   #---------------------------------------------------------------------------------------#
}else{
   #---------------------------------------------------------------------------------------#
   #     Find the gap-filled time series using other methods.                              #
   #---------------------------------------------------------------------------------------#

   #----- Initialise the data frame that will contain the gap-filled dataset. -------------#
   eft     = fspline
   ntimes  = length(eft$when)
   #---------------------------------------------------------------------------------------#






   #=======================================================================================#
   #=======================================================================================#
   #    METHOD 3: Fill the remaining variables with some linear model.  Although this is   #
   #              similar to 1, this may use gap-eft variables which will likely        #
   #              increase the uncertainty (and we must take that into account too when    #
   #              estimating the error).                                                   #
   #---------------------------------------------------------------------------------------#
   cat0( " + Run a simple parametric model for some variables.")
   for (vv in sequence(nvargap)){
      #----- Copy the structure values to some scratch variables. -------------------------#
      this.gap         = vargap
      this.vnam        = vargap$vari      [vv]
      this.desc        = vargap$desc      [vv]
      this.unit        = vargap$unit      [vv]
      this.harm        = vargap$harm      [vv]
      this.harm.tback  = vargap$harm.tback[vv]
      this.oban        = vargap$oban      [vv]
      this.form        = vargap$form      [vv]
      #------------------------------------------------------------------------------------#




      #----- Make the names of the gap-filling flag and error. ----------------------------#
      this.gfflg = paste0("gfflg.",this.vnam)
      this.gferr = paste0("gferr.",this.vnam)
      #------------------------------------------------------------------------------------#


      if (! this.harm && (! is.na(this.form))){
         cat0("  - Fill ",this.desc," with parametric filling.")

         #----- Create a method 3 list in case it hasn't been created yet. ----------------#
         if (! "m3" %in% names(filllist)) filllist$m3 = list()
         #---------------------------------------------------------------------------------#



         #----- Apply the model then predict with the eft time series. --------------------#
         fit.in      = rlm(formula = as.formula(this.form),data = clean,maxit=100)
         if (! fit.in$converged){
             stop(paste0("Parametric fit didn't converge for ",this.desc,"."))
         }#end if
         summ.in     = summary(fit.in)
         this.pred   = predict(object=fit.in,newdata=eft)
         coeff       = coefficients(fit.in)
         ncoeff      = length(coeff)
         name.coeff  = names(coeff)
         pvalue.in   = 2. * pt(-abs(summ.in$coefficients[,3]),df=summ.in$df[2])
         #---------------------------------------------------------------------------------#



         #----- Find the adjusted R2. -----------------------------------------------------#
         ss.err                   = sum(summ.in$residuals^2)
         df.err                   = summ.in$df[2]
         mean.y                   = mean(eft[[this.vnam]],na.rm=TRUE)
         ss.tot                   = sum((eft[[this.vnam]]-mean.y)^2,na.rm=TRUE)
         df.tot                   = sum(is.finite(eft[[this.vnam]]))-1
         r2.in                    = 1.0 - ss.err * df.tot / (ss.tot * df.err)
         filllist$m3[[this.vnam]] = list( fit       = fit.in
                                        , summ      = summ.in
                                        , pvalue.in = pvalue.in
                                        , r2.in     = r2.in
                                        )#end list
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #      Fill in missing data with the model prediction.                            #
         #---------------------------------------------------------------------------------#
         miss                       = is.na(eft[[this.vnam]])
         eft[[this.vnam ]][miss] = this.pred[miss]
         eft[[this.gfflg]][miss] = 3
         #---------------------------------------------------------------------------------#


         #---------------------------------------------------------------------------------#
         #     Estimate the error.  The error is going to be the RMSE of the model fit and #
         # the mean error associated with previous gap filling of the explaining vari-     #
         # ables.                                                                          #
         #---------------------------------------------------------------------------------#
         #----- Add contribution from the gap filling uncertainty. ------------------------#
         err2       = (summ.in$sigma)^2
         for (cc in sequence(ncoeff)[-1]){
            #---- Propagate the error due to the other variables. -------------------------#
            exp.vnam = paste("gferr",name.coeff[cc],sep=".")
            err2.var = (coeff[cc] * eft[[exp.vnam]])^2
            err2     = err2 + err2.var
         }#end for
         #----- The is the square root of sum of the errors. ------------------------------#
         eft[[this.gferr]][miss] = sqrt(err2[miss])
         #---------------------------------------------------------------------------------#
      }#end if (! this.harm && (! is.na(this.form)))
      #------------------------------------------------------------------------------------#
   }#end for (vv in sequence(nvargap))
   #=======================================================================================#
   #=======================================================================================#






   #=======================================================================================#
   #=======================================================================================#
   #     METHOD 4.  Not a single method... Each of the remaining variables will be eft  #
   #                using some alternative model or known relationships between the        #
   #                models.                                                                #
   #---------------------------------------------------------------------------------------#

      #----- Initialise the list that will hold the model fits. ---------------------------#
      filllist$m4 = list()
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #     Find specific humidity using the partial pressure of water vapour and          #
      # pressure.                                                                          #
      #------------------------------------------------------------------------------------#
      cat0(" - Find relative and specific humidity, and ensure that they are bounded.")
      eft$atm.pvap       = pmin(eft$atm.pvap,eslif(eft$atm.tmp))
      eft$atm.shv        = ( ep * eft$atm.pvap
                           / (eft$atm.prss + (1. - ep) * eft$atm.pvap) )
      eft$atm.rhv        = eft$atm.pvap / eslif(eft$atm.tmp)
      eft$atm.tdew       = tslif(pres=eft$atm.prss,hum=eft$atm.pvap,type.hum="pvap")
      eft$atm.tamb       = rep(NA,times=ntimes)
      eft$atm.tson       = rep(NA,times=ntimes)
      eft$gfflg.atm.tdew = eft$gfflg.atm.pvap
      eft$gfflg.atm.tamb = rep(NA,times=ntimes)
      eft$gfflg.atm.tson = rep(NA,times=ntimes)
      #------------------------------------------------------------------------------------#




      #====================================================================================#
      #====================================================================================#
      #     Load the rainfall from a daily precipitation file from INMET (or NCDC, but     #
      # previously converted to INMET format.  This is a lot less accurate than the        #
      # objective analysis and should be used sparingly.                                   #
      #------------------------------------------------------------------------------------#
      cat0( " + Fill in rainfall with daily estimate of rainfall.")

         #----- Load the target daily rainfall. -------------------------------------------#
         target       = read.table( file=raindayfile
                                  , nrows        = -1
                                  , header       = FALSE
                                  , skip         = 17
                                  , sep          = ";"
                                  , na.strings   = c("NA","-9999")
                                  , comment.char = ""
                                  ,colClasses    = rep("numeric",times=8)
                                  )#end read.table
         names(target) = c("wmo","day","month","year","hour","minu","seco","rain")
         target$today  = chron(paste(target$month,target$day,target$year,sep="/"))
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #      Flag days when rainfall is missing.  Because weather stations report rain  #
         # at 12 UTC, shift one day in case it is past the rainfall time.                  #
         #---------------------------------------------------------------------------------#
         miss         = is.na(eft$rain)
         today.rain   = chron(eft$today + as.numeric(eft$hour > 12))
         gap.days     = unique(today.rain[miss])
         table.days   = gap.days[   gap.days %in% target$today ]
         sample.days  = gap.days[ ! gap.days %in% target$today ]
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #     Find the probability of rain and mean rainfall rates.  The probability is   #
         # scaled so the sum of probabilities is 1.  This will preserve the ratios amongst #
         # different hours whilst reducing the number of failed "dice tossing".            #
         #---------------------------------------------------------------------------------#
         nvalid    = aggregate( x     = is.finite(eft$rain)
                              , by    = list(eft$month,eft$hour)
                              , FUN   = sum
                              , na.rm = TRUE)
         nreve     = aggregate( x     = eft$rain > 0
                              , by    = list(eft$month,eft$hour)
                              , FUN   = sum
                              , na.rm = TRUE)
         mon       = unique(nvalid[,1])
         nmon      = length(mon)
         hr        = unique(nvalid[,2])
         nhr       = length(hr)
         prob.rain = matrix(nreve[,3]/nvalid[,3],ncol=nmon,nrow=nhr,byrow=TRUE)
         tot.prob  = matrix(colSums(prob.rain,na.rm=TRUE)
                           ,ncol=nmon,nrow=nhr,byrow=TRUE)
         prob.rain = prob.rain/tot.prob
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #     Create weights that are relative to the mean precipitation rate.  This will #
         # make most of the rain fall when rainfall is the strongest in case more than one #
         # hour is selected to have rain.                                                  #
         #---------------------------------------------------------------------------------#
         mean.rain    = aggregate(x=eft$rain,by=list(eft$month,eft$hour)
                                 ,FUN=mean,na.rm=TRUE)
         weight       = matrix(mean.rain[,3],ncol=nmon,nrow=nhr,byrow=TRUE)
         #---------------------------------------------------------------------------------#


         #---------------------------------------------------------------------------------#
         #     Name the table rows and columns to make it easier to read and access.       #
         #---------------------------------------------------------------------------------#
         dimnames(prob.rain) = dimnames(weight) = list(hr,mon)
         #---------------------------------------------------------------------------------#


         #---------------------------------------------------------------------------------#
         #     Initialise the gap-filled rain with zeroes.  We will only fill in missing   #
         # days for which we have additional data.                                         #
         #---------------------------------------------------------------------------------#
         rough.rain          = rep(0,times=ntimes)
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #     Fill in the precipitation with the daily rainfall for days we have outside  #
         # information.  This must be done day by day, so we loop over each day that has   #
         # gaps.                                                                           #
         #---------------------------------------------------------------------------------#
         if (length(table.days) > 0){
            cat0("   - Fill in with external source of information:")
            for (this.day in table.days){
               cat0("     * Day: ",paste(chron(this.day)),".")


               #----- Select gap day. -----------------------------------------------------#
               eft.day  = today.rain == this.day
               miss.today  = miss & eft.day
               #---------------------------------------------------------------------------#


               #----- Find the total rain for today and the target total rainfall. --------#
               eft.rain = mean(eft$rain[eft.day]*day.sec,na.rm=TRUE)
               if (! is.finite(eft.rain)) eft.rain = 0.0
               target.day  = target$today == this.day
               target.rain = target$rain[target.day]
               #---------------------------------------------------------------------------#


               #----- Select the probabilities and weights for this day. ------------------#
               miss.month  = as.character(nummonths(this.day))
               miss.hour   = as.character(hours(eft$when[miss.today]))
               this.prob   = prob.rain[miss.hour,miss.month]
               this.weight = weight   [miss.hour,miss.month]
               nfill       = length(this.weight)
               #---------------------------------------------------------------------------#


               #---------------------------------------------------------------------------#
               #      Start the rain with zeroes.                                          #
               #---------------------------------------------------------------------------#
               rain.fill   = rep(0,times=nfill)
               #---------------------------------------------------------------------------#


               #---------------------------------------------------------------------------#
               #      We add rain only if the partial amount of precipitation is less than #
               # the target rainfall.                                                      #
               #---------------------------------------------------------------------------#
               if (eft.rain < target.rain){
                  #----- Find how much rain is "missing". ---------------------------------#
                  add.rain = (target.rain - eft.rain) / hr.sec
                  #------------------------------------------------------------------------#

                  #------------------------------------------------------------------------#
                  #    Toss the dice to pick up the hours that will contain rain.  We must #
                  # toss the dice until at least one of the hours is flagged as rain hour. #
                  #------------------------------------------------------------------------#
                  iterate = TRUE
                  while (iterate){
                     epsil   = .Machine$double.neg.eps
                     dice    = runif(n=nfill,min=epsil,max=1.-epsil)
                     iterate = ! any(dice <= this.prob)
                  }#end while
                  #------------------------------------------------------------------------#


                  #------------------------------------------------------------------------#
                  #    Spread the rain across the rain hours, scaling the amount by their  #
                  # mean rain rates.                                                       #
                  #------------------------------------------------------------------------#
                  yes.rain            = dice <= this.prob
                  rain.fill[yes.rain] = ( add.rain * this.weight[yes.rain]
                                        / sum(this.weight[yes.rain])       )
                  #------------------------------------------------------------------------#
               }#end if (eft.rain < target.rain)
               #---------------------------------------------------------------------------#

               #---------------------------------------------------------------------------#
               #      Copy the rainfall to the gap-eft vector.                             #
               #---------------------------------------------------------------------------#
               rough.rain[miss.today] = rain.fill
               #---------------------------------------------------------------------------#
            }#end for (this.day in table.days)
            #------------------------------------------------------------------------------#
         }#end if (length(table.days) > 0)
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #     Complete the time series with information sampled from the time series it-  #
         # self.                                                                           #
         #---------------------------------------------------------------------------------#
         if (length(sample.days) > 0 ){
            cat0("   - Fill with days randomly sampled from the time series.")
            for (this.day in sample.days){
               cat0("     * Day: ",paste(chron(this.day)),".")

               #---------------------------------------------------------------------------#
               #    List candidates, but exclude days that are not in the time series or   #
               # that have missing data themselves.                                        #
               #---------------------------------------------------------------------------#
               candidates = chron(this.day + c(-1,1)*seq(from=1,to=5,by=1))
               valid.days = candidates %in% today.rain & (! candidates %in% gap.days)
               candidates = candidates[valid.days]
               #---------------------------------------------------------------------------#



               #---------------------------------------------------------------------------#
               #    If we run out of candidates, then sample from any day from the same    #
               # month (including other years).                                            #
               #---------------------------------------------------------------------------#
               if (length(candidates) == 0){
                  sel        = nummonths(today.rain) == nummonths(this.day)
                  candidates = unique(today.rain[sel])
                  valid.days = candidates %in% today.rain & (! candidates %in% gap.days)
                  candidates = candidates[valid.days]
                  if (length(candidates) == 0){
                     stop ("Long period with no rainfall data!")
                  }#end if
               }#end if
               #---------------------------------------------------------------------------#


               #---------------------------------------------------------------------------#
               #    Pick up one day, and replace missing rainfall by the rainfall on this  #
               # day.                                                                      #
               #---------------------------------------------------------------------------#
               if (length(candidates) == 1){
                  replace.day = candidates
               }else{
                  replace.day = chron(sample(as.numeric(candidates),size=1))
               }#end if
               cat0("       > Replace rain from ",paste(chron(this.day))
                   ," by ",paste(chron(replace.day)),".")
               idx.miss           = sort(which(today.rain == this.day   ))
               idx.replace        = sort(which(today.rain == replace.day))
               rough.rain[idx.miss] = eft$rain[idx.replace]
               #---------------------------------------------------------------------------#
            }#end for
            #------------------------------------------------------------------------------#
         }#end if
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #      We replace only the times that we have no data...                          #
         #---------------------------------------------------------------------------------#
         eft$rain      [miss] = rough.rain[miss]
         eft$gfflg.rain[miss] = 4
         #---------------------------------------------------------------------------------#
      #====================================================================================#
      #====================================================================================#






      #====================================================================================#
      #====================================================================================#
      #      Fill in incoming long wave radiation using a very simple parametrisation,     #
      # very loosely based on the general equation from:                                   #
      #                                                                                    #
      # Monteith, J. L., Unsworth, M. H., 2008: Principles of environmental physics,       #
      #     Academic Press, 3rd. edition, London.                                          #
      #------------------------------------------------------------------------------------#
      cat0( " + Fill in incoming longwave radiation.")
      mmi.fit   = rlong.in.mmi.optim   ( datum     = eft
                                       , run.optim = FALSE
                                       , keep.day  = TRUE
                                       , keep.ngt  = TRUE
                                       , verbose   = FALSE
                                       )#end rlong.in.mmi.optim
      mmi.mat   = sapply(X=mmi.fit,FUN=sapply,c,simplify=TRUE)
      mmi.names = rownames(mmi.mat)
      summ.mmi  = data.frame( r.square   = unlist(mmi.mat[match("r.square"  ,mmi.names),])
                            , sigma      = unlist(mmi.mat[match("sigma"     ,mmi.names),])
                            , support    = unlist(mmi.mat[match("support"   ,mmi.names),])
                            , emiss.csky = sapply(mmi.mat[match("emiss.csky",mmi.names),]
                                                 ,mean,na.rm=TRUE)
                            , emiss.eff  = sapply(mmi.mat[match("emiss.eff" ,mmi.names),]
                                                 ,mean,na.rm=TRUE)
                            , f.cloud    = sapply(mmi.mat[match("f.cloud"   ,mmi.names),]
                                                 ,mean,na.rm=TRUE)
                            , dtemp      = ( mean(eft$atm.tmp,na.rm=TRUE)
                                           - sapply(mmi.mat[match("cld.tmp" ,mmi.names),]
                                                   ,mean,na.rm=TRUE)
                                           )#end
                            )#end summ.mmi
      m24.fit   = rlong.in.mmi.optim   ( datum     = eft
                                       , run.optim = TRUE
                                       , keep.day  = TRUE
                                       , keep.ngt  = TRUE
                                       , verbose   = FALSE
                                       )#end rlong.in.mmi.optim
      m24.mat   = sapply(X=m24.fit,FUN=sapply,c,simplify=TRUE)
      m24.names = rownames(m24.mat)
      summ.m24  = data.frame( r.square   = unlist(m24.mat[match("r.square"  ,m24.names),])
                            , sigma      = unlist(m24.mat[match("sigma"     ,m24.names),])
                            , support    = unlist(m24.mat[match("support"   ,m24.names),])
                            , emiss.csky = sapply(m24.mat[match("emiss.csky",m24.names),]
                                                 ,mean,na.rm=TRUE)
                            , emiss.eff  = sapply(m24.mat[match("emiss.eff" ,m24.names),]
                                                 ,mean,na.rm=TRUE)
                            , f.cloud    = sapply(m24.mat[match("f.cloud"   ,m24.names),]
                                                 ,mean,na.rm=TRUE)
                            , dtemp      = ( mean(eft$atm.tmp,na.rm=TRUE)
                                           - sapply(m24.mat[match("cld.tmp" ,m24.names),]
                                                   ,mean,na.rm=TRUE)
                                           )#end
                            )#end summ.m24

      use.fit                       = m24.fit$aml
      filllist$m4$rlong.in$default  = mmi.fit
      filllist$m4$rlong.in$summ.mmi = summ.mmi
      filllist$m4$rlong.in$tested   = m24.fit
      filllist$m4$rlong.in$summ.m24 = summ.m24
      filllist$m4$rlong.in$used     = use.fit
      fill                          = is.na(eft$rlong.in) & is.finite(use.fit$fitted.values)
      eft$rlong.in        [fill]    = use.fit$fitted.values[fill]
      eft$gfflg.rlong.in  [fill]    = 4
      eft$gferr.rlong.in  [fill]    = use.fit$sigma
      #====================================================================================#
      #====================================================================================#
   #=======================================================================================#
   #=======================================================================================#






   #=======================================================================================#
   #=======================================================================================#
   #     Additional variables.   These variables are found only after the gap-filling is   #
   # applied to all "input" variables.                                                     #
   #---------------------------------------------------------------------------------------#

      #------------------------------------------------------------------------------------#
      #     Albedo.  We need upwelling and downwelling shortwave radiation.                #
      #------------------------------------------------------------------------------------#
      cat( " + Fill in albedo.")
      eft$albedo                = NA * eft$rshort.in
      eft$albedo[eft$daytime]   = ( eft$rshort.out[eft$daytime] 
                                      / eft$rshort.in [eft$daytime] )
      malb                            = mean(eft$albedo,na.rm=TRUE)
      eft$albedo[eft$nighttime] = malb
      #------------------------------------------------------------------------------------#



      #----- Find the shortwave radiation breakdown. --------------------------------------#
      cat( " + Find the shortwave radiation partition into its components.")
      rsbreak                 = rshort.bdown( rad.in   = eft$rshort.in
                                            , atm.prss = ref.prss + 0. * eft$cosz
                                            , cosz     = eft$cosz
                                            , rad.type = "rshort"
                                            )#end rshort.bdown
      miss                    = is.na(eft$par.in)
      par.off                 = eft$par.in > eft$rshort.in
      par.off[is.na(par.off)] = FALSE
      miss                    = miss | par.off
      #----- Fill in missing values with the result from the model. -----------------------#
      eft$par.in      [miss] = rsbreak$par.full[miss]
      eft$gfflg.par.in[miss] = 4
      erte                   = ! miss & eft$daytime 
      eft$gferr.par.in[miss] = sqrt(mean((eft$par.in[erte] - rsbreak$par.full[erte])^2))
      #----- We find the break down one more time, using PAR as the reference value. ------#
      rsbreak                 = rshort.bdown( rad.in   = eft$par.in
                                            , atm.prss = ref.prss + 0. * eft$cosz
                                            , cosz     = eft$cosz
                                            , rad.type = "par")
      eft$par.beam = rsbreak$par.beam
      eft$par.diff = rsbreak$par.diff
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #     We fill in near infrared as the difference between SW and PAR.  For the break  #
      # down between direct and diffuse we should scale between the model and the          #
      # difference so everything adds up, even if the time series doesn't look that great. #
      #------------------------------------------------------------------------------------#
      eft$nir.in        = eft$rshort.in - eft$par.in
      eft$nir.beam      = 0 * eft$rshort.in
      eft$nir.diff      = 0 * eft$rshort.in
      sel               = eft$nir.in > 0 & rsbreak$nir.full > 0
      eft$nir.beam[sel] = rsbreak$nir.beam[sel] * eft$nir.in[sel] / rsbreak$nir.full[sel]
      eft$nir.diff[sel] = rsbreak$nir.diff[sel] * eft$nir.in[sel] / rsbreak$nir.full[sel]
      #------------------------------------------------------------------------------------#



      #----- Remove negative shortwave radiation ------------------------------------------#
      cat0(" - Remove negative values from shortwave.")
      eft$rshort.in [eft$rshort.in  < 0.0] = 0.0
      eft$rshort.out[eft$rshort.out < 0.0] = 0.0
      eft$par.in    [eft$par.in     < 0.0] = 0.0
      eft$par.beam  [eft$par.beam   < 0.0] = 0.0
      eft$par.diff  [eft$par.diff   < 0.0] = 0.0
      eft$nir.in    [eft$nir.in     < 0.0] = 0.0
      eft$nir.beam  [eft$nir.beam   < 0.0] = 0.0
      eft$nir.diff  [eft$nir.diff   < 0.0] = 0.0
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #     Recalculate the wind components.  The magnitude is better represented by the   #
      # interpolation of magnitude rather than by the wind components, which are used just #
      # to determine the direction (and we don't really care about the direction in ED).   #
      #------------------------------------------------------------------------------------#
      cat0(" - Make wind speed positive definite and determine components.")
      eft$ustar    = pmax(eft$ustar   ,ustar.min   )
      eft$atm.vels = pmax(eft$atm.vels,atm.vels.min)
      eft$trigo    = atan2(eft$atm.vspd,eft$atm.uspd)
      eft$atm.uspd = eft$atm.vels * cos(eft$trigo)
      eft$atm.vspd = eft$atm.vels * sin(eft$trigo)
      eft$atm.vdir = (270. - eft$trigo * 180. / pi) %% 360.
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #     Find air density.                                                              #
      #------------------------------------------------------------------------------------#
      cat0(" - Find air density.")
      eft$atm.rhos = idealdenssh(eft$atm.prss,eft$atm.tmp,eft$atm.shv)
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #     Find the vapour pressure deficit.                                              #
      #------------------------------------------------------------------------------------#
      cat0(" - Find vapour pressure deficit.")
      eft$atm.vpd      = eslif(eft$atm.tmp) - eft$atm.pvap
      sat              = eft$atm.vpd < 0
      eft$atm.vpd[sat] = 0.
      #------------------------------------------------------------------------------------#
   #=======================================================================================#
   #=======================================================================================#



   #---------------------------------------------------------------------------------------#
   #     Save the gap-filled dataset (except NEE) so we don't load them again.             #
   #---------------------------------------------------------------------------------------#
   cat0(" + Save gap-filled (except NEE) data to file ",basename(session.metfill),".")
   save(eft,filllist,file=session.metfill)
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #    Copy eft to metfill.  From this point on, the data will not be called eft, but     #
   # metfill because we always use eft as the structure name that is saved to the RData    #
   # regardless of the level.                                                              #
   #---------------------------------------------------------------------------------------#
   metfill = eft
   #---------------------------------------------------------------------------------------#
}#end if
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
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#    Zoom in on the gap.                                                                   #
#------------------------------------------------------------------------------------------#
if (plotgapfill){
   cat0( " + Plot some zoomed-in time series of gap-filled data:")
   for (vv in sequence(nvargap)){
      #----- Load the variable properties. ------------------------------------------------#
      thisvar   = vargap$vari[vv]
      thisdesc  = vargap$desc[vv]
      thisunit  = vargap$unit[vv]
      thisadd0  = vargap$add0[vv]
      thismult  = vargap$mult[vv]
      thisharm  = vargap$harm[vv]
      thisform  = vargap$form[vv]
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #       Apply transformation if sought.                                              #
      #------------------------------------------------------------------------------------#
      if (is.na(thisadd0)) add0=0. else add0=eval(parse(text=thisadd0))
      if (is.na(thismult)) mult=1. else mult=eval(parse(text=thismult))
      #------------------------------------------------------------------------------------#



      cat0("   - Variable:",thisdesc,".")


      #----- Find every missing value. ----------------------------------------------------#
      miss          = is.na(clean[[thisvar]])
      nmiss         = length(clean[[thisvar]][miss])
      nfill         = is.finite(metfill[[thisvar]]) & miss
      whenmiss      = as.numeric(metfill$when[miss])

      #----- If there was no gap or if the time series is empty, skip this variable. ------#
      if (nmiss > 0 & nmiss < ntimes & nfill > 0){
         gapidx        = cumsum(miss & (! c(FALSE,miss[-ntimes])))
         gapidx[!miss] = NA
         gapa          = unique( unlist( tapply( X     = gapidx
                                               , INDEX = gapidx
                                               , FUN   = match
                                               , table = gapidx
                                               )  )  )
         gapz          = ( ntimes + 1
                         - unique( unlist( tapply( X     = gapidx
                                                 , INDEX = gapidx
                                                 , FUN   = match
                                                 , table = rev(gapidx)
                                                 ) ) ) )


         #---------------------------------------------------------------------------------#
         #      Check whether there is a directory or not. In case it doesn't, make one.   #
         #---------------------------------------------------------------------------------#
         outgap = file.path(outroot,"gaps")
         outdir = file.path(outgap,thisvar)
         if (! file.exists(outgap)) dir.create(outgap)
         if (! file.exists(outdir)) dir.create(outdir)
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #    Make the list of gaps, by combining continuous blocks of missing data.  In   #
         # case they are too frequent, merge gaps.                                         #
         #---------------------------------------------------------------------------------#
         ngaps    = length(gapa)
         mergegap = ngaps > maxgap
         if (mergegap) cat0("   * Combine multiple gaps into groups:")
         itgap = 0
         while (mergegap){
            #----- Update the iterative loop counter. -------------------------------------#
            itgap    = itgap + 1
            deltagap = c(Inf,gapa[-1]-gapz[-ngaps])
            #------------------------------------------------------------------------------#


            #------------------------------------------------------------------------------#
            #     Determine which gaps can be combined.  We chose gaps that are next to    #
            # each other.                                                                  #
            #------------------------------------------------------------------------------#
            merge.prev  = deltagap <= 2^(itgap-1)
            which.merge = which(merge.prev)
            which.stay  = which(! merge.prev)
            #----- If we found gaps that we can merge, merge them, ------------------------#
            if (any(merge.prev)){
               #----- Replace the gap index. ----------------------------------------------#
               for (w in which.merge){
                  which.bind  = max(which.stay[which.stay < w])
                  sel         = is.finite(gapidx) & gapidx == w
                  gapidx[sel] = which.bind
               }#end for
               #---------------------------------------------------------------------------#


               #---------------------------------------------------------------------------#
               #      Collapse the gap indices, so they are in order.                      #
               #---------------------------------------------------------------------------#
               gap         = ! is.na(gapidx)
               gapidx[gap] = match(gapidx[gap],unique(gapidx[gap]))
               #----- Update the beginning and the end of the gaps. -----------------------#
               gapa        = unique( unlist( tapply( X     = gapidx
                                                   , INDEX = gapidx
                                                   , FUN   = match
                                                   , table = gapidx
                                                   )  )  )
               gapz        = ( ntimes + 1
                             - unique( unlist( tapply( X     = gapidx
                                                     , INDEX = gapidx
                                                     , FUN   = match
                                                     , table = rev(gapidx)
                                                     ) ) ) )
               #---------------------------------------------------------------------------#
            }#end if
            #------------------------------------------------------------------------------#



            #------------------------------------------------------------------------------#
            #      Count the number of gaps and check whether we must collapse more.       #
            #------------------------------------------------------------------------------#
            ngaps       = length(gapa)
            mergegap    = ngaps > maxgap
            cat0("     > IT = ",itgap,"  NGAPS = ",ngaps,".")
            #------------------------------------------------------------------------------#
         }#end while
         #---------------------------------------------------------------------------------#


         #---------------------------------------------------------------------------------#
         #     Loop over every gap.                                                        #
         #---------------------------------------------------------------------------------#
         cat0("   * Plot time series.")
         for (g in sequence(ngaps)){
            cgap  = sprintf("%4.4i",g)

            cat0("     > Gap number: ",g,".")
            whena    = clean$when[gapa[g]]
            whenz    = clean$when[gapz[g]]
            whena    = chron(paste(nummonths(whena),numdays(whena),numyears(whena),sep="/")
                            ,paste(hours(whena),minutes(whena),seconds(whena),sep=":"))
            whenz    = chron(paste(nummonths(whenz),numdays(whenz),numyears(whenz),sep="/")
                            ,paste(hours(whenz),minutes(whenz),seconds(whenz),sep=":"))

            nw       = length(metfill$when)
            whenplta = chron(max(metfill$when[1],whena-2.))
            whenpltz = chron(min(metfill$when[nw],whenz+2.))
            whenplta = chron(paste(nummonths(whenplta)
                                  ,numdays(whenplta)
                                  ,numyears(whenplta)
                                  ,sep="/")
                            ,paste(hours(whenplta)
                                  ,minutes(whenplta)
                                  ,seconds(whenplta)
                                  ,sep=":"))
            whenpltz = chron(paste(nummonths(whenpltz)
                                  ,numdays(whenpltz)
                                  ,numyears(whenpltz)
                                  ,sep="/")
                            ,paste(hours(whenpltz)
                                  ,minutes(whenpltz)
                                  ,seconds(whenpltz)
                                  ,sep=":"))
            #------------------------------------------------------------------------------#



            #----- Define a suitable scale for those time series that uses when... --------#
            whenplot = pretty.time(c(whenplta,whenpltz),n=8)
            xlimit   = chron(range(whenplot$levels))
            #------------------------------------------------------------------------------#

            #----- Select the period to plot. ---------------------------------------------#
            sel        = metfill$when %wr% c(xlimit,xlimit[2])
            whenvar    = metfill$when[sel]
            cleanvar   = add0 + mult *clean  [[thisvar]][sel]
            metfillvar = add0 + mult *metfill[[thisvar]][sel]
            #------------------------------------------------------------------------------#



            #----- Define a suitable scale for the y axis... ------------------------------#
            ylimit = pretty.xylim(u=c(cleanvar,metfillvar),fracexp=0.0)
            #------------------------------------------------------------------------------#

            #----- Define plot annotation. ------------------------------------------------#
            letitre = paste0("Gap filling - ",longname,"\n","Gap period number: ",g,".")
            lex     = desc.unit(desc="Time"  ,unit=untab$utc        )
            ley     = desc.unit(desc=thisdesc,unit=untab[[thisunit]])
            #------------------------------------------------------------------------------#



            #------------------------------------------------------------------------------#
            #     Loop over formats.                                                       #
            #------------------------------------------------------------------------------#
            for (o in sequence(nout)){
               fichier = file.path(outdir,paste0(thisvar,"-",cgap,".",outform[o]))
               dummy   = open.file( fichier = fichier
                                  , outform = outform[o]
                                  , size    = ex.size
                                  , ptsz    = ptsz
                                  , depth   = depth
                                  )#end open.file
               #---------------------------------------------------------------------------#

               #----- Split device into two (legend, main plot). --------------------------#
               par(par.user)
               layout(mat=rbind(2,1),heights=c(1.-f.leg,f.leg))
               #---------------------------------------------------------------------------#


               #----- Legend. -------------------------------------------------------------#
               par(mar=c(0.1,0.1,0.1,0.1))
               plot.new()
               plot.window(c(0,1),c(0,1))
               legend( x      = "bottom"
                     , inset  = 0.0
                     , legend = c("Raw data","Gap filled")
                     , bg     = "transparent"
                     , col    = c("dodgerblue3","orangered")
                     , lwd    = c(2.0,2.0)
                     , lty    = c("solid","solid")
                     , xpd    = TRUE
                     , bty    = "n"
                     )#end legend
               #---------------------------------------------------------------------------#




               #----- Make the plot window. -----------------------------------------------#
               par(par.user)
               par(mar=c(4.1,4.6,3.1,0.6))
               plot.new()
               plot.window(xlim=xlimit,ylim=ylimit)
               title(main=letitre,xlab=lex,ylab=ley)
               #----- Plot time axis. -----------------------------------------------------#
               axis(side=1,at=whenplot$levels,labels=whenplot$labels,padj=whenplot$padj)
               axis(side=2,las=1)

               #---------------------------------------------------------------------------#
               #     Make a grid that matches the time, and add beginning and end of gap.  #
               #---------------------------------------------------------------------------#
               if (plotgrid){ 
                  abline(v=whenplot$levels,h=axTicks(side=2)
                        ,col="lightgrey",lty="dotted",lwd=0.25)
                  abline(v=c(whena,whenz),col="chartreuse",lty="dotdash",lwd=2)
               }#end if
               #----- Plot the original and the gapfilled time series. --------------------#
               lines(x = whenvar , y = metfillvar, col="orangered"  ,lty="solid"  ,lwd=2.0)
               lines(x = whenvar , y = cleanvar  , col="dodgerblue3",lty="solid"  ,lwd=2.0)
               #----- Lastly, add the box. ------------------------------------------------#
               box()
               #---------------------------------------------------------------------------# 




               #---------------------------------------------------------------------------# 
               #   Close the device, and if this is a X11 session, wait until the user has #
               # clicked on the figure.                                                    #
               #---------------------------------------------------------------------------# 
               dummy = close.plot(outform=outform[o])
               #---------------------------------------------------------------------------# 
            }#end for (o in sequence(nout))
            #------------------------------------------------------------------------------#
         }#end for (g in sequence(gap))
         #---------------------------------------------------------------------------------#
      }#end if (miss > 0)
      #------------------------------------------------------------------------------------#
   }#end for (vv in sequence(ngapvar))
   #---------------------------------------------------------------------------------------#
}#end if (plotgapfill)
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#

