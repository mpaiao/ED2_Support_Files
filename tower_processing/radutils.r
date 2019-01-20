#------------------------------------------------------------------------------------------#
# General Radiation constants.                                                             #
#------------------------------------------------------------------------------------------#
pio180        <<- pi/ 180.
capri         <<- -23.44 * pio180 # Tropic of Capricornium latitude             [      rad]
shsummer      <<- -10             # Day of year of S.Hemisphere summer solstice [      day]
solar         <<- 1.3533e3        # Solar constant                              [     W/m2]
prefsea       <<- 101325.         # Reference sea level pressure                [       Pa]
twothirds     <<- 2./3.          # 2/3                                          [      ---]
cosz.min      <<- 0.03            # Minimum cosine of zenith angle
cosz.highsun  <<- cos(84*pi/180)  # Zenith angle to not be called sunrise or sunset
cosz.twilight <<- cos(96*pi/180)  # Cosine of the end of civil twilight
fvis.beam.def <<- 0.43
fnir.beam.def <<- 1.0 - fvis.beam.def
fvis.diff.def <<- 0.52
fnir.diff.def <<- 1.0 - fvis.diff.def
#------------------------------------------------------------------------------------------#



#==========================================================================================#
#==========================================================================================#
#    This function finds the cosine of the zenith angle either for the right instant, or   #
# to the interval between two consecutive times.                                           #
#                                                                                          #
# Input variables:                                                                         #
#    - lon       - longitude of the point.  Mandatory, one value only.                     #
#    - lat       - latitude of the point.  Mandatory, one value only.                      #
#    - when      - time.  Mandatory, one point only or a vector.                           #
#    - ed21      - I shall use ED-2.1 method (TRUE/FALSE).  Default is TRUE                #
#    - zeronight - The cosine of zenith angle shall be set to zero at night.               #
#                  Default is FALSE                                                        #
#    - meanval   - I shall find the mean cosine of the integration time.  The beginning    #
#                  and the end are given by variable imetavg.  Default is FALSE.  In case  #
#                  it is TRUE but "when" has just one point, this flag will be ignored and #
#                  it will be solved as instantaneous.                                     #
#    - imetavg   - Which kind of time average was used?                                    #
#                  1 - averages ending at the reference time;                              #
#                  2 - averages beginning at the reference time;                           #
#                  3 - averages centred at the reference time.                             #
#    - nmean     - Number of intermediate points for the average                           #
#                                                                                          #
# The output is going to be a list with the following values:                              #
#    - cosz      - Cosine of zenith angle                                                  #
#    - zen       - The zenith angle in degrees                                             #
#    - height    - The sun height in degrees                                               #
#    - declin    - Declination in degrees                                                  #
#    - day       - Daytime (Sun above horizon)                                             #
#    - night     - Night time (Sun below 6 degrees below the horizon                       #
#    (N.B. When both day and night are false, we consider it twilight.                     #
#------------------------------------------------------------------------------------------#
ed.zen = function (lon,lat,when,ed21=TRUE,zeronight=FALSE,meanval=FALSE,imetavg=1
                  ,nmean=120,...){
   #------ Constants. ---------------------------------------------------------------------#
   dcoeff   = c( 0.006918, -0.399912,  0.070257, -0.006758,  0.000907, -0.002697,  0.001480)
   #---------------------------------------------------------------------------------------#


   #------ Find the number of elements. ---------------------------------------------------#
   ntimes  = length(when)
   if ((! meanval) | ntimes == 1) nmean = 1
   #---------------------------------------------------------------------------------------#


   #------ Make matrix of times to make the results time averages if needed be. -----------#
   if (nmean > 1){
      #------------------------------------------------------------------------------------#
      #     The minimum difference is safer than the mean in case the time series has      #
      # gaps.                                                                              #
      #------------------------------------------------------------------------------------#
      dwhen = diff(as.numeric(when))
      sel   = is.finite(dwhen)
      dwhen = dwhen[sel]
      dwhen = min(dwhen[dwhen > 0])
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #    Decide the beginning and ending times depending on imetavg.                     #
      #------------------------------------------------------------------------------------#
      if (imetavg == 1){
         #----- Averages ending at the reference time. ------------------------------------#
         na = 1 - nmean
         nz = 0
      }else if (imetavg == 2){
         #----- Averages starting at the reference time. ----------------------------------#
         na = 0
         nz = nmean - 1
      }else if (imetavg == 3){
         #---------------------------------------------------------------------------------#
         #     Averages centered at the reference time. The initial and ending times na    #
         # and nz will be slightly different depending on whether the number of mean       #
         # points is odd or even.                                                          #
         #---------------------------------------------------------------------------------#
         nz = floor(nmean/2) + 0.5 * ((nmean %% 2) - 1.0)
         na = - nz
      }else{
         print(paste(" ---> In function ed.zen: imetavg =",imetavg,".",sep=""))
         stop ("Invalid imetavg, it must be 1, 2, or 3!")
      }#end if
      #------------------------------------------------------------------------------------#



      #----- Averages ending at the reference time. ---------------------------------------#
      dtidx = seq(from=na,to=nz,by=1) / (nz - na + 1)
      WHEN  = chron( matrix(as.numeric(when),ncol=nmean,nrow=ntimes)
                   + matrix(dtidx,ncol=nmean,nrow=ntimes,byrow=TRUE) * dwhen)
      #------------------------------------------------------------------------------------#
   }else{
      #----- Single time, use only the instantaneous value. -------------------------------#
      WHEN  = matrix(as.numeric(when),ncol=nmean,nrow=ntimes)
   }#end if
   empty = as.numeric(WHEN) * NA
   #---------------------------------------------------------------------------------------#



   #------ Find the day of year, list of leap year times, and sun hour. -------------------#
   doy     = matrix(dayofyear(when)           ,ncol=nmean,nrow=ntimes)
   leap    = matrix(is.leap  (when)           ,ncol=nmean,nrow=ntimes)
   fracday = matrix(hms2frac (as.vector(WHEN)),ncol=nmean,nrow=ntimes)
   sunhr   = (fracday * day.hr + lon / 15. + day.hr) %% day.hr
   #---------------------------------------------------------------------------------------#



   #------ Find the hour angle and its cosine. --------------------------------------------#
   hrangle = 15 * (sunhr - 12) * pio180
   chra    = cos(hrangle)
   #---------------------------------------------------------------------------------------#



   #------ Find the declination
   if (ed21){
      doyfun = empty
      doyfun[!leap] = 2 * pi * (doy[!leap] - shsummer) / 365.
      doyfun[ leap] = 2 * pi * (doy[ leap] - shsummer) / 366.

      declin = capri * cos(doyfun)
   }else{
      doyfun = empty
      doyfun[!leap] = 2 * pi * (doy[!leap] - 1) / 365.
      doyfun[ leap] = 2 * pi * (doy[ leap] - 1) / 366.

      declin = ( dcoeff[1]
               + dcoeff[2] * cos(1.*doyfun) + dcoeff[3] * sin(1.*doyfun)
               + dcoeff[4] * cos(2.*doyfun) + dcoeff[5] * sin(2.*doyfun)
               + dcoeff[6] * cos(3.*doyfun) + dcoeff[7] * sin(3.*doyfun) )
   }#end if
   #---------------------------------------------------------------------------------------#

   #------ Find the cosine and sine of latitude and declination. --------------------------#
   clat = cos(pio180*lat)
   slat = sin(pio180*lat)
   cdec = matrix(cos(declin),ncol=nmean,nrow=ntimes)
   sdec = matrix(sin(declin),ncol=nmean,nrow=ntimes)

   #------ Find the cosine of the zenith angle, the zenith angle, and day/night flag. -----#
   cosz   = rowMeans(slat * sdec + clat * cdec * chra,...)
   zen    = acos(cosz) / pio180
   hgt    = 90. - zen
   declin = declin / pio180
   night  = cosz <  cosz.twilight
   day    = cosz >= cosz.min

   if (zeronight){
      cosz[night] =  0.
      hgt [night] =  0.
      zen [night] = 90.
   }#end if

   ans = list(cosz=cosz,zen=zen,hgt=hgt,declin=declin,day=day,night=night)
   return(ans)
}#end function ed.zen
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#      This subroutine computes the split between direct and diffuse radiation, and        #
# between visible and near-infrared radiation using the method suggested by:               #
#                                                                                          #
# Weiss, A., J. M. Norman, 1985: Partitioning solar radiation into direct and diffuse,     #
#     visible and near-infrared components.  Agric. For. Meteorol., 34, 205-213. (WN85)    #
#                                                                                          #
# Input variables:                                                                         #
#                                                                                          #
#    * rad.in   - The incoming radiation at surface, in W/m2.  This can be either PAR,     #
#                 NIR, or the total shortwave radiation, but it must be in W/m2 in any of  #
#                 the cases.                                                               #
#    * atm.prss - The atmospheric pressure at the surface, in Pa.  An actual measurement   #
#                 is better, but if you don't have one just use some typical value given   #
#                 the place altitude (higher elevation sites get more radiation).          #
#    * cosz     - The cosine of zenith angle.  This can be estimated using function ed.zen #
#                 in file zen.r
#    * rad.type - The type of radiation provided in rad.in.  Default is total shortwave    #
#                 radiation, but the function also accepts PAR or NIR.  The value is case  #
#                 insensitive and only the first letter is checked.  "p" means PAR, "n"    #
#                 means NIR, and any other letter will be assumed shortwave.               #
#------------------------------------------------------------------------------------------#
rshort.bdown = function(rad.in,atm.prss,cosz,rad.type="rshort"){
   #---------------------------------------------------------------------------------------#
   #    Local constants.                                                                   #
   #---------------------------------------------------------------------------------------#
   #----- Extinction coefficient. (equations 1 and 4 of WN85) -----------------------------#
   par.beam.expext  = -0.185
   nir.beam.expext  = -0.060
   #----- This is the typical conversion of diffuse radiation in sunny days. --------------#
   par2diff.sun = 0.400
   nir2diff.sun = 0.600
   #----- Coefficients for various equations in WN85. -------------------------------------#
   wn85.06 = c( -1.1950, 0.4459, -0.0345 )
   wn85.11 = c(    0.90, 0.70  )
   wn85.12 = c(    0.88, 0.68  )
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Make rad type case insensitive, and retain only the first letter.                 #
   #---------------------------------------------------------------------------------------#
   rty = substring(tolower(rad.type),1,1)
   #---------------------------------------------------------------------------------------#


   #------ Initialise the radiation with NAs. ---------------------------------------------#
   par.beam    = NA * rad.in
   nir.beam    = NA * rad.in
   par.diff    = NA * rad.in
   nir.diff    = NA * rad.in
   par.full    = NA * rad.in
   nir.full    = NA * rad.in
   rshort.beam = NA * rad.in
   rshort.diff = NA * rad.in
   rshort.full = NA * rad.in
   par.max     = NA * rad.in
   nir.max     = NA * rad.in
   rshort.max  = NA * rad.in


   #------ Make day and night flags. ------------------------------------------------------#
   ntimes = length(cosz)
   night  = cosz <= cosz.min
   day    = ! night 
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     First thing to check is whether this is daytime or "night-time".  If the zenith   #
   # angle is too close to horizon, we assume it's dawn/dusk and all radiation goes to     #
   # diffuse.                                                                              #
   #---------------------------------------------------------------------------------------#
   par.beam    [night] = 0.0
   nir.beam    [night] = 0.0
   if (rty == "p"){
      par.diff    [night] = rad.in[night]
      nir.diff    [night] = fnir.diff.def * rad.in[night] / fvis.diff.def
   }else if(rty == "n"){
      par.diff    [night] = fvis.diff.def * rad.in[night] / fnir.diff.def
      nir.diff    [night] = rad.in[night]
   }else{
      par.diff    [night] = fvis.diff.def * rad.in[night]
      nir.diff    [night] = fnir.diff.def * rad.in[night]
   }#end if
   par.full    [night] = par.beam   [night] + par.diff   [night]
   nir.full    [night] = nir.beam   [night] + nir.diff   [night]
   rshort.beam [night] = par.beam   [night] + nir.beam   [night]
   rshort.diff [night] = par.diff   [night] + nir.diff   [night]
   rshort.full [night] = rshort.beam[night] + rshort.diff[night]
   par.max     [night] = 0.0
   nir.max     [night] = 0.0
   rshort.max  [night] = 0.0
   #---------------------------------------------------------------------------------------#



   #----- Save 1/cos(zen), which is the secant.  We will use this several times. ----------#
   secz      = 1. / cosz[day]
   log10secz = log10(secz)
   #---------------------------------------------------------------------------------------#


   #----- Total radiation at the top [  W/m2], using ED defaults. -------------------------#
   par.beam.top = fvis.beam.def * solar
   nir.beam.top = fnir.beam.def * solar
   #---------------------------------------------------------------------------------------#

   #---------------------------------------------------------------------------------------#
   #    Find the potential PAR components (beam, diffuse, total), using equations 1, 3,    #
   # and 9 of WN85.                                                                        #
   #---------------------------------------------------------------------------------------#
   par.beam.pot = ( par.beam.top
                  * exp ( par.beam.expext * (atm.prss[day] / prefsea) * secz) * cosz[day])
   par.diff.pot = par2diff.sun * (par.beam.top - par.beam.pot) * cosz[day]
   par.full.pot = par.beam.pot + par.diff.pot
   #------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Find the NIR absorption of 10 mm of precipitable water, using WN85 equation 6.    #
   #---------------------------------------------------------------------------------------#
   w10 = solar * 10 ** ((wn85.06[1]) + log10secz * (wn85.06[2] + wn85.06[3] * log10secz))
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Find the potential direct and diffuse near-infrared radiation, using equations    #
   # 4, 5, and 10 of WN85.                                                                 #
   #---------------------------------------------------------------------------------------#
   nir.beam.pot = ( ( nir.beam.top
                    * exp ( nir.beam.expext * (atm.prss[day] / prefsea) * secz) - w10 )
                  * cosz[day] )
   nir.diff.pot = nir2diff.sun * ( nir.beam.top - nir.beam.pot - w10 ) * cosz[day]
   nir.full.pot = nir.beam.pot + nir.diff.pot
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Total maximum radiation.                                                          #
   #---------------------------------------------------------------------------------------#
   par.max   [day] = par.full.pot
   nir.max   [day] = nir.full.pot
   rshort.max[day] = par.full.pot + nir.full.pot
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Find the actual total for PAR and NIR, using equations 7 and 8.                   #
   #---------------------------------------------------------------------------------------#
   if (rty == "p"){
      ratio      = rad.in[day] / par.full.pot
   }else if (rty == "n"){
      ratio      = rad.in[day] / nir.full.pot
   }else{
      ratio      = rad.in[day] / (par.full.pot + nir.full.pot)
   }#end if
   par.full[day] = ratio * par.full.pot
   nir.full[day] = ratio * nir.full.pot
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Find the fraction of PAR and NIR that stays as beam, using equations 11 and 12    #
   # of WN85.                                                                              #
   #---------------------------------------------------------------------------------------#
   #----- Make sure that the ratio is bounded. --------------------------------------------#
   aux.par  = pmin(wn85.11[1],pmax(0.,ratio))
   aux.nir  = pmin(wn85.12[1],pmax(0.,ratio))

   fvis.beam.act = ( par.beam.pot 
                   * (1. - ((wn85.11[1] - aux.par)/wn85.11[2]) ^ twothirds)
                   / par.full.pot )
   fvis.beam.act = pmin(1.,pmax(0.,fvis.beam.act))

   fnir.beam.act = ( nir.beam.pot 
                   * (1. - ((wn85.12[1] - aux.nir)/wn85.12[2]) ^ twothirds)
                   / nir.full.pot )
   fnir.beam.act = pmin(1.,pmax(0.,fvis.beam.act))

   fvis.diff.act = 1. - fvis.beam.act
   fnir.diff.act = 1. - fnir.beam.act
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Find the radiation components.                                                    #
   #---------------------------------------------------------------------------------------#
   par.beam    [day] = fvis.beam.act * par.full[day]
   par.diff    [day] = fvis.diff.act * par.full[day]
   nir.beam    [day] = fnir.beam.act * nir.full[day]
   nir.diff    [day] = fnir.diff.act * nir.full[day]
   rshort.beam [day] = par.beam   [day] + nir.beam   [day]
   rshort.diff [day] = par.diff   [day] + nir.diff   [day]
   rshort.full [day] = rshort.beam[day] + rshort.diff[day]
   #---------------------------------------------------------------------------------------#
   rshort.bdown = list( par.beam    = par.beam
                      , par.diff    = par.diff
                      , par.full    = par.full
                      , nir.beam    = nir.beam
                      , nir.diff    = nir.diff
                      , nir.full    = nir.full
                      , rshort.beam = rshort.beam
                      , rshort.diff = rshort.diff
                      , rshort.full = rshort.full
                      , par.max     = par.max
                      , nir.max     = nir.max
                      , rshort.max  = rshort.max
                      )#end list
   return(rshort.bdown)
}#end function rshort.bdown
#==========================================================================================#
#==========================================================================================#
