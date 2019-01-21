#==========================================================================================#
#==========================================================================================#
#       Function that converts DBH to height.                                              #
#------------------------------------------------------------------------------------------#
dbh2h <<- function(dbh,ipft,use.crit=TRUE){
   #----- Make sure that the PFT variable has the same length as dbh. ---------------------#
   if (length(ipft) == 1){
     zpft = rep(ipft,times=length(dbh))
   }else{
     zpft = ipft
   }#end if
   #---------------------------------------------------------------------------------------#



   #----- Limit dbh to dbh.crit. ----------------------------------------------------------#
   if (use.crit){
      dbhuse = pmin(pft$dbh.crit[zpft],dbh) + 0. * dbh
   }else{
      dbhuse = dbh
   }#end if (use.crit)
   #---------------------------------------------------------------------------------------#


   #----- Find out which functional form to use. ------------------------------------------#
   tropo         = pft$tropical[zpft] & iallom %in% c(0,1,3)
   tropn         = pft$tropical[zpft] & iallom %in% c(2)
   tempe         = ! pft$tropical[zpft]
   #---------------------------------------------------------------------------------------#


   #----- Handy aliases. ------------------------------------------------------------------#
   hgt.ref = pft$hgt.ref[zpft]
   b1Ht    = pft$b1Ht   [zpft]
   b2Ht    = pft$b2Ht   [zpft]
   #---------------------------------------------------------------------------------------#



   #----- Find height. --------------------------------------------------------------------#
   h = ifelse( test = tempe
             , yes  = hgt.ref + b1Ht * (1. - exp(b2Ht * dbhuse) )
             , no   = ifelse( test = tropn
                            , yes  = hgt.ref * (1.-exp(-b1Ht * dbhuse^b2Ht ) )
                            , no   = exp(b1Ht + b2Ht * log(dbhuse) )
                            )#end ifelse
             )#end ifelse
   #---------------------------------------------------------------------------------------#

   return(h)
}#end function dbh2h
#==========================================================================================!
#==========================================================================================!






#==========================================================================================#
#==========================================================================================#
#       Function that converts size (DBH and Height) to biomass of living tissues.         #
#------------------------------------------------------------------------------------------#
size2ba <<- function(dbh,hgt,ipft,use.crit=TRUE){
   #----- Make sure that the PFT variable has the same length as dbh. ---------------------#
   if (length(ipft) == 1){
     zpft = rep(ipft,times=length(dbh))
   }else{
     zpft = ipft
   }#end if
   #---------------------------------------------------------------------------------------#


   #----- Limit dbh to dbh.crit. ----------------------------------------------------------#
   if (use.crit){
      dbhuse  = pmin(dbh,pft$dbh.crit[zpft]) + 0. * dbh
   }else{
      dbhuse  = dbh
   }#end if (use.crit)
   #---------------------------------------------------------------------------------------#


   #----- Decide which variable to use as dependent variable (DBH or DBH^2*Hgt). ----------#
   size     = ifelse( test = pft$tropical[zpft] & (! pft$liana[zpft]) & (iallom %in% 3)
                    , yes  = dbhuse * dbhuse * hgt
                    , no   = dbhuse
                    )#end ifelse
   #---------------------------------------------------------------------------------------#



   #----- Find leaf biomass. --------------------------------------------------------------#
   bleaf  = pft$b1Bl[zpft] / C2B * size ^ pft$b2Bl[zpft]
   broot  = pft$qroot[zpft] * bleaf
   bsapw  = pft$qsw  [zpft] * hgt * bleaf
   bbark  = pft$qbark[zpft] * hgt * bleaf
   balive = bleaf + broot + bsapw + bbark
   #---------------------------------------------------------------------------------------#


   return(balive)
}# end function size2ba
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#       Function that converts size (DBH and Height) to biomass of structural tissues.     #
#------------------------------------------------------------------------------------------#
size2bd <<- function(dbh,hgt,ipft){
   #----- Make sure that the PFT variable has the same length as dbh. ---------------------#
   if (length(ipft) == 1){
     zpft = rep(ipft,times=length(dbh))
   }else{
     zpft = ipft
   }#end if
   #---------------------------------------------------------------------------------------#


   #----- Decide which variable to use as dependent variable (DBH or DBH^2*Hgt). ----------#
   size     = ifelse( test = pft$tropical[zpft] & (! pft$liana[zpft]) & (iallom %in% 3)
                    , yes  = dbh * dbh * dbh2h(dbh,ipft=zpft)
                    , no   = dbh
                    )#end ifelse
   #---------------------------------------------------------------------------------------#



   #----- Select allometric parameters based on the size. ---------------------------------#
   bdead = ifelse( test = dbh < pft$dbh.crit[zpft]
                 , yes  = pft$b1Bs.small[zpft] / C2B * size ^ pft$b2Bs.small[zpft]
                 , no   = pft$b1Bs.large[zpft] / C2B * size ^ pft$b2Bs.large[zpft]
                 )#end ifelse
   #---------------------------------------------------------------------------------------#

   return(bdead)
}# end function size2bd
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#       Function that converts size (DBH and Height) to individual leaf area index.        #
#------------------------------------------------------------------------------------------#
size2lai <<- function(dbh,hgt,nplant,ipft,use.crit=TRUE){
   #----- Make sure that the PFT variable has the same length as dbh. ---------------------#
   if (length(ipft) == 1){
     zpft = rep(ipft,times=length(dbh))
   }else{
     zpft = ipft
   }#end if
   #---------------------------------------------------------------------------------------#


   #----- Limit dbh to dbh.crit. ----------------------------------------------------------#
   if (use.crit){
      dbhuse  = pmin(dbh,pft$dbh.crit[zpft]) + 0. * dbh
   }else{
      dbhuse  = dbh
   }#end if (use.crit)
   #---------------------------------------------------------------------------------------#


   #----- Decide which variable to use as dependent variable (DBH or DBH^2*Hgt). ----------#
   size     = ifelse( test = pft$tropical[zpft] & (! pft$liana[zpft]) & (iallom %in% 3)
                    , yes  = dbhuse * dbhuse * hgt
                    , no   = dbhuse
                    )#end ifelse
   #---------------------------------------------------------------------------------------#



   #----- Find leaf biomass. --------------------------------------------------------------#
   bleaf  = pft$b1Bl[zpft] / C2B * size ^ pft$b2Bl[zpft]
   lai    = nplant * pft$SLA[zpft] * bleaf
   #---------------------------------------------------------------------------------------#


   return(lai)
}# end function size2lai
#==========================================================================================#
#==========================================================================================#
