Script tower_gapfill.r processes data from eddy covariance towers to fill gaps in 
   meteorological variables.  The current example uses data from nearby weather stations, 
   and requires that script gen_othersites.r (in sub-directory other_sites) is run first.

The script itself has quite a few steps that are specific to the data sets being read, so
  expect to edit it extensively when running for other locations.

Because I do now own the original data, I am leaving just the first few lines of the 
  input file to show how the input data should be formatted.  
  
The example given had been pre-processed with EddyPro and the variables follow the format
  and units from EddyPro. For any meaningful use within ED2, it is a good idea to provide 
  a few years of data.  In addition, the script takes daily precipitation data from the
  nearby weather station to fill in any remaining missing data for rainfall.  This follows
  the format provided by the historic records from INMET, and I left the first few lines
  as an example of the data set format.

Script model_fill_nee.r fills gaps for NEE to compute GPP and ecosystem respiration.
   Definitely check the settings at the beginning of the script.

Script plot_ts.r creates simple plots to show the time series, which can be useful for 
  checking the results.
