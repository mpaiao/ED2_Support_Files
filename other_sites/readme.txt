Script gen_othersites.r creates a pool of nearby stations to help filling gaps in the
main time series.  It should be run prior of running the main script in the parent
directory.

The script itself has quite a few steps that are specific to the data sets being read, so
expect to edit it extensively when running for other locations.

Because I do now own the original data, I am leaving just the first few lines to show 
how the input data should be formatted.  For any meaningful use within ED2, it is a good
idea to provide a few years of data.

-------------------------------
Current format of input data:
-------------------------------
year      -- Calendar year
month     -- Calendar month
day       -- Calendar date
hour      -- Calendar hour  (both local time and UTC work, but you may need to adjust 
             variable "sites" in gen_othersites.r accordingly).
atm.tmp   -- Air temperature              [    degC]
atm.rhv   -- Relative humidity            [       %]
atm.tdew  -- Dewpoint temperature         [    degC]
atm.prss  -- Atmospheric pressure         [     hPa]
atm.vdir  -- Wind direction               [     deg]
atm.vels  -- Wind speed                   [     m/s]
rain      -- Precipitation                [      mm]
rshort.in -- Incoming shortwave radiation [kJ/m2/hr]
-------------------------------
Note: If your input data has different units, you may also adjust variable "inmet" in 
      gen_othersites.r
-------------------------------
