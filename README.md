# pv_calc

Calculates potential vorticity off of netcdf files based on relative vorticity, u component of wind, v component of wind, and temperature.

Also required is logarthmic surface pressure levels, and other pressure levels are calculated through Ak and Bk and exponentiation. 

levels.py includes constants required for calculation. Ak and Bk are used for pressure calculation, and are sourced from ECMWF documentation. These should not be changed. const\_pressures can be changed, but should not exceed 59 indices, and should be increasing monotonically.

