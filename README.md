# Historical Snow Analysis

This repository contains MATLAB codes for analyzing, downloading, and creating insights into global snow regimes using snow cover datasets. These scripts have been developed by Jeremy Johnston (jeremy.johnston[at]unh.edu) at the University of New Hampshire Earth Systems Research Center (ESRC). Please reach out to Jeremy with questions on usage or troubleshooting.
A companion paper 'Global Snow Seasonality Regimes from Satellite Records of Snow Cover' is undergoing review, (WHEN PUBLISHED WILL ADD LINK/DOI HERE).

![alt text](https://github.com/jjohns60/HistoricalSnowAnalysis/blob/main/SCAheader_image.jpg?raw=true)

## Folders & Files

* `Regrid_Global_SCA`
<br> MATLAB codes for converting tiled MODIS or VIIRS snow cover data (NDSI) to high-resolution, gap-filled global snow cover rasters.
* `Analysis_Functions`
<br> MATLAB functions for deriving snow metrics from the produced snow cover rasters.
* `Subset_SCA`
<br> MATLAB codes for extracting, interpolating, and re-gridding region specific NDSI observations from VIIRS or MODIS.
* `Download_Historical_Point_Observations`
<br> MATLAB function for extracting Global Historical Climatology Network (GHCN) observations from any site (or region) globally.


## Relevant Data Links and Sources

* MODIS Cloud-gap-filled (CGF) Snow Products (https://modis-snow-ice.gsfc.nasa.gov/?c=MOD10A1F)
* VIIRS Daily CGF Snow Product, VNP10A1F (https://nsidc.org/data/vnp10a1f/versions/1)
* Article on product accuracy, Hall et al., 2019 (https://hess.copernicus.org/articles/23/5227/2019/)
* Landing page for daily observational data from GHCN (https://www.ncei.noaa.gov/products/land-based-station/global-historical-climatology-network-daily)
* Rutgers University Global Snow Lab, Northern Hemisphere Snow Cover Portal (https://climate.rutgers.edu/snowcover/)
