%% Example Application of accessing and subsetting NSIDC data
% Date: April 5, 2023
% Author: Jeremy Johnston
%
% This script provides example usage of functions to access, download,
% re-project, and interpolate NSIDC datasets (specifically VIIRS and MODIS 
% NDSI) and save as GeoTiff(s) to a local directory. Example 1 shows how to
% download a two individual granules from the VIIRS VNP10A1F product. 
% Example 2 extracts and re-grids NDSI data for a bounding region in 
% southern New Jersey and Delaware, USA
%

%% Example 1: use NSIDC_HTTPS_ACCESS function to download a individual VIIRS NDSI granules
HTTPS_PATH = 'https://n5eil01u.ecs.nsidc.org/VIIRS/VNP10A1F.001/'; %HTTPS path to dataset on NSIDC
username = ''; % ADD YOUR NSIDC USERNAME HERE
password = ''; % ADD YOUR NSIDC PASSWORD HERE
wget_path = ''; %%PATH TO WGET EXECUTABLE ON YOUR SYSTEM
date = datetime(2015,1,27); %identify day of interest as datetime object
str_array = {'VNP10A1F','h5','h12v04','h12v05'}; %specify the dataset ID, filetype, then tiles of interest (can add multiple)
savepath = [pwd '/VIIRS_NDSI_data/']; %set save location in local directory

%make local directory to store outputs
mkdir(savepath);

%run function
NSIDC_HTTPS_ACCESS(HTTPS_PATH,username,password,wget_path,date,str_array,savepath)


%% Example 2: downloadSubsetNDSI which will download NDSI and interpolate for a specific region
savepath = [pwd '/MODIS_Regridded_NDSI_data/']; %set save location in local directory
dataset_id = 'MOD10A1F'; %specify dataset type (currently, only VIIRS and MODIS CGF NDSI are supported)
bounding_region = [39 40 -76 -75 0.001]; %set region and desired spatial resolution [lat_min lat_max lon_min lon_max res_degree]
start_date = datetime(2015,1,27); %start date of period of interest, as datetime object
end_date = datetime(2015,1,29); %end date of period of interest, as datetime object
interpolation_method = 'natural'; %specify interpolation method, see downloadSubsetNDSI for supported interpolation methods

%make local directory to store outputs
mkdir(savepath);

%run function, outputs will be located in created folder at 'savepath'
downloadSubsetNDSI(savepath,dataset_id,username,password,wget_path,...
    bounding_region,start_date,end_date,interpolation_method);
%% Notes: 
% These functions use MATLAB parallelization to speed the download of 
% multiple files, thus the download numbers may not be sequential.
% During times of high NSIDC traffic the download function may hang