%% Example of globalRegrid() function to map NDSI granules to global grid
% Date: April 5, 2023
% Author: Jeremy Johnston
%
% This script provides an example usage of the globalRegrid() function by
% first downloading all MODIS NDSI granules for a single date then mapping 
% them to a global lat/lon grid with a resolution of 0.05-degrees, provides
% an example doing the same for VIIRS NDSI.
%
% Note: globalRegrid() runs the fastest and most reliable when downloading
% files using wget. Thus, it requires wget be installed on your local
% machine and the path to the executable indicated. MATLAB built in
% functions websave() and urlwrite() can be used with slight code
% modifications, but have proven unreliable and can be very slow.

%% Example 1: Map NDSI globally for a single date using MODIS data
raw_savepath = [pwd '/MODIS_raw/']; %create location to download files locally
processed_savepath = [pwd '/MODIS_processed/']; %create location to store processed global grids lcoally
username = 'jjohns60'; %%INPUT YOUR NSIDC USERNAME
password = '2080Bent!'; %%INPUT YOUR NSIDC PASSWORD
wget_path = '/Users/jjohns/Desktop/wget'; %%PATH TO WGET EXECUTABLE ON YOUR SYSTEM
HTTPS_PATH = 'https://n5eil01u.ecs.nsidc.org/MOST/MOD10A1F.061/'; %path to dataset at NSIDC server
DATASET_ID = 'MOD10A1F'; %dataset id
FILE_EXT = 'hdf'; %file extension of MOD10A1F data
VAR_PATH = 'CGF_NDSI_Snow_Cover'; %path to variable in downloaded file
START_DATE = datetime(2020,10,31); %start date to process
END_DATE = datetime(2020,10,31); %last date to process
TARGET_RES = 0.05; %fit to 0.05-degree global lat/lon grid

%create directories to store datasets
mkdir(raw_savepath);
mkdir(processed_savepath);

%run regridding function (outputs will be located in the specified folder 'processed_savepath')
globalRegrid(raw_savepath,processed_savepath,username,password,wget_path,HTTPS_PATH,DATASET_ID,FILE_EXT,VAR_PATH,START_DATE,END_DATE,TARGET_RES);
%% Note:
% This function takes ~3-5 minutes to process a single day (globally). The
% function may hang if NSIDC is receiving too many data requests. If quit, 
% restarting processing should resume only downloading the still missing 
% files for the specified date.


%% Example 2: Map NDSI globally for a single date using VIIRS data
raw_savepath = [pwd '/VIIRS_raw/']; %create location to download files locally
processed_savepath = [pwd '/VIIRS_processed/']; %create location to store processed global grids lcoally
HTTPS_PATH = 'https://n5eil01u.ecs.nsidc.org/VIIRS/VNP10A1F.001/'; %path to dataset at NSIDC server
DATASET_ID = 'VNP10A1F'; %dataset id
FILE_EXT = 'h5'; %file extension of MOD10A1F data
VAR_PATH = '/HDFEOS/GRIDS/NPP_Grid_IMG_2D/Data Fields/CGF_NDSI_Snow_Cover'; %path to variable in downloaded file
START_DATE = datetime(2020,10,31); %start date to process
END_DATE = datetime(2020,10,31); %last date to process
TARGET_RES = 0.05; %re-grid to 0.05-degree global grid

%create directories to store datasets
mkdir(raw_savepath);
mkdir(processed_savepath);

%run regridding function (outputs will be located in the specified folder 'processed_savepath')
globalRegrid(raw_savepath,processed_savepath,username,password,wget_path,HTTPS_PATH,DATASET_ID,FILE_EXT,VAR_PATH,START_DATE,END_DATE,TARGET_RES);
%% Note:
% This function takes ~3-5 minutes to process a single day (globally). The
% function may hang if NSIDC is receiving too many data requests. If quit, 
% restarting processing should resume only downloading the still missing 
% files for the specified date.