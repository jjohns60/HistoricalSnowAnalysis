%% Example of globalRegrid() function to map NDSI granules to global grid
% Date: April 5, 2023
% Author: Jeremy Johnston
%
% This script provides an example usage of the globalRegrid() function by
% first downloading all MODIS NDSI granules for a single date then mapping 
% them to a global lat/lon grid with a resolution of 0.05-degrees (Example 
% 1), provides an example doing the same for VIIRS NDSI (Example 2). The
% final example (Example 3), applies functions to convert NDSI to SCA, then
% produce completely gap-filled SCA rasters (0 = land, 1 = snow, 2 = water)
%
% Note: globalRegrid() runs the fastest and most reliable when downloading
% files using wget. Thus, it requires wget be installed on your local
% machine and the path to the executable indicated. MATLAB built in
% functions websave() and urlwrite() can be used with slight code
% modifications, but have proven unreliable and can be very slow.

%% Example 1: Map NDSI globally for a single date using MODIS data
raw_savepath = [pwd '/MODIS_raw/']; %create location to download files locally
processed_savepath = [pwd '/MODIS_processed/']; %create location to store processed global grids lcoally
username = ''; %%INPUT YOUR NSIDC USERNAME
password = ''; %%INPUT YOUR NSIDC PASSWORD
wget_path = ''; %%PATH TO WGET EXECUTABLE ON YOUR SYSTEM
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


%% Example 2: Map NDSI globally using VIIRS data for a 6-day period
raw_savepath = [pwd '/VIIRS_raw/']; %create location to download files locally
processed_savepath = [pwd '/VIIRS_processed/']; %create location to store processed global grids lcoally
HTTPS_PATH = 'https://n5eil01u.ecs.nsidc.org/VIIRS/VNP10A1F.001/'; %path to dataset at NSIDC server
DATASET_ID = 'VNP10A1F'; %dataset id
FILE_EXT = 'h5'; %file extension of MOD10A1F data
VAR_PATH = '/HDFEOS/GRIDS/NPP_Grid_IMG_2D/Data Fields/CGF_NDSI_Snow_Cover'; %path to variable in downloaded file
START_DATE = datetime(2020,10,31); %start date to process
END_DATE = datetime(2020,11,5); %last date to process
TARGET_RES = 0.05; %re-grid to 0.05-degree global grid

%create directories to store datasets
mkdir(raw_savepath);
mkdir(processed_savepath);

%run regridding function (outputs will be located in the specified folder 'processed_savepath')
globalRegrid(raw_savepath,processed_savepath,username,password,wget_path,HTTPS_PATH,DATASET_ID,FILE_EXT,VAR_PATH,START_DATE,END_DATE,TARGET_RES);
%% Note:
% This function takes ~3-5 minutes to process a single day (globally). The
% function may hang if NSIDC is receiving too many data requests. If quit, 
% restarting processing should resume by only downloading the still missing 
% files for the specified date.


%% Example 3: Using regridded VIIRS NDSI datasets, convert to binary SCA, and apply gap filling procedure
NDSI_path = [pwd '/VIIRS_processed/']; %provide location of the downloaded global VIIRS  files
SCA_savepath = [pwd '/VIIRS_SCA/']; %will create new path storing files converted to binary SCA
SCA_GF_savepath = [pwd '/VIIRS_SCA_GapFilled/']; %create new path storing the completely gap-filled SCA files
method = 1; %method applies simple NDSI > 0.1 threshold for snow cover detection (see function documentation for other methods)

%create directories to store output SCA datasets
mkdir(SCA_savepath);
mkdir(SCA_GF_savepath);

%loop through all files on the NDSI_path
files = dir([NDSI_path '*.tif']);
for i = 1:length(files)
    %create full path to NDSI file
    file = files(i).name;
    filepath = [NDSI_path file]; disp(file); %show file being processed
    %convert to SCA using a simple NDSI threshold of 0.1 (10 in rescaled GeoTiffs)
    %in addition to saving the files, also returns outputs including the 
    %full data grid, and corresponding datetime
    [SCA,dt] = convertNDSItoSCA(filepath,SCA_savepath,method);
end

%using binarized SCA as inputs, apply gap filling procedure to remove no
%data cells (=255)


%% Note:
% For methods 4 (based on land cover) and method 6 (based on NDVI) there is
% a third argument required by the convertNDSItoSCA() function, signifying
% the path to the corresponding dataset. See function description for
% details. Also, the gap-filling procedure is described in detail in the
% accompanying manuscript, 'Global Snow Seasonality Regimes from Satellite
% Records of Snow Cover', under review in the Journal of Hydrometeorology,
% 2023 (DOI PENDING - WILL ADD HERE ONCE PUBLISHED)


