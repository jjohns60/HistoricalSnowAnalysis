%% Example Application of GHCN_HTTPS() MATLAB function
% Date: April 5, 2023
% Author: Jeremy Johnston
%
% This script provides example usage of the GHCN_HTTPS() function. First, 
% we return a list of sites containing daily precipitation and snowfall in
% south eastern New Hampshire having at least 50-year long records. We then
% provide examples of downloading all of the data meeting the criteria as
% well as using only a subset of the identified sites.
%

%% Example 1: for identifying GHCN sites meeting a given criteria in a specific region of interest (ROI)
filepath = [pwd '/GHCN_data/']; %will create a folder in the current working directory to store GHCN data
variables = {'PRCP','SNOW'}; %download precipitation (mm or in) and snowfall (mm or in)
search_region = [42.92 43.34 -71.40 -70.84]; %latitude and longitude coordinates of bounding region in SE New Hampshire [lat_min,lat_max,lon_min,lon_max]
min_length = 50; %minumum record length in years to extract (== 50 years)
download = 0;

%create directory to store outputs
mkdir(filepath);

%get site information tables (returns a MATLAB table object with site IDs,
%locations, variables of interest, and record length)
site_list1 = GHCN_HTTPS(filepath,variables,search_region,min_length,download);
%% Note on Example 1: 
% In addition to creating a MATLAB table object (site_list), this function
% also saves a .csv version to the 'filepath' location


%% Example 2: for downloading data from the sites in the ROI
%returns exactly the same outputs as Example 1, however will also download 
%the data to the created directory
download = 1;
site_list2 = GHCN_HTTPS(filepath,variables,search_region,min_length,download);


%% Example 3: download only a subset of sites, based on user input table or path to .csv file
%modify site_list2 to only contain the first site
%also creates a new .csv table, titled '_stationInfoMod.csv'
site_list2 = site_list2(1,:);
GHCN_HTTPS(filepath,site_list2); %download indicated files only
