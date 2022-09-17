function downloadSubsetNDSI(out_path,dataset_id,bounds,start_date,end_date,interp)
%downloadNDSIsubset Downloads and subsets Normalized Difference Snow Index
%observations from MODIS or VIIRS to a specified region
%
%   This function exctracts gap-filled snow cover products from the NSIDC
%   database, crops to the region of interest over the dates of interest, 
%   matches the data to an input grid either specified by an input GeoTiff 
%   or coordinates (must also then include the target resolution), and 
%   saves the output as a GeoTiff file in the specified location. A
%   dependent download function (NSIDC_HTTPS_ACCESS) uses parallel
%   computing to speed tile downloads (can be removed by changing parfor
%   loops to standard fot loops within that function). It is advised to use
%   relatively small regions (< than a few square lat/lon degrees) as
%   larger regions can be very slow during interpolation to fine grids (<0.01)
%
%   Example usage: downloadSubsetNDSI('your_path','MOD10A1F',[39 40 -76 -75 0.001],datetime(2015,1,27),datetime(2019,12,31),'natural')
%
%   INPUTS:
%
%   out_path: The location to download data, and store output files
%   
%   dataset_id: Specify the snow cover dataset to pull from, currently,
%       this can be (1) 'MOD10A1F' (MODIS Terra 500m Cloud Gap Filled Snow
%       Cover; 2000-present) or (2) 'VNP10A1F' (VIIRS 350m Cloud Gap Filled
%       Snow Cover; 2014-present)
%
%   bounds: Either an input grid in the form [lat_min lat_max lon_min
%       lon_max target_resolution] where all inputs are in degrees, an 
%       input raster to pull a referencing object from, or a matlab 
%       referencing object
%
%   start_date: The start date of which to process data from as a MATLAB
%       datetime object
%
%   end_date: The end date of which to process data from as a MATLAB
%       datetime object
%
%   interp: A string specifying the interpolation technique. Possible
%   inputs are 'linear' (linear), 'nearest' (nearest neighbor), 'natural'
%   (default, combines linear and cubic for naturally realistc output),
%   'cubic' (cubic), and 'v4' (very slow, MATLAB Biharmonic spline). This 
%   input is optional, with the default being 'natural'
%
%   OUTPUT:
%
%   Files will be downloaded to the specified path (then deleted) and
%   replaced with a subset .tif, programmatic output is null

%default setting
if nargin < 5
    error('Missing inputs')
elseif nargin == 5
    interp = 'natural';
end


%Prepare For File Downloads
%determine which datapath and variables to use
if strcmp(dataset_id,'MOD10A1F')
    HTTPS_PATH = 'https://n5eil01u.ecs.nsidc.org/MOST/MOD10A1F.061/';
    FILE_EXT = 'hdf'; %file suffix, or file type extension
    VAR_PATH = 'CGF_NDSI_Snow_Cover';
elseif strcmp(dataset_id,'VNP10A1F')
    HTTPS_PATH = 'https://n5eil01u.ecs.nsidc.org/VIIRS/VNP10A1F.001/';
    FILE_EXT = 'h5'; %file suffix, or file type extension
    VAR_PATH = '/HDFEOS/GRIDS/NPP_Grid_IMG_2D/Data Fields/CGF_NDSI_Snow_Cover';
else
    error('Invalid dataset ID')
end

%get referencing object from raster and bounding coordinates (if is string)
if ischar(bounds)
    [~,R] = readgeoraster(bounds);
    bounds = [R.LatitudeLimits R.LongitudeLimits];

    %get grid cell size in degrees
    sLon = R.CellExtentInLongitude;
    sLat = R.CellExtentInLatitude;

    %create grid matching input raster over AOI (use grid centers)
    LON = bounds(3) + 1/2*sLon:sLon:bounds(4) - 1/2*sLon;
    LAT = bounds(1) + 1/2*sLat:sLat:bounds(2) - 1/2*sLat;
    LON = repmat(LON,[length(LAT) 1]);
    LAT = repmat(flipud(LAT'),[1 size(LON,2)]);

%if inputs are numeric
elseif isnumeric(bounds)
    %get grid cell size in degrees
    sLon = bounds(5);
    sLat = bounds(5);

    %create grid (use grid centers at indicated resolution)
    LON = bounds(3) + 1/2*sLon:sLon:bounds(4) - 1/2*sLon;
    LAT = bounds(1) + 1/2*sLat:sLat:bounds(2) - 1/2*sLat;
    LON = repmat(LON,[length(LAT) 1]);
    LAT = repmat(flipud(LAT'),[1 size(LON,2)]);

    %create referencing object
    R = georefcells([min(LAT(:)) max(LAT(:))],[min(LON(:)) max(LON(:))],size(LAT));
    R.ColumnsStartFrom = 'North';

%assume the data input is a referencing object
else
    R = bounds;
    bounds = [R.LatitudeLimits R.LongitudeLimits];

    %get grid cell size in degrees
    sLon = R.CellExtentInLongitude;
    sLat = R.CellExtentInLatitude;

    %create grid matching input raster over AOI (use grid centers)
    LON = bounds(3) + 1/2*sLon:sLon:bounds(4) - 1/2*sLon;
    LAT = bounds(1) + 1/2*sLat:sLat:bounds(2) - 1/2*sLat;
    LON = repmat(LON,[length(LAT) 1]);
    LAT = repmat(flipud(LAT'),[1 size(LON,2)]);
end

%convert date range to list of daily dates (inclusive)
DATES = start_date:end_date;

%txt file stores the extent of each tile on MODIS and VIIRS sinusoidal
%global grid. File must be in same path as this function (or, the full path
%must be specified)
BOUND_10deg_PATH = 'sn_bound_10deg.txt'; 

%load in bounding coordinates of each tile
BOUNDING_COORDS = readtable(BOUND_10deg_PATH);
BOUNDING_COORDS = BOUNDING_COORDS(1:end-1,1:6);
BOUNDING_COORDS = BOUNDING_COORDS(BOUNDING_COORDS.lon_min ~= -999,1:6);

%identify sinusoidal grid tiles that cover the region of interest
idx = zeros(size(BOUNDING_COORDS.iv));
UL = [bounds(2) bounds(3)];
UR = [bounds(2) bounds(4)];
LL = [bounds(1) bounds(3)];
LR = [bounds(1) bounds(4)];
%loop through all tile grids and check if they contain any corner points
% for the region of interest
for i = 1:length(idx)
    lat_min_i = BOUNDING_COORDS.lat_min(i);
    lat_max_i = BOUNDING_COORDS.lat_max(i);
    lon_min_i = BOUNDING_COORDS.lon_min(i);
    lon_max_i = BOUNDING_COORDS.lon_max(i);
    UL_ind = lat_min_i <= UL(1) & lat_max_i >= UL(1)...
        & lon_min_i <= UL(2) & lon_max_i >= UL(2);
    UR_ind = lat_min_i <= UR(1) & lat_max_i >= UR(1)...
        & lon_min_i <= UR(2) & lon_max_i >= UR(2);
    LL_ind = lat_min_i <= LL(1) & lat_max_i >= LL(1)...
        & lon_min_i <= LL(2) & lon_max_i >= LL(2);
    LR_ind = lat_min_i <= LR(1) & lat_max_i >= LR(1)...
        & lon_min_i <= LR(2) & lon_max_i >= LR(2);
    %if any corner point is included within the tile, will mark as a tile
    %which covers to AOI
    if sum(UL_ind + UR_ind + LL_ind + LR_ind) > 0
        idx(i) = 1;
    end
end

%convert to logical and trim tile list, identify tiles of interest
idx = logical(idx);
BOUNDING_COORDS = BOUNDING_COORDS(idx,:);
tiles = cell([height(BOUNDING_COORDS) 1]);
for i = 1:height(BOUNDING_COORDS)
    tiles{i} = strcat('h',num2str(BOUNDING_COORDS.ih(i),'%02.f'),'v',num2str(BOUNDING_COORDS.iv(i),'%02.f'));
end

%show which tiles cover the AOI
disp('Tiles covering AOI:')
for i = 1:length(tiles)
    disp(tiles{i})
end

%% (3) Loop through each date, download, re-grid, and save
for i = 1:length(DATES)
    skip_i = '';

    %determine target date
    date = DATES(i);
    disp(date)
    
    
    %Download relevant datafiles from NSIDC
    try
        NSIDC_HTTPS_ACCESS(HTTPS_PATH,date,{dataset_id FILE_EXT tiles{:}},out_path);
        
    catch %if server returns error (is likely from being busy), wait and retry
        skip_i = '*'; %keep track of if there were errors during download 
        pause(rand*3)
        NSIDC_HTTPS_ACCESS(HTTPS_PATH,date,{dataset_id FILE_EXT tiles{:}},out_path);
    end
    

    %Loop through temporary downloaded files in data folder
    if contains(dataset_id,'VNP')
        files = dir([out_path '*.h5']);
    elseif contains(dataset_id,'MOD')
        files = dir([out_path '*.hdf']);
    end    
    files = {files.name}';

    %ensure files are from correct date
    d_str = strcat("A",num2str(year(date)),num2str(day(date,'dayofyear'),'%03.f'));
    idx = contains(files,d_str);
    files = files(idx);

    %loop through each relevant tile
    L = length(files);
    all_lats = [];
    all_lons = [];
    all_D = [];
    for ii = 1:L
        %read file
        file = files{ii};

        %determine MODIS/VIIRS sinusoidal grid tile info using filename
        hv = split(file,'.');
        hv = hv{3};
        H = str2double(hv(2:3));
        V = str2double(hv(5:6));

        %load in gap filled snow cover data for each tile
        try %files may be corrupted or unreadable (rare)
            if contains(dataset_id,'VNP')
                D = h5read([out_path file],VAR_PATH);
            elseif contains(dataset_id,'MOD')
                D = hdfread([out_path file],VAR_PATH);
            end    
        catch %if file is unreadable, use previous file and replace with all fill values
            skip_i = '*';
            disp("File unreadable, replacing with fill")
            D = zeros(size(D),'uint8') + 255;
        end

        %reorient so north is the top and west is to the left
        if contains(dataset_id,'VNP')
            D = flipud(rot90(D));
        elseif contains(dataset_id,'MOD')
            %files are already oriented properly
        end

        %determine corresponding geographic coordinates
        [lat_ii,lon_ii] = inverseMappingSINgrid(size(D,1),H,V);

        %trim to only retain data in the area of interest (retain a bit extra, 0.1 degree border, to ensure good interpolation)
        idx = lat_ii >= bounds(1)-0.1 & lat_ii <= bounds(2)+0.1 & lon_ii >= bounds(3)-0.1 & lon_ii <= bounds(4)+0.1;
        all_lats = [all_lats; lat_ii(idx)];
        all_lons = [all_lons; lon_ii(idx)];
        all_D = [all_D; D(idx)];

        %delete raw downloaded file (comment out to keep downloaded files)
        delete([out_path file])
        
    end
    
    %use the vector data to create a new interpolated gridded product at 1 arc-second 
    D = double(all_D);
    D(D > 100) = NaN; %set non-NDSI values to NaN
    D = griddata(all_lons,all_lats,D,LON,LAT,interp);

    %write to geotiff file
    name = strcat(dataset_id,"_",datestr(date,'yyyy_mm_dd'),'_',num2str(bounds(5)),skip_i,'.tif');
    geotiffwrite(strcat(out_path,name),D,R)

end

end