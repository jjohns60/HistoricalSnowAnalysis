function globalRegrid(raw_savepath,processed_savepath,username,password,wget_path,HTTPS_PATH,DATASET_ID,FILE_EXT,VAR_PATH,START_DATE,END_DATE,TARGET_RES)
%globalRegrid Summary of this function goes here
%   Effectively downloads, re-grids, and resamples global NDSI observations
%   to a grid with a specified input resolution. Functionality has been 
%   tested effectively for VIIRS and MODIS cloud gap filled NDSI snow cover
%   products located at NSIDC.
%
%   INPUTS:
%   
%   raw_savepath - local directory to save downloaded granules
%   processed_savepath - local directory to save re-gridded GeoTIFFs
%   username - NSIDC logon username as string
%   password - NSIDC logon password as string
%   wget_path - path to 'wget' executable on your local system
%   HTTPS_PATH - HTTPS address to dataset on NSIDC
%   DATASET_ID - the dataset id (i.e., 'MOD10A1F')
%   FILE_EXT - extension of the dataset files to download (i.e., 'hdf')
%   VAR_PATH - path within downloaded datafile to the dataset of interest
%   START_DATE - datetime object specifying the first date to process
%   END_DATE - datetime object specifying the last date to process
%   TARGET_RES - resolution of global target grid (should be larger than
%   0.005-degree)
%
%
%   For working with NSIDC Gap Filled Snow Cover Datasets
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   HTTPS_PATH should be set as either:
%   for MODIS: 'https://n5eil01u.ecs.nsidc.org/MOST/MOD10A1F.061/'
%   for VIIRS: 'https://n5eil01u.ecs.nsidc.org/VIIRS/VNP10A1F.001/'
%
%   DATASET_ID should be set as either:
%   for MODIS: 'MOD10A1F'
%   for VIIRS: 'VNP10A1F'
%
%   FILE_EXT should be set as either:
%   for MODIS: 'hdf'
%   for VIIRS: 'h5'
%
%   VAR_PATH should be set as either:
%   for MODIS: 'CGF_NDSI_Snow_Cover'
%   for VIIRS: '/HDFEOS/GRIDS/NPP_Grid_IMG_2D/Data Fields/CGF_NDSI_Snow_Cover'
%
%   NOTE: As is, requires MATLAB Parallel Computing Toolbox (but can be easily revised to use normal MATLAB loops)

%Set local directories to store files
TEMP_FILES_PATH = raw_savepath; %this is where daily data will be stored for processing
OUT_PATH = processed_savepath; %this is where the processed, re-gridded files will be stored

%tile (10-deg grid) coordinates file must be located in same working
%directoy
BOUND_10deg_PATH = [pwd '/sn_bound_10deg.txt']; %indicates the bounding coordinates of each tile in the global sinusoidal grid

%convert date range to list of daily dates (inclusive)
DATES = START_DATE:END_DATE;

%create lat/lon coordinate grid at target_res & save to OUT_PATH
%   Note: Both uses grid centers
LAT = single(90 - 0.5*TARGET_RES:-TARGET_RES:-90 + 0.5*TARGET_RES);
LON = single(-180 + 0.5*TARGET_RES:TARGET_RES:180 - 0.5*TARGET_RES);
LAT = repmat(LAT',[1 length(LON)]);
LON = repmat(LON,[size(LAT,1) 1]);
%get all bounding coordinates that fall on grid
BOUNDING_COORDS = readtable(BOUND_10deg_PATH);
BOUNDING_COORDS = BOUNDING_COORDS(1:end-1,1:6);
BOUNDING_COORDS = BOUNDING_COORDS(BOUNDING_COORDS.lon_min ~= -999,1:6);

f = waitbar(0,'Saving grid spec file and starting parallel pool...');
f.Children.Title.Interpreter = 'none';
%uncomment to save the created grid
%save([OUT_PATH 'grid_' num2str(TARGET_RES) '.mat'],'LAT','LON','-v7.3');
S = load('handel.mat'); %loads sound to play after each iteration
%start up parallelpool for parallel processing
workers = 8; %10 workers shown to be too many, 4-8 work well
poolobj = parpool(workers);

%Loop through all target dates
t = []; %used to store elapsed times
est_rem_msg = strcat("Estimated total time remaining (0/",num2str(length(DATES))," days complete): TBD"); %time remaining estimate in minutes
for i = 1:length(DATES)
    skip_i = '';
    t1 = tic;

    %determine target date
    date = DATES(i);
    msg1 = ['Downloading files for ' datestr(date) '...'];
    msg2 = est_rem_msg;
    waitbar(0.1,f,{msg1,msg2});
    
    %Download relevant datafiles from NSIDC
    try
        NSIDC_HTTPS_ACCESS(HTTPS_PATH,username,password,wget_path,date,{DATASET_ID FILE_EXT},TEMP_FILES_PATH);
        
    catch %if server returns error (is likely from being busy), wait and retry
        skip_i = '*'; %keep track of if there were errors during download 
        pause(rand*3)
        NSIDC_HTTPS_ACCESS(HTTPS_PATH,username,password,wget_path,date,{DATASET_ID FILE_EXT},TEMP_FILES_PATH);
    end
    msg1 = 'File downloads complete, starting processing...';
    msg2 = est_rem_msg;
    waitbar(1,f,{msg1,msg2});

    %Loop through downloaded files in temporary data folder
    if contains(DATASET_ID,'VNP')
        files = dir([TEMP_FILES_PATH '*.h5']);
    elseif contains(DATASET_ID,'MOD')
        files = dir([TEMP_FILES_PATH '*.hdf']);
    end    
    files = {files.name};
    %ensure files are from correct date
    d_str = strcat("A",num2str(year(date)),num2str(day(date,'dayofyear'),'%03.f'));
    idx = contains(files,d_str);
    files = files(idx);

    %create grids at target resolution (each loop)
    GRID_OUT = zeros([180/TARGET_RES 360/TARGET_RES]);
    GRID_OUT = uint8(GRID_OUT); %stores final output
    MASK_OUT = GRID_OUT; %stores locations of valid data
    %N_OUT = GRID_OUT; %to store the number of data points contributing to each cell

    %loop through and place files data on rescaled global grid (~1km)
    L = length(files);
    for ii = 1:L
        
        p = ii/L;
        msg1 = strcat("Processing tile ",num2str(ii)," of ",num2str(L)," (",datestr(date),")");
        msg2 = est_rem_msg;
        waitbar(p,f,{msg1,msg2});
        
        %read file to target and display filename
        file = files{ii};
        %disp(file)

        %determine MODIS/VIIRS sinusoidal grid tile info using filename
        hv = split(file,'.');
        hv = hv{3};
        H = str2double(hv(2:3));
        V = str2double(hv(5:6));

        %load in gap filled snow cover data for each tile
        try %files may be corrupted or unreadable (rare)
            if contains(DATASET_ID,'VNP')
                D = h5read([TEMP_FILES_PATH file],VAR_PATH);
            elseif contains(DATASET_ID,'MOD')
                D = hdfread([TEMP_FILES_PATH file],VAR_PATH);
            end    
        catch %if file is unreadable, use previous file and replace with all fill values
            skip_i = '*';
            disp(['File unreadable (' file '), replacing with fill'])
            D = zeros(size(D),'uint8') + 255;
        end

        %reorient so north is the top and west is to the left
        if contains(DATASET_ID,'VNP')
            D = flipud(rot90(D));
        elseif contains(DATASET_ID,'MOD')
            %no reorientation currently required for MODIS files
        end
        
        %determine corresponding geographic coordinates
        [lat_ii,lon_ii] = inverseMappingSINgrid(size(D,1),H,V);

        %use H and V + bounding coordinates to isolate the general lat/lon area
        idx = BOUNDING_COORDS.iv == V & BOUNDING_COORDS.ih == H;
        bounds = BOUNDING_COORDS(idx,3:end);
        %convert bounds to rounded coordinates
        bounds = matchBounds2grid(bounds,TARGET_RES);

        %find indices of grid from within large grid using bounding coordinates
        [~,col_min] = min(abs(LON(1,:) - bounds.lon_min));
        [~,col_max] = min(abs(LON(1,:) - bounds.lon_max));
        [~,row_max] = min(abs(LAT(:,1) - bounds.lat_min));
        [~,row_min] = min(abs(LAT(:,1) - bounds.lat_max));

        %get correct portion of grid
        lat_grid = LAT(row_min:row_max,col_min:col_max);
        lon_grid = LON(row_min:row_max,col_min:col_max);

        %takes ~1 sec prior to this line
        %identify nearest pixels in small target grid, place in cell array
        [DATA,~,MASK_cnt] = gridvaluesearch(single(D),lat_ii,lon_ii,lat_grid,lon_grid,'SCA');

        %store locations with valid data on global mask
        %Mgrid = uint8((DATA_cnt >= MASK_cnt) & (DATA_cnt > 0)); %valid where there are the same number or more valid data pixels than fill
        
        %store number of valid data points going into each cell
        %Ngrid = uint8(DATA_cnt);
        %Ngrid(Mgrid == 0) = 0;

        %add to final output grids
        % GRID_OUT contains sums of valid data,
        % MASK_OUT contains number of individual valid data points
        GRID_OUT(row_min:row_max,col_min:col_max) = GRID_OUT(row_min:row_max,col_min:col_max) + uint8(DATA);
        MASK_OUT(row_min:row_max,col_min:col_max) = MASK_OUT(row_min:row_max,col_min:col_max) + uint8(MASK_cnt);
        %N_OUT(row_min:row_max,col_min:col_max) = N_OUT(row_min:row_max,col_min:col_max) + Ngrid;
        
        %clear contents of temporary data folder
        delete([TEMP_FILES_PATH file])
    end
    
    %get average of valid values within each cell by dividing sum by N
    GRID_OUT = GRID_OUT./MASK_OUT;
    
    %fill all cells that did not have any valid data
    GRID_OUT(MASK_OUT == 0) = 255;

    %store the durations in a list (can use to get time remaining
    %estimates)
    t = [t, toc(t1)]; %will store the processing time for each loop (in seconds)
   
    %save data grid and update progress bar
    name = strcat(DATASET_ID,"_",datestr(date,'yyyy_mm_dd'),"_",num2str(TARGET_RES),skip_i,'.tif');
    msg1 = strcat("Processed. Saving file ",name," ...");
    msg2 = strcat(num2str(t(end)/60,3)," min elapsed");
    waitbar(1,f,{msg1,msg2});
    %save(strcat(OUT_PATH,name),'GRID_OUT','-v7.3'); %for saving as .mat
    R = georefcells([LAT(end) LAT(1)],[LON(1) LON(end)],size(GRID_OUT));
    R.ColumnsStartFrom = 'north';
    geotiffwrite(strcat(OUT_PATH,name),GRID_OUT,R)
     
    %show estimate of time remaining
    N = (length(DATES) - i);
    est_rem = mean(t)*N;
    est_rem_msg = strcat("Estimated total time remaining (",num2str(i),"/",num2str(length(DATES))," days complete): ",num2str(est_rem/(60*60),3)," hrs"); %time remaining estimate in minutes
    
    %play sound each time a day is completely processed
    sound(S.y(1000:18000),S.Fs)
end

%close parallel processing pool
delete(poolobj);
%close progress bar
delete(f);
end