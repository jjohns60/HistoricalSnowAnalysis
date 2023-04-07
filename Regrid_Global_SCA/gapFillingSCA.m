function gapFillingSCA(savepath,datapath,method)
%gapFillingSCA Using input binary SCA dataset (with missing values set to
%255) replaces all missing data using gap-filling procedure
%   DETAILED DESCRIPTION - Gap Filling Missing Land/Ice Pixels, water
%   filled based on masks developed by PAPER 2015---
%   
%   
%   
%
%   INPUTS:
%   savepath - path to store gap filled files
%   datapath - path to SCA data GeoTiffs
%   method = if using previous timestep gap filling approach set == 1, if 
%   using nearest valid data approach (will do this automatically on the
%   first iteration if there is no prior timestep), set == 2
%

%set default
if nargin == 2
    method = 1;
end

%identify files for gap-filling procedure
files = dir([datapath '*.tif']);

%% Create paths to supporting datasets
%create path to dataset including the absolute number of days from the 
%climatological maximum mean air temperature
%Created from ERA5-Land 2m air temperature data (2000-2022)
D_Tmax_path = [pwd '/Data Inputs/DaysFromMaxTair/'];

%path to never snow mask folder, contains indication of areas that have
%less than a 0.1% probability of snow cover for the given D_Tmax value
nosnow_mask_path = [pwd '/Data Inputs/NeverSnowMasks/'];

%load in the hyperplanes defining 99% snow probability (snow) vs. 1% snow probability (no-snow) (only applied over locations with missing data)
%these were learned by summarizing the probability of snow by latitude,
%elevation, and D_Tmax (days since max temp, seasonality) bands across all
%MODIS SCA observations over the 2000-2022 period of record
SCAp_99_function = load([pwd '/Data Inputs/GeneralSCApModel/SCAp_99.mat']);
SCAp_1_function = load([pwd '/Data Inputs/GeneralSCApModel/SCAp_1.mat']);
SCAp_99_function = SCAp_99_function.SCAp_99_function;
SCAp_1_function = SCAp_1_function.SCAp_1_function;

%persistent variables to be created only on first function call (to speed processing)
%create a grid containing absolute latitudes
%get image resolution from first file
f1 = files(1).name;
f1 = split(f1,'_');
res = str2double(f1{5});
persistent AbsLat
%grid center absolute latitude (0.01-degree)
AbsLat = abs(gridcreate(res));

%create path to elevation (in meters above sea level, 0.01-degree)
persistent E
E = readgeoraster([pwd '/Data Inputs/Elevation/elevation_SRTMplus_0_01.tif']);
E = imresize(E,[180/res 360/res],'nearest');

%load in 0.01 degree water mask, based on Mikelsons 2021 paper
water_mask = load([pwd '/Data Inputs/WaterMask/watermask_01_Mikelsons2021.mat']);
water_mask = water_mask.watermask;
water_mask = imresize(water_mask,[180/res 360/res],'nearest');

% Iterate Across all SCA Files/Days & perform gap filling steps
%loop through all days
iter = 0;
for i = 1:length(files)   

    tic
    %get file name
    file = files(i).name; disp(file);
    %get month and day of month from file (assumes file naming convention produced by prior processing steps)
    f = split(file,'_');
    mm = f{3};
    dd = f{4};
    %load in SCA
    [SCA_i,R] = readgeoraster([datapath file]);
    
    %apply water mask from Mikelsons 2021 manuscript (at 0.01 resolution)
    SCA_i(water_mask == 1) = 2;

    %report initial % of missing data
    disp(['Missing Data: ' num2str((sum(SCA_i(:) == 255)/(size(SCA_i,1)*size(SCA_i,2))) * 100,3) '%'])
    
    %set below 85S to always snow, south pole (speeds processing)
    SCA_i(size(SCA_i,1) - size(SCA_i,1)*5/180:size(SCA_i,1),:) = 1;
  
    %report % of missing data after Antarctica filling (step 1)
    disp(['Post Step 1: ' num2str((sum(SCA_i(:) == 255)/(size(SCA_i,1)*size(SCA_i,2))) * 100,3) '%'])

    %select and use pre-defined 'never snow' mask for the given date (snow line based on the probability of snow being <0.1%)
    NoSnowMask = load([nosnow_mask_path 'NoSnowMask_SCAp0_01Less_' mm '_' dd '.mat']);
    NoSnowMask = NoSnowMask.no_snow_mask;
    %resample nosnow mask to the size of the input grid
    NoSnowMask = imresize(NoSnowMask,size(SCA_i),'nearest');
    SCA_i(NoSnowMask == 1) = 0; %apply mask, ensure masked areas are set as land pixels

    %report % of missing data after applying no snow mask (step 2)
    disp(['Post Step 2: ' num2str((sum(SCA_i(:) == 255)/(size(SCA_i,1)*size(SCA_i,2))) * 100,3) '%'])

    % (1) Using snow line to gap-fill
    % (1a) identify and load in appropriate DaysfromTmax grid
    D_Tmax_i = load([D_Tmax_path 'DaysFromTmax_' mm '_' dd '.mat']);
    D_Tmax_i = D_Tmax_i.A;
    %resample D_Tmax to size of input SCA grid
    D_Tmax_i = imresize(D_Tmax_i,size(SCA_i),'nearest');

    % (1b) identify no data cells and apply snow line masking (>99% snow probability = snow & <1% = no snow)
    %resample E and AbsLat to the size of the input SCA grid
    ind = (SCA_i == 255) & (D_Tmax_i ~= 999); %identify no data locations with valid predictor data
    x = double(D_Tmax_i(ind));
    y = AbsLat(ind);
    E_i = E(ind);
    fill = uint8(zeros(size(x)) + 255); %create array to store gap-filled values
    
    % (1c) apply masks
    %for SCAp 99% if elevation predicted snow line is below actual elevation, set to snow
    E_snowline = feval(SCAp_99_function,x,y);
    idx = E_snowline < E_i;
    fill(idx) = 1; %update gap fill array with snow  
    
    %for SCAp 1% if elevation predicted snow line is below above actual elevation, set to no snow
    fill2 = fill(~idx); %perform only on remaining no-data pixels
    E_snowline = feval(SCAp_1_function,x(~idx),y(~idx));
    idx2 = E_snowline > E_i(~idx); %identify no snow locations
    fill2(idx2) = 0;
    fill(~idx) = fill2; %update fill array with no snow locations

    %update missing locations with gap filled data (where applicable)
    SCA_i(ind) = fill;
    
    %report % of missing data after applying snow likely and snow unlikely masks (step 3)
    disp(['Post Step 3: ' num2str((sum(SCA_i(:) == 255)/(size(SCA_i,1)*size(SCA_i,2))) * 100,3) '%'])
    
    %final gap filling approach
    if method == 1 && iter > 0 %must occur as second timestep and beyond

        % (1d) fill remaining missing cells using previous timestep
        % (dependent on having a fully gap-filled previous date)
        % It is important to ensure that this previous date had minimal
        % missing data - APPLY IN POST PROCESSING TO CORRECT ERRANT FILES
        savefiles = dir([savepath '*.tif']);
        SCA_GF = readgeoraster([savepath savefiles(i - 1).name]);
        ind = (SCA_i == 255);
        SCA_i(ind) = SCA_GF(ind);

    %this method will be used on the first loop regardless of the selected
    %method
    elseif method == 2 || iter == 0 %slower, especially in cases of a lot of missing data

        % (1d) fill remaining missing cells using 33% rule (if >33% of valid surrounding pixels are snow, is considered snow)
        % similar sliding window approach as: Riggs, Hall, Vuyovich, and GiGirolamo 2022
        %NOTE: should be edited if using grid size larger than 0.1
        w = round(0.1/res); %window size is 0.1 x 0.1-degree (they used 0.25 x 0.25-degree)
        [row,col] = find(SCA_i == 255); %locate all no data cells by row and col numbers
        for ii = 1:length(row)

            %identify individual pixel
            row_ii = row(ii);
            col_ii = col(ii);

            %create lists of rows/columns to extract (ensure only valid indices)
            row_list = (row_ii - w/2):(row_ii + w/2); row_list = row_list(row_list > 0 & row_list <= size(SCA_i,1));
            col_list = (col_ii - w/2):(col_ii + w/2); col_list = col_list(col_list > 0 & col_list <= size(SCA_i,2));

            %get surrounding data
            d = SCA_i(row_list,col_list);

            %remove fill or water cells
            ind = (d == 255) | (d == 2);
            d = d(~ind);

            %apply 33% rule for snow and set pixel to new value
            if sum(d(:) == 1)/length(d) > 0.33

                SCA_i(row_ii,col_ii) = 1;

            else %otherwise is set as land (is also the default if there is missing data)

                SCA_i(row_ii,col_ii) = 0;

            end

        end

    end

    %verify that all missing data was removed
    if (sum(SCA_i(:) == 255) == 0)
        disp('Gap filling completed')
        savefile = [file(1:9) 'GF_' file(10:end)];
    else
        disp("warning, all no data values NOT removed, saving with leading '_'")
        savefile = ['_' file(1:9) 'GF_' file(10:end)]; %will save file with leading underscore
    end

    %save data
    geotiffwrite([savepath savefile],SCA_i,R);
    toc

    iter = iter + 1; %track number of completed iterations

end