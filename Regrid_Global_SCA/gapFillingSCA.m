%% Gap Filling Missing Land/Ice Pixels
% For MODIS Terra CGF record WYs 2000-01 to 2021-22

%% Inputs

%if using previous timestep gap filling approach set == 1, if using nearest
%valid data approach, set to 2
method = 1;

%save path for new gap filled data
savepath = '/Volumes/GRA_Data_Backup/POSTDOC_UNH/SCA_REGRID/MODIS_CGF_NDSI10_GapFilled2/';

%file path for original MODIS data
filepath = '/Volumes/GRA_Data_Backup/POSTDOC_UNH/SCA_REGRID/MODIS_CGF_NDSI10/';
files = dir([filepath '*.tif']);

%path to file which includes the absolute number of days from the
%day with the climatological max mean daily air temperature (based on ERA5-Land 2m air temperature data from 2000-2022)
D_Tmax_path = '/Volumes/GRA_Data_Backup/POSTDOC_UNH/SCA_CLIMATOLOGY/GeneralSCApModel/Inputs/DaysFromMaxTair_expandedCoverage/';

%load in elevation data (in meters above sea level, 0.01-degree)
E = readgeoraster('/Volumes/GRA_Data_Backup/POSTDOC_UNH/SCA_CLIMATOLOGY/GeneralSCApModel/Inputs/elevation_SRTMplus_0_01.tif');

%% MAKE PERSISTENT IN THE FUNCTION (so it only computes once)
%create a grid containing absolute latitudes
AbsLat = abs(gridcreate(0.01)); %grid center absolute latitude (0.01-degree)

%path to never snow mask folder, contains indication of areas that have
%less than a 0.1% probability of snow cover for the given D_Tmax value
nosnow_mask_path = '/Volumes/GRA_Data_Backup/POSTDOC_UNH/SCA_CLIMATOLOGY/GeneralSCApModel/Inputs/NeverSnowMasks/';

%load in the hyperplanes defining 99% snow probability (snow) vs. 1% snow probability (no-snow) (only applied over locations with missing data)
%these were learned by summarizing the probability of snow by latitude,
%elevation, and D_Tmax (days since max temp, seasonality) bands across all
%MODIS SCA observations over the 2000-2022 period of record
SCAp_99_function = load('/Volumes/GRA_Data_Backup/POSTDOC_UNH/SCA_CLIMATOLOGY/GeneralSCApModel/Models/SnowMaskingHyperplanes/SCAp_99.mat');
SCAp_1_function = load('/Volumes/GRA_Data_Backup/POSTDOC_UNH/SCA_CLIMATOLOGY/GeneralSCApModel/Models/SnowMaskingHyperplanes/SCAp_1.mat');
SCAp_99_function = SCAp_99_function.SCAp_99_function;
SCAp_1_function = SCAp_1_function.SCAp_1_function;


%% Iterate Across all SCA Files/Days & perform gap filling steps
%loop through all days
iter = 0;
for i = 2:98%length(files)   

    tic
    %get file name and load in SCA
    file = files(i).name; disp(file)
    [SCA_i,R] = readgeoraster([filepath file]);

    %report initial % of missing data
    disp('Current missing data')
    disp([num2str((sum(SCA_i(:) == 255)/648000000) * 100,3) '%'])

    %set below 85S to always snow (speeds processing)
    SCA_i(17500:18000,:) = 1;
  
    %{
    %set land values between 70S and 80S to no data (will fill with probabilistic snow vs. no snow -> should be largely SNOW)
    D = SCA_i(16000:17000,:);
    idx = (D == 0);
    D(idx) = 255;
    SCA_i(16000:17000,:) = D;
    %}
    
    %report % of missing data after Antarctica filling (step 1)
    disp([num2str((sum(SCA_i(:) == 255)/648000000) * 100,3) '%'])

    %select and use pre-defined 'never snow' mask for the given date (snow line based on the probability of snow being <0.1%)
    NoSnowMask = load([nosnow_mask_path 'NoSnowMask_SCAp0_01Less_' file(26:30) '.mat']);
    NoSnowMask = NoSnowMask.no_snow_mask;
    SCA_i(NoSnowMask == 1) = 0; %apply mask, ensure masked areas are set as land pixels

    %report % of missing data after applying no snow mask (step 2)
    disp([num2str((sum(SCA_i(:) == 255)/648000000) * 100,3) '%'])

    % (1) Using snow line to gap-fill
    % (1a) identify and load in appropriate DaysfromTmax grid
    D_Tmax_i = load([D_Tmax_path 'DaysFromTmax_' file(26:30) '.mat']);
    D_Tmax_i = D_Tmax_i.A;

    % (1b) identify no data cells and apply snow line masking (>99% snow probability = snow & <1% = no snow)
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
    disp([num2str((sum(SCA_i(:) == 255)/648000000) * 100,3) '%'])
    
    %final gap filling approach
    if method == 1 && iter > 0 %must occur as second timestep and beyond

        % (1d) fill remaining missing cells using previous timestep
        % (dependent on having a fully gap-filled previous date)
        % It is important to ensure that this previous date had minimal
        % missing data - APPLY IN POST PROCESSING TO CORRECT ERRANT FILES
        SCA_GF = readgeoraster([savepath savefile]);
        ind = (SCA_i == 255);
        SCA_i(ind) = SCA_GF(ind);

    %this method will be used on the first loop regardless of the selected
    %method
    elseif method == 2 || iter == 0 %slower, especially in cases of a lot of missing data

        % (1d) fill remaining missing cells using 33% rule (if >33% of valid surrounding pixels are snow, is considered snow)
        % similar sliding window approach as: Riggs, Hall, Vuyovich, and GiGirolamo 2022
        w = 10; %window size is 0.1 x 0.1-degree (they used 0.25 x 0.25-degree)
        [row,col] = find(SCA_i == 255); %locate all no data cells by row and col numbers
        for ii = 1:length(row)

            %identify individual pixel
            row_ii = row(ii);
            col_ii = col(ii);

            %create lists of rows/columns to extract (ensure only valid indices)
            row_list = (row_ii - w/2):(row_ii + w/2); row_list = row_list(row_list > 0 & row_list <= 18000);
            col_list = (col_ii - w/2):(col_ii + w/2); col_list = col_list(col_list > 0 & col_list <= 36000);

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
        disp('error')
        savefile = ['_' file(1:9) 'GF_' file(10:end)]; %will save file with leading underscore
    end

    %save data
    geotiffwrite([savepath savefile],SCA_i,R);
    toc

    iter = iter + 1; %track number of completed iterations

end
