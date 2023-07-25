function SCD = getSCD(path,savepath,date_range,valid_snow,valid_land,min_count,calc_SP)
%getSCD Returns the snow cover duration (SCD, defined by input variable 
% resolution) considering an input date range. Requires an input folder of 
% gridded snow cover data. Currently, only works for Tiff/GeoTiff files.
%
% INPUTS:
% path - path to folder containing binary snow cover files
% savepath - location to store output data
% date_range - specify 1 for snow year (default), 2 for calendar year, or 
%   input a datetime array in the format: [start_date end_date]
% valid_snow - specify the variable in the .tif file representing snow
%   covered pixels (default == 1)
% valid_land - specify the variable in the .tif file representing bare land
%   pixels, all others will be considered as no data (default == 0)
% min_count - set the minimum number of files (as a proportion) in a given 
%   timeseries to still return outputs. Also must include a value indicative
%   of a full data period (Default = [0.95 365])
% calc_SP - setting to determine if to also compute the snow persistence
%   (SP), which is a unitless proportion of snow cover during the input
%   period. Default is 1 (save SP), set to 0 to only store SCD
%
% OUTPUTS:
% SCD - for 2D spatial inputs, will return .tif files with snow cover
%   duration counts, saved to the specified 'savepath.' Will also store the
%   most recently calculated SCD as a variable. Note: If calc_SP is set to 
%   1, will also save the SP to the 'savepath.' -1 is assigned as the no 
%   data values, while SP is rescaled from 0 - 1 to 0 - 1000
%
%
% Note: the larger the dataset, the slower processing (e.g., takes ~30 
% minutes to process a single year at ~1km)


%error handling
if nargin < 2
    error('Missing function inputs')
elseif nargin == 2
    date_range = 1;
    valid_snow = 1;
    valid_land = 0;
    min_count = [0.95 365];
    calc_SP = 1;
elseif nargin == 3
    valid_snow = 1;
    valid_land = 0;
    min_count = [0.95 365];
    calc_SP = 1;
elseif nargin == 4
    valid_land = 0;
    min_count = [0.95 365];
    calc_SP = 1;
elseif nargin == 5
    min_count = [0.95 365];
    calc_SP = 1;
elseif nargin == 6
    calc_SP = 1;
end

%extract timestamps from each file (assumes consistent naming structure)
files = dir([path '*.tif']);
dt = NaT(size(files));
for i = 1:length(dt)
    file = files(i).name;
    f = split(file,'_');
    dt(i) = datetime(str2double(f{end-5}),str2double(f{end-4}),str2double(f{end-3}));
end

%check for the date_range input type and build index specifying file ranges
%to count snow cover durations
if isa(date_range,"double") && date_range == 1 %use the hemisphere specific snow year and re-combine
    %NH: Aug 1, Y -> Sep 30, Y + 1
    %SH: Mar 1, Y -> Feb 28/29, Y + 1
    
    % identify each snow year within the dt list
    year_list = unique(year(dt));
    for i = 1:(length(year_list) - 1)

        %identify files meeting condition
        SY_idx_NH = dt >= datetime(year_list(i),8,1) & dt < datetime(year_list(i+1),8,1);
        SY_idx_SH = dt >= datetime(year_list(i),3,1) & dt < datetime(year_list(i+1),3,1);
        file_idx = zeros(size(files));
        file_idx(SY_idx_NH & ~SY_idx_SH) = 1; %indicate files to only add to NH counts
        file_idx(SY_idx_SH & ~SY_idx_NH) = 2; %indicate files to only add to SH counts
        file_idx(SY_idx_NH & SY_idx_SH) = 3; %indicate shared files to count for both

        %check if both indices meet the minimum criteria for valid
        %timesteps
        cond1 = sum(SY_idx_NH)/min_count(2) >= min_count(1);
        cond2 = sum(SY_idx_SH)/min_count(2) >= min_count(1);
        if cond1 && cond2
            
            %identify files to read
            files_i = files(file_idx > 0);
            file_idx = file_idx(file_idx > 0);
            %loop through all files appending snow cover counts to the
            %appropriate hemisphere
            for ii = 1:length(files_i)

                tic

                %load in file data and georeferencing information
                file_ii = files_i(ii).name;
                var = file_idx(ii);
                [D,R] = readgeoraster([path file_ii]);

                %identify locations of snow and land
                snow_idx = int16(D == valid_snow);
                land_idx = int16(D == valid_land);

                %initialize counting matrices, and increment counts
                if ii == 1

                    %identify indices of southern vs. northern hemisphere
                    lats = R.LatitudeLimits;
                    lats = lats(2):- (lats(2) - lats(1))/size(D,1):lats(1);
                    NH_rows = (lats > 0);
                    NH_rows = [find(NH_rows,1,'first') find(NH_rows,1,'last')];
                    SH_rows = (lats <= 0);
                    SH_rows = [find(SH_rows,1,'first') find(SH_rows,1,'last') - 1];

                    %check if there is data in the Northern Hemisphere
                    NH_exist = 1;
                    if isempty(NH_rows)
                        NH_exist = 0;
                    end

                    %check if there is data in the Southern Hemisphere
                    SH_exist = 1;
                    if isempty(SH_rows)
                        SH_exist = 0;
                    end
      
                    %stores total number of snow classifications by
                    %hemisphere
                    SCD = zeros([size(D,1) size(D,2)],'int16');

                    %stores total number of valid land or snow classifications
                    VALID = zeros([size(D,1) size(D,2)],'int16');

                    %check file_idx to determine which values to increment
                    if (var == 1) && (NH_exist == 1) % == 1, only increment NH counts
                        
                        if SH_exist == 1
                            %set all values below equator to 0
                            snow_idx(SH_rows(1):SH_rows(2),:) = 0;
                            land_idx(SH_rows(1):SH_rows(2),:) = 0;
                        else
                            SCD = SCD + snow_idx;
                            VALID = VALID + snow_idx + land_idx;
                        end
                                             
                    elseif (var == 2) && (SH_exist == 1) % == 2, only increment SH counts
                        
                        if NH_exist == 1
                            %set all values above equator to 0
                            snow_idx(NH_rows(1):NH_rows(2),:) = 0;
                            land_idx(NH_rows(1):NH_rows(2),:) = 0;
                        else
                            SCD = SCD + snow_idx;
                            VALID = VALID + snow_idx + land_idx;
                        end

                    elseif var == 3 % == 3, increment both

                        %increment all locations
                        SCD = SCD + snow_idx;
                        VALID = VALID + snow_idx + land_idx;
                    
                    else    
                        continue;
                    end
              
                %increment counts
                else

                    %check file_idx to determine which values to increment
                    if (var == 1) && (NH_exist == 1) % == 1, only increment NH counts

                        if SH_exist == 1
                            %set all values below equator to 0
                            snow_idx(SH_rows(1):SH_rows(2),:) = 0;
                            land_idx(SH_rows(1):SH_rows(2),:) = 0;
                        else
                            SCD = SCD + snow_idx;
                            VALID = VALID + snow_idx + land_idx;
                        end

                    elseif (var == 2) && (SH_exist == 1) % == 2, only increment SH counts

                        if NH_exist == 1
                            %set all values above equator to 0
                            snow_idx(NH_rows(1):NH_rows(2),:) = 0;
                            land_idx(NH_rows(1):NH_rows(2),:) = 0;
                        else
                            SCD = SCD + snow_idx;
                            VALID = VALID + snow_idx + land_idx;
                        end
                        
                    elseif var == 3 % == 3, increment both

                        %increment all locations
                        SCD = SCD + snow_idx;
                        VALID = VALID + snow_idx + land_idx;
                    
                    else    
                        continue;    
                    end
                       
                end

            disp(file_ii);
            toc
            
            end
            
        else
            %if there are not enough valid data frames, will skip to to the
            % next iteration without saving
            continue
        end

        %save outputs
        f = split(file_ii,'_');
        fname = ['SCD_' f{1} '_' f{2} '_SY' num2str(year_list(i)) '-' num2str(year_list(i + 1)) '_' f{end}];
        SCD(VALID == 0) = -1;
        geotiffwrite([savepath fname],SCD,R);
            
        %compute SP and save, if user selected
        if calc_SP == 1
            SP = single(SCD)./single(VALID);
            %rescale 0 (0%) - 1000 (100%)
            SP = int16(SP .* 1000);
            SP(VALID == 0) = -1;
            fname = ['SP_' f{1} '_' f{2} '_SY' num2str(year_list(i)) '-' num2str(year_list(i + 1)) '_' f{end}];
            geotiffwrite([savepath fname],SP,R);
        end

    end    

elseif isa(date_range,"double") && date_range == 2 %use the calendar year
    
    % identify all the unique calendar years (January 1 - December 31)
    year_list = unique(year(dt));
    for i = 1:length(year_list)
        
        %identify files falling within the specified calendar year
        Y_idx = dt >= datetime(year_list(i),1,1) & dt <= datetime(year_list(i),12,31);

        %check if the minimum criteria for number of valid timesteps is met
        cond1 = sum(Y_idx)/min_count(2) >= min_count(1);
        if cond1

            %identify files to read
            files_i = files(Y_idx);
            %loop through all files appending snow cover counts to the
            %appropriate hemisphere
            for ii = 1:length(files_i)

                tic

                %load in file data and georeferencing information
                file_ii = files_i(ii).name; disp(file_ii);
                [D,R] = readgeoraster([path file_ii]);

                %identify locations of snow and land
                snow_idx = int16(D == valid_snow);
                land_idx = int16(D == valid_land);

                %initialize counting matrices, and increment counts
                if ii == 1

                    %stores total number of snow classifications
                    SCD = zeros([size(D,1) size(D,2)],'int16');

                    %stores total number of valid land or snow classifications
                    VALID = zeros([size(D,1) size(D,2)],'int16');

                    %increment all locations
                    SCD = SCD + snow_idx;
                    VALID = VALID + snow_idx + land_idx;

                else

                    %increment all locations
                    SCD = SCD + snow_idx;
                    VALID = VALID + snow_idx + land_idx;

                end

                toc
            end
        else
            %if there are not enough valid data frames, will skip to to the
            % next iteration without saving
            continue
        end

        %save outputs
        f = split(file_ii,'_');
        fname = ['SCD_' f{1} '_' f{2} '_' num2str(year_list(i)) '_' f{end}];
        SCD(VALID == 0) = -1;
        geotiffwrite([savepath fname],SCD,R);
            
        %compute SP and save, if user selected
        if calc_SP == 1
            SP = single(SCD)./single(VALID);
            %rescale 0 (0%) - 1000 (100%)
            SP = int16(SP .* 1000);
            SP(VALID == 0) = -1;
            fname = ['SP_' f{1} '_' f{2} '_' num2str(year_list(i)) '_' f{end}];
            geotiffwrite([savepath fname],SP,R);
        end

    end
    
elseif isa(date_range,'datetime') %uses the specified date range
    
    %identify files meeting the condition
    idx = dt >= date_range(1) & dt <= date_range(2);

    %check if there is at least one valid file in the input range
    if sum(idx) > 0

        %identify all files and loop
        files_i = files(idx);
        for i = 1:length(files_i)

            tic

            %load in file data and georeferencing information
            file_i = files_i(i).name; disp(file_i);
            [D,R] = readgeoraster([path file_i]);

            %identify locations of snow and land
            snow_idx = int16(D == valid_snow);
            land_idx = int16(D == valid_land);

            %initialize counting matrices, and increment counts
            if i == 1

                %stores total number of snow classifications
                SCD = zeros([size(D,1) size(D,2)],'int16');

                %stores total number of valid land or snow classifications
                VALID = zeros([size(D,1) size(D,2)],'int16');

                %increment all locations
                SCD = SCD + snow_idx;
                VALID = VALID + snow_idx + land_idx;

            else

                %increment all locations
                SCD = SCD + snow_idx;
                VALID = VALID + snow_idx + land_idx;

            end

            toc

        end

    else
        error('No valid SCA files found within the input date range')
    end

    %save outputs
    d1 = date_range(1); d1.Format = 'uuuu-MM-dd';
    d2 = date_range(2); d2.Format = 'uuuu-MM-dd';
    f = split(file_i,'_');
    fname = ['SCD_' f{1} '_' f{2} '_' char(d1) '_' char(d2) '_' f{end}];
    SCD(VALID == 0) = -1;
    geotiffwrite([savepath fname],SCD,R);

    %compute SP and save, if user selected
    if calc_SP == 1
        SP = single(SCD)./single(VALID);
        %rescale 0 (0%) - 1000 (100%)
        SP = int16(SP .* 1000);
        SP(VALID == 0) = -1;
        fname = ['SP_' f{1} '_' f{2} '_' char(d1) '_' char(d2) '_' f{end}];
        geotiffwrite([savepath fname],SP,R);
    end

%catch errors
else
    error("Invalid 'date_range' input")
end

end