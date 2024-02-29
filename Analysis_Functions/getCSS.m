function CSS = getCSS(path,savepath,date_range,valid_snow,valid_land,min_count)
%getCSS Returns the length of the core snow season (CSS) considering an 
% input date range (see 'date_range'). Requires an input folder to gridded 
% snow cover data. Currently, works for Tiff/GeoTiff files.
%
% INPUTS:
% path - path to folder containing binary snow cover files
% savepath - location to store output data
% date_range - specify 1 for snow year (default), 2 for calendar year 
%   (not recommended), or input a datetime array in the format: [start_date end_date]
% valid_snow - specify the variable in the .tif file representing snow
%   covered pixels (default == 1)
% valid_land - specify the variable in the .tif file representing bare land
%   pixels, all others will be considered as no data (default == 0)
% min_count - set the minimum number of files (as a proportion) in a given 
%   timeseries to still return outputs. Also must include a value 
%   indicative of a full data period (default = [0.95 365 (days)])
%
% OUTPUTS:
% CSS - for 2D spatial inputs, will return .tif files with length of the
%   core snow season (FSS) saved to the specified 'savepath.' CSS is
%   defined as the duration (depending on input temporal resolution, i.e.,
%   weekly or daily), of the longest unbroken snow covered period within a
%   specified date range (snow year, annual, or user input). This function
%   will also store the most recently calculated CSS as a variable. -1 is 
%   assigned as the no data value
%
% Note: this is the most computationally expensive function. The larger the
% dataset, the slower processing (e.g., takes ~100 minutes to process a 
% single global snow year at ~1km)


%error handling
if nargin < 2
    error('Missing function inputs')
elseif nargin == 2
    date_range = 1;
    valid_snow = 1;
    valid_land = 0;
    min_count = [0.95 365];
elseif nargin == 3
    valid_snow = 1;
    valid_land = 0;
    min_count = [0.95 365];
elseif nargin == 4
    valid_land = 0;
    min_count = [0.95 365];
elseif nargin == 5
    min_count = [0.95 365];
end

%extract timestamps from each file (assumes consistent naming structure)
files = dir([path '*.tif']);
dt = NaT(size(files));
for i = 1:length(dt)
    file = files(i).name;
    f = split(file,'_');
    dt(i) = datetime(str2double(f{end-3}),str2double(f{end-2}),str2double(f{end-1}));
end

%ensure files are sorted in chronological order
[dt,I] = sort(dt);
files = files(I);

%check for the date_range input type and build index specifying file ranges
%to count snow cover durations
if isa(date_range,"double") && date_range == 1 %use the hemisphere specific snow year and re-combine
    %NH: Aug 1, Y -> July 31, Y + 1
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
            %identify index of last NH and SH file
            last_NH = find(file_idx == 1 | file_idx == 3,1,'last'); %NH
            last_SH = find(file_idx == 2 | file_idx == 3,1,'last'); %SH
            %loop through all files appending snow cover counts to the
            %appropriate hemisphere
            for ii = 1:length(files_i)

                tic

                %load in file data and georeferencing information
                file_ii = files_i(ii).name;
                var = file_idx(ii);
                [D,R] = readgeoraster([path file_ii]);

                %identify locations with snow
                snow_idx = (D == valid_snow);

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
                    
                    %to keep track of consecutive snow lengths
                    SLcount = zeros(size(D),'int16');
                    CSS = zeros(size(D),'int16');

                    %store as previous timestep grid for comparison
                    SCA_prev = D;

                    %create placeholder logical indices (only for first loop)
                    idxSnow2NoSnow = (SLcount ~= 0); %all false in first iteration
                    idxNoSnow2Snow = snow_idx; %true in locations starting with snow
                    SLcount_idx = idxNoSnow2Snow; %locations to increment counter

                    %create mask of locations without valid observations
                    mask = ~(snow_idx | D == valid_land);

                    %check file_idx to determine which values to update
                    if (var == 1) && (NH_exist == 1) % == 1, only increment NH counts

                        %only update NH values

                        %keep track of longest run of snow observations between no snow observations
                        %increment counts at all sites during snow on period
                        SLcount_i = SLcount(NH_rows(1):NH_rows(2),:);
                        SLcount_idx_i = SLcount_idx(NH_rows(1):NH_rows(2),:);
                        SLcount_i(SLcount_idx_i) = SLcount_i(SLcount_idx_i) + 1;
                        SLcount(NH_rows(1):NH_rows(2),:) = SLcount_i;

                        %set snow to no snow index to false in the southern
                        %hemisphere
                        if SH_exist == 1
                            idxSnow2NoSnow(SH_rows(1):SH_rows(2),:) = false;
                        end

                        %check if the counts are greater than the current max run, but only
                        % in locations that are at the end of a snow on period
                        if ii == last_NH %special case if CSS runs through end of the time series
                            idx = ((SLcount > CSS) & idxSnow2NoSnow) | ((SLcount > CSS) & D == valid_snow);
                        else %all other cases, update when changing from snow to no snow if new SLcount is longer than existing longest snow on period
                            idx = (SLcount > CSS) & idxSnow2NoSnow;
                        end
                        %if so, update the max no snow run (between snow obs)
                        CSS(idx) = SLcount(idx);
                        %reset counter at all locations that switched back to no snow
                        SLcount(idxSnow2NoSnow) = 0;
                             
                    elseif (var == 2) && (SH_exist == 1) % == 2, only update SH counts

                        %only update SH values

                        %keep track of longest run of snow observations between no snow observations
                        %increment counts at all sites during snow on period
                        SLcount_i = SLcount(SH_rows(1):SH_rows(2),:);
                        SLcount_idx_i = SLcount_idx(SH_rows(1):SH_rows(2),:);
                        SLcount_i(SLcount_idx_i) = SLcount_i(SLcount_idx_i) + 1;
                        SLcount(SH_rows(1):SH_rows(2),:) = SLcount_i;

                        %set snow to no snow index to false in the southern
                        %hemisphere
                        if NH_exist == 1
                            idxSnow2NoSnow(NH_rows(1):NH_rows(2),:) = false;
                        end

                        %check if the counts are greater than the current max run, but only
                        % in locations that are at the end of a snow on period
                        if ii == last_SH %special case if CSS runs through end of the time series
                            idx = ((SLcount > CSS) & idxSnow2NoSnow) | ((SLcount > CSS) & D == valid_snow);
                        else %all other cases, update when changing from snow to no snow if new SLcount is longer than existing longest snow on period
                            idx = (SLcount > CSS) & idxSnow2NoSnow;
                        end
                        %if so, update the max no snow run (between snow obs)
                        CSS(idx) = SLcount(idx);
                        %reset counter at all locations that switched back to no snow
                        SLcount(idxSnow2NoSnow) = 0;

                    elseif var == 3 % == 3, update both

                        %keep track of longest run of snow observations between no snow observations
                        %increment counts at all sites during snow on period
                        SLcount(SLcount_idx) = SLcount(SLcount_idx) + 1;
                        %check if the counts are greater than the current max run, but only
                        % in locations that are at the end of a snow on period
                        if (ii == last_NH || ii == last_SH) %special case if CSS runs through end of the time series
                            idx = ((SLcount > CSS) & idxSnow2NoSnow) | ((SLcount > CSS) & D == valid_snow);
                        else %all other cases, update when changing from snow to no snow if new SLcount is longer than existing longest snow on period
                            idx = (SLcount > CSS) & idxSnow2NoSnow;
                        end
                        %if so, update the max no snow run (between snow obs)
                        CSS(idx) = SLcount(idx);
                        %reset counter at all locations that switched back to no snow
                        SLcount(idxSnow2NoSnow) = 0;
                    
                    else %if no conditions are met
                        continue
                    end
              
                %increment counts
                else

                    %identify locations going from snow to no snow
                    idxSnow2NoSnow = (SCA_prev == 1) & (D == 0);
                    %identify locations going from no snow to snow
                    idxNoSnow2Snow = (SCA_prev == 0) & (snow_idx);

                    %start counting at a location just switching from no snow
                    %to snow (includes prior iterations that switched) and stop
                    %counts at locations switching from snow to no snow
                    SLcount_idx(idxNoSnow2Snow) = true; %update locations to start counting
                    SLcount_idx(idxSnow2NoSnow) = false; %update locations to stop counting

                    %update previous timestep array
                    SCA_prev = D;           
                    
                    %check file_idx to determine which values to increment
                    if (var == 1) && (NH_exist == 1) % == 1, only increment NH counts

                        %only update NH values

                        %keep track of longest run of snow observations between no snow observations
                        %increment counts at all sites during snow on period
                        SLcount_i = SLcount(NH_rows(1):NH_rows(2),:);
                        SLcount_idx_i = SLcount_idx(NH_rows(1):NH_rows(2),:);
                        SLcount_i(SLcount_idx_i) = SLcount_i(SLcount_idx_i) + 1;
                        SLcount(NH_rows(1):NH_rows(2),:) = SLcount_i;

                        %set snow to no snow index to false in the southern
                        %hemisphere
                        if SH_exist == 1
                            idxSnow2NoSnow(SH_rows(1):SH_rows(2),:) = false;
                        end

                        %check if the counts are greater than the current max run, but only
                        % in locations that are at the end of a snow on period
                        if ii == last_NH %special case if CSS runs through end of the time series
                            idx = ((SLcount > CSS) & idxSnow2NoSnow) | ((SLcount > CSS) & D == valid_snow);
                        else %all other cases, update when changing from snow to no snow if new SLcount is longer than existing longest snow on period
                            idx = (SLcount > CSS) & idxSnow2NoSnow;
                        end
                        %if so, update the max no snow run (between snow obs)
                        CSS(idx) = SLcount(idx);
                        %reset counter at all locations that switched back to no snow
                        SLcount(idxSnow2NoSnow) = 0;

                                               
                    elseif (var == 2) && (SH_exist == 1) % == 2, only increment SH counts

                        %only update SH values

                        %keep track of longest run of snow observations between no snow observations
                        %increment counts at all sites during snow on period
                        SLcount_i = SLcount(SH_rows(1):SH_rows(2),:);
                        SLcount_idx_i = SLcount_idx(SH_rows(1):SH_rows(2),:);
                        SLcount_i(SLcount_idx_i) = SLcount_i(SLcount_idx_i) + 1;
                        SLcount(SH_rows(1):SH_rows(2),:) = SLcount_i;

                        %set snow to no snow index to false in the northern               
                        %hemisphere
                        if NH_exist == 1
                            idxSnow2NoSnow(NH_rows(1):NH_rows(2),:) = false;
                        end

                        %check if the counts are greater than the current max run, but only
                        % in locations that are at the end of a snow on period
                        if ii == last_SH %special case if CSS runs through end of the time series
                            idx = ((SLcount > CSS) & idxSnow2NoSnow) | ((SLcount > CSS) & D == valid_snow);
                        else %all other cases, update when changing from snow to no snow if new SLcount is longer than existing longest snow on period
                            idx = (SLcount > CSS) & idxSnow2NoSnow;
                        end
                        %if so, update the max no snow run (between snow obs)
                        CSS(idx) = SLcount(idx);
                        %reset counter at all locations that switched back to no snow
                        SLcount(idxSnow2NoSnow) = 0;
                 
                    elseif var == 3 % == 3, increment both

                        %keep track of longest run of snow observations between no snow observations
                        %increment counts at all sites during snow on period
                        SLcount(SLcount_idx) = SLcount(SLcount_idx) + 1;
                        %check if the counts are greater than the current max run, but only
                        % in locations that are at the end of a snow on period
                        if (ii == last_NH || ii == last_SH) %special case if CSS runs through end of the time series
                            idx = ((SLcount > CSS) & idxSnow2NoSnow) | ((SLcount > CSS) & D == valid_snow);
                        else %all other cases, update when changing from snow to no snow if new SLcount is longer than existing longest snow on period
                            idx = (SLcount > CSS) & idxSnow2NoSnow;
                        end
                        %if so, update the max no snow run (between snow obs)
                        CSS(idx) = SLcount(idx);
                        %reset counter at all locations that switched back to no snow
                        SLcount(idxSnow2NoSnow) = 0;
                    
                    else %if no conditions are met
                        continue
                    end
                       
                end

                disp(file_ii);
                toc

            end
            
        else
            %create empty array for CSS (to ensure an output argument)
            CSS = NaN;
            %if there are not enough valid data frames, will skip to to the
            % next iteration without saving
            disp(['Not processed, too few timesteps for year: ' num2str(year_list(i))])
            continue
        end

        %save outputs
        f = split(file_ii,'_');
        fname = ['CSS_' f{1} '_' f{2} '_SY' num2str(year_list(i)) '-' num2str(year_list(i + 1)) '_' f{end}];
        CSS(mask) = -1;
        geotiffwrite([savepath fname],CSS,R);

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

                %identify locations with snow
                snow_idx = (D == valid_snow);

                %initialize counting matrices, and increment counts
                if ii == 1

                    %to keep track of consecutive snow lengths
                    SLcount = zeros(size(D),'int16');
                    CSS = zeros(size(D),'int16');

                    %store as previous timestep grid for comparison
                    SCA_prev = D;

                    %create placeholder logical indices (only for first loop)
                    idxSnow2NoSnow = (SLcount ~= 0); %all false in first iteration
                    idxNoSnow2Snow = snow_idx; %true in locations starting with snow
                    SLcount_idx = idxNoSnow2Snow; %locations to increment counter

                    %create NH mask
                    mask = ~(snow_idx | D == valid_land);

                else

                    %identify locations going from snow to no snow
                    idxSnow2NoSnow = (SCA_prev == valid_snow) & (D == valid_land);
                    %identify locations going from no snow to snow
                    idxNoSnow2Snow = (SCA_prev == valid_land) & (snow_idx);

                    %start counting at a location just switching from no snow
                    %to snow (includes prior iterations that switched) and stop
                    %counts at locations switching from snow to no snow
                    SLcount_idx(idxNoSnow2Snow) = true; %update locations to start counting
                    SLcount_idx(idxSnow2NoSnow) = false; %update locations to stop counting

                    %update previous timestep array
                    SCA_prev = D;

                end

                %keep track of longest run of snow observations between no snow observations
                %increment counts at all sites during snow on period
                SLcount(SLcount_idx) = SLcount(SLcount_idx) + 1;
                %check if the counts are greater than the current max run, but only
                % in locations that are at the end of a snow on period
                if ii == length(files_i) %special case if CSS runs through end of the time series
                    idx = ((SLcount > CSS) & idxSnow2NoSnow) | ((SLcount > CSS) & D == valid_snow);
                else %all other cases, update when changing from snow to no snow if new SLcount is longer than existing longest snow on period
                    idx = (SLcount > CSS) & idxSnow2NoSnow;
                end
                %if so, update the max no snow run (between snow obs)
                CSS(idx) = SLcount(idx);
                %reset counter at all locations that switched back to no snow
                SLcount(idxSnow2NoSnow) = 0;

                toc

            end
        else
            %create empty array for CSS (to ensure an output argument)
            CSS = NaN;
            %if there are not enough valid data frames, will skip to to the
            % next iteration without saving
            disp(['Not processed, too few timesteps for year: ' num2str(year_list(i))])
            continue
        end

        %save outputs
        f = split(file_ii,'_');
        fname = ['CSS_' f{1} '_' f{2} '_' num2str(year_list(i)) '_' f{end}];
        CSS(mask) = -1;
        geotiffwrite([savepath fname],CSS,R);

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

            %identify locations with snow
            snow_idx = (D == valid_snow);

            %initialize counting matrices, and increment counts
            if i == 1

                %to keep track of consecutive snow lengths
                SLcount = zeros(size(D),'int16');
                CSS = zeros(size(D),'int16');

                %store as previous timestep grid for comparison
                SCA_prev = D;

                %create placeholder logical indices (only for first loop) 
                idxSnow2NoSnow = (SLcount ~= 0); %all false in first iteration
                idxNoSnow2Snow = snow_idx; %true in locations starting with snow
                SLcount_idx = idxNoSnow2Snow; %locations to increment counter

                %create NH mask
                mask = ~(snow_idx | D == valid_land);

            else
                
                %identify locations going from snow to no snow
                idxSnow2NoSnow = (SCA_prev == 1) & (D == 0);
                %identify locations going from no snow to snow
                idxNoSnow2Snow = (SCA_prev == 0) & (snow_idx);

                %start counting at a location just switching from no snow
                %to snow (includes prior iterations that switched) and stop
                %counts at locations switching from snow to no snow
                SLcount_idx(idxNoSnow2Snow) = true; %update locations to start counting
                SLcount_idx(idxSnow2NoSnow) = false; %update locations to stop counting

                %update previous timestep array
                SCA_prev = D;

            end

            %keep track of longest run of snow observations between no snow observations
            %increment counts at all sites during snow on period
            SLcount(SLcount_idx) = SLcount(SLcount_idx) + 1;
            %check if the counts are greater than the current max run, but only
            % in locations that are at the end of a snow on period
            if i == length(files_i) %special case if SLmax runs through end of the time series
                idx = ((SLcount > CSS) & idxSnow2NoSnow) | ((SLcount > CSS) & D == valid_snow);
            else %all other cases, update when changing from snow to no snow if new SLcount is longer than existing longest snow on period
                idx = (SLcount > CSS) & idxSnow2NoSnow;
            end
            %if so, update the max no snow run (between snow obs)
            CSS(idx) = SLcount(idx);
            %reset counter at all locations that switched back to no snow
            SLcount(idxSnow2NoSnow) = 0;

            toc

        end

    else
        error('No valid SCA files found within the input date range')
    end

    %save outputs
    d1 = date_range(1); d1.Format = 'uuuu-MM-dd';
    d2 = date_range(2); d2.Format = 'uuuu-MM-dd';
    f = split(file_i,'_');
    fname = ['CSS_' f{1} '_' f{2} '_' char(d1) '_' char(d2) '_' f{end}];
    CSS(mask) = -1;
    geotiffwrite([savepath fname],CSS,R);

%catch errors
else
    error("Invalid 'date_range' input")
end

end