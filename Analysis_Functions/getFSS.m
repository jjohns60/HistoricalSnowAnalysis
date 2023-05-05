function FSS = getFSS(path,savepath,date_range,valid_snow,valid_land,min_count)
%getFSS Returns the length of the full snow season (FSS) considering an 
% input date range (see 'date_range'). Requires an input folder to gridded 
% snow cover data. Currently, works for Tiff/GeoTiff files.
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
%   timeseries to still return outputs. Also must include a value 
%   indicative of a full data period (default = [0.95 365 (days)])
%
% OUTPUTS:
% FSS - for 2D spatial inputs, will return .tif files with length of the
%   full snow season (FSS) saved to the specified 'savepath.' FSS is
%   defined as the duration (depending on input temporal resolution, i.e.,
%   weekly or daily), from the first occurance of snow cover to the last.
%   Will also store the most recently calculated FSS as a variable. -1 is 
%   assigned as the no data value
%
% Note: the larger the dataset, the slower processing (e.g., takes ~60 
% minutes to process a single global snow year at ~1km)


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
                    
                    %stores the iteration of the first occurance of snow cover
                    first_SCA = zeros([size(D,1) size(D,2)],'int16');
                    %stores the iteration of the last occurance of snow cover
                    last_SCA = zeros([size(D,1) size(D,2)],'int16');

                    %create mask identifying locations of no data pixels
                    snow_idx = (D == valid_snow);
                    land_idx = (D == valid_land);
                    first_SCA_idx = (snow_idx & first_SCA == 0);
                    mask = ~(snow_idx | land_idx);

                    %check file_idx to determine which values to update
                    if (var == 1 && NH_exist == 1) % == 1, only increment NH counts

                        %only update NH values (trim to NH before assigning number)
                        first_SCA_i = first_SCA(NH_rows(1):NH_rows(2),:); %extract only the NH from the matrix storing first snow occurances
                        first_SCA_idx = first_SCA_idx(NH_rows(1):NH_rows(2),:); %extract only the NH from the matrix storing the locations of places with their first snow cover
                        first_SCA_i(first_SCA_idx) = ii; %update the locations with the new number representing the timestep with first snow occurance
                        first_SCA(NH_rows(1):NH_rows(2),:) = first_SCA_i; %apply these back to the original/global data

                        last_SCA_i = last_SCA(NH_rows(1):NH_rows(2),:); %extract only the NH from the matrix storing first snow occurances
                        snow_idx = snow_idx(NH_rows(1):NH_rows(2),:); %extract only the NH from the matrix storing the locations of places with their first snow cover
                        last_SCA_i(snow_idx) = ii; %update the locations with the new number representing the timestep with first snow occurance
                        last_SCA(NH_rows(1):NH_rows(2),:) = last_SCA_i; %apply these back to the original/global data
                                
                    elseif (var == 2 && SH_exist) % == 2, only update SH counts

                        %only update SH values (trim to NH before assigning number)
                        first_SCA_i = first_SCA(SH_rows(1):SH_rows(2),:); %extract only the NH from the matrix storing first snow occurances
                        first_SCA_idx = first_SCA_idx(SH_rows(1):SH_rows(2),:); %extract only the NH from the matrix storing the locations of places with their first snow cover
                        first_SCA_i(first_SCA_idx) = ii; %update the locations with the new number representing the timestep with first snow occurance
                        first_SCA(SH_rows(1):SH_rows(2),:) = first_SCA_i; %apply these back to the original/global data

                        last_SCA_i = last_SCA(SH_rows(1):SH_rows(2),:); %extract only the NH from the matrix storing first snow occurances
                        snow_idx = snow_idx(SH_rows(1):SH_rows(2),:); %extract only the NH from the matrix storing the locations of places with their first snow cover
                        last_SCA_i(snow_idx) = ii; %update the locations with the new number representing the timestep with first snow occurance
                        last_SCA(SH_rows(1):SH_rows(2),:) = last_SCA_i; %apply these back to the original/global data
                         
                        
                    elseif var == 3 % == 3, update both

                        %check for first occurances of snow cover, update
                        first_SCA(first_SCA_idx) = ii;

                        %keep running count of most recent snow cover period
                        last_SCA(snow_idx) = ii;

                    else %if no conditions are met
                        continue  
                    end
              
                %increment counts
                else

                    %identify locations with snow cover
                    snow_idx = (D == valid_snow);
                    first_SCA_idx = (snow_idx & first_SCA == 0);          
                    
                    %check file_idx to determine which values to increment
                    if (var == 1 && NH_exist == 1) % == 1, only increment NH counts

                        %only update NH values (trim to NH before assigning number)
                        first_SCA_i = first_SCA(NH_rows(1):NH_rows(2),:); %extract only the NH from the matrix storing first snow occurances
                        first_SCA_idx = first_SCA_idx(NH_rows(1):NH_rows(2),:); %extract only the NH from the matrix storing the locations of places with their first snow cover
                        first_SCA_i(first_SCA_idx) = ii; %update the locations with the new number representing the timestep with first snow occurance
                        first_SCA(NH_rows(1):NH_rows(2),:) = first_SCA_i; %apply these back to the original/global data

                        last_SCA_i = last_SCA(NH_rows(1):NH_rows(2),:); %extract only the NH from the matrix storing first snow occurances
                        snow_idx = snow_idx(NH_rows(1):NH_rows(2),:); %extract only the NH from the matrix storing the locations of places with their first snow cover
                        last_SCA_i(snow_idx) = ii; %update the locations with the new number representing the timestep with first snow occurance
                        last_SCA(NH_rows(1):NH_rows(2),:) = last_SCA_i; %apply these back to the original/global data         
                                               
                    elseif (var == 2 && SH_exist) % == 2, only increment SH counts

                        %only update SH values (trim to NH before assigning number)
                        first_SCA_i = first_SCA(SH_rows(1):SH_rows(2),:); %extract only the NH from the matrix storing first snow occurances
                        first_SCA_idx = first_SCA_idx(SH_rows(1):SH_rows(2),:); %extract only the NH from the matrix storing the locations of places with their first snow cover
                        first_SCA_i(first_SCA_idx) = ii; %update the locations with the new number representing the timestep with first snow occurance
                        first_SCA(SH_rows(1):SH_rows(2),:) = first_SCA_i; %apply these back to the original/global data

                        last_SCA_i = last_SCA(SH_rows(1):SH_rows(2),:); %extract only the NH from the matrix storing first snow occurances
                        snow_idx = snow_idx(SH_rows(1):SH_rows(2),:); %extract only the NH from the matrix storing the locations of places with their first snow cover
                        last_SCA_i(snow_idx) = ii; %update the locations with the new number representing the timestep with first snow occurance
                        last_SCA(SH_rows(1):SH_rows(2),:) = last_SCA_i; %apply these back to the original/global data
                         
                    elseif var == 3 % == 3, increment both

                        %check for first occurances of snow cover, update
                        first_SCA(first_SCA_idx) = ii;

                        %keep running count of most recent snow cover period
                        last_SCA(snow_idx) = ii;
                        
                    else %if no conditions are met
                        continue
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

        %compute FSS length, inclusive
        FSS = (last_SCA - first_SCA) + 1;

        %save outputs
        f = split(file_ii,'_');
        fname = ['FSS_' f{1} '_' f{2} '_SY' num2str(year_list(i)) '-' num2str(year_list(i + 1)) '_' f{end}];
        FSS(mask) = -1;
        geotiffwrite([savepath fname],FSS,R);

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

                %initialize counting matrices, and increment counts
                if ii == 1

                    %stores the iteration of the first occurance of snow cover
                    first_SCA = zeros([size(D,1) size(D,2)],'int16');
                    %stores the iteration of the last occurance of snow cover
                    last_SCA = zeros([size(D,1) size(D,2)],'int16');

                    %create mask identifying locations of no data pixels
                    snow_idx = (D == valid_snow);
                    land_idx = (D == valid_land);
                    mask = ~(snow_idx | land_idx);

                    %identify locations with snow cover in the first
                    %timestep, update the matrix
                    first_SCA(snow_idx) = ii;
                    %keep running count of most recent snow cover period
                    last_SCA(snow_idx) = ii;

                else

                    %identify locations with snow cover
                    snow_idx = (D == valid_snow);
                    first_SCA_idx = (snow_idx & first_SCA == 0);

                    %check for first occurances of snow cover, update
                    first_SCA(first_SCA_idx) = ii;

                    %keep running count of most recent snow cover period
                    last_SCA(snow_idx) = ii;

                end

                toc

            end
        else
            %if there are not enough valid data frames, will skip to to the
            % next iteration without saving
            continue
        end

        %compute FSS length, inclusive
        FSS = (last_SCA - first_SCA) + 1;

        %save outputs
        f = split(file_ii,'_');
        fname = ['FSS_' f{1} '_' f{2} '_' num2str(year_list(i)) '_' f{end}];
        FSS(mask) = -1;
        geotiffwrite([savepath fname],FSS,R);

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

            %initialize counting matrices, and increment counts
            if i == 1

                %stores the iteration of the first occurance of snow cover
                first_SCA = zeros([size(D,1) size(D,2)],'int16');
                %stores the iteration of the last occurance of snow cover
                last_SCA = zeros([size(D,1) size(D,2)],'int16');

                %create mask identifying locations of no data pixels
                snow_idx = (D == valid_snow);
                land_idx = (D == valid_land);
                mask = ~(snow_idx | land_idx);

                %identify locations with snow cover in the first
                %timestep, update the matrix
                first_SCA(snow_idx) = i;
                %keep running count of most recent snow cover period
                last_SCA(snow_idx) = i;

            else
                
                %identify locations with snow cover    
                snow_idx = (D == valid_snow);
                first_SCA_idx = (snow_idx & first_SCA == 0);

                %check for first occurances of snow cover, update
                first_SCA(first_SCA_idx) = i;

                %keep running count of most recent snow cover period
                last_SCA(snow_idx) = i;

            end

            toc

        end

    else
        error('No valid SCA files found within the input date range')
    end

    %compute FSS length, inclusive
    FSS = (last_SCA - first_SCA) + 1;

    %save outputs
    d1 = date_range(1); d1.Format = 'uuuu-MM-dd';
    d2 = date_range(2); d2.Format = 'uuuu-MM-dd';
    f = split(file_i,'_');
    fname = ['FSS_' f{1} '_' f{2} '_' char(d1) '_' char(d2) '_' f{end}];
    FSS(mask) = -1;
    geotiffwrite([savepath fname],FSS,R);

%catch errors
else
    error("Invalid 'date_range' input")
end

end
