function Ncyc = getNcyc(path,savepath,date_range,valid_snow,valid_land,min_count)
%getNcyc Returns the number of snow covered (on) to non-snow covered (off) 
% cycles considering an input date range (see 'date_range'). Requires an 
% input folder of gridded snow cover data. Currently, works for 
% Tiff/GeoTiff files.
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
% Ncyc - for 2D spatial inputs, will return .tif files with snow cover
%   cycle counts, saved to the specified 'savepath.' Will also store the
%   most recently calculated Ncyc as a variable. 255 is assigned as the no 
%   data value
%
%
% Note: the larger the dataset, the slower processing (e.g., takes ~15 
% minutes to process a single snow year at ~1km)


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
            %loop through all files appending snow cycle counts to
            %appropriate matrix
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
                    
                    %stores number of snow on off cycles
                    Ncyc = zeros([size(D,1) size(D,2)],'uint8');

                    %create mask identifying locations of no data pixels
                    snow_idx = (D == valid_snow);
                    land_idx = (D == valid_land);
                    mask = ~(snow_idx | land_idx);

                    %identify locations with snow cover in the first
                    %timestep
                    SCA_prev = snow_idx;
              
                %increment counts
                else

                    %checks locations that previously had snow and currently do not
                    snow_on_off = (SCA_prev) & (D == valid_land);

                    %reset previous to current time step and continue
                    SCA_prev = (D == valid_snow);

                    %check file_idx to determine which values to increment
                    if (var == 1) && (NH_exist == 1) % == 1, only increment NH counts

                        if SH_exist == 1
                            %set all values below equator to 0
                            snow_on_off(SH_rows(1):SH_rows(2),:) = 0;
                        end

                        %keep running sum of cycles, only increment NH
                        Ncyc = Ncyc + uint8(snow_on_off);
                                               
                    elseif (var == 2) && (SH_exist == 1) % == 2, only increment SH counts

                        if NH_exist == 1
                            %set all values above equator to 0
                            snow_on_off(NH_rows(1):NH_rows(2),:) = 0;
                        end

                        %keep running sum of cycles, only increment SH
                        Ncyc = Ncyc + uint8(snow_on_off);
                        
                    elseif var == 3 % == 3, increment both

                        %increment all locations
                        Ncyc = Ncyc + uint8(snow_on_off);

                    else %if none of the conditions are met, skip to next iteration
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

        %save outputs
        f = split(file_ii,'_');
        fname = ['Ncyc_' f{1} '_' f{2} '_SY' num2str(year_list(i)) '-' num2str(year_list(i + 1)) '_' f{end}];
        Ncyc(mask) = 255;
        geotiffwrite([savepath fname],Ncyc,R);

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
            %loop through all files appending snow cycle counts to
            %appropriate matrix
            for ii = 1:length(files_i)

                tic

                %load in file data and georeferencing information
                file_ii = files_i(ii).name; disp(file_ii);
                [D,R] = readgeoraster([path file_ii]);

                %initialize counting matrices, and increment counts
                if ii == 1

                    %stores number of snow on off cycles
                    Ncyc = zeros([size(D,1) size(D,2)],'uint8');

                    %create mask identifying locations of no data pixels
                    snow_idx = (D == valid_snow);
                    land_idx = (D == valid_land);
                    mask = ~(snow_idx | land_idx);

                    %identify locations with snow cover in the first
                    %timestep
                    SCA_prev = snow_idx;

                else

                    %checks locations that previously had snow and currently do not
                    snow_on_off = (SCA_prev) & (D == valid_land);

                    %reset previous to current time step and continue
                    SCA_prev = (D == valid_snow);

                    %keep running sum of cycles, only increment NH
                    Ncyc = Ncyc + uint8(snow_on_off);

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
        fname = ['Ncyc_' f{1} '_' f{2} '_' num2str(year_list(i)) '_' f{end}];
        Ncyc(mask) = 255;
        geotiffwrite([savepath fname],Ncyc,R);

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

                %stores total number of snow on off cycles
                Ncyc = zeros([size(D,1) size(D,2)],'uint8');

                %create mask identifying locations of no data pixels
                snow_idx = (D == valid_snow);
                land_idx = (D == valid_land);
                mask = ~(snow_idx | land_idx);

                %identify locations with snow cover in the first
                %timestep
                SCA_prev = snow_idx;

            else

                %checks locations that previously had snow and currently do not
                snow_on_off = (SCA_prev) & (D == valid_land);

                %reset previous to current time step and continue
                SCA_prev = (D == valid_snow);

                %keep running sum of cycles, only increment NH
                Ncyc = Ncyc + uint8(snow_on_off);

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
    fname = ['Ncyc_' f{1} '_' f{2} '_' char(d1) '_' char(d2) '_' f{end}];
    Ncyc(mask) = 255;
    geotiffwrite([savepath fname],Ncyc,R);

%catch errors
else
    error("Invalid 'date_range' input")
end

end
