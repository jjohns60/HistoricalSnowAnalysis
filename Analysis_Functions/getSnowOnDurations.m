function [snow_on_durations,Pvalid] = getSnowOnDurations(SCE,t,valid_snow,valid_land,type)
%getSnowOnDurations Takes an input 3D matrix of snow cover information
%(SCE), a list of times (t), the values corresponding to a snow
%classification (valid_snow), and the value corresponding to a non-snow
%land classification (valid_land). The function returns a list of snow on
%durations for each location in the SCE grid, converting it from 3D into a
%2D matrix.
%
% OUTPUTS include a cell array of lists of snow durations and a grid
% showing the proportion of valid periods in the input data by location


if type == 1

    %(1) create arrays to store list of consecutive snow periods
    valid_data_periods = zeros(size(SCE,1),size(SCE,2)); %store total number of timesteps with valid (snow or no-snow data)

    %to keep running counts of consecutive snow or snow off days
    snow_days_sum = zeros(size(SCE,1),size(SCE,2));

    %to store lists of snow on or off durations
    snow_on_durations = cell([size(SCE,1) size(SCE,2)]);

    %number of timesteps
    L = size(SCE,3);
    disp(['Processing ' num2str(L) ' timesteps....'])
    for i = 1:L
        disp([datestr(t(i)) ': ' num2str(i) '/' num2str(L)])
        %   (2) convert data to binary snow, no-snow grid
        SCE_i = SCE(:,:,i);
        valid_data_periods = valid_data_periods + double(SCE_i == valid_land | SCE_i == valid_snow);
        SCE_i = (SCE_i == valid_snow);

        %       (3) check if is snow and previous day was snow (case 1)
        %           check if is snow and previous day was no-snow (case 2)
        %           check if is no-snow and previous day was snow (case 3)
        %           check if is no-snow and previous day was no-snow (case 4)
        if i == 1
            %count day of snow/no snow for first iteration
            snow_days_sum = snow_days_sum + double(SCE_i == 1);

            SCE_prev = SCE_i;
            continue
        else
            case_1 = SCE_i == 1 & SCE_prev == 1;
            case_2 = SCE_i == 1 & SCE_prev == 0;
            case_3 = SCE_i == 0 & SCE_prev == 1;
            %case_4 = SCE_i == 0 & SCE_prev == 0;

            %append to summary matrices, case 1 indicates back to back
            %snow days, case 2 indicates the start of a snow on period
            snow_days_sum = snow_days_sum + double(case_1) + double(case_2);

            %if case 3 (was snow on, now snow off), get snow day sum values and append to
            %snow_on_durations list
            [row_case_3,col_case_3] = find(case_3);
            idx = sub2ind(size(case_3),row_case_3,col_case_3);

            %loop through all locations in index and append lists at
            %locations in which case 3 occurs (is pretty fast, even with loop)
            for ii = 1:length(idx)
                snow_on_durations{idx(ii)} = [snow_on_durations{idx(ii)} snow_days_sum(idx(ii))];
            end

            %reset locations of case_3
            snow_days_sum(case_3) = 0;

            %set previous day
            SCE_prev = SCE_i;

            %on last iteration, append all remaining snow day sums to the snow
            %durations cell array
            if i == L
                [row,col] = find(snow_days_sum > 0);
                idx = sub2ind(size(snow_days_sum),row,col);

                %loop through all locations in index and append lists at
                %locations in which case 3 occurs (is pretty fast, even with loop)
                for ii = 1:length(idx)
                    snow_on_durations{idx(ii)} = [snow_on_durations{idx(ii)}, snow_days_sum(idx(ii))];
                end
            end
        end
    end

    %compute map of proportions of valid data (can mask regions with little to no valid data)
    Pvalid = valid_data_periods./L;

%if input is a folder
elseif type == 2

    %number of timesteps
    L = length(t);
    disp(['Processing ' num2str(L) ' timesteps....'])
    for i = 1:L
        disp([datestr(t(i)) ': ' num2str(i) '/' num2str(L)])
        
        d = ncread(SCE{i},'IMS_Surface_Values');
        if i == 1
            %(1) create arrays to store list of consecutive snow periods
            valid_data_periods = zeros(size(d)); %store total number of timesteps with valid (snow or no-snow data)
            %to keep running counts of consecutive snow or snow off days
            snow_days_sum = zeros(size(d));
            %to store lists of snow on or off durations
            snow_on_durations = cell(size(d));
        end

        %(2) convert data to binary snow, no-snow grid
        SCE_i = d;
        valid_data_periods = valid_data_periods + double(SCE_i == valid_land | SCE_i == valid_snow);
        SCE_i = (SCE_i == valid_snow);

        %       (3) check if is snow and previous day was snow (case 1)
        %           check if is snow and previous day was no-snow (case 2)
        %           check if is no-snow and previous day was snow (case 3)
        %           check if is no-snow and previous day was no-snow (case 4)
        if i == 1

            %count day of snow/no snow for first iteration
            snow_days_sum = snow_days_sum + double(SCE_i == 1);

            SCE_prev = SCE_i;
            continue

        else
            case_1 = SCE_i == 1 & SCE_prev == 1;
            case_2 = SCE_i == 1 & SCE_prev == 0;
            case_3 = SCE_i == 0 & SCE_prev == 1;
            %case_4 = SCE_i == 0 & SCE_prev == 0;

            %append to summary matrices, case 1 indicates back to back
            %snow days, case 2 indicates the start of a snow on period
            snow_days_sum = snow_days_sum + double(case_1) + double(case_2);

            %if case 3 (was snow on, now snow off), get snow day sum values and append to
            %snow_on_durations list
            [row_case_3,col_case_3] = find(case_3);
            idx = sub2ind(size(case_3),row_case_3,col_case_3);

            %loop through all locations in index and append lists at
            %locations in which case 3 occurs (is pretty fast, even with loop)
            for ii = 1:length(idx)
                snow_on_durations{idx(ii)} = [snow_on_durations{idx(ii)} snow_days_sum(idx(ii))];
            end

            %reset locations of case_3
            snow_days_sum(case_3) = 0;

            %set previous day
            SCE_prev = SCE_i;

            %on last iteration, append all remaining snow day sums to the snow
            %durations cell array
            if i == L
                [row,col] = find(snow_days_sum > 0);
                idx = sub2ind(size(snow_days_sum),row,col);

                %loop through all locations in index and append lists at
                %locations in which case 3 occurs (is pretty fast, even with loop)
                for ii = 1:length(idx)
                    snow_on_durations{idx(ii)} = [snow_on_durations{idx(ii)}, snow_days_sum(idx(ii))];
                end
            end
        end
    end

    %compute map of proportions of valid data (can mask regions with little to no valid data)
    Pvalid = valid_data_periods./L;

%for extra large arrays, split into parts    
elseif type == 3
    %number to parts to partition the data into
    N = 36;
    %number of timesteps
    L = length(t);
    disp(['Processing ' num2str(L) ' timesteps....'])
    I = ncinfo(SCE{1});
    dims = I.Variables(4).Size;
    row_w = dims(1)./(N.^0.5);
    col_w = dims(2)./(N.^0.5);
    row_starts = 1:row_w:dims(2);
    col_starts = 1:col_w:dims(1);
    row_starts = [8193];
    col_starts = [4097];
    for r = 1:length(row_starts)
        row_start = row_starts(r);
        for c = 1:length(col_starts)
            col_start = col_starts(c);
            tic
            for i = 1:L
                disp([datestr(t(i)) ': ' num2str(i) '/' num2str(L)])

                %read data in by section
                d = ncread(SCE{i},'IMS_Surface_Values',[row_start col_start 1],[row_w col_w 1]);
                %figure; imagesc(d)
                
                if i == 1
                    %(1) create arrays to store list of consecutive snow periods
                    valid_data_periods = zeros(size(d)); %store total number of timesteps with valid (snow or no-snow data)
                    %to keep running counts of consecutive snow or snow off days
                    snow_days_sum = zeros(size(d));
                    %to store lists of snow on or off durations
                    snow_on_durations = cell(size(d));
                end

                %(2) convert data to binary snow, no-snow grid
                SCE_i = d;
                valid_data_periods = valid_data_periods + double(SCE_i == valid_land | SCE_i == valid_snow);
                SCE_i = (SCE_i == valid_snow);

                %       (3) check if is snow and previous day was snow (case 1)
                %           check if is snow and previous day was no-snow (case 2)
                %           check if is no-snow and previous day was snow (case 3)
                %           check if is no-snow and previous day was no-snow (case 4)
                if i == 1

                    %count day of snow/no snow for first iteration
                    snow_days_sum = snow_days_sum + double(SCE_i == 1);

                    SCE_prev = SCE_i;
                    continue

                else
                    case_1 = SCE_i == 1 & SCE_prev == 1;
                    case_2 = SCE_i == 1 & SCE_prev == 0;
                    case_3 = SCE_i == 0 & SCE_prev == 1;
                    %case_4 = SCE_i == 0 & SCE_prev == 0;

                    %append to summary matrices, case 1 indicates back to back
                    %snow days, case 2 indicates the start of a snow on period
                    snow_days_sum = snow_days_sum + double(case_1) + double(case_2);

                    %if case 3 (was snow on, now snow off), get snow day sum values and append to
                    %snow_on_durations list
                    [row_case_3,col_case_3] = find(case_3);
                    idx = sub2ind(size(case_3),row_case_3,col_case_3);

                    %loop through all locations in index and append lists at
                    %locations in which case 3 occurs (is pretty fast, even with loop)
                    for ii = 1:length(idx)
                        snow_on_durations{idx(ii)} = [snow_on_durations{idx(ii)} snow_days_sum(idx(ii))];
                    end

                    %reset locations of case_3
                    snow_days_sum(case_3) = 0;

                    %set previous day
                    SCE_prev = SCE_i;

                    %on last iteration, append all remaining snow day sums to the snow
                    %durations cell array
                    if i == L
                        [row,col] = find(snow_days_sum > 0);
                        idx = sub2ind(size(snow_days_sum),row,col);

                        %loop through all locations in index and append lists at
                        %locations in which case 3 occurs (is pretty fast, even with loop)
                        for ii = 1:length(idx)
                            snow_on_durations{idx(ii)} = [snow_on_durations{idx(ii)}, snow_days_sum(idx(ii))];
                        end
                    end
                end
            end

            %compute map of proportions of valid data (can mask regions with little to no valid data)
            Pvalid = valid_data_periods./L;
            save(['SnowDurations_r' num2str(row_start) '_c' num2str(col_start) '.mat'],'snow_on_durations','Pvalid','-v7.3')
            toc
        end
    end

%for 1D timeseries with binary snow on/off values
elseif type == 4

    %to keep running counts of consecutive snow or snow off days
    snow_days_sum = 0;

    %to store lists of snow on or off durations
    snow_on_durations = [];

    %number of timesteps
    L = length(SCE);
    for i = 1:L
        %get individual timestep snow value (already binary)
        SCE_i = SCE(i);

        %       (3) check if is snow and previous day was snow (case 1)
        %           check if is snow and previous day was no-snow (case 2)
        %           check if is no-snow and previous day was snow (case 3)
        %           check if is no-snow and previous day was no-snow (case 4)
        if i == 1
            %count day of snow/no snow for first iteration
            snow_days_sum = snow_days_sum + double(SCE_i);

            SCE_prev = SCE_i;
            continue
        else
            
            if SCE_i == 1 && SCE_prev == 1 %case 1, consecutive snow days
                %keep running total of the sum of consecutive snow
                snow_days_sum = snow_days_sum + 1;
            
            elseif SCE_i == 1 && SCE_prev == 0 %case 2, current day is snow, but previous was not
                %keep running total of the sum of consecutive snow
                snow_days_sum = snow_days_sum + 1;

            elseif SCE_i == 0 && SCE_prev == 1 %case 3, current day is no snow, but previous was
                %append snow duration to durations list
                snow_on_durations = [snow_on_durations snow_days_sum];
                %reset the running sum of consecutive snow days
                snow_days_sum = 0;
            end



            %set previous day
            SCE_prev = SCE_i;

            %on last iteration, append all remaining snow day sums to the snow
            %durations cell array
            if i == L  
               snow_on_durations = [snow_on_durations snow_days_sum];
            end
        end
    end

    %compute map of proportions of valid data (can mask regions with little to no valid data)
    Pvalid = NaN;

end

end