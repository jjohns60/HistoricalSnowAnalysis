function GHCN_stations = GHCN_HTTPS(filepath,variables,search_region,min_length,download)
%GHCN_HTTPS() Accesses climate data daily summaries from the NOAA Global
%Historical Climatology Network (GHCN) database and saves them to a local
%directory
%
%   This function can search for variables included in the GHCN-Daily
%   summary dataset within a specified region:
%   (https://www.ncei.noaa.gov/data/global-historical-climatology-network-daily/doc/GHCND_documentation.pdf)
%
%   INPUTS
%       filepath: path to directory in which files will be saved [required]
%
%       variables: a cell array of variable names, corresponding to those
%           in the GHCN documentation (e.g., {'PRCP','SNOW'}) [required]
%           Note: 'STATION','DATE','LATITUDE','LONGITUDE','ELEVATION',
%           'NAME' fields are automatically included in the output
%           Also note that the variables must be surrounded with brackets {}
%           OR
%           the path to a table or .csv with station names of interest. In 
%           this case, all other inputs should be left null. Requires a
%           column with the name 'STATION' and 'VARIABLE' if input as .csv 
%           or Table
%
%       search_region: a list defining the geographic boundary of which to
%           download stations, in the form [lat_min lat_max lon_min lon_max]
%           [required]
%
%       min_length: defines the minimum period of the climate record
%           required to be downloaded in years (default: 30)
%
%       download: either 1 (yes) or 0 (no), determines if GHCN files are
%           downloaded onto the specified path. OR a special case (2), in
%           which the variables input serves as a list or table containing
%           information of GHCN sites to download
%
%   OUTPUT
%       'stationID'.csv files stored at the indicated filepath and
%       a _stationInfo.csv file which contains a summary of all extracted
%       stations. This is also returned directly to workspace as a table. 
%       Note: The fill value of 9999 may be present within data ouputs,
%       the user is responsible for replacing this value with NULL, NaN, or
%       other before performing analysis with the dataset

%captures special case, in which filepath is the save location, but
%variables is a table including list of files for download
if nargin == 2
    download = 2;
elseif nargin == 3
    min_length = 30;
    download = 0;
end    

%{
%testing
filepath = '/Users/jjohns/Documents/Research/Projects Ongoing/Snow Climatology and Classification Project/GHCN_DATA/';
variables = '/Users/jjohns/Documents/Research/Projects Ongoing/Snow Climatology and Classification Project/GHCN_DATA/_stationInfo.csv';
search_region = [0 90 -180 180];
min_length = 10;
download = 0;
%}

% data is pulled from following path
GHCN_path = 'https://www.ncei.noaa.gov/data/global-historical-climatology-network-daily/access/';

%options for reading data from the web
options = weboptions('ContentReader',@readtable,'Timeout',60);

%list of stations with variable type information located here
GHCN_stations_path = 'https://noaa-ghcn-pds.s3.amazonaws.com/ghcnd-inventory.txt';

%if special case, will just download files in table provided then exit
if download == 2

    %determine data type
    if isa(variables,'table')
        %table with station names of interset
        GHCN_stations = variables;
    elseif ischar(variables)
        %load in as table
        GHCN_stations = readtable(variables);
    end
    
    %stores number of files
    n = height(GHCN_stations);

    %arrays to store long name versions of files
    NAMES = cell(n,1);
    ELEVATIONS = NaN(n,1);

    %display number of stations meeting criteria
    f = waitbar(0,['Starting download of ' num2str(n) ' files']);
    pause(1.5)
    waitbar(0,f,'Starting download.');
    pause(0.5)
    waitbar(0,f,'Starting download..');
    pause(0.5)
    waitbar(0,f,'Starting download...');
    pause(0.5)

    %to correct special case when variable does not exist in the file (even though it did in the station inventory)
    missing_variable = 0;
    flag_file = 0;
    for i = 1:n
        prog = i/n;

        %identify file to download
        file = [GHCN_stations.STATION{i} '.csv'];

        waitbar(prog,f,['Downloading file: ' file ' (' num2str(i) '/' num2str(n) ')']);

        %read file data into MATLAB
        data = webread([GHCN_path file],options);

        %store long NAME information into 'NAMES' array
        NAMES(i) = data.NAME(1);
        ELEVATIONS(i) = data.ELEVATION(1);

        %crop down to include only variables of interest
        T = table(data.STATION,data.DATE,data.LATITUDE,data.LONGITUDE,...
            data.ELEVATION,data.NAME,'VariableNames',{'STATION','DATE',...
            'LATITUDE','LONGITUDE','ELEVATION_m','NAME'});
        
        %determine search variables from input table or .csv
        vars = split(GHCN_stations.VARIABLE{1},',');
      
        %append search variables to table
        for ii = 1:length(vars)

            %identify variable to search for in data table
            var = vars{ii};

            %get all variable names from data table
            data_names = data.Properties.VariableNames;

            %identify location of match
            idx = cellfun(@(x) strcmp(x,var),data_names);
            col = find(idx == 1);

            if col > 0
                %get data from column
                d = data.(col);

                if strcmp(var,'PRCP')
                    d = d/10; %to give mm (stored in database w/o decimals)
                elseif strcmp(var,'SNWD')
                    d = d/10;
                elseif strcmp(var,'TMAX')
                    d = d/10;
                elseif strcmp(var,'TMIN')
                    d = d/10;
                end
                
                %if data includes over 90% NaNs, add flag to file name
                if sum(isnan(d))/length(d) > 0.9
                    flag_file = flag_file + 1;
                    disp(['High number of NaNs for ' var ' in ' file])
                end
                %append to table using identified location
                T = addvars(T,d,'NewVariableNames',var);
            else
                %variable was not found, will not download
                missing_variable = 1;
                disp(['Missing variable in ' file ' will not download'])
            end
        end

        %save file data to download location
        if missing_variable == 0
            %add indicators of high amounts of missing data. Each star
            %indicates a variables with less than 90% complete daily data
            if flag_file > 0
                flag_str = repmat('*',[1 flag_file]);
                FILE = strcat(filepath,flag_str,file);
                writetable(T,FILE);
            else
                writetable(T,[filepath file]);
            end
        end
        %reset missing variable identifier
        missing_variable = 0;
        flag_file = 0;

    end

    %close progress bar
    waitbar(1,f,'Downloads complete');
    pause(1)
    close(f)
    %save .csv file list of all files
    try
        GHCN_stations = addvars(GHCN_stations,NAMES,'NewVariableNames','NAME','Before','LATITUDE');
    catch
        disp('Station inventory not updated with long site names')
    end

    try
        GHCN_stations = addvars(GHCN_stations,ELEVATIONS,'NewVariableNames','ELEVATION_m','After','LONGITUDE');
    catch
        disp('Station inventory not updated with elevations')
    end
    writetable(GHCN_stations,[filepath '_stationInfoMod.csv'])


%if more than just savepath and download list is provided
else

    % (1) Identify list of stations that meet the input criteria
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %load in station data (will timeout after a minute, depends on connection)
    GHCN_stations = webread(GHCN_stations_path,options);

    %add variable names, create new variable storing the record duration
    GHCN_stations.Properties.VariableNames = {'STATION','LATITUDE','LONGITUDE','VARIABLE','START','END'};
    GHCN_stations = addvars(GHCN_stations,GHCN_stations.END - GHCN_stations.START,'NewVariableNames','DURATION');

    %subset to stations with at least 40 year long records
    GHCN_stations = GHCN_stations(GHCN_stations.DURATION >= min_length,:);

    %subset to stations within bounding coordinates
    lat_idx = (GHCN_stations.LATITUDE >= search_region(1) & GHCN_stations.LATITUDE <= search_region(2)); %within latitude range
    lon_idx = (GHCN_stations.LONGITUDE >= search_region(3) & GHCN_stations.LONGITUDE <= search_region(4)); %within longitude range
    GHCN_stations = GHCN_stations((lat_idx & lon_idx),:);

    %subset to stations recording the indicated variables in 'variables'
    v_names = GHCN_stations.VARIABLE;
    var_idx = zeros(length(v_names),1);
    %loop through all search variables in 'variables'
    for i = 1:length(variables)
        %identify variable to search for from variable list
        var = variables{i};

        %create single string containing all variable names
        if i == 1
            v_string = var;
        else
            v_string = strcat(v_string,',',var);
        end

        %combine logicals for different variables
        var_idx = var_idx + cellfun(@(x) contains(x,var),v_names);

        %convert to logical array and crop station list
        if i == length(variables)

            %identify index of stations that contain the variables of interest
            var_idx = var_idx > 0;

            %crops here, to get final station info table
            GHCN_stations = GHCN_stations(var_idx,:);

            %only store stations having all search variables
            station_IDs = GHCN_stations.STATION; %return array of all stations

            %identify unique stations
            [A,idx] = unique(station_IDs,'stable');
            GHCN_stations = GHCN_stations(idx,:);

            %count the occurances of each station
            n_occ = cellfun(@(x) sum(ismember(station_IDs,x)),A,'un',0);

            %check that the file contains all 'variables'
            idx = cellfun(@(x) isequal(x,length(variables)),n_occ);

            %crop to final
            GHCN_stations = GHCN_stations(idx,:);

            %count number of stations meeting criteria
            n = height(GHCN_stations.STATION);

            %update station list with all variables included in file
            C = cell([n 1]);
            C(:) = {v_string};
            GHCN_stations.VARIABLE = C;
        end
    end


    % (2) Download identified files
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if download == 1
        %array to store long name versions of files
        NAMES = cell(n,1);
        ELEVATIONS = cell(n,1);
        %display number of stations meeting criteria
        f = waitbar(0,['Found ' num2str(n) ' GHCN stations meeting search criteria']);
        pause(1.5)
        waitbar(0,f,'Starting download.');
        pause(0.5)
        waitbar(0,f,'Starting download..');
        pause(0.5)
        waitbar(0,f,'Starting download...');
        pause(0.5)
        %to correct special case when variable does not exist in the file (even though it did in the station inventory)
        missing_variable = 0;
        for i = 1:n
            prog = i/n;
            %identify file to download
            file = [GHCN_stations.STATION{i} '.csv'];

            waitbar(prog,f,['Downloading file: ' file ' (' num2str(i) '/' num2str(n) ')']);

            %read file data into MATLAB
            data = webread([GHCN_path file],options);

            %store long NAME information into 'NAMES' array
            NAMES(i) = data.NAME(1);
            %disp(data.ELEVATION(1))
            ELEVATIONS{i} = data.ELEVATION(1);

            %crop down to include only variables of interest
            T = table(data.STATION,data.DATE,data.LATITUDE,data.LONGITUDE,...
                data.ELEVATION,data.NAME,'VariableNames',{'STATION','DATE',...
                'LATITUDE','LONGITUDE','ELEVATION_m','NAME'});

            %append search variables to table
            for ii = 1:length(variables)

                %identify variable to search for in data table
                var = variables{ii};

                %get all variable names from data table
                data_names = data.Properties.VariableNames;

                %identify location of match
                idx = cellfun(@(x) strcmp(x,var),data_names);
                col = find(idx == 1);

                if col > 0
                    d = data.(col);

                    if strcmp(var,'PRCP')
                        d = d/10; %to give mm (stored in database w/o decimals)
                    end

                    %append to table using identified location
                    T = addvars(T,d,'NewVariableNames',var);
                else
                    %variable was not found, will not download
                    missing_variable = 1;
                    disp(['Missing variable in ' file ' will not download'])
                end
            end

            %save file data to download location
            if missing_variable == 0
                writetable(T,[filepath file]);
            end
            %reset missing variable identifier
            missing_variable = 0;

        end

        %close progress bar
        waitbar(1,f,'Downloads complete');
        pause(1)
        close(f)
    end

    %if downloading, add additional descriptive variables
    if download == 1
        %save .csv file list of all files
        GHCN_stations = addvars(GHCN_stations,NAMES,'NewVariableNames','NAME','Before','LATITUDE');
        GHCN_stations = addvars(GHCN_stations,ELEVATIONS,'NewVariableNames','ELEVATION_m','After','LONGITUDE');
    end

    writetable(GHCN_stations,[filepath '_stationInfo.csv'])

end

end