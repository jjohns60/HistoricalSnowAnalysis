%% INPUTS
% data_hires: high resolution data input, of which values are extracted from
% lat_hires: higher resolution latitude values corresponding to input data
% lon_hires: higher resolution longitude values corresponding to input data
% lat_lowres: center grid cell latitudes to fit data to
% lon_lowres: center grid cell longitudes to fit data to
% method: 'list' returns cell array with lists containing all data values falling within each cell
%         'average' returns array with single average value, 'SCA' is
%         special function for processing snow cover data
%% OUTPUTS
% DATA: a cell array with a list of all values that fell within the indicated grid
% DATA_cnt: for averaging or SCA function, will return the total number of
% valid values per cell
% MASK_cnt: for SCA function, will return the total number non-nan mask
% values per cell
function [DATA,DATA_cnt,MASK_cnt] = gridvaluesearch(data_hires,lat_hires,lon_hires,lat_lowres,lon_lowres,method)
    
    if nargin == 5
        method = 'list';
    end
   
    %speeds up processing time by only considering in range SCA values
    if strcmp(method,'SCA')
    %replace all locations with out of range data with NaN in lat/lon grids
        idx = data_hires > 100;
        lon_hires(idx) = NaN;
        %only need to convert 1 to NaN since they are considered as data pairs
        %lat_hires(idx) = NaN;
    end
    

    %identify bounds of all target resolution cells
    Xedges = lon_lowres(~isnan(lon_lowres));
    Xedges = unique(Xedges);
    Xedges_i = zeros(length(Xedges)+1,1);
    Xedges_i(1) = Xedges(1) + ((Xedges(1) - Xedges(2))/2);
    Xedges_i(end) = Xedges(end) + ((Xedges(end) - Xedges(end-1))/2);
    Xedges_i(2:end-1) = .5*(Xedges(1:end-1) + Xedges(2:end));
    
    Yedges = lat_lowres(~isnan(lat_lowres));
    Yedges = unique(Yedges);
    Yedges_i = zeros(length(Yedges)+1,1);
    Yedges_i(1) = Yedges(1) + ((Yedges(1) - Yedges(2))/2);
    Yedges_i(end) = Yedges(end) + ((Yedges(end) - Yedges(end-1))/2);
    Yedges_i(2:end-1) = .5*(Yedges(1:end-1) + Yedges(2:end));
    
    %identify all cells that include points
    [~,~,~,Xidx,Yidx] = histcounts2(lon_hires,lat_hires,Xedges_i,Yedges_i);
    
    
    %pre-allocate data array
    if strcmp(method,'list') || strcmp(method,'mode')
        DATA = cell(size(lon_lowres));
    elseif strcmp(method,'average') || strcmp(method,'SCA')
        DATA = zeros(size(lon_lowres));
        DATA_cnt = zeros(size(lon_lowres));
        MASK_cnt = zeros(size(lon_lowres));
    end
    
    %prep data inputs to column vectors
    Xidx = repmat(Xidx(:),size(data_hires,3),1);
    Yidx = repmat(Yidx(:),size(data_hires,3),1);
    Aidx = (1:(length(Xidx)))'; %all positions in data input as col vector
    data_hires = data_hires(:); %convert to column vector
    
    %trim to locations with data
    idx = (Xidx > 0);
    Xidx = Xidx(idx);
    Yidx = Yidx(idx);
    Aidx = Aidx(idx);
        
    %use the identified index of cells with data to extract all of the higher resolution data in the specified cells
    %append these data lists to a cell array (DATA)
    if strcmp(method,'list')
        
        for i = 1:(length(Xidx))
            DATA{Yidx(i),Xidx(i)} = [DATA{Yidx(i),Xidx(i)} data_hires(Aidx(i))];
        end
        
        %correct orientation
        DATA = flipud(DATA);

    elseif strcmp(method,'mode')

        for i = 1:(length(Xidx))
            DATA{Yidx(i),Xidx(i)} = [DATA{Yidx(i),Xidx(i)} data_hires(Aidx(i))];
        end
        
        DATA = cellfun(@mode,DATA);
        DATA = flipud(DATA);

    elseif strcmp(method,'average')
        
        for i = 1:(length(Xidx))
            %get data value at each index
            d = data_hires(Aidx(i));
            if ~isnan(d)
                DATA(Yidx(i),Xidx(i)) = DATA(Yidx(i),Xidx(i)) + d; %sum
                DATA_cnt(Yidx(i),Xidx(i)) = DATA_cnt(Yidx(i),Xidx(i)) + 1; %total number of values at each index
            end
        end

        %correct orientation
        DATA = flipud(DATA./DATA_cnt);
        DATA_cnt = flipud(DATA_cnt);
        MASK_cnt = flipud(MASK_cnt);

    elseif strcmp(method,'SCA')
        
        %loop through every valid data point
        for i = 1:length(Xidx)
            
            %as long as there is a valid data point in the grid cell, it
            %will be counted as valid in the output
            %add each valid value to the appropriate location in the data grid
            DATA(Yidx(i),Xidx(i)) = DATA(Yidx(i),Xidx(i)) + data_hires(Aidx(i)); %sum
            DATA_cnt(Yidx(i),Xidx(i)) = DATA_cnt(Yidx(i),Xidx(i)) + 1; %total number of values at each index
            
        end

        %correct orientation
        DATA = flipud(DATA./DATA_cnt);
        DATA_cnt = flipud(DATA_cnt);
        %stores grid specifying where valid data is located
        MASK_cnt = DATA_cnt > 0; 

    end

end