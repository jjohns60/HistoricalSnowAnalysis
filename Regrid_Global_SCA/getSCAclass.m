function D = getSCAclass(D_cell)
%getSCAclass Determines the dominant classification based on an input list
%of SCA data (see below).
%
%   Uses an input grid from VIIRS or MODIS. In which:
%   0 - 100: are NDSI values
%   237, 239: indicate inland and ocean water, respectively
%   250: cloud cover
%   all else: fill, invalid, missing or unusable data
%
%   Returns a single grids:
%   D = the dominant classification
%       (1) Returns either the mean of all valid NDSI values within the 
%           list if >50% of the data is NDSI
%       (2) If less than 

%store number of data points
L = length(D_cell);

%if over 50% data included are NDSI values
idx = D_cell >= 0 & D_cell <= 100;
if sum(idx)/L > 0.5

    %return mean value
    D = uint8(round(mean(D_cell(idx))));

else
    D = uint8(0);
end

%{
%if the list is empty    
elseif isempty(D_cell)
    D = uint8(0); %will ensure that cells with no observations remain at zero

%if less than 50% of the data are valid NDSI    
else

    %returns mode value (if multiple modes, returns the first one)
    d = mode(D_cell(~idx));
    
    if d == 250 %is cloud cover
        D = uint8(201);
    elseif d == 237 || d == 239 %is water
        D = uint8(240);
    else %all other fill values (no data, masked, etc..)
        D = uint8(255);
    end

end
%}
end