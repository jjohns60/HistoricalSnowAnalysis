function [NsnowCycles,Pvalid] = getSnowOnOff(SCE,t,valid_snow,valid_land,type)
%getSnowOnOff Calculates the number of snow on/off cycles using snow cover
%   data. Pvalid returns the proportion of periods in which there was valid
%   data while NsnowCycles stores the total number of on -> off cycles. Counts
%   every time there is a change from snow on to snow off, NOT off -> on. This
%   provides a better calculation of the total number of unique snow on
%   periods

%for using single 3D data structure
if type == 1
    %loop through and return number of valid obs per location, and sum of snow
    %on observations
    NsnowCycles = zeros(size(SCE(:,:,1)));
    SC_valid = zeros(size(SCE(:,:,1)));
    L = length(t);
    disp(['Processing ' num2str(L) ' timesteps....'])
    for i = 1:L
        disp([datestr(t(i)) ': ' num2str(i) '/' num2str(L)])
        SCE_i = double(SCE(:,:,i));

        if i == 1
            %stores locations with snow cover
            SCE_prev = (SCE_i == valid_snow);
            
        else

            %checks locations that previously had snow and currently do not
            snow_on_off = (SCE_prev) & (SCE_i == valid_land);

            %keep running sum of cycles
            NsnowCycles = NsnowCycles + double(snow_on_off);

            %reset previous to current time step and continue
            SCE_prev = (SCE_i == valid_snow);
        end

        %keep count of valid values
        SC_valid = SC_valid + double(SCE_i == valid_snow | SCE_i == valid_land);

    end

%for using paths
elseif type == 2
%loop through and return number of valid obs per location, and sum of snow
    %on observations
    d = ncread(SCE{1},'IMS_Surface_Values');
    NsnowCycles = zeros(size(d));
    SC_valid = zeros(size(d));
    L = length(t);
    disp(['Processing ' num2str(L) ' timesteps....'])
    for i = 1:L
        disp([datestr(t(i)) ': ' num2str(i) '/' num2str(L)])
        SCE_i = ncread(SCE{i},'IMS_Surface_Values');

        if i == 1
            %stores locations with snow cover
            SCE_prev = (SCE_i == valid_snow);
            
        else

            %checks locations that previously had snow and currently do not
            snow_on_off = (SCE_prev) & (SCE_i == valid_land);

            %keep running sum of cycles
            NsnowCycles = NsnowCycles + uint8(snow_on_off);

            %reset previous to current time step and continue
            SCE_prev = (SCE_i == valid_snow);
        end

        %keep count of valid values
        SC_valid = SC_valid + double(SCE_i == valid_snow | SCE_i == valid_land);

    end
end
Pvalid = SC_valid ./ L;



end