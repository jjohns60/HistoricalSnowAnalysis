function [SCAp,Pvalid] = getSCAp(SCE,t,valid_snow,valid_land,type)
%getSCAp Determines the snow cover proportion gridwise, given an input 3D
%grid of snow cover, corresponding times, and definitions of the variables
%defining snow covered land and non snow covered land

if type == 1
    %loop through and return number of valid obs per location, and sum of snow
    %on observations
    SC_sum = zeros(size(SCE(:,:,1)));
    SC_valid = zeros(size(SCE(:,:,1)));
    L = length(t);
    disp(['Processing ' num2str(L) ' timesteps....'])
    for i = 1:L
        disp([datestr(t(i)) ': ' num2str(i) '/' num2str(L)])
        SCE_i = double(SCE(:,:,i));

        SC_valid = SC_valid + double(SCE_i == valid_snow | SCE_i == valid_land);
        SC_sum = SC_sum + double(SCE_i == valid_snow);

    end

    SCAp = SC_sum ./ SC_valid;
    Pvalid = SC_valid ./ L;

elseif type == 2
    %loop through and return number of valid obs per location, and sum of snow
    %on observations
    d = ncread(SCE{1},'IMS_Surface_Values');
    SC_sum = zeros(size(d));
    SC_valid = zeros(size(d));
    L = length(t);
    disp(['Processing ' num2str(L) ' timesteps....'])
    for i = 1:L
        disp([datestr(t(i)) ': ' num2str(i) '/' num2str(L)])
        SCE_i = ncread(SCE{i},'IMS_Surface_Values');

        SC_valid = SC_valid + double(SCE_i == valid_snow | SCE_i == valid_land);
        SC_sum = SC_sum + double(SCE_i == valid_snow);

    end

    SCAp = SC_sum ./ SC_valid;
    Pvalid = SC_valid ./ L;

%special case for larger data sets to speed processing, return integer
%values only
elseif type == 3
    d = ncread(SCE{1},'IMS_Surface_Values');
    SC_sum = zeros(size(d),'single');
    SC_valid = zeros(size(d),'single');
    L = length(t);
    disp(['Processing ' num2str(L) ' timesteps....'])
    for i = 1:L
        disp([datestr(t(i)) ': ' num2str(i) '/' num2str(L)])
        SCE_i = ncread(SCE{i},'IMS_Surface_Values');

        SC_valid = SC_valid + single(SCE_i == valid_snow | SCE_i == valid_land);
        SC_sum = SC_sum + single(SCE_i == valid_snow);

    end

    SCAp = SC_sum; %count of number of snow on periods
    Pvalid = SC_valid; %count of the number of valid periods

end

end