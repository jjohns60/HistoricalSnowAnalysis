function [SCA,t] = convertNDSItoSCA(filepath,savepath,method,datapath)
%convertNDSItoSCA Takes an input path of MODIS NDSI (must be 0 - 100) and 
% converts it to a binary snow/no snow determination
%
% INPUTS
% filepath - the full path to the file of interest (MODIS Terra SCA)
% savepath - local directory of which to save binarized SCA outputs
% method - 1: universal low NDSI threshold (0.1 or 10)
%          2: universal high NDSI threshold (0.4 or 40)
%          3: universal moderate NDSI threshold (0.3 or 30)
%          4: use land cover class variant universal thresholds, determined using Otsu's method
%          5: using the legacy MODIS C5 fSCA formula, and a fSCA threshold of 0.5, for MODIS Terra 
%          6: uses the Klein 1998 method which combines NDVI & NDSI       
% datapath - for method 4, this must include the path on your local machine
%            to a folder containing MODIS land cover (MCD12C1_v061)
%            for method 6, this must include the path on your local machine
%            to MODIS NDVI data (MOD13C2_v061)
%
% Note: The literature supporting each of these approaches is included
% within each relevant section of this script. Areas with no data are
% returned using the fill value of 255 in the 8-bit GeoTiffs

if nargin == 0
    error('missing arguments')
elseif nargin == 2 %default to simple application of NDSI 0.1 threshold
    method = 1;
elseif nargin == 3
    if method == 4 || method == 6
        error('requires additional data input, see function description')
    end    
elseif nargin > 4
    error('too many arguments')
end

%load in file NDSI (valid range 0 - 100)
[NDSI,R] = readgeoraster(filepath);
no_data_idx = NDSI > 100;

%get date from file
t = split(filepath,'/');
t = t{end};
t = datetime(str2double(t(10:13)),str2double(t(15:16)),str2double(t(18:19)));


%method 1 using the simple low threshold at 10 (or NDSI >= 0.1)
%Yin et al., 2013; Zhang et al., 2019; Riggs, Hall, Vuyovich, and GiGirolamo 2022
if method == 1
    
    SCA = uint8(NDSI >= 10);
    SCA(no_data_idx) = 255;
    
    %to visualize binary result (comment in)
    %figure; imagesc(SCA,[0 2])

%method 2 using a simple high threshold at 40 (NDSI >= 0.4); NO LONGER RECOMMENDED
%Hall et al., 2002; and other early MODIS product release/ATBD (based on Dozier 1989)
elseif method == 2

    SCA = uint8(NDSI >= 40);
    SCA(no_data_idx) = 255;  

%method 3 using moderate threshold of 30 (NDSI >= 0.3)
%Riggs, Hall, and Roman 2015, 2019 (VIIRS ATBD)
elseif method == 3

    SCA = uint8(NDSI >= 30);
    SCA(no_data_idx) = 255;

%method 4 using OTSU's method thresholds for each land cover type
%Yin et al., 2013; Otsu 1979
elseif method == 4

    %(1) Access IGBP land cover for corresponding year (MCD12C1 v006 -- for corresponding year, 0.05 -> resampled to 0.01)
    LC_files = dir([datapath '*.hdf']);
    LC_files = {LC_files.name};

    %loop through and return list of corresponding years
    LC_year = zeros(size(LC_files));
    for i = 1:length(LC_year)
        f = LC_files{i};
        LC_year(i) = str2double(f(10:13));
    end

    %compare to year of file and load in closest year
    dif = abs(year(t) - LC_year);
    idx = dif == min(dif);
    LC_file = LC_files{idx};

    %IGBP land cover (17-classes)
    LC = hdfread([datapath LC_file],'Majority_Land_Cover_Type_1');
    LC = imresize(LC,5,"nearest"); %resample using nearest interpolate to convert to 0.01-degree grid

    %create binary SCA placeholder
    SCA = zeros(size(NDSI),"uint8");

    %(2) Apply (pre-learned via Otsu) land cover based NDSI thresholds
    LC_class = 2; %Evergreen broadleaf
    NDSI_thresh = 0.10; %threshold for the above class
    idx = (LC == LC_class) & (NDSI >= NDSI_thresh*100);
    SCA(idx) = 1;

    LC_class = 13; %Urban
    NDSI_thresh = 0.24; %threshold for the above class
    idx = (LC == LC_class) & (NDSI >= NDSI_thresh*100);
    SCA(idx) = 1;

    LC_class = [1 3 4 5]; %Needleaf/Mixed Forests & decidious broadleaf
    NDSI_thresh = 0.25; %threshold for the above class
    idx = (LC == LC_class(1) | LC == LC_class(2) | LC == LC_class(3) | LC == LC_class(4)) & (NDSI >= NDSI_thresh*100);
    SCA(idx) = 1;

    LC_class = [8 12 14]; %Woody savannas and croplands
    NDSI_thresh = 0.26; %threshold for the above class
    idx = (LC == LC_class(1) | LC == LC_class(2) | LC == LC_class(3)) & (NDSI >= NDSI_thresh*100);
    SCA(idx) = 1;
    
    LC_class = 9; %Savannas
    NDSI_thresh = 0.29; %threshold for the above class
    idx = (LC == LC_class) & (NDSI >= NDSI_thresh*100);
    SCA(idx) = 1;

    LC_class = 6; %Closed shrublands
    NDSI_thresh = 0.30; %threshold for the above class
    idx = (LC == LC_class) & (NDSI >= NDSI_thresh*100);
    SCA(idx) = 1;

    LC_class = [7 10 11 16]; %Sparsely vegetated classes
    NDSI_thresh = 0.31; %threshold for the above class
    idx = (LC == LC_class(1) | LC == LC_class(2) | LC == LC_class(3) | LC == LC_class(4)) & (NDSI >= NDSI_thresh*100);
    SCA(idx) = 1;

    LC_class = 15; %Snow and ice
    NDSI_thresh = 0.40; %threshold for the above class
    idx = (LC == LC_class) & (NDSI >= NDSI_thresh*100);
    SCA(idx) = 1;
    
    %masks locations with water or no valid data
    SCA(no_data_idx) = 255;


%method 5 using the legacy MODIS C5 fSCA formula, for MODIS Terra   
%Riggs, Hall, and Roman 2016
elseif method == 5

    SCA = (-0.01  + 1.45 * single(NDSI)/100); %modified to take input NDSI 0-100
    SCA = uint8(SCA >= 0.5); %if fSCA is >=50% consider as snow covered
    SCA(no_data_idx) = 255;

    %figure; imagesc(SCA,[0 2])

%method 6 the Klein method and NDVI
%Klein et al., 1998; Riggs et al., 2006
elseif method == 6

    %GET NDVI (5km) -> access from folder, nearest month to input date
    NDVI_files = dir([datapath '*.hdf']);
    NDVI_files = {NDVI_files.name};

    %loop through and return list of corresponding datetimes
    NDVI_t = NaT(size(NDVI_files));
    for i = 1:length(NDVI_t)
        f = NDVI_files{i};
        NDVI_t(i) = datetime(f(10:16),"InputFormat","uuuuDDD");
    end

    %identify file to download, nearest NDVI value to input date
    dif = t - NDVI_t;
    idx = dif == min(abs(dif));
    NDVI_file = NDVI_files{idx};
    NDVI = single(hdfread([datapath NDVI_file],'CMG 0.05 Deg Monthly NDVI'));
    NDVI(NDVI > 10000 | NDVI < -2000) = NaN; %constrain to valid range
    NDVI = NDVI./10000; %apply scale factor
    NDVI = imresize(NDVI,5,"nearest"); %resample using nearest interpolate to convert to 0.01-degree grid

    %classification based on Klein 1998 approach
    idx1 = NDSI >= 40; %snow w/ confidence
    idx2 = (NDSI >= 10 & NDSI < 40) & ...
           (NDVI >= -0.5*(single(NDSI)/100) + 0.3) & ...
           (NDVI <= -4.5 *(single(NDSI)/100).^2 + 4.75 * (single(NDSI)/100) - 0.18); %NDVI based algorithm, additional snow pixels
    
    SCA = uint8((idx1 | idx2) & ~no_data_idx);
    SCA(no_data_idx) = 255;

%PENDING/FUTURE: method 7 - suggested by Rittger (and Hao), using spectral unmixing
%is better than empirical thresholds
%Rittger, Painter, and Dozier 2013; Hao et al., 2018
elseif method == 7

    %spectral unmixing approach, would require more than just NDSI

end

%save binarized SCA data to the indicated savepath
file = split(filepath,'/');
file = file{end};
file = [file(1:end-4) '_SCA_method' num2str(method) '.tif'];
geotiffwrite([savepath file],SCA,R);

end