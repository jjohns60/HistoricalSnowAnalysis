function SSP = getSSP(SCD_path,FSS_path,savepath,fill)
%getSSP Computes the snow season persistence (SSP) given input folders
% containing snow cover duration (SCD) and full snow season (FSS)
% rasters. SSP is calculated as the total snow cover duration for a given
% period divided by the time between the first and last occurance of snow
% for the same period. Works for GeoTiff files.
%
% INPUTS:
% SCD_path - path to directory containing SCD rasters
% FSS_path - path to directory containing FSS rasters
% savepath - location to save output files
% fill - specifies fill value within the SCD & FSS dataset
%
% OUTPUTS:
% SSP - is saved locally to savepath for each period with corresponding SCD
%   and FSS files. The most recently calculated SSP can also be returned as
%   a MATLAB variable. Is output as a 16-bit raster, in which 0 - 1 is 
%   rescaled 0 - 1000, the no data/fill value is set as -1
%

%error handling
if nargin < 4
    error('Not enough inputs')
end

SCD_path = '/Users/jjohns/Desktop/SCA/SCD/';
FSS_path = '/Users/jjohns/Desktop/SCA/FSS/';
savepath = '/Users/jjohns/Desktop/SCA/SSP/';
fill = -1;

%get all file names
SCD_files = dir([SCD_path '*.tif']); SCD_files = {SCD_files.name};
FSS_files = dir([FSS_path '*.tif']); FSS_files = {FSS_files.name};

% determine the shorter list of files
SCD_N = length(SCD_files); 
FSS_N = length(FSS_files);
if SCD_N < FSS_N
    N = SCD_N;
    type = 1;
else
    N = FSS_N;
    type = 2;
end

% identify corresponding files (assumes file prefix is only difference)
for i = 1:N
    if type == 1 %if the SCD file list is shorter

        %get search string from SCD file
        SCD_file = SCD_files{i};
        f_search = SCD_file;
        idx = regexp(f_search,'_','once');
        f_search = f_search(idx:end);

        %check in the other list of files
        idx = contains(FSS_files,f_search);

        %get corresponding FSS file
        FSS_file = FSS_files{idx};

    elseif type == 2 %if the FSS list is shorter, or they are the same

        %get search string from FSS file
        FSS_file = FSS_files{i};
        f_search = FSS_file;
        idx = regexp(f_search,'_','once');
        f_search = f_search(idx:end);

        %check in the other list of files
        idx = contains(SCD_files,f_search);

        %get corresponding FSS file
        SCD_file = SCD_files{idx};
    end

    %load in proper files
    [SCD,R] = readgeoraster([SCD_path SCD_file]);
    FSS = readgeoraster([FSS_path FSS_file]);

    %calculate SSP
    mask = (SCD == fill) | (FSS == fill);
    SSP = double(SCD)./double(FSS);

    % Save SSP raster, rescaled to a 16-bit image
    SSP = int16(SSP * 1000);
    SSP(mask) = -1;
    geotiffwrite([savepath 'SSP' f_search],SSP,R);

end
      
end