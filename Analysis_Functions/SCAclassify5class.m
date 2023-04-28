function C = SCAclassify5class(savepath,CSS,SSP,SP,mask,split_points)
%SCAclassify5class Applys the simplified decision-tree based 5-class snow 
% seasonality classification (see Johnston et al., 2023)
%
% INPUTS:
% savepath - path to save the output raster (if rasters are input),
%   otherwise, will only store output as a MATLAB variable ('C')
% CSS - core snow season length (path to raster or workspace variable),
%   should be normalized (scaled, 0 - 1). Will be done automatically for
%   rasters that contain CSS as total days
% SP - snow persistence (path to raster or workspace variable)
% SSP - snow season persistence (path to raster or workspace variable)
% mask - a logical array specifying which cells to NOT classify. Input is
%   optional
% split_points - decision tree split values, to partition classes. Input as
%   a numeric list: [CSS split 1, CSS split 2, SSP split, SP split]. The
%   complete tree structure is presented in Johnston et al., 2023. 
%   Default values are [0.25 0.50 0.80 0.20], based on Johnston et al.,
%   2023
%
% OUTPUTS:
% C - MATLAB variable containing the snow cover regime classifications. If
%   inputs are rasters, will also save a GeoTiff to 'savepath' titled
%   'SCAclassification.tif', if not, will save as 'SCAclassification.mat'
%

if nargin < 5
    error('Not enough input arguments');
elseif nargin == 5  
    split_points = [0.25 0.50 0.80 0.20];
end

%set split points
CSS1 = split_points(1);
CSS2 = split_points(2);
SSP1 = split_points(3);
SP1 = split_points(4);

%check if the inputs are paths to rasters or MATLAB variables. A mixture of
%each type is also valid as long as they produce the same size data grids
%load in and format data
raster = 0;
if isa(CSS,"char")
    raster = 1;
    [CSS,R] = readgeoraster(CSS);
    CSS = double(CSS);
else
    CSS = double(CSS);
end

if isa(SSP,"char")
    raster = 1;
    [SSP,R] = readgeoraster(SSP);
    SSP = double(SSP);
else
    SSP = double(SSP);
end

if isa(SP,"char")
    raster = 1;
    [SP,R] = readgeoraster(SP);
    SP = double(SP);
else
    SP = double(SP);
end

%Check if neccessary, then convert variables to normalized values (0 - 1)
%Assumes the mask input will mask all fill values, may not work propertly
%if that is not the case and the fill value is larger than the data
M = max(CSS(~mask),[],"all");
if  M > 1
    CSS = CSS./M; %assumed CSS is calculated over a snow year
end

M = max(SSP(~mask),[],"all");
if M > 1
    SSP = SSP./M;
end

M = max(SP(~mask),[],"all");
if M > 1
    SP = SP./M;
end

%verify inputs are same size
idx = isequal(size(CSS),size(SSP)) & isequal(size(SSP),size(SP));
if ~idx
    error('Inputs grid dimensions differ');
end

%set the file savename
if raster == 1
    savepath = [savepath 'SCAclassification.tif'];
else
    savepath = [savepath 'SCAclassification.mat'];
end

%create empty classification array
C = zeros(size(CSS),'int8') - 1;
%ephemeral 1
C(CSS < CSS1 & SSP < SSP1) = 1;
C(CSS < CSS1 & SSP >= SSP1 & SP < SP1) = 1;
%ephemeral 2
%C(SLmax < SLmax1 & SSP >= SSP1 & SP < SP1) = 2;
%transitional
C(CSS < CSS1 & SSP >= SSP1 & SP >= SP1) = 2;
C(CSS >= CSS1 & CSS < CSS2) = 2;
%seasonal
C(CSS >= CSS2) = 3;

%apply masks
C(SP == 0) = 0; %no snow mask
C(SP > 11/12) = 4; %perennial snow mask
%check if an input mask exists, specifying fill areas
if exist('mask','var')
    C(mask) = -1;
end

%save file to savepath
if raster == 1
    geotiffwrite(savepath,C,R); %save as .tif
else
    save(savepath,'C','-v7.3'); %save as .mat
end    

end