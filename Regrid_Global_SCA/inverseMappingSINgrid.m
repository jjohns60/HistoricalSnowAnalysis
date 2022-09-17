function [lat,lon] = inverseMappingSINgrid(DIM,H,V)
%inverseMappingSINgrid Computes geographic coordinates associated with
%MODIS and VIIRS global SINusoidal grid products
%   
%   Inputs are DIM (defines the width/height in pixels of the input file),
%   H (which defines the horizontal position of the tile), and V (which
%   defines the vertical position of the tile
%   
%   NOTE: Assumes that tiles are square, number of rows == number of cols

%static variables are correct for the 0 -> 35 (cols) and 0 - 17 (rows) projection
R = 6371007.181; %radius of idealized sphere representing Earth (m)
T = 1111950; %height and width of each tile in the projection plane (m)
xmin = -20015109; %western limit of the projection plane (m)
ymax = 10007555; %northern limit of the projection plane (m)

%specify file specific parameters
w = T/DIM; %actual height/width of each grid cell (m)
row_i = repmat((0:DIM-1)',[1 DIM]); %grid with all row numbers
col_j = repmat((0:DIM-1),[DIM 1]); %grid with all column numbers

%compute center of each grid cell on the global sinusoidal grid
x = (col_j + 0.5).*w + H.*T + xmin;
y = ymax - (row_i + 0.5).*w - V.*T;

%get latitude and longitude, in radians, then convert to degrees
lat = y./R;
lon = x./(R.*cos(lat));
lat = lat .* 180./pi;
lon = lon .* 180./pi;

%if out of valid range, make NaN
idx = lat > 90 | lat < -90 | lon > 180 | lon < -180;
lat(idx) = NaN;
lon(idx) = NaN;

end