function bounds_new = matchBounds2grid(bounds,TARGET_RES)
%matchBounds2grid Ensures that coordinate boundaries match grid resolution. 
% Assumes the grid is in degrees. 'bounds' should be input table with the
% variable names: "lon_min, lon_max, lat_min, lat_max"
% Method also assumes target grid uses the TARGET_RES (resolution), but
% each input signifies a grid cell center. Thus, this approach finds the
% appropriate grid cell centers.

%determine decimal place value of TARGET_RES input
n = 0;
while (floor(TARGET_RES*10^n)~=TARGET_RES*10^n)
    n = n + 1;
end

%for minimums
lon_min = bounds.lon_min;
lat_min = bounds.lat_min;

%find target decimals it falls between
lon_min1 = round(lon_min,n);
if lon_min > lon_min1 %if the data is rounded to a lower value
    lon_min2 = lon_min1 + TARGET_RES;
    lon_minx = (lon_min1 + lon_min2)/2;
elseif lon_min < lon_min1 %if the data is rounded to a higher value
    lon_min2 = lon_min1 - TARGET_RES;
    lon_minx = (lon_min1 + lon_min2)/2;
else %if rounding does not change the border values, round down
    lon_minx = lon_min1 + TARGET_RES/2;
end

lat_min1 = round(lat_min,n);
if lat_min > lat_min1 %if the data is rounded to a lower value
    lat_min2 = lat_min1 + TARGET_RES;
    lat_minx = (lat_min1 + lat_min2)/2;
elseif lat_min < lat_min1 %if the data is rounded to a higher value
    lat_min2 = lat_min1 - TARGET_RES;
    lat_minx = (lat_min1 + lat_min2)/2;
else %if rounding does not change the border values, round down
    lat_minx = lat_min1 + TARGET_RES/2;
end

%for maximums, round up if on cell border
lon_max = bounds.lon_max;
lat_max = bounds.lat_max;

%find target decimals it falls between
lon_max1 = round(lon_max,n);
if lon_max > lon_max1 %if the data is rounded to a lower value
    lon_max2 = lon_max1 + TARGET_RES;
    lon_maxx = (lon_max1 + lon_max2)/2;
elseif lon_max < lon_max1 %if the data is rounded to a higher value
    lon_max2 = lon_max1 - TARGET_RES;
    lon_maxx = (lon_max1 + lon_max2)/2;
else %if rounding does not change the border values, round up
    lon_maxx = lon_max1 - TARGET_RES/2;
end

lat_max1 = round(lat_max,n);
if lat_max > lat_max1 %if the data is rounded to a lower value
    lat_max2 = lat_max1 + TARGET_RES;
    lat_maxx = (lat_max1 + lat_max2)/2;
elseif lat_max < lat_max1 %if the data is rounded to a higher value
    lat_max2 = lat_max1 - TARGET_RES;
    lat_maxx = (lat_max1 + lat_max2)/2;
else %if rounding does not change the border values, round up
    lat_maxx = lat_max1 - TARGET_RES/2;
end

%update bounds table with new values
bounds_new = bounds;
bounds_new.lon_min = single(round(lon_minx,n+1));
bounds_new.lon_max = single(round(lon_maxx,n+1));
bounds_new.lat_min = single(round(lat_minx,n+1));
bounds_new.lat_max = single(round(lat_maxx,n+1));

end