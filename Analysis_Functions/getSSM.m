function D = getSSM(snow_durations,d)
%getSSM Computes the snow seasonality metric depending on the duration
%threshold (d, days) for a list of snow cover durations

days_seasonal = sum(snow_durations(snow_durations >= d));

days_ephemeral = sum(snow_durations(snow_durations < d));

days_snow = sum(snow_durations(~isnan(snow_durations)));

D = (days_seasonal - days_ephemeral)./(days_snow);

end