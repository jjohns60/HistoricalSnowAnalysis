function [SCI,CCI] = getSCI(SCE,t,valid_snow,valid_land,valid_cloud,type)
%getSCI computes the Snow Condition Index (SCI) which was presented in
%Richer et al., 2013, Moore et al., 2015, and Saavedra et al., 2017. This
%metric computes the proportion of which a given day had snow cover through
%a multi-year time series. Also reports the Cloud Cover Index (CCI) which
%keeps a count of the proportion of days that are cloud covered through the
%years, or had missing data. Output is a sorted 3D matrix in which the 3rd
%dimension represents DOY (1 = Jan 1, 365/366 = Dec 31).

%% NOTE: type 1 will not work for very large data arrays (>> 1000 x 1000)

%get each unique day included in record
[~,ind] = unique(datestr(t,'mmdd'),'rows','stable');
t_days = t(ind);
[~,ind] = sort(day(t_days,'dayofyear'));
t_days = t_days(ind);
[~,ind] = sort(month(t_days));
t_days = t_days(ind);

%for full data array
if type == 1
  
    %prepare arrays to store data
    SCI = NaN([size(SCE(:,:,1),1) size(SCE(:,:,1),2) length(t_days)]);
    CCI = NaN([size(SCE(:,:,1),1) size(SCE(:,:,1),2) length(t_days)]);


    %loop through each unique date in the time list t
    L = length(t_days);
    disp(['Processing ' num2str(L) ' timesteps....'])
    for i = 1:L
        disp([num2str(i) '/' num2str(L)])

        %get all timesteps matching the specific date (for all years)
        idx = day(t) == day(t_days(i)) & month(t) == month(t_days(i));

        %extract all relevant data grids
        SCE_i = double(SCE(:,:,idx));

        %return SCI for the given data and append to array
        SC_total = sum(SCE_i == valid_snow,3);
        SC_valid_total = sum(SCE_i == valid_snow,3) + sum(SCE_i == valid_land,3);
        SCI(:,:,i) = SC_total./SC_valid_total;


        %return CCI or missing data proportion
        CC_total = sum(SCE_i ~= valid_snow & SCE_i ~=  valid_land | SCE_i == valid_cloud,3);
        CCI(:,:,i) = CC_total./size(SCE_i,3);

    end

%to handle larger arrays
elseif type == 2
    
    %loop through all unique days
    L = length(t_days);
    disp(['Processing ' num2str(L) ' timesteps....'])
    for i = 1:L
        tic
        disp([num2str(i) '/' num2str(L)])

        %identify files from the given path that match the date
        idx = month(t) == month(t_days(i)) & day(t) == day(t_days(i));
        SCE_files_i = SCE(idx);

        %loop through all years
        for ii = 1:length(SCE_files_i)
            d = ncread(SCE_files_i{ii},'IMS_Surface_Values');
            
            %get dimensions once
            if i == 1 && ii == 1
                SCI_sum = zeros(size(d));
                CCI_sum = zeros(size(d));
            end

            if ii == 1

                SCI_i = zeros(size(d)) + double(d == valid_snow);
                SCI_sum_valid = zeros(size(d)) + double(d == valid_snow) + double(d == valid_land); 
                CCI_i = zeros(size(d)) + double(d == valid_cloud);

            elseif ii == length(SCE_files_i)

                SCI_i = SCI_i + double(d == valid_snow);
                SCI_sum_valid = SCI_sum_valid + double(d == valid_snow) + double(d == valid_land); 
                CCI_i = CCI_i + double(d == valid_cloud);
                
                %get final value for the given date
                SCI_i = SCI_i./SCI_sum_valid;
                CCI_i = CCI_i./length(SCE_files_i);

            else
                SCI_i = SCI_i + double(d == valid_snow);
                SCI_sum_valid = SCI_sum_valid + double(d == valid_snow) + double(d == valid_land); 
                CCI_i = CCI_i + double(d == valid_cloud);
            end
            toc
        end
        
        %update arrays for computing SP (average SCI across all days)
        SCI_sum = SCI_sum + SCI_i;
        CCI_sum = CCI_sum + CCI_i;
    end
    %compute SP (but assign to SCI for this type)
    SCI = SCI_sum ./ length(t_days);
    CCI = CCI_sum ./ length(t_days);
end
end