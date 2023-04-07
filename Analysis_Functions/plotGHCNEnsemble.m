function plotGHCNEnsemble(GHCN_path,start_date,end_date)
%plotGHCNEnsemble Creates an ensemble visualization of GHCN observations
%   Input is a folder containing GHCN site data in .csv form. Currently,
%   only configured to plot min and max temperature and snow depth

files = dir([GHCN_path '*.csv']);
files = {files.name};
files = files(~startsWith(files,'_'));

%loop through each file and aggregate observations into large vectors
t = [];
Tmax = [];
Tmin = [];
SD = [];
I = [];
for i = 1:length(files)

    %read in data
    T = readtable([GHCN_path files{i}]);

    %append to list with all timesteps
    t = [t;datetime(T.DATE)];
    %append to list with all max temp data
    Tmax = [Tmax;T.TMAX];
    %append to list with all min temp data
    Tmin = [Tmin;T.TMIN];
    %append to list with all snow depth data
    SD = [SD;T.SNWD];
    %create index for each unique site
    I = [I;repmat(i,size(T.DATE))];

    %display progress
    disp(strcat("Completed ",num2str(i),"/",num2str(length(files))))
end

%convert to proper units (celsius and cm)
Tmax = Tmax/10; %C
Tmin = Tmin/10; %C
SD = SD/10; %cm

%quality control of data
Tmax(Tmax > 100) = NaN;
Tmax(Tmax < -100) = NaN;
Tmin(Tmin > 100) = NaN;
Tmin(Tmin < -100) = NaN;
SD(SD > 1000) = NaN;
SD(SD < 0) = NaN;

%loop through data for each unique day in the record
U = unique(t);
%apply limitation on dates
U = U(U >= start_date & U <= end_date);
Tmax_avg = NaN(size(U));
Tmin_avg = NaN(size(U));
SD_avg = NaN(size(U));
Tmax_Std = NaN(size(U));
Tmin_Std = NaN(size(U));
SD_Std = NaN(size(U));
for i = 1:length(U)

    %specify day of interest and index
    d = U(i);
    idx = (t == d);

    %compute average daily values
    Tmax_avg(i) = nanmean(Tmax(idx));
    Tmin_avg(i) = nanmean(Tmin(idx));
    SD_avg(i) = nanmean(SD(idx));
    
    %if at least 5 values, compute statistics
    if sum(idx) >= 5
        %compute standard deviation, store 2 above and below the mean
        Tmax_Std(i) = nanstd(Tmax(idx),[],'all');
        Tmin_Std(i) = nanstd(Tmin(idx),[],'all');
        SD_Std(i) = nanstd(SD(idx),[],'all');
    end

    %display progress
    disp(strcat("Processed ",datestr(d)))
end

%create plot
figure;
subplot(2,1,1);
%min temps
plot(U,Tmin_avg,'b','LineWidth',1.5); hold on;
X = [U;flip(U)];
Y = [Tmin_avg + (2*Tmin_Std);flip(Tmin_avg - (2*Tmin_Std))];
idx = isnan(Y);
fill(X(~idx),Y(~idx),'b','FaceAlpha',0.5,'EdgeAlpha',0);
%max temps
plot(U,Tmax_avg,'r','LineWidth',1.5); hold on;
X = [U;flip(U)];
Y = [Tmax_avg + (2*Tmax_Std);flip(Tmax_avg - (2*Tmax_Std))];
idx = isnan(Y);
fill(X(~idx),Y(~idx),'r','FaceAlpha',0.5,'EdgeAlpha',0);

%formatting
xlim([start_date end_date])
grid on;
ylabel(['Temperature (' char(176) 'C)'])
set(gca,'Fontsize',22)

%add snow depth
subplot(2,1,2);
plot(U,SD_avg,'Color',[0.8500 0.3250 0.0980],'LineWidth',1.5); hold on;
X = [U;flip(U)];
Y = [SD_avg + (2*SD_Std);flip(SD_avg - (2*SD_Std))];
idx = isnan(Y);
fill(X(~idx),Y(~idx),[0.8500 0.3250 0.0980],'FaceAlpha',0.5,'EdgeAlpha',0);

%formatting
xlim([start_date end_date])
grid on;
ylabel(['Temperature (' char(176) 'C)'])
set(gca,'Fontsize',22)
grid on;
ylabel('Snow Depth (cm)')
ylim([0 300])
set(gca,'Fontsize',22)





