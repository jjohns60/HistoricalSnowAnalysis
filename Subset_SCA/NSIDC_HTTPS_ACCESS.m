function NSIDC_HTTPS_ACCESS(HTTPS_PATH,date,str_array,savepath)
%NSIDC_HTTPS_ACCESS Takes an input path to a particular dataset on the
%NSIDC server,a target datetime day, and a file search string array used to
%indicate the dataset type and the file type of interest. Files will be 
%stored at the directory indicated by 'savepath'

%{
%testing
HTTPS_PATH = 'https://n5eil01u.ecs.nsidc.org/VIIRS/VNP10A1F.001/';
date = datetime(2012,2,15);
str_array = {'VNP10A1F','h5'};
savepath = '/Volumes/GRA_Data_Backup/UNH Postdoc/VIIRS Cloud Filled/';
%}


%set permissions by updating options and headerfields
username = 'jjohns60';
password = '2080Bent!';
options = weboptions('HeaderFields',{'Authorization',...
    ['Basic ' matlab.net.base64encode([username ':' password])]});
%options = weboptions;
%options.Username = username;
%options.Password = password;
options.Timeout = 30;

%Get list of all filenames meeting conditions (1st two string arguments are
%required)
date_str = datestr(date,'yyyy.mm.dd');
html_text = webread([HTTPS_PATH date_str],options);
links = regexp(html_text,[str_array{1} '.\w*.\w*.\w*.\w*.' str_array{2}],'match');
links = unique(links);

%if there are more than 2 arguments, use to further trim the list of links
if length(str_array) > 2
    idx = zeros(size(links));
    for i = 3:length(str_array)
        match_str = str_array{i};
        idx = idx + double(contains(links,match_str));
    end
    %trim to files of interest
    links = links(idx > 0);
end

%loop through files and save to path
existing_files = dir(savepath);
existing_files = {existing_files.name};
L = length(links);

%loop through with parallelization to speed up processing
parfor i = 1:L

    %keep track of iteration and get individual file name
    file = links{i};
    %check if file already exists
    if sum(strcmp(file,existing_files)) > 0
        continue
        %if it does not, then download
    else
        try
            websave([savepath file],[HTTPS_PATH date_str '/' file],options);
            %show progress of download
            disp(strcat("Downloading ",num2str(i)," of ",num2str(L)))
        catch
            pause(rand*1) %to prevent 503 errors
            websave([savepath file],[HTTPS_PATH date_str '/' file],options);
            %show progress of download
            disp(strcat("Downloading ",num2str(i)," of ",num2str(L)))
        end
    end

end

end