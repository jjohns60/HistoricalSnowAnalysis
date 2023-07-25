function NSIDC_HTTPS_ACCESS(HTTPS_PATH,username,password,wget_path,date,str_array,savepath)
%NSIDC_HTTPS_ACCESS Takes an input path to a particular dataset on the
%NSIDC server,a target datetime day, and a file search string array used to
%indicate the dataset type and the file type of interest. Files will be 
%stored at the directory indicated by 'savepath'
%
%
% INPUTS:
%   HTTPS_PATH - path to dataset (i.e., 'https://n5eil01u.ecs.nsidc.org/VIIRS/VNP10A1F.001/')
%   username - NSIDC login username as string, 'username'
%   password - NSIDC login password as string, 'password'
%   wget_path - path to wget executable on your system (code can be easily modified to use websave or urlwrite, but these are unreliable)
%   date - date of interest (to download)
%   str_array - cell array of strings to match to file names, allows download of only file types of interest (i.e., {'VNP10A1F','h5'})
%   savepath - location of which to save the downloaded files locally



%set permissions by updating options and headerfields (websave method only)
options = weboptions('HeaderFields',{'Authorization',...
    ['Basic ' matlab.net.base64encode([username ':' password])]});
%options = weboptions;
%options.Username = username;
%options.Password = password;
options.Timeout = 5;
options.UserAgent = 'Mozilla/5.0'; %simulates requests coming from web browser
options.ContentType = 'text';

%Get list of all filenames meeting conditions (1st two string arguments are
%required)
date_str = datestr(date,'yyyy.mm.dd');
%% THIS LINE HANGS IF SCRIPT MAKES TOO MANY REQUESTS QUICKLY
disp('1a')
html_text = webread([HTTPS_PATH date_str],options);
disp('1b')
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

%{
%initiate parallel computing pool (if not already)
poolobj = gcp('nocreate'); % If no pool, do not create new one.
if isempty(poolobj)
    parpool(4);
end
%}

%loop through with parallelization to speed up processing
%parfor i = 1:L
for i = 1:L

    %keep track of iteration and get individual file name
    file = links{i};
    %check if file already exists
    if sum(strcmp(file,existing_files)) > 0
        continue
        %if it does not, then download
    else
        
        %create wget command
        f_path = [HTTPS_PATH date_str '/' file];
        command = [wget_path ' --http-user=' username ' --http-password=' ...
            password [' -nH -nd -q --no-check-certificate --auth-no-challenge=on ' ...
            '-r --reject "index.html*" -np -e robots=off --directory-prefix=']...
            savepath ' ' f_path];

        try
            %MATLAB builtin methods websave and urlwrite are comparably slow
            %websave([savepath file],[HTTPS_PATH date_str '/' file],options);
            %urlwrite([HTTPS_PATH date_str '/' file],[savepath file],'Authentication','Basic','Username',username,'Password',password);
            
            %run the wget command
            system(command);

            %show progress of download
            disp(strcat("Downloading ",num2str(i)," of ",num2str(L)))
        catch
            pause(rand*1) %to prevent 503 errors
            
            %urlwrite([HTTPS_PATH date_str '/' file],[savepath file],'Authentication','Basic','Username',username,'Password',password);
            %websave([savepath file],[HTTPS_PATH date_str '/' file],options);
            
            %run the wget command
            system(command);

            %show progress of download
            disp(strcat("Downloading ",num2str(i)," of ",num2str(L)))
        end
    end

end

end