function Motive_Spine_Marker_Sorting (filePath, varargin)

% This function takes either a file or a folder of Motive exported .csv
% files and sorts the markers first along the Superior Inferior axis and
% then along the Mediolateral axis. The Markers are then re-exported with
% the same name to a new file called sorted with appropriate headers. This
% function can be used if a full grid of markers isn't used by specifiying
% gridSize and missingMrk parameters. 

% filePath is a string specifiying the path to either a directory or file.
%          If filePath is a directory it will perform Marker_Sort on every
%          file containting the string '.csv'
% nameString (optional) will only read in files containing the nameString 
%          (default: '*.csv')
% sortOrder (optional) specifies the order of the rows that the program
%          should sort along. Default is [2,1] because most likely the
%          superior/inferior axis is along 'y' and the mediolateral axis is
%          along 'x'. The program does not sort along the third axis
%          (AP/z). For headers to work +x should be Left, +y = Up, +z = Ant
% nColumns (optional) number of columns in the grid (default = 3)
% startLvl (optional) type cell: {'S1'} default. First spine level of
%          markers
% missingMrk (optional) is an array of markers that are missing in the
%          grid. Each row is a missing value. Column1 is the Column# (Left,
%          right, center) and column2 is the spine level. default is [];
% gridSize (optional) gives the number of columns and spine levels. If
%          default [], grid size will solve based on number of markers.
% visualize (optional) should be set to 'on' or 'off' (default) in order to
%          see a labelled figure of the back. If visualized, data will not 
%          write to file unless the user accepts the sorted data.
%
% Written by Derek Zwambag April 12, 2018
% Version 1.0

%% Set up Input Parser
p = inputParser;

% filePath must be a file or a directory
validationFcn = @(s) assert( exist(s,'dir') | exist(s,'file') ,'filePath is not valid');
addRequired(p,'filePath',validationFcn); 
% nameString 
default = '*.csv';
addParameter(p,'nameString',default);
% sortOrder must be an array of two values 
validationFcn = @(d) assert(numel(d)==2 & isnumeric(d) & any(d)<3 ,'sortOrder should be two integers');    
defaultVal = [2,1];
addOptional(p,'sortOrder',defaultVal,validationFcn); 
% nColumns must be an integer
validationFcn = @(i) assert(isnumeric(i) && isscalar(i), 'nColumns should be a single integer');
defaultVal = 3;
addParameter(p,'nColumns',defaultVal,validationFcn);
% startLevel must match a cell in spineLvl
spineLvl = {'S2','S1','L5','L4','L3','L2','L1','T12','T11','T10','T9','T8','T7','T6','T5','T4','T3','T2','T1','C7','C6'}; 
validationFcn = @(s) assert(any(strcmp(s,spineLvl)),'startLvl is not valid');
default = {'S1'};
addParameter(p,'startLvl',default,validationFcn);
% missingMrk must be a [nx2] matrix of integers
validationFcn = @(d) assert(isnumeric(d) && size(d,2) == 2,'missingMrk must be [nx2] matrix');
default = [];
addParameter(p,'missingMrk',default,validationFcn);
% grid size must be two integers
validationFcn = @(d) assert(numel(d)==2 & isnumeric(d),'sortOrder should be two integers');    
defaultVal = [];
addOptional(p,'gridSize',defaultVal,validationFcn); 
% visualize must be 'on' or 'off'
validationFcn = @(s) assert(any(strcmp(s,{'off','on'})),'visualize should be ''on'' or ''off''');
default = 'off';
addParameter(p,'visualize',default,validationFcn);

% Parse Inputs
parse(p,filePath,varargin{:});
filePath  = p.Results.filePath;
nameString = p.Results.nameString;
sortOrder = p.Results.sortOrder;
nColumns  = p.Results.nColumns;
startLvl  = p.Results.startLvl;
missingMrk = p.Results.missingMrk;
gridSize = p.Results.gridSize;
visualize = p.Results.visualize;

%% Read in all files contained in Folder that include nameString
if exist(filePath,'dir')
    files = dir( fullfile(filePath,nameString) );
    for k = numel(files):-1:1
        file{k} = fullfile(filePath,files(k).name);
    end
else
    file = {filePath};
end

%% Create Sorted Folder if it doesn't exist
[folder,~,ext] = fileparts( file{1} ); 
if exist([folder,'/Sorted'],'dir')==0
    status = mkdir([folder,'/Sorted']);
    if status == 0
        disp('Problem creating folder')
    end
end

%% Perform marker sort by looping through each file
for kk = 1:numel(file)
    
    % Get file header lines
    fid  = fopen( file{kk} );
    header1 = fgetl( fid ); fgetl( fid ); fgetl( fid ); fgetl( fid ); fgetl( fid ); fgetl( fid );
    header6 = fgetl( fid );
    fclose( fid );
    
    % Read in data (Motive has 7 lines of headers)
    data = csvread( file{kk}, 7, 0);
    n = (size(data,2) - 2) / 3; % Number of markers (+Frame and Time)
    f = size(data,1);           % Number of frames
    

   %% NEW FOR DEVON TO USE WITH VICON
    
    % Vicon uses a +X = Left; +Y = Post; +Z = Up;
    data_new = data;
    data_new(:, 4:3:end) =  data(:, 5:3:end);
    data_new(:, 5:3:end) = -data(:, 4:3:end); 
    %OR
    % Vicon uses a +X = Right; +Y = Ant; +Z = Up;
    % Change to Motive coordinate system where
    % +X = Left; +Y = Up; +Z = Ant
    
%     data_new = data;
%     data_new(:, 3:3:end) =  -data(:, 3:3:end);
%     data_new(:, 4:3:end) =  data(:, 5:3:end);
%     data_new(:, 5:3:end) = data(:, 4:3:end); 
    
    data = data_new; clear data_new

    %% Create a grid of markers similar to the collection
    if ~isempty(gridSize)
        m = ones(gridSize);
    else
        m = ones(nColumns,ceil(n/nColumns));
    end
    
    % Delete missing markers
    if ~isempty(missingMrk)
        for k = 1:size(missingMrk,1);
            m(missingMrk(k,1),missingMrk(k,2)) = 0;
        end
    end
   
    [i,j] = find(m);
    % i represents the Left, Middle, and Right Columns
    % j represents the Spine Level
    
    if numel(i) ~= n && exist(filePath,'dir')
        warning([file{kk},': Skipped']) % Display Status to User
        continue
    elseif numel(i) ~= n
        dataFirst = data(1,3:end);
        dataFirst = reshape(dataFirst, 3, n );
        scatter( dataFirst(1,:), dataFirst(2,:));
        text(dataFirst(1,:) + 0.005,dataFirst(2,:),num2str((1:n)'))
        error('Number of markers doesn''t match grid in %s',file{kk})
    else
  
    % Rearrange first line of data into a [3xn] column of [x;y;z] values
    dataFirst = data(1,3:end);
    dataFirst = reshape(dataFirst, 3, n );
       
    % Sort along 1st value in sortOrder
    [~,mOrder] = sort( dataFirst(sortOrder(1),:) );
    
    % For each spine level
    for k = 1:max(j);
        mLvl = mOrder(j==k); % Identify Marker numbers for this spine level
        [~,sortIdx] = sort( dataFirst( sortOrder(2),mLvl) ); % Sort these markers along axis identified by sortOrder(2)
        mOrder(j==k) = mLvl(sortIdx); % Correct marker numbers along ML axis
    end 
    
    % Old method
    % Each column of [i] represents the nColumns of markers at each spine level
    %i = reshape(i,nColumns,n/nColumns);
    
    % Resort again along 2nd value in sortOrder. 
    %for j = 1:size(i,2)
    %   [~,ii] = sort( dataFirst(sortOrder(2),i(:,j)) );
    %   i(:,j) = i(ii,j);
    %end
    
    dataSorted = reshape(data(:,3:end),f,3,n); % Reshape to make it easy to sort using mOrder
    dataSorted = dataSorted(:,:,mOrder);
    dataSorted = reshape(dataSorted,f,3*n);    % Reshape back to [f,3*n] matrix
    
    % Figure if 'visualize' is 'on'
    if strcmp(visualize,'on')
        figure('Position',[450,200,800,500]);
        %scatter( dataSorted(1,sortOrder(2):3:end), dataSorted(1,sortOrder(1):3:end) )
        scatter3( dataSorted(1,1:3:end), dataSorted(1,3:3:end), dataSorted(1,2:3:end) )
        hold on; axis equal
        %text(dataSorted(1,sortOrder(2):3:end) + 0.005,dataSorted(1,sortOrder(1):3:end),num2str((1:n)'))
        text( dataSorted(1,1:3:end)+ 0.005, dataSorted(1,3:3:end)+ 0.005, dataSorted(1,2:3:end)+0.005,num2str((1:n)'))
        
        % Users must confirm everything is good if 'visualize' is 'on'
        button = menu('Good to write to file and continue?','Yes','Quit');
        
        if button ~= 1
            return
        end
        close all
    end
    
    %% Write to File.
    dataOut = [data(:,1:2),dataSorted]; % Output data
    
    [~,name,~] = fileparts(file{kk});    % Output file
    fileOut = [folder,'/Sorted/',name,ext];
    
    % Identify starting Spine Lvl
    [~,idx] = find(strcmp(spineLvl,startLvl));
    spineLvl = spineLvl(idx:end);
    
    % set up headers
    header2 = ['Data sorted with Motive_Spine_Marker_Sorting version 1.0 written by Derek Zwambag on ',date];
    header3 = repmat({','},1,n*3+2); header3{2} = 'Spine Level,';
    header4 = repmat({','},1,n*3+2); header4{2} = 'Row,';
    header5 = repmat({','},1,n*3+2); header5{2} = 'Column,';
    for k = 1:n
        header3{k*3} = [spineLvl{j(k)},','];
        header4{k*3} = [num2str(i(k)),','];
        header5{k*3} = [num2str(j(k)),','];
    end
     
    fid = fopen(fileOut,'w');
    fprintf(fid,'\n');
    fprintf(fid,[header2,'\n']);
    fprintf(fid,[cell2mat(header3),'\n']);
    fprintf(fid,[cell2mat(header4),'\n']);
    fprintf(fid,[cell2mat(header5),'\n']);
    fprintf(fid,[header6,'\n']);
    fclose(fid);
    dlmwrite(fileOut,dataOut,'-append');
    
    disp([file{kk},': Complete']) % Display Status to User
    end % End if number of markers matches expected
end