function Intersegmental_Spine_Model(filePath,varargin)
% This function performs the Intersegmental_Spine_Model on data that have
% been collectd in Motive and sorted using Motive_Spine_Marker_Sorting.m
% function. 

% Input: 
% filePath (required): can be a file or folder
% nameString (optional - default: '*csv'): Only performs the function on 
%      files containing the 'nameString'
% outFreq (optional - default: 30): Target frequency of the data written to
%      file
% lcsLevels (optional - default: []): LCS that are calculated (S1 == 1
%      C7 == 19). If default is chosen then a LCS will be calcuated at each
%      level that a marker is provided.
% includeLevels (optional - default: [1,19]: This option will drop marker
%      data from spine levels outside of 'includeLevels'. Can be used to
%      calculate lower T9-S1 model for males. This is done before padding
% nPadding (optional - default: 3): Number of padded rows of markers
%      included at the top and bottom of each column.
% visualize (optional - default: 'off'): Can be used to look at the
%      scatterplot of markers from the first frame. LCS are not calculated 
%      if this is 'on'.
% KnotsOption (optional - default: []): specifies the number of knots used 
%      in the piecewise polynomial where pieces = nKnots - 1;
%      The original ABME paper used 6 knots across 21 levels yeilding 4
%      spine levels per piece (nSpineLvls - 1) / (k-1) ).
%      When padding, there were 4.8 spine levels / piece
%      I don't think it's a good idea to go below 4 spine levels/peice 
%      because then the spline is too flexible (can go through every point.
%      If default is selected it will choose the number of knots to ensure 
%      at least 4 markers/segment. Alternatively specify a scalar (the
%      number of knots to be used) or a vector of desired knot placements.
% writeToFolder (optional - default '/IntersegmentalLCS/'): name of the
%      folder created to output data. This is always a child folder of the
%      input folder. Don't forget to have starting and ending '/'!
% notify (optional - default: 'off'): If 'on' plays a train whistle at the
%      end to notify the program is done. 

% Written by Derek Zwambag April 15, 2018

%% Set up Input Parser
p = inputParser;

% filePath must be a file or a directory
validationFcn = @(s) assert( exist(s,'dir') | exist(s,'file') ,'filePath is not valid');
addRequired(p,'filePath',validationFcn); 

% nameString 
default = '*.csv';
addParameter(p,'nameString',default);

% outFreq should be the output frequency
validationFcn = @(d) assert( isnumeric(d) & isscalar(d),'outFreq must be a scalar number');
default = 30;
addParameter(p,'outFreq',default,validationFcn);

% lcsLevels
validationFcn = @(d) assert( isnumeric(d) & all(size(d(:)) == [2,1]), 'lcsLevelsOption should be two numbers');
default = [];
addParameter(p,'lcsLevelsOption',default,validationFcn);

% includeLevels
validationFcn = @(d) assert( isnumeric(d) & all(size(d(:)) == [2,1]), 'includeLevels should be two numbers');
default = [1,19];
addParameter(p,'includeLevels',default,validationFcn);

% nPadding
validationFcn = @(d) assert( isnumeric(d) & isscalar(d),'nPadding must be a scalar number');
default = 3;
addParameter(p,'nPadding',default,validationFcn);

% nKnots
validationFcn = @(d) assert( isnumeric(d),'KnotsOption must be a numeric number');
default = [];
addParameter(p,'KnotsOption',default,validationFcn);

% visualize must be 'on' or 'off'
validationFcn = @(s) assert(any(strcmp(s,{'off','on'})),'visualize should be ''on'' or ''off''');
default = 'off';
addParameter(p,'visualize',default,validationFcn);

% writeToFolder 
default = '/IntersegmentalLCS/';
addParameter(p,'writeToFolder',default);

% sound 
default = 'off';
addParameter(p,'notify',default);

% Parse Inputs
parse(p,filePath,varargin{:});
filePath = p.Results.filePath;
nameString = p.Results.nameString;
outFreq = p.Results.outFreq;
lcsLevelsOption = p.Results.lcsLevelsOption;
nPadding = p.Results.nPadding;
KnotsOption = p.Results.KnotsOption;
includeLevels = p.Results.includeLevels;
visualize = p.Results.visualize;
writeToFolder = p.Results.writeToFolder;
notify = p.Results.notify;

%% Read in all files contained in Folder that include nameString

if exist(filePath,'dir')
    files = dir( fullfile(filePath,nameString) );
    for k = numel(files):-1:1
        file{k} = fullfile(filePath,files(k).name);
    end
else
    file = {filePath};
end

%% Perform Intersegmental Model on each file
for kk = numel(file):-1:1
    
    % Read in data (Motive has 7 lines of headers)
    data = csvread( file{kk}, 7, 0);
    t = data(:,1:2);   % Frames and Time
    m = data(:,3:end); % Marker data
    [f,n] = size(m);   % Number of frames and columns
    global Outlength %Allow for variable to be called in other functions
    Outlength = f %Don't want f as a constant
    clear data
    
    %% Downsample if outFreq is less than input frequency
    
    dt = t(2,2) - t(1,2); % sample period
    if round(1/dt) > outFreq 
        dsFactor = round( round(1/dt) / outFreq );
        t = downsample(t,dsFactor);
        m = downsample(m,dsFactor);
        [f,~] = size(m);   % Update number of frames
        
        dt = t(2,2) - t(1,2); % Update sample period
        
        if round(1/dt) ~= outFreq
            warning('Output Frequency is %d Hz',round(1/dt))
        end
    elseif round(1/dt) < outFreq
        warning('Input frequency of %d Hz is less than target Output Frequency',round(1/dt))
    end
    clear dt dsFactor
    
    %% Convert data to preferred GCS
    
    % Input GCS:      +x Left     ; +y Up; +z Anterior
    % Preferred GCS:  +x Posterior; +y Up; +z Left
    
    X = -m(:,3:3:end);
    Y =  m(:,2:3:end);
    Z =  m(:,1:3:end);
    clear m
    %% Obtain Marker Columns and Spine Levels
    
    fid = fopen(file{kk});
    fgetl(fid); fgetl(fid);
    SPINE = fgetl(fid);
    fclose(fid);
    clear fid
    
    % SPINE is a cell of vertebral body names
    SPINE = strsplit(SPINE,',');
    SPINEstr = SPINE{3};
    
    SPINE = {'S2','S1','L5','L4','L3','L2','L1','T12','T11','T10','T9','T8','T7','T6','T5','T4','T3','T2','T1','C7'};
    SPINEstr = find(cellfun(@any,regexp(SPINE,SPINEstr)));
    
    SPINE = SPINE(SPINEstr:end);
    
    % Obtain indicies for the column number and the spineLvl
    mrkGrid = csvread(file{kk},3,2,[3,2,4,2+n-1]);
    column = mrkGrid(1,1:3:end);
    spineLvl = mrkGrid(2,1:3:end);
    clear mrkGrid n
    
    % number of columns of markers
    nColumn = max(column);
    
    %% Option to only use specific spine levels as input
    k = (spineLvl < min(includeLevels) | spineLvl >  max(includeLevels)) ;
    column(k) = [];
    spineLvl(k) = [];
    X(:,k) = [];
    Y(:,k) = [];
    Z(:,k) = [];
    clear k
    
    %% Pad markers here
    lmin = min(spineLvl);
    lmax = max(spineLvl);
    if nPadding > 0
        
        % Pad markers at top and bottom
        BotPadX = mean(X(:,spineLvl==lmin),2) - mean(X(:,spineLvl==lmin+1),2); % Row 1 (Mrks 1-3) - Row 2 (Mrks 4-6)
        BotPadY = mean(Y(:,spineLvl==lmin),2) - mean(Y(:,spineLvl==lmin+1),2);
        BotPadZ = mean(Z(:,spineLvl==lmin),2) - mean(Z(:,spineLvl==lmin+1),2);
        TopPadX = mean(X(:,spineLvl==lmax),2) - mean(X(:,spineLvl==lmax-1),2); % Row 19 (Mrks 55-57) - Row 18 (Mrks 52-44)
        TopPadY = mean(Y(:,spineLvl==lmax),2) - mean(Y(:,spineLvl==lmax-1),2);
        TopPadZ = mean(Z(:,spineLvl==lmax),2) - mean(Z(:,spineLvl==lmax-1),2);
        
        ii = nColumn*2;
        Xpad = zeros(f,nColumn*2*nPadding); Ypad = Xpad; Zpad = Xpad;
        for k = 1:nPadding
             Xpad(:,(k-1)*ii+1:k*ii) = [ X(:,spineLvl==lmin) + k*repmat(BotPadX,1,nColumn) , X(:,spineLvl==lmax) + k*repmat(TopPadX,1,nColumn) ];
             Ypad(:,(k-1)*ii+1:k*ii) = [ Y(:,spineLvl==lmin) + k*repmat(BotPadY,1,nColumn) , Y(:,spineLvl==lmax) + k*repmat(TopPadY,1,nColumn) ];
             Zpad(:,(k-1)*ii+1:k*ii) = [ Z(:,spineLvl==lmin) + k*repmat(BotPadZ,1,nColumn) , Z(:,spineLvl==lmax) + k*repmat(TopPadZ,1,nColumn) ];
        end
        
        X = [X,Xpad];
        Y = [Y,Ypad];
        Z = [Z,Zpad];
        
        columnPad = repmat(1:nColumn,1,nPadding*2);
        spineLvlPad1 = repmat([lmin lmax],nColumn,nPadding);  
        spineLvlPad2 = repmat(1:nPadding,nColumn*2,1).*repmat([-ones(nColumn,1);ones(nColumn,1)],1,nPadding);
        spineLvlPad = (spineLvlPad1(:) + spineLvlPad2(:))';
        
        column = [column,columnPad];
        spineLvl = [spineLvl,spineLvlPad];
    end % end of padding
    clear BotPadX BotPadY BotPadZ TopPadZ TopPadY TopPadX Xpad Ypad Zpad
    clear columnPad spineLvlPad spineLvlPad1 spineLvlPad2
    
    %% Determine number of Knots to use
    % If lcsLevels is empty (default) solve Local Coordinate Systems for
    % all spine levels between the top and bottom marker.
    if isempty(lcsLevelsOption)
        lcsLevels = lmin:lmax;
    else
        lcsLevels = lcsLevelsOption;
    end
    
    % Determine the number of knots in the piecewise spline
    % number of Spine Levels is the max - min + 1
    % #Knots is (nSpineLvls - 1)/ levelsPerSegment + 1; 
    % Need to ciel to ensure it is an integer
    nSpineLvls = max(spineLvl) - min(spineLvl) + 1;
    if isempty(KnotsOption)
        levelsPerSegment = 5;
        Knots_input = 1 + ceil( (nSpineLvls -1) / levelsPerSegment );
        nKnots = Knots_input;
    else
        Knots_input = KnotsOption;
        if length(Knots_input) > 1
            nKnots = length(Knots_input);
        else
            nKnots = Knots_input;
        end
    end
    trueMarkerPerLevel =  (nSpineLvls -1) / (nKnots - 1);
    clear nKnots nSpineLvls levelsPerSegment 
    
    %% Fit Piecewise Linear Polynomial to each frame
    out = zeros(12,length(lcsLevels),f);
    for k = f:-1:1;  
        % Right column
        RslmX = slmengine_DZ(spineLvl(column==1),X(k,column==1),'knots',Knots_input,'plot',visualize);
        RslmY = slmengine_DZ(spineLvl(column==1),Y(k,column==1),'knots',Knots_input,'plot',visualize);
        RslmZ = slmengine_DZ(spineLvl(column==1),Z(k,column==1),'knots',Knots_input,'plot',visualize);
        
        % Middle column
        MslmX = slmengine_DZ(spineLvl(column==2),X(k,column==2),'knots',Knots_input,'plot',visualize);
        MslmY = slmengine_DZ(spineLvl(column==2),Y(k,column==2),'knots',Knots_input,'plot',visualize);
        MslmZ = slmengine_DZ(spineLvl(column==2),Z(k,column==2),'knots',Knots_input,'plot',visualize);
        
        % Left column
        LslmX = slmengine_DZ(spineLvl(column==3),X(k,column==3),'knots',Knots_input,'plot',visualize);
        LslmY = slmengine_DZ(spineLvl(column==3),Y(k,column==3),'knots',Knots_input,'plot',visualize);
        LslmZ = slmengine_DZ(spineLvl(column==3),Z(k,column==3),'knots',Knots_input,'plot',visualize);
        
        if strcmp(visualize,'on')
            
                %% plot if visualize is 'on'
   
            figure; hold on
            scatter3(X(k,column==1),Y(k,column==1),Z(k,column==1),'r','filled')
            scatter3(X(k,column==2),Y(k,column==2),Z(k,column==2),'g','filled')
            scatter3(X(k,column==3),Y(k,column==3),Z(k,column==3),'b','filled')
            scatter3(X(k,spineLvl<lmin),Y(k,spineLvl<lmin),Z(k,spineLvl<lmin),'markeredgecolor','k')
            scatter3(X(k,spineLvl>lmax),Y(k,spineLvl>lmax),Z(k,spineLvl>lmax),'markeredgecolor','k') 
            plot3(slmeval(lmin:0.1:lmax,RslmX,0),slmeval(lmin:0.1:lmax,RslmY,0),slmeval(lmin:0.1:lmax,RslmZ,0),'r')
            plot3(slmeval(lmin:0.1:lmax,MslmX,0),slmeval(lmin:0.1:lmax,MslmY,0),slmeval(lmin:0.1:lmax,MslmZ,0),'g')
            plot3(slmeval(lmin:0.1:lmax,LslmX,0),slmeval(lmin:0.1:lmax,LslmY,0),slmeval(lmin:0.1:lmax,LslmZ,0),'b')
            axis equal
            
            disp('visualize stops progress before LCS are calculated')
            return
        end
        
        %% Evaluate Splines to get LCS
        % ML = Left - Right splines
        MLvecs = [slmeval(lcsLevels,LslmX,0);slmeval(lcsLevels,LslmY,0);slmeval(lcsLevels,LslmZ,0)]-...
                 [slmeval(lcsLevels,RslmX,0);slmeval(lcsLevels,RslmY,0);slmeval(lcsLevels,RslmZ,0)];
        MLdir = MLvecs./repmat(sqrt(sum(MLvecs.^2)),3,1); % Unit vector

        SIvecs = [slmeval(lcsLevels,MslmX,1);slmeval(lcsLevels,MslmY,1);slmeval(lcsLevels,MslmZ,1)];
        SIdir = SIvecs./repmat(sqrt(sum(SIvecs.^2)),3,1); % Unit vector

        APvecs = cross(SIdir,MLdir);
        APdir = APvecs./repmat(sqrt(sum(APvecs.^2)),3,1); % Unit vector

        % Recross ML vector 
        MLdir = cross(APdir,SIdir);
        
        %% Ouput structure
        out(:,:,k) = [slmeval(lcsLevels,MslmX,0);slmeval(lcsLevels,MslmY,0);slmeval(lcsLevels,MslmZ,0);APdir;SIdir;MLdir];
        
    end % end loop through each frame
       
    %% Setup output for writing to file.
    
    out = [t,reshape(permute(out,[3,1,2]),f,12*length(lcsLevels))];
    
    [folder,fileName,ext] = fileparts( file{kk} ); 
    if exist([folder,writeToFolder],'dir')==0
        status = mkdir([folder,writeToFolder]);
        if status == 0
            disp('Problem creating folder')
        end
    end

    %% Write to File
    
    fileOut = fullfile([folder,writeToFolder,fileName,ext]);
    fid = fopen(fileOut,'w');
    fprintf(fid,'File created with: Intersegmental_Spine_Model function by Derek Zwambag\n');
    fprintf(fid,['File: ',fileName,'\n']);
    fprintf(fid,['Date:,',num2str(date),',#Knots:,',num2str(Knots_input),',#markers/level:,',num2str(trueMarkerPerLevel),'\n,,']);    
    for i = 1:length(lcsLevels)
        fprintf(fid,[SPINE{i},',,,,,,,,,,,,']);
    end
    fprintf(fid,'\nFrame,Time(s),');
    Head = repmat({'X,Y,Z,APi,APj,APk,SIi,SIj,SIk,MLi,MLj,MLk,'},1,length(lcsLevels));
    for i = 1:numel(Head)
        fprintf(fid,Head{i});
    end
    fprintf(fid,'\n');
    fclose(fid);
    dlmwrite(fileOut,out,'-append');
    
    disp([file{kk},': Complete']) % Display Status to User

end % end of processing loop through files.
disp(['Finished at:',datestr(now)]);

if strcmp(notify,'on')
    load train
    sound(y,Fs)
end

end % end of function