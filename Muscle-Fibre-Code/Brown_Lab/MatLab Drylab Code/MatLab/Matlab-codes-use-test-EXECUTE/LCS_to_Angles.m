function LCS_to_Angles(filePath,varargin)

% This function reads in a file or folder of data output by
% Intersegmental_Spine_Model.m function and performs the following steps.

% 1) The Intersegmental, Total (C7/S1), Thoracic (C7/T12), and Lumbar 
%    (T12/S1) Cardan Angles are calculated.
% 2) Positions and Angles are low pass filtered
% 3) Whole time series are time normalized to include the same number of
%    frames. 

% Input: 
% filePath (required): can be a file or folder
% nameString (optional - default: '*csv'): Only performs the function on 
%      files containing the 'nameString'
% outFrames (optional - default: 101): Number of frames written in output
% writeToFolder (optional - default '\Processed\'): name of the
%      folder created to output data. This is always a child folder of the
%      input folder. 
% regionalAngle (optional - default: {'S1',T12','C7'}; This can be used to
%      specify the LCS which define the lumbar, thoracic, and total angle.

% Written by Derek Zwambag April 16, 2018

%% Set up Input Parser
p = inputParser;

% filePath must be a file or a directory
validationFcn = @(s) assert( exist(s,'dir') | exist(s,'file') ,'filePath is not valid');
addRequired(p,'filePath',validationFcn); 

% nameString 
default = '*.csv';
addParameter(p,'nameString',default);

% outFrames is the number of frames in output
validationFcn = @(d) assert( isnumeric(d) & isscalar(d),'outFrames must be a scalar number');
default = 101;
addParameter(p,'outFrames',default,validationFcn);

% regionalAngle must be a cell of strings defining the borders of the
% regions
validationFcn = @(d) assert( numel(d) == 3 & iscell(d), 'regionalAngle should be a 1x3 cell');
default = {'S1','T12','C7'};
addParameter(p,'regionalAngle',default, validationFcn);

% writeToFolder 
default = '/Processed/';
addParameter(p,'writeToFolder',default);

% Parse Inputs
parse(p,filePath,varargin{:});
filePath   = p.Results.filePath;
nameString = p.Results.nameString;
outFrames  = p.Results.outFrames;
regionalAngle = p.Results.regionalAngle;
writeToFolder = p.Results.writeToFolder;

%% Read in all files contained in Folder that include nameString

if exist(filePath,'dir')
    files = dir( fullfile(filePath,nameString) );
    for k = numel(files):-1:1
        file{k} = fullfile(filePath,files(k).name);
    end
else
    file = {filePath};
end

%% Perform function on each file

for kk = numel(file):-1:1
    %% Parse into postion and angular data
    
    % Read in data (Intersegmental_Spine_Model files have 5 lines of headers)
    data = csvread( file{kk}, 5, 0);
    t = data(:,1:2);     % Frame #s and Time
    [f,nc] = size(data);   % Number of frames and columns
    
    lcs = reshape(data(:,3:end),f,12,(nc-2)/12); % lcs
    n = size(lcs,3);
    
    position = lcs(:,1:3,:); % position of each LCS
    position = reshape(position,f,3*n);
    dcm = lcs(:,4:12,:);     % direction cosine matrix

    % Read in Spine Levels
    fid = fopen(file{kk});
    fgetl(fid); fgetl(fid); fgetl(fid);
    SPINE = fgetl(fid);
    fclose(fid);
    SPINE = strsplit(SPINE,',');
    SPINE = SPINE(2:end); % First spot is blank
    
    %% Direction Cosine Matricies
    
    % Relative direction cosine matrix between adjacent levels
    % This has been vectorized to speed it up (5.5x faster)
    rel_dcm = zeros(f,9,n-1);
    for k = n-1:-1:1
     rel_dcm(:,:,k) = [dot(dcm(:,1:3,k),dcm(:,1:3,k+1),2),... APi
                       dot(dcm(:,4:6,k),dcm(:,1:3,k+1),2),... APj
                       dot(dcm(:,7:9,k),dcm(:,1:3,k+1),2),... APk
                       dot(dcm(:,1:3,k),dcm(:,4:6,k+1),2),... SIi
                       dot(dcm(:,4:6,k),dcm(:,4:6,k+1),2),... SIj
                       dot(dcm(:,7:9,k),dcm(:,4:6,k+1),2),... SIk
                       dot(dcm(:,1:3,k),dcm(:,7:9,k+1),2),... MLi
                       dot(dcm(:,4:6,k),dcm(:,7:9,k+1),2),... MLj
                       dot(dcm(:,7:9,k),dcm(:,7:9,k+1),2)]; % MLk
    end

    % Regional direction cosines
    s1  = find(strcmp(SPINE,regionalAngle(1)));
    t12 = find(strcmp(SPINE,regionalAngle(2)));
    c7  = find(strcmp(SPINE,regionalAngle(3)));
    
    nReg = zeros(2,3);
    if ~isempty(s1 & c7)
        nReg(:,1) = [s1;c7];
    end
    if ~isempty(t12 & c7)
        nReg(:,2) = [t12;c7];
    end 
    if ~isempty(s1 & t12)
        nReg(:,3) = [s1;t12];
    end
    
    reg_dcm = zeros(f,9,3);
    for k = 1:3
        if all(nReg(:,k)) % Don't try to compute if spine level is missing; in that case export zeros
     reg_dcm(:,:,k) = [dot( dcm(:,1:3,nReg(1,k)) , dcm(:,1:3,nReg(2,k)) ,2),... APi
                       dot( dcm(:,4:6,nReg(1,k)) , dcm(:,1:3,nReg(2,k)) ,2),... APj
                       dot( dcm(:,7:9,nReg(1,k)) , dcm(:,1:3,nReg(2,k)) ,2),... APk
                       dot( dcm(:,1:3,nReg(1,k)) , dcm(:,4:6,nReg(2,k)) ,2),... SIi
                       dot( dcm(:,4:6,nReg(1,k)) , dcm(:,4:6,nReg(2,k)) ,2),... SIj
                       dot( dcm(:,7:9,nReg(1,k)) , dcm(:,4:6,nReg(2,k)) ,2),... SIk
                       dot( dcm(:,1:3,nReg(1,k)) , dcm(:,7:9,nReg(2,k)) ,2),... MLi
                       dot( dcm(:,4:6,nReg(1,k)) , dcm(:,7:9,nReg(2,k)) ,2),... MLj
                       dot( dcm(:,7:9,nReg(1,k)) , dcm(:,7:9,nReg(2,k)) ,2)]; % MLk
        end % end of if statement
    end % end of loop through Regional Levels
    
    %% Calculate Angles
    
    % Absolute Angles
    angABS = zeros(f,3,n);
    for kvb = 1:n % Loop through spine levels
        %r1 is about Y (AT); r2 is about X (LB); r3 is about Z (FE);
        [r1,r2,r3] = dcm2angle(reshape(dcm(:,:,kvb)',3,3,f),'YXZ','ZeroR3');
        angABS(:,:,kvb) = [r1,r2,r3];
    end
    angABS = reshape(angABS,f,3*n);
    
    % Intersegmental Angles
    angREL = zeros(f,3,n-1);
    for kvb = 1:n-1 % Loop through spine levels
        %r1 is about Y (AT); r2 is about X (LB); r3 is about Z (FE);
        [r1,r2,r3] = dcm2angle(reshape(rel_dcm(:,:,kvb)',3,3,f),'YXZ','ZeroR3');
        angREL(:,:,kvb) = [r1,r2,r3];
    end
    angREL = reshape(angREL, f, 3*(n-1) );
    
    % Regional Angles (Total, Thoracic, Lumbar)
    angREG = zeros(f,3,size(nReg,2));
    for kvb = 1:size(nReg,2) % Loop through spine levels
        %r1 is about Y (AT); r2 is about X (LB); r3 is about Z (FE);
        [r1,r2,r3] = dcm2angle(reshape(reg_dcm(:,:,kvb)',3,3,f),'YXZ','ZeroR3');
        angREG(:,:,kvb) = [r1,r2,r3];
    end
    angREG = reshape(angREG, f, 3*size(nReg,2) );  
    
    %% Filter data
    
    out = [position,angABS,angREL,angREG];
    
    Fs = round(1/(t(2,2)-t(1,2))); % Sample Frequency
    order = 2;
    fc = 6; % MAKE SURE YOU SET THIS FOR YOUR STUDY
    [b,a] = butter(order,(fc/0.802)/(Fs/2));
    out = filtfilt(b,a,out);
    
    %% Rubberband REINSERTED
    
    out = rubberband(outFrames,[t,out]);
    
    %% Set up for Writing to File
    
    % Create folder if it doesn't exist
    [folder,fileName,ext] = fileparts( file{kk} ); 
    if exist([folder,writeToFolder],'dir')==0
        status = mkdir([folder,writeToFolder]);
        if status == 0
            disp('Problem creating folder')
        end
    end
    
    % Headers
    Header = cell(1,size(out,2));
    Header(1) = {'Frames,'}; Header(2) = {'Time (s),'};
    % Postion and Absolute Headers
    i = (1:3:3*n) + 2;
    for k = 1:n
        Header(i(k)+0) = {[SPINE{k},': Post,']};
        Header(i(k)+1) = {[SPINE{k},': Up,']};
        Header(i(k)+2) = {[SPINE{k},': Left,']};
        Header(i(k)+(n*3+0)) = {[SPINE{k},': AT,']};
        Header(i(k)+(n*3+1)) = {[SPINE{k},': LB,']};
        Header(i(k)+(n*3+2)) = {[SPINE{k},': FE,']};
    end
    % Relative Angle Headers
    i = (1:3:3 * (n-1) ) + n*6+2;
    for k = 1:n-1
        Header(i(k)+0) = {[SPINE{k},'/',SPINE{k+1},': AT,']};
        Header(i(k)+1) = {[SPINE{k},'/',SPINE{k+1},': LB,']};
        Header(i(k)+2) = {[SPINE{k},'/',SPINE{k+1},': FE,']};
    end
    % Regional Angles
    i = 1 + (3*n-1)*3+2;
    Header(i:end) = {'Total: AT,','Total: LB,','Total: FE,','Thoracic: AT,','Thoracic: LB,','Thoracic: FE,','Lumbar: AT,','Lumbar: LB,','Lumbar: FE\n'};
        
    
    %% Write to File
    
    fileOut = fullfile([folder,writeToFolder,fileName,ext]);
    fid = fopen(fileOut,'w');
    fprintf(fid,'File created with: LCS to Angles function by Derek Zwambag\n');
    fprintf(fid,['File: ',fileName,'\n']);
    fprintf(fid,['Date:,',date,'\n']); 
    fprintf(fid,cell2mat(Header));
    fclose(fid);
    dlmwrite(fileOut,out,'-append');
    
    disp([file{kk},': Complete']) % Display Status to User
end % end of loop through files

end % end of function 