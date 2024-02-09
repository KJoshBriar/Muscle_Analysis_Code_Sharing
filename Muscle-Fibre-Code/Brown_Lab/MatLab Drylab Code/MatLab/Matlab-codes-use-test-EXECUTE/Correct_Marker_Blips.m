function Correct_Marker_Blips(filePath, varargin)

% This function should be run after data have been exported from Motive.
% It is used to make sure that there are no marker blips of switches in the
% data. It computes the instantaneous marker acceleration (total of all
% three vector components). If acceleration is greater than 'thresh'
% then a figure will pop up and data will be rewritten to file.

% filePath is a string specifiying the path to either a directory or file.
%          If filePath is a directory it will perform Correct_Marker_Blips
%          on every file containting the nameString
% nameString (optional) will only read in files containing the nameString 
%          (default: '*.csv')
% thresh (optional) is the threshold for determining if user input is
%          required. If maximum aL is < threshold, data do not need
%          cleaning. Default is 10. Set to -1 to force visualization of all


% Written by Derek Zwambag April 14, 2018
% Version 2.0

%% Set up Input Parser
p = inputParser;

% filePath must be a file or a directory
validationFcn = @(s) assert( exist(s,'dir') | exist(s,'file') ,'filePath is not valid');
addRequired(p,'filePath',validationFcn); 
% nameString 
default = '*.csv';
addParameter(p,'nameString',default);
% thresh
validationFcn = @(d) assert( isnumeric(d) & isscalar(d),'thresh must be a scalar number');
default = 10;
addParameter(p,'thresh',default,validationFcn);

% Parse Inputs
parse(p,filePath,varargin{:});
filePath   = p.Results.filePath;
nameString = p.Results.nameString;
thresh     = p.Results.thresh;

%% Read in all files contained in Folder that include nameString

if exist(filePath,'dir')
    files = dir( fullfile(filePath,nameString) );
    for k = numel(files):-1:1
        file{k} = fullfile(filePath,files(k).name);
    end
else
    file = {filePath};
end

%% Correct Marker Blips by looping through each file

% Create Figure once
figure('position',[40,140,1100,550]);
subplot(4,1,1:2); ylabel('Marker Displacement (m)')
subplot(4,1,3); ylabel('Marker Acceleration (m/s^2')
subplot(4,1,4); xlabel('Marker Number'); ylabel('Max Acceleration')
for kk = numel(file):-1:1

    % Read in data (Motive has 7 lines of headers)
    data = csvread( file{kk}, 7, 0);
    % V3D has 5 lines of header
    %data = csvread( file{kk}, 5, 0);
    t = data(:,1:2);   % Frames and Time
    m = data(:,3:end); % Marker data
    [~,n] = size(m);   % Number of columns
    
    dt = t(2,2) - t(1,2); % sample period
    
    d = [zeros(1,n);diff(m)]; % component length between frames
    dL = sqrt(d(:,1:3:end).^2 + d(:,2:3:end).^2 + d(:,3:3:end).^2);       % marker displacement
    vL = [zeros(1,n/3);dL(3:end,:)-dL(1:end-2,:);zeros(1,n/3)] ./ (2*dt); % marker velocity
    aL = [zeros(1,n/3);vL(3:end,:)-vL(1:end-2,:);zeros(1,n/3)] ./ (2*dt); % marker acceleration
     
    % Only show figure and ask for input if max(aL) is greater than
    % threshold or the user asks for it.
    if max(aL(:)) > thresh 
        
        counter = 0;
        choice = 0;
        while choice ~= 1
            
            % plot data
            subplot(4,1,1:2); plot(dL)
            subplot(4,1,3);   plot(abs(aL)) 
            subplot(4,1,4);   plot(max(abs(aL)))
            set(gca,'xtick',(1:n/3),'xaxislocation','top')
            [~,fName,~] = fileparts(file{kk});
            fName(fName=='_') = ' '; % Drop underscores which cause subscript
            title(fName)

            choice = menu('Do you want to interpolate?','No','Yes','Quit Run');
            if choice == 3;
                disp(['Quit by User on file: ',fName]);
                return
            elseif choice == 2;
                [range,~] = ginput(2);
                range = floor(range(1)):floor(range(2));
                channel = cell2mat(inputdlg('What channel?'));
                channel = str2double(channel);

                idxin = [range(1),range(end)];
                in = m(idxin,channel*3-2:channel*3);
                xout = interp1(idxin,in(:,1),range,'pchip');
                yout = interp1(idxin,in(:,2),range,'pchip');
                zout = interp1(idxin,in(:,3),range,'pchip');

                m(range,channel*3-2:channel*3) = [xout;yout;zout]';
                
                % recalculate accelerations for next iteration through
                % while loop
                d = [zeros(1,n);diff(m)]; % component length between frames
                dL = sqrt(d(:,1:3:end).^2 + d(:,2:3:end).^2 + d(:,3:3:end).^2);     % marker displacement
                vL = [zeros(1,n/3);dL(3:end,:)-dL(1:end-2,:);zeros(1,n/3)] ./ (2*dt); % marker velocity
                aL = [zeros(1,n/3);vL(3:end,:)-vL(1:end-2,:);zeros(1,n/3)] ./ (2*dt); % marker acceleration
                
                counter = counter+1; % Used to determine if any changes were made.
            end % choice to interpolate
        end % end while loop  

        %% Write to file if any changes were made
        if any(counter) 
            
            % Get File Headers...
            fid = fopen(file{kk});
            header = cell(7,1);
            for i = 1:7
                header{i} = fgetl(fid);
            end
            fclose(fid);
            
            % Append header
            header{1} =[header{1},',Cleaned with Correct_Marker_Blips.m,By Derek Zwambag'];

            fid = fopen(file{kk},'w');
            for i = [1,3:7]
                fprintf(fid,[header{i},'\n']);
            end
            fclose(fid);
            dlmwrite(file{kk},[t,m],'-append');
        end % End write to file loop. 
        
    end % end if statement
    disp([file{kk},': Complete']) % Display Status to User
end % end loop through files

close all % Close function

end % end function