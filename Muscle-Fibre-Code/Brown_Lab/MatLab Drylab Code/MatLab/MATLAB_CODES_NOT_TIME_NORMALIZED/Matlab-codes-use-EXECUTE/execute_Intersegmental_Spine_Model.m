% Run the Intersegmental Analysis Program

% You should create an M-file for each study. This is where you will set up
% your data and call all the required functions.  

% Don't run using the F5 or run command. You'll want to do each line
% separately in order to check for error messages. 

%% Indicate your raw files stored on the computer.
%folder = 'C:/Users/Devon/Documents/Devon_School/UofG/Research/ROWING_2018-2019/analysis/V3D/';

folder = 'C:\Cathrine\FLEXION_EXTENSION\TEST_CODES\FE_P04_FOR_MATLAB';


% If you want to run a single file you can replace 'folder' with 'file' in
% all the function calls below.
%file = 'C:/Users/Devon/Documents/Devon_School/UofG/Research/ROWING_2018-2019/analysis/V3D/TEST_EXPORT_REGULAR.csv';
%file= 'C:/Users/Devon/Documents/Devon_School/UofG/Research/ROWING_2018-2019/analysis/V3D/SHIFTED_EXPORT.csv';

% If you have lots of subjects in the same folder and you only want to run
% a few files in the folder you can use 'nameString' which will read in all
% the files in 'folder' that contain 'nameString'. The '*' is a wildcard
% and can be replaced by anything. 

% e.g. to run all the .csv files in a folder simply use 
% nameString = % '*.csv';
nameString = '*.csv'; % This will run all the trials for subject 'FT_13'

% Check if your folder exists.

% Step 1: Clean the data
Correct_Marker_Blips(folder,'nameString',nameString);
%Correct_Marker_Blips(file);

% For this step you'll be shown marker displacements and accelerations. If
% there are large blips then you can choose to interpolate. 
% Using ginput click on either side of the blip (as close as you can) and
% then input which channel needs to be cleaned. You can only interpolate
% one region per channel at a time, but you can click 'interpolate' as many
% times in a row as you'd like. 

%%
% Step 2: Sort the Markers
Motive_Spine_Marker_Sorting(folder,'sortOrder',[2,1],'nameString',nameString);

% You don't have to really do anything at this step, except watch to see if
% there were any problems with the data. This program is expecting 3xN data
% points where N is the number of spine levels collected from. If there are
% not 3xN levels (maybe you had a marker fall off or is missing) then an
% warninng message will pop up and the program will continue. You will need
% to go back and fix that trial. 

% If you have extra markers, go into Motive, delete the extra marker and
% then reexport. You'll need to re-correct Marker blips.

% If you have missing markers you can tell the program this using the
% following example line. 

% file = 'C:\Users\Devon\Documents\Devon_School\UofG\Research\ROWING_2018-2019\analysis\INTERSEGMENTAL_SPINE\P08_R1_INTERSEG.csv';
% Motive_Spine_Marker_Sorting(file,'missingMrk',[3,12],'visualize','on');
%Motive_Spine_Marker_Sorting(file,'sortOrder',[2,1],'visualize','on'); %Sort data when starting in flexed posture.

% In this example there is a missing marker from the 3rd column (most
% positive ML) and the 12th row from the bottom. If you had to missing
% markers you could specify [3,12; 2,15];
% 'visualize','on' means that a figure will pop up and data will not be
% written to file until you accept that it is correct. 
%%
% Step 3: Run the Spine Model.
folder = fullfile(folder, 'Sorted');
%folder = [folder,'Sorted/'];
Intersegmental_Spine_Model(folder,'outFreq',120,'nameString', nameString);
%%
% Step 4: Export Angles
folder = fullfile(folder, 'IntersegmentalLCS');
LCS_to_Angles(folder,'nameString',nameString);

% 'regionalAngle' is used to determine which regions define Total,
% Thoracic, and Lumbar. 
