clear

MVC1 = csvread('c:\Data\LatissimusDorsiStudy\P02_Exports\Pilot2_MVC_003_tracked_NIDAQ.csv',15,3); 
MVC2 = csvread('c:\Data\LatissimusDorsiStudy\P02_Exports\Pilot2_MVC_003_tracked_NIDAQ.csv',15,3); 
MVC3 = csvread('c:\Data\LatissimusDorsiStudy\P02_Exports\Pilot2_MVC_003_tracked_NIDAQ.csv',15,3); 

order = 2;
samplerate = 1920; %in Hz

%channels: UpperLats = 1 LowerLats = 2
%MVC1 = MVC1(:,1); MVC2 = MVC2(:,1); MVC3 = MVC3(:,1);
%MVC4 = MVC4(:,2); MVC5 = MVC5(:,2); MVC6 = MVC6(:,2);

[r,c] = size(MVC1); [r2,c2] = size(MVC2); [r3,c3] = size(MVC3);

%MVCtriallength = r / samplerate; %(time in seconds)

%MVC1 = MVC1(16:r,4:5); MVC2 = MVC2(16:r2,4:5); MVC3 = MVC3(16:r3,4:5);
%plots both EMG channels for all three MVC trials for visualization
figure(1);
plot(MVC1(:,1))
figure(2);
plot(MVC1(:,2))
figure(3);
plot(MVC2(:,1))
figure(4);
plot(MVC2(:,2))
figure(5);
plot(MVC3(:,1))
figure(6);
plot(MVC3(:,2))

%remove bias

BMVC1 = MVC1 - ones(r,1)*mean(MVC1);
BMVC2 = MVC2 - ones(r2,1)*mean(MVC2);
BMVC3 = MVC3 - ones(r3,1)*mean(MVC3);

%rectify and low-pass filter at 2.5 Hz
RMVC1 = abs(BMVC1); 
RMVC2 = abs(BMVC2); 
RMVC3 = abs(BMVC3);  

[B,A] = butter(order,2.5/(samplerate/2));

LMVC1 = filter(B,A,RMVC1);
LMVC2 = filter(B,A,RMVC2);
LMVC3 = filter(B,A,RMVC3);
 
%find max point in each trial
MMVC1 = max(LMVC1);  
MMVC2 = max(LMVC2);  
MMVC3 = max(LMVC3);  

MVCArray = [MMVC1; MMVC2; MMVC3];
FinalMVCs = max(MVCArray);

%plot(LMVC4(:,4))

%ensure that the MVCs come from the correct MVC trials
%MAXarrayRightNeck = [MMVC1(:,1); MMVC2(:,1); MMVC3(:,1)];
%[MaxvalRneck,indexMAXvalRneck] = max(MAXarrayRightNeck);   %index allows us to determine what MVC trial the max values came from

%MAXarrayLeftNeck = [MMVC1(:,2); MMVC2(:,2); MMVC3(:,2)];
%[MaxvalLneck,indexMAXvalLneck] = max(MAXarrayLeftNeck);

%MAXarrayRightTES = [MMVC4(:,3); MMVC5(:,3); MMVC6(:,3)];
%[MaxvalRTES,indexMAXvalRTES] = max(MAXarrayRightTES);

%MAXarrayLeftTES = [MMVC4(:,4); MMVC5(:,4); MMVC6(:,4)];
%[MaxvalLTES,indexMAXvalLTES] = max(MAXarrayLeftTES);

%MAXarrayRLES = [MMVC4(:,5); MMVC5(:,5); MMVC6(:,5)];
%[MaxvalRLES,indexMAXvalRLES] = max(MAXarrayRLES);

%MAXarrayLLES = [MMVC4(:,6); MMVC5(:,6); MMVC6(:,6)];
%[MaxvalLLES,indexMAXvalLLES] = max(MAXarrayLLES);


%AllMaxes = [MaxvalRneck MaxvalLneck MaxvalRTES MaxvalLTES MaxvalRLES MaxvalLLES];

%dlmwrite('c:\Data\Kristin_deskheightstudy\MVCValues.csv',AllMaxes);




