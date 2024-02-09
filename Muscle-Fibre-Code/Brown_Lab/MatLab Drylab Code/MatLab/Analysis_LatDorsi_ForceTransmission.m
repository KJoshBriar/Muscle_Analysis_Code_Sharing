clear

EMG = csvread('c:\Data\LatissimusDorsiStudy\P02_Exports\Pilot2_t2hip_001_tracked_NIDAQ.csv',15,3); 
Motion = csvread('c:\Data\LatissimusDorsiStudy\P02_Exports\Pilot2_t2hip_001_tracked.csv',4,176);

%Need to enter these two values based on MVC trials
MVCvalue_Channel1 = 0.8392;
MVCvalue_Channel2 = 2.6329;

MVCvalue = [MVCvalue_Channel1 MVCvalue_Channel2];

order = 2;
samplerate = 1920;
motionsamplerate = 120;

%EMGchannels: UpperLats = 1 LowerLats = 2

[r,c] = size(EMG);
[r2,c2] = size(Motion);

TrialLength = r2 / motionsamplerate; % length of trial in seconds

%MVCtriallength = r / samplerate; %(time in seconds)

%MVC1 = MVC1(16:r,4:5); MVC2 = MVC2(16:r2,4:5); MVC3 = MVC3(16:r3,4:5);
figure(1);
plot(EMG(:,1))
%hold on
figure(2);
plot(EMG(:,2))
%hold off

%remove bias

BEMG = EMG - ones(r,1)*mean(EMG);

%rectify and low-pass filter at 2.5 Hz
REMG = abs(BEMG);   

[B,A] = butter(order,2.5/(samplerate/2));

LEMG = filter(B,A,REMG);
 
%normalize to MVC
NEMG = (LEMG ./ MVCvalue) * 100;  

figure(3);
plot(NEMG(:,1))
figure(4);
plot(NEMG(:,2))

%finding mean EMG and spine angles over each of the first and second phases

EMG_mean_phase1 = mean(NEMG((samplerate+1):(samplerate*4+1),:));
EMG_mean_phase2 = mean(NEMG((r-(samplerate*4+1)):(r-(samplerate+1)),:));

Motion_mean_phase1 = mean(Motion((motionsamplerate+1):(motionsamplerate*4+1),:));
Motion_mean_phase2 = mean(Motion((r2-(motionsamplerate*4+1)):(r2-(motionsamplerate+1)),:));

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




