% Microsaccade_detection_example

% The example script below shows how you can use the functions of the
% EYE-EEG toolbox to detect (micro)saccades. I used srate = 250 Hz, which
% seems to work well. The authors of the detection algorithm used 500 Hz. 
% In the example data file, large eyemovements have not been interpolated 
% or removed, this should be done first if you only want to detect 
% microsaccades.
% The toolbox van be found here: www2.hu-berlin.de/eyetracking-eeg/
%
% Thomas Meindertsma, December 2016

%% step 1 Import data

load('Example_data_Microsaccade_detection.mat');
srate = data.fsample;

%% step 2 Calculate threshold for microsaccades in x- and y-signalls

% - smoothlevel 0: no smoothing, simple diff()
% - smoothlevel 1: 3-point window
% - smoothlevel 2: 5-point window
smoothlevel = 2;

for itrl = 1:length(data.trial)
    l = data.trial{itrl}(1:2,:)';
    vl = vecvel(l,srate,smoothlevel);
    [l_msdx(itrl) l_msdy(itrl)] = velthresh(vl);
end

%% step 3 Detect saccades

vfac = 4; % velocity factor (scales with microsaccade detection threshold)
mindur = 4; % Engbert & Mergenthaler (2006) use 6 samples with SR = 500.

for itrl = 1:length(data.trial)
    l = data.trial{itrl}(1:2,:)';
    vl = vecvel(l,srate,smoothlevel); % transpose time x chan
    sac{itrl} = microsacc_plugin(l,vl,vfac,mindur,nanmean(l_msdx),nanmean(l_msdy));
    
    
    %% merge nearby saccades
    
    clusterdist= 5; %value in sampling points that defines the minimum allowed fixation duration between two saccades.
    clustermode = 4; %4= clustered saccades are merged into one longer saccade.
    
    sac{itrl} = saccpar([sac{itrl} sac{itrl}]); % calculates binocular microsaccades (necassary for function mergesacc)
    sac{itrl} = mergesacc(sac{itrl},l,clusterdist,clustermode); % merges microsaccade clusters.
end

% columns of [sac]:
% 1: saccade onset (sample)
% 2: saccade offset (sample)
% 3: duration (samples)
% 4: delay between eyes (samples)
% 5: vpeak (peak velocity)
% 6: saccade "distance" (eucly. dist. between start and endpoint)
% 7: saccade angle (based on saccade "distance")
% 8: saccade "amplitude" (eucly. dist. of min/max in full saccade trajectory)
% 9: saccade angle (based on saccade "amplitude")



