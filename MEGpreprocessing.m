function [] = MEGpreprocessing()
% default pipeline for rejecting bad trials, which takes epoched raw data
% in Fieldtrip-format as input. 
%
% See also: Good practice for conducting and reporting MEG research
% http://www.sciencedirect.com/science/article/pii/S1053811912009895
%
% Anne Urai, 2016

addpath('~/Documents/fieldtrip');
ft_defaults;

% LOAD IN SUBJECT SPECIFICS AND READ DATA

% this file should be epoched
load('yourdata.mat');

% ==================================================================
% 1. REMOVE TRIALS WITH EXCESSIVE HEAD MOTION
% see http://www.fieldtriptoolbox.org/example/how_to_incorporate_head_movements_in_meg_analysis
% Head motion in reference to the start of your measurement can be computed 
% based on the location of the 3 fiducial coils. Trials with excessive head 
% movement can be detected using an outlier test (Grubb?s test). This will not 
% detect gradual changes in head position (e.g. sliding down in the helmet).
% ==================================================================

% compute head rotation wrt first trial
cc_rel = computeHeadRotation(data);

% plot the rotation of the head
subplot(4,4,cnt); cnt = cnt + 1;
plot(cc_rel); ylabel('HeadM');
axis tight; box off;

% find outliers
[~, idx] = deleteoutliers(cc_rel);
[t,~]    = ind2sub(size(cc_rel),idx);

% only take those where the deviation is more than 6 mm
t = t(any(abs(cc_rel(t, :)) > 6, 2));

% show those on the plot
hold on;
for thist = 1:length(t),
    plot([t(thist) t(thist)], [max(get(gca, 'ylim')) max(get(gca, 'ylim'))], 'k.');
end

% remove those trials
cfg                     = [];
cfg.trials              = true(1, length(data.trial));
cfg.trials(unique(t))   = false; % remove these trials
data                    = ft_selectdata(cfg, data);
fprintf('removing %d excessive head motion trials \n', length(find(cfg.trials == 0)));

subplot(4,4,cnt); cnt = cnt + 1;
if isempty(t),
    title('No motion'); axis off;
else
    % show head motion without those removed
    cc_rel = computeHeadRotation(data);
    
    % plot the rotation of the head
    plot(cc_rel); ylabel('Motion resid');
    axis tight; box off;
end

% plot a quick power spectrum
% save those cfgs for later plotting
cfgfreq             = [];
cfgfreq.method      = 'mtmfft';
cfgfreq.output      = 'pow';
cfgfreq.taper       = 'hanning';
cfgfreq.channel     = 'MEG';
cfgfreq.foi         = 1:130;
cfgfreq.keeptrials  = 'no';
freq                = ft_freqanalysis(cfgfreq, data);

% plot those data and save for visual inspection
subplot(4,4,cnt); cnt = cnt + 1;
loglog(freq.freq, freq.powspctrm, 'linewidth', 0.1); hold on;
loglog(freq.freq, mean(freq.powspctrm), 'k', 'linewidth', 1);
axis tight; axis square; box off;
set(gca, 'xtick', [10 50 100], 'tickdir', 'out', 'xticklabel', []);

% ==================================================================
% 2. REMOVE TRIALS WITH EYEBLINKS (only during beginning of trial)
% Bandpass filter the vertical EOG channel between 1-15 Hz and z-transform 
% this filtered time course. Select complete trials that exceed a threshold of  
% z =4 (alternatively you can set the z-threshold per data file or per subject 
% with the ?interactive? mode in ft_artifact_zvalue function). Reject trials 
% that contain blink artifacts before going on to the next step. For monitoring 
% purposes, plot the time courses of your trials before and after blink rejection.
% ==================================================================

cfg                              = [];
cfg.continuous                   = 'no'; % data has been epoched

% channel selection, cutoff and padding
cfg.artfctdef.zvalue.channel     = {'EOGV'};

% 001, 006, 0012 and 0018 are the vertical and horizontal eog chans
cfg.artfctdef.zvalue.trlpadding  = 0; % avoid filter edge artefacts by setting to negative
cfg.artfctdef.zvalue.fltpadding  = 0;
cfg.artfctdef.zvalue.artpadding  = 0.05; % go a bit to the sides of blinks

% algorithmic parameters
cfg.artfctdef.zvalue.bpfilter   = 'yes';
cfg.artfctdef.zvalue.bpfilttype = 'but';
cfg.artfctdef.zvalue.bpfreq     = [1 15];
cfg.artfctdef.zvalue.bpfiltord  = 4;
cfg.artfctdef.zvalue.hilbert    = 'yes';

% set cutoff
cfg.artfctdef.zvalue.cutoff     = 2; % to detect all blinks, be strict
% cfg.artfctdef.zvalue.interactive = 'yes';
[~, artifact_eog]               = ft_artifact_zvalue(cfg, data);

cfg                             = [];
cfg.artfctdef.reject            = 'complete';
cfg.artfctdef.eog.artifact      = artifact_eog;

% reject blinks only when they occur between fix and stim offset
crittoilim = [ data.trialinfo(:,2) - data.trialinfo(:,1) - 0.4*data.fsample ...
    data.trialinfo(:,5) - data.trialinfo(:,1) + 0.8*data.fsample]  / data.fsample;
cfg.artfctdef.crittoilim        = crittoilim;
data                            = ft_rejectartifact(cfg, data);

% ==================================================================
% 3. REMOVE TRIALS WITH SACCADES (only during beginning of trial)
% Remove trials with (horizontal) saccades (EOGH). Use the same settings as 
% for the EOGV-based blinks detection. The z-threshold can be set a bit higher 
% (z = [4 6]). Reject all trials that contain saccades before going further.
% ==================================================================

cfg                              = [];
cfg.continuous                   = 'no'; % data has been epoched

% channel selection, cutoff and padding
cfg.artfctdef.zvalue.channel     = {'EOGH'};

% 001, 006, 0012 and 0018 are the vertical and horizontal eog chans
cfg.artfctdef.zvalue.trlpadding  = 0; % padding doesnt work for data thats already on disk
cfg.artfctdef.zvalue.fltpadding  = 0; % 0.2; this crashes the artifact func!
cfg.artfctdef.zvalue.artpadding  = 0.05; % go a bit to the sides of blinks

% algorithmic parameters
cfg.artfctdef.zvalue.bpfilter   = 'yes';
cfg.artfctdef.zvalue.bpfilttype = 'but';
cfg.artfctdef.zvalue.bpfreq     = [1 15];
cfg.artfctdef.zvalue.bpfiltord  = 4;
cfg.artfctdef.zvalue.hilbert    = 'yes';

% set cutoff
cfg.artfctdef.zvalue.cutoff     = 4;
% cfg.artfctdef.zvalue.interactive = 'yes';
[~, artifact_eog]               = ft_artifact_zvalue(cfg, data);

cfg                             = [];
% cfg.artfctdef.reject            = 'complete';
cfg.artfctdef.eog.artifact      = artifact_eog;

% reject blinks when they occur between ref and the offset of the stim
crittoilim = [data.trialinfo(:,2) - data.trialinfo(:,1) - 0.4*data.fsample ...
    data.trialinfo(:,5) - data.trialinfo(:,1) + 0.8*data.fsample] / data.fsample;
cfg.artfctdef.crittoilim        = crittoilim;
data                            = ft_rejectartifact(cfg, data);

% ==================================================================
% 4. REMOVE TRIALS WITH JUMPS
% Compute the power spectrum of all trials and a linear line on the loglog-
% transformed power spectrum. Jumps cause broad range increase in the power 
% spectrum so trials containing jumps can be selected by detecting outliers 
% in the intercepts of the fitted lines (using Grubb?s test for outliers).
% ==================================================================

% detrend and demean
cfg             = [];
cfg.detrend     = 'yes';
cfg.demean      = 'yes';
data            = ft_preprocessing(cfg, data);

% get the fourier spectrum per trial and sensor
cfgfreq.keeptrials  = 'yes';
freq                = ft_freqanalysis(cfgfreq, data);

% compute the intercept of the loglog fourier spectrum on each trial
disp('searching for trials with squid jumps...');
intercept       = nan(size(freq.powspctrm, 1), size(freq.powspctrm, 2));
x = [ones(size(freq.freq))' log(freq.freq)'];

for t = 1:size(freq.powspctrm, 1),
    for c = 1:size(freq.powspctrm, 2),
        b = x\log(squeeze(freq.powspctrm(t,c,:)));
        intercept(t,c) = b(1);
    end
end

% detect jumps as outliers
[~, idx] = deleteoutliers(intercept(:));
subplot(4,4,cnt); cnt = cnt + 1;
if isempty(idx),
    fprintf('no squid jump trials found \n');
    title('No jumps'); axis off;
else
    fprintf('removing %d squid jump trials \n', length(unique(t)));
    [t,~] = ind2sub(size(intercept),idx);
    
    % remove those trials
    cfg                 = [];
    cfg.trials          = true(1, length(data.trial));
    cfg.trials(unique(t)) = false; % remove these trials
    data                = ft_selectdata(cfg, data);
    
    % plot the spectrum again
    cfgfreq.keeptrials = 'no';
    freq            = ft_freqanalysis(cfgfreq, data);
    loglog(freq.freq, freq.powspctrm, 'linewidth', 0.1); hold on;
    loglog(freq.freq, mean(freq.powspctrm), 'k', 'linewidth', 1);
    axis tight; axis square; box off;
    set(gca, 'xtick', [10 50 100], 'tickdir', 'out', 'xticklabel', []);
    title(sprintf('%d jumps removed', length(unique(t))));
end

% ==================================================================
% 5. REMOVE LINE NOISE
% ==================================================================

cfg             = [];
cfg.bsfilter    = 'yes';
cfg.bsfreq      = [49 51; 99 101; 149 151];
data            = ft_preprocessing(cfg, data);

% plot power spectrum
freq            = ft_freqanalysis(cfgfreq, data);
subplot(4,4,cnt); cnt = cnt + 1;
loglog(freq.freq, freq.powspctrm, 'linewidth', 0.5); hold on;
loglog(freq.freq, mean(freq.powspctrm), 'k', 'linewidth', 1);
axis tight; ylim(ylims); axis square; box off;
title('After bandstop');
set(gca, 'xtick', [10 50 100], 'tickdir', 'out', 'xticklabel', []);

% ==================================================================
% 6. REMOVE CARS BASED ON THRESHOLD
% Cars moving past the MEG lab cause big slow signal changes. Trials 
% containing these artifacts can be selected and removed by computing 
% the maximum range of the data for every trial. Trials with a larger 
% range than a threshold (standard = 0.75e-11) can be rejected (the standard 
% threshold might be low if you have long trials).
% ==================================================================

disp('Looking for CAR artifacts...');
cfg = [];
cfg.trials = true(1, length(data.trial));
worstChanRange = nan(1, length(data.trial));
for t = 1:length(data.trial),
    % compute the range as the maximum of the peak-to-peak values within each channel
    ptpval = max(data.trial{t}, [], 2) - min(data.trial{t}, [], 2);
    % determine range and index of 'worst' channel
    worstChanRange(t) = max(ptpval);
end

% default range for peak-to-peak
artfctdef.range           = 0.75e-11;

% decide whether to reject this trial
cfg.trials = (worstChanRange < artfctdef.range);
fprintf('removing %d CAR trials \n', length(find(cfg.trials == 0)));
data = ft_selectdata(cfg, data);

% ==================================================================
% 7. REMOVE TRIALS WITH MUSCLE BURSTS BEFORE RESPONSE
% Remove muscle using the same z-value-based approach as for the eye 
% channels. Filter the data between 110-140 Hz and use a z-value threshold of 10.
% ==================================================================

cfg                              = [];
cfg.continuous                   = 'no'; % data has been epoched

% channel selection, cutoff and padding
cfg.artfctdef.zvalue.channel     = {'MEG'}; % make sure there are no NaNs
cfg.artfctdef.zvalue.trlpadding  = 0;
cfg.artfctdef.zvalue.fltpadding  = 0; % 0.2; - this crashes ft_artifact_zvalue!
cfg.artfctdef.zvalue.artpadding  = 0.1;

% algorithmic parameters
cfg.artfctdef.zvalue.bpfilter    = 'yes';
cfg.artfctdef.zvalue.bpfreq      = [110 140];
cfg.artfctdef.zvalue.bpfiltord   = 9;
cfg.artfctdef.zvalue.bpfilttype  = 'but';
cfg.artfctdef.zvalue.hilbert     = 'yes';
cfg.artfctdef.zvalue.boxcar      = 0.2;

% set cutoff
cfg.artfctdef.zvalue.cutoff      = 20;
[~, artifact_muscle]             = ft_artifact_zvalue(cfg, data);

cfg                              = [];
cfg.artfctdef.reject             = 'complete';
cfg.artfctdef.muscle.artifact    = artifact_muscle;

% only remove muscle bursts before the response
crittoilim = [data.trialinfo(:,1) - data.trialinfo(:,1) ...
    data.trialinfo(:,9) - data.trialinfo(:,1)]  ./ data.fsample;
cfg.artfctdef.crittoilim        = crittoilim;
data                            = ft_rejectartifact(cfg, data);

% plot final power spectrum
freq            = ft_freqanalysis(cfgfreq, data);
subplot(4,4,cnt); 
loglog(freq.freq, freq.powspctrm, 'linewidth', 0.5); hold on;
loglog(freq.freq, mean(freq.powspctrm), 'k', 'linewidth', 1);
axis tight; axis square; box off; ylim(ylims);
set(gca, 'xtick', [10 50 100], 'tickdir', 'out');

end % rec


function cc_rel = computeHeadRotation(data)

% take only head position channels
cfg         = [];
cfg.channel = {'HLC0011','HLC0012','HLC0013', ...
    'HLC0021','HLC0022','HLC0023', ...
    'HLC0031','HLC0032','HLC0033'};
hpos        = ft_selectdata(cfg, data);

% calculate the mean coil position per trial
coil1 = nan(3, length(hpos.trial));
coil2 = nan(3, length(hpos.trial));
coil3 = nan(3, length(hpos.trial));

for t = 1:length(hpos.trial),
    coil1(:,t) = [mean(hpos.trial{1,t}(1,:)); mean(hpos.trial{1,t}(2,:)); mean(hpos.trial{1,t}(3,:))];
    coil2(:,t) = [mean(hpos.trial{1,t}(4,:)); mean(hpos.trial{1,t}(5,:)); mean(hpos.trial{1,t}(6,:))];
    coil3(:,t) = [mean(hpos.trial{1,t}(7,:)); mean(hpos.trial{1,t}(8,:)); mean(hpos.trial{1,t}(9,:))];
end

% calculate the headposition and orientation per trial (function at the
% bottom of this script)
cc     = circumcenter(coil1, coil2, coil3);

% compute relative to the first trial
cc_rel = [cc - repmat(cc(:,1),1,size(cc,2))]';
cc_rel = 1000*cc_rel(:, 1:3); % translation in mm

end
