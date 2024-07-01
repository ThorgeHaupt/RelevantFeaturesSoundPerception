function [EEG,PATH] = OT_preprocessing(s,k,sbj,lp)
%% Does the preprocessing for Onset tracking project
%
%Input:
%s = subject number
%k = condition number
%sbj = array of participant indicies
%
%
%
%Output:
%EEG = EEG structure (completely preprocessed)
%
%
%Notes:
%requires all the other previous pre-processing scripts to have run

%set up the path to the data
OT_setup

%start eeglab for no reason
PATH = fullfile(DATAPATH, sprintf(sbj{s}),filesep);
file = {sprintf('%s_ica_narrow_game_added_trigger.set',sbj{s}),sprintf('%s_ica_wide_game_added_trigger.set',sbj{s})};
    


[ALLEEG EEG CURRENTSET ALLCOM] = eeglab('nogui');
%load the ica set
EEG = pop_loadset(fullfile(PATH,file{k}));

%reject the absolute shit out of this data set
EEG = pop_icflag(EEG, [NaN NaN;0.7 1;0.7 1;0.6 1;0.7 1;0.7 1;NaN NaN]);
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
EEG = pop_subcomp(EEG, [])
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);


%% start the preprocessing

EEG = pop_firws(EEG, 'fcutoff', lp, 'ftype', 'lowpass', 'wtype', 'rectangular', 'forder',100)%,'plotfresp',1);
EEG = pop_resample(EEG,100);
EEG = pop_firws(EEG, 'fcutoff', 0.3, 'ftype', 'highpass', 'wtype', 'hann',...
    'forder', 518)%,'plotfresp',1);
EEG = pop_interp(EEG,EEG.urchanlocs,'spherical');
EEG = pop_reref(EEG,{'TP9','TP10'},'keepref','on');
EEG = pop_select(EEG,'nochannel',{'TP9','TP10'});
