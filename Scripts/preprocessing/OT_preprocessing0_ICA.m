%% Preprocessing - ICA weights
OT_setup

%% the loop was only done to derive the ICA weights
for s=1:length(sbj)
    sb = sbj{s};
    %setup path file
    PATH = fullfile(DATAPATH, sb,filesep);
    cd(PATH)
    [ALLEEG EEG CURRENTSET ALLCOM] = eeglab('nogui');
    file = {'narrow_game_added_trigger.set','wide_game_added_trigger.set'};  %set up name
    badchan = [];

    for i=1:2
        EEG = pop_loadset(fullfile(PATH,file{i}));
        if i ==1
            urchanlocs = EEG.chanlocs;
        end
        
        ALLEEG{i} = EEG;
    end
    
    %merge the data set
    EEG = pop_mergeset( ALLEEG{1}, ALLEEG{2} );
    
    %downsample
    EEG.data = double(EEG.data);
    EEG = pop_resample(EEG,250);
    EEG.data = double(EEG.data);
    
    %highpass
    EEG = pop_firws(EEG, 'fcutoff', 1, 'ftype', 'highpass', 'wtype', 'hann',...
        'forder', 568);
    
    %lowpass
    EEG = pop_firws(EEG, 'fcutoff', 42, 'ftype', 'lowpass', 'wtype', 'hann',...
        'forder', 128);
    
    %clean the channels 
    EEG = clean_channels(EEG);
    
    if isfield(EEG.chaninfo,'removedchans')
        badchan = EEG.chaninfo.removedchans;
        chan_rej(s,1) = size(badchan,2);
    else
        chan_rej(s,1) = 0;
        
    end
    
%     %reference to average
%     EEG = pop_reref(EEG,[]);
   
    
    %run ICA
    EEG = eeg_regepochs(EEG, 'recurrence', 1);
    EEG = eeg_checkset(EEG, 'eventconsistency');
    EEG.data = double(EEG.data);
    
    % remove epochs with artefacts to improve ICA training
    PRUNE = 3;
    EEG = pop_jointprob(EEG, 1, [1:size(EEG.data,1)], PRUNE, PRUNE, 0, 1, 0);
    EEG = eeg_checkset(EEG, 'eventconsistency');

    EEG.data = double(EEG.data);

    % compute ICA weights
    EEG = pop_runica(EEG, 'icatype', 'runica','extended',1,'interrupt','on','concatenate','on');
    
    
    icawinv = EEG.icawinv;
    icasphere = EEG.icasphere;
    icaweights = EEG.icaweights;
    icachansind = EEG.icachansind;

    
    clear global 
    [ALLEEG EEG CURRENTSET ALLCOM] = eeglab('nogui');

    %load the two data sets and add the ICA weights
    for i=1:2
        EEG = pop_loadset(fullfile(PATH,file{i}));
        if ~isempty(badchan)
            badch = {badchan.labels};
            EEG = pop_select(EEG,'nochannel',badch);
        end
        EEG = eeg_checkset(EEG, 'eventconsistency');
        EEG.urchanlocs = urchanlocs;
        
        
        EEG.icawinv = icawinv;
        EEG.icasphere = icasphere;
        EEG.icaweights = icaweights;
        EEG.icachansind = icachansind;
        
        EEG = eeg_checkset( EEG );
        EEG = pop_iclabel(EEG, 'default');
        [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
        pop_saveset(EEG,'filename',sprintf('%s_ica_%s',sbj{s},file{i}),'filepath',PATH);
    end
    clear badchan badch
end