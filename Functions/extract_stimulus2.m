 function stim = extract_stimulus2(EEG, PATH, label, k, sbj,task)
% This function extract the stimulus for the variance partitioning model 
%Input:
%EEG    STRUCT the whole EEG structure
%PATH   STR the path where the feature stimulus matrix is saved
%label  STR labeling the feature of interest
%   label options:
%   onset               onsets derived from onset detection master
%   envelope            envelope
%   envelope onset      focus on the gain of the envelope
%   mTRF envelope       envelope based on the mTRFenvelope function
%   onsenv              onsets and envelope combined
%   onsmenv             onsets and mTRF envelope
%   ensenvmenv          onsets, envelope, and mTRF envelope
%   alarm               alarm sound extracted from experimental labels
%   irregular           irregular sounds extracted from experimental labels
%   odd                 every unmasked odd tone
%   oddinalarm          odd tones masked by alarm tones
%   oddinirr            odd tones masked by irregular tones
%   alarmenv            both alarm and envelope information 
%   random              returns an empty array of the length of the data
%   alarmons            combines both OnsetDetection onsets and alarm
%                       onsets
%   irrons              irregular and normal onsets(cleaned)
%   oddons              odd and automatic onsets(cleaned)
%   irralarmons         irregular alarm and automatic onsets(cleaned)
%   irroddons           irregular odd and automatic onsets(cleaned)
%   alarmoddons         odd alarm and automatic onsets(cleaned)
%   irralarmoddons      everything
%
%k      the condition 
%
%
%Output:
%stim   vector of the stimulus that you want

TriggerSoundDelay = 0.019;

%sets the selection threshold off onsets to be removed
thresh = 25;
%load the stimulus matrix
cd(PATH)
if k == 1
    filename = sprintf('stims_narrow_%s.mat',sbj);
    load(filename)
    [audioIn,fs] =audioread(sprintf('narrow_audio_game_%s.wav',sbj));
    maudio = miraudio(sprintf('narrow_audio_game_%s.wav',sbj),'Sampling',24000);
else
    filename = sprintf('stims_wide_%s.mat',sbj);
    load(filename)
    [audioIn,fs] =audioread(sprintf('wide_audio_game_%s.wav',sbj));
    maudio = miraudio(sprintf('wide_audio_game_%s.wav',sbj),'Sampling',24000);

end

switch label
    
    case 'onset'
        stim = stims(1,:)';
    case 'envelope'
        temp_dat= stims(2,:)';
        stim = (temp_dat - min(temp_dat,[],'all')) / ( max(temp_dat,[],'all') - min(temp_dat,[],'all') );
    case 'envelope onset'
        stim = stims(3,:)';
    case 'mTRF envelope'
        temp_dat= stims(4,:)';
        stim = (temp_dat - min(temp_dat,[],'all')) / ( max(temp_dat,[],'all') - min(temp_dat,[],'all') );
    case 'onsenv'
        stim_ons = stims(1,:)';
        temp_dat = stims(2,:)';
        stim_env = (temp_dat - min(temp_dat,[],'all')) / ( max(temp_dat,[],'all') - min(temp_dat,[],'all') );
        stim = [stim_ons stim_env];
        
    case 'onsmenv'
        stim_ons = stims(1,:)';
        temp_dat = stims(4,:)';
        stim_env = (temp_dat - min(temp_dat,[],'all')) / ( max(temp_dat,[],'all') - min(temp_dat,[],'all') );
        stim = [stim_ons stim_env];
    case 'onsenvmenv'
        stim_ons = stims(1,:)';
        temp_dat = stims(2,:)';
        stim_env = (temp_dat - min(temp_dat,[],'all')) / ( max(temp_dat,[],'all') - min(temp_dat,[],'all') );
        temp_dat = stims(4,:)';
        stim_menv = (temp_dat - min(temp_dat,[],'all')) / ( max(temp_dat,[],'all') - min(temp_dat,[],'all') );
        
        stim = [stim_ons stim_env stim_menv];
        
    case 'alarm'
        stim = zeros(1,EEG.pnts);
        
        stim_label = 'alarm_tone_post';
        %find the indicies of events
        aoidx = strncmp(stim_label,{EEG.event.type} ,strlength(stim_label));
        aoidxx = find(double(aoidx)~=0);
        
        %get the latencies of the events
        train_lat = zeros(size(aoidxx,2),1);
        for lt = 1:length(aoidxx)
            train_lat(lt,:) = EEG.event(aoidxx(lt)).latency + TriggerSoundDelay*EEG.srate;
            %generate the feature vector
            stim(1,round(train_lat(lt,:))) = 1;
        end
        stim = stim';
        
    case 'irregular'
        stim = zeros(1,EEG.pnts);
        
        stim_label = 'irr_sound';
        
        %find the indicies of events
        aoidx = strncmp(stim_label,{EEG.event.type} ,strlength(stim_label));
        aoidxx_f = find(double(aoidx)~=0);
        
        %take only every second value
        aoidxx = aoidxx_f(2:2:end);
        
        %get the latencies of the events
        train_lat = zeros(size(aoidxx,2),1);
        for lt = 1:length(aoidxx)
            train_lat(lt,:) = EEG.event(aoidxx(lt)).latency + TriggerSoundDelay*EEG.srate;
            %generate the feature vector
            stim(1,round(train_lat(lt,:))) = 1;
        end
        stim = stim';
        
    case 'odd'
        stim = zeros(1,EEG.pnts);
        
        stim_label = 'odd_tone';
        
        %find the indicies of events
        aoidx = strncmp(stim_label,{EEG.event.type} ,strlength(stim_label));
        aoidxx_f = find(double(aoidx)~=0);
        
        %remove double values
        aoidxx = aoidxx_f(2:2:end);

        %get the latencies of the events
        train_lat = zeros(size(aoidxx,2),1);
        for lt = 1:length(aoidxx)
            train_lat(lt,:) = EEG.event(aoidxx(lt)).latency + TriggerSoundDelay*EEG.srate;
            %generate the feature vector
            stim(1,round(train_lat(lt,:))) = 1;
        end
        stim = stim';
        
        
    case 'odd_in_alarm'
        stim = zeros(1,EEG.pnts);
        
        stim_label = 'odd_in_alarm';
        
        %find the indicies of events
        aoidx = strncmp(stim_label,{EEG.event.type} ,strlength(stim_label));
        aoidxx = find(double(aoidx)~=0);
        
        %get the latencies of the events
        train_lat = zeros(size(aoidxx,2),1);
        for lt = 1:length(aoidxx)
            train_lat(lt,:) = EEG.event(aoidxx(lt)).latency + TriggerSoundDelay*EEG.srate;
            %generate the feature vector
            stim(1,round(train_lat(lt,:))) = 1;
        end
        stim = stim';
        
        
    case 'odd_in_irr'
        stim = zeros(1,EEG.pnts);
        
        stim_label = 'odd_in_irr';
        
        %find the indicies of events
        aoidx = strncmp(stim_label,{EEG.event.type} ,strlength(stim_label));
        aoidxx = find(double(aoidx)~=0);
        
        %get the latencies of the events
        train_lat = zeros(size(aoidxx,2),1);
        for lt = 1:length(aoidxx)
            train_lat(lt,:) = EEG.event(aoidxx(lt)).latency + TriggerSoundDelay*EEG.srate;
            %generate the feature vector
            stim(1,round(train_lat(lt,:))) = 1;
        end
        stim = stim';
        
        
    case 'alarmenv'
        stim_env = stims(2,:)';
        stim_env = (stim_env - min(stim_env,[],'all')) / ( max(stim_env,[],'all') - min(stim_env,[],'all') );
        stim = zeros(1,EEG.pnts);
        
        stim_label = 'alarm_tone_post';
        %find the indicies of events
        aoidx = strncmp(stim_label,{EEG.event.type} ,strlength(stim_label));
        aoidxx = find(double(aoidx)~=0);
        
        %get the latencies of the events
        train_lat = zeros(size(aoidxx,2),1);
        for lt = 1:length(aoidxx)
            train_lat(lt,:) = EEG.event(aoidxx(lt)).latency + TriggerSoundDelay*EEG.srate;
            %generate the feature vector
            stim(1,round(train_lat(lt,:))) = 1;
        end
        stim = stim';
        stim = cat(2,stim,stim_env);
%% onset area               
    case 'random'
        nr_ons = sum(stims(1,:));
        rand_idx = randperm(EEG.pnts,nr_ons);
        stim = zeros(EEG.pnts,1);
        stim(rand_idx) = 1;
        
    case 'alarmons'
        stim_ons = stims(1,:)';
        stim_alarm = zeros(1,EEG.pnts);
        
        stim_label = 'alarm_tone_post';
        %find the indicies of events
        aoidx = strncmp(stim_label,{EEG.event.type} ,strlength(stim_label));
        aoidxx = find(double(aoidx)~=0);
        
        %get the latencies of the events
        train_lat = zeros(size(aoidxx,2),1);
        for lt = 1:length(aoidxx)
            train_lat(lt,:) = EEG.event(aoidxx(lt)).latency + TriggerSoundDelay*EEG.srate;
            %generate the feature vector
            stim_alarm(1,round(train_lat(lt,:))) = 1;
        end

         
        stim = cat(2,stim_alarm',stim_ons);
        
    case 'irrons'
        stim_ons = stims(1,:)';
        stim_irr = zeros(1,EEG.pnts);
        
        stim_label = 'irr_sound';
        %find the indicies of events
        aoidx = strncmp(stim_label,{EEG.event.type} ,strlength(stim_label));
        aoidxx_f = find(double(aoidx)~=0);
        
        aoidxx = aoidxx_f(2:2:end);
        %get the latencies of the events
        train_lat = zeros(size(aoidxx,2),1);
        for lt = 1:length(aoidxx)
            train_lat(lt,:) = EEG.event(aoidxx(lt)).latency + TriggerSoundDelay*EEG.srate;
            %generate the feature vector
            stim_irr(1,round(train_lat(lt,:))) = 1;
        end
        
       
        
        stim = cat(2,stim_irr',stim_ons);
        
        
    case 'oddons'
        stim_ons = stims(1,:)';
        stim_odd = zeros(1,EEG.pnts);
        
        stim_label = 'odd_tone';
        %find the indicies of events
        aoidx = strncmp(stim_label,{EEG.event.type} ,strlength(stim_label));
        aoidxx_f = find(double(aoidx)~=0);
        
        aoidxx = aoidxx_f(2:2:end);
        %get the latencies of the events
        train_lat = zeros(size(aoidxx,2),1);
        for lt = 1:length(aoidxx)
            train_lat(lt,:) = EEG.event(aoidxx(lt)).latency + TriggerSoundDelay*EEG.srate;
            %generate the feature vector
            stim_odd(1,round(train_lat(lt,:))) = 1;
        end
        
        
        
        stim = cat(2,stim_odd',stim_ons);
        
    case 'oddirr'
        stim_oddirr = zeros(2,EEG.pnts);
        
        stim_label = {'odd_tone','irr_sound'};
        for stl = 1:length(stim_label)
            %find the indicies of events
            aoidx = strncmp(stim_label{stl},{EEG.event.type} ,strlength(stim_label{stl}));
            aoidxx_f = find(double(aoidx)~=0);
            
            aoidxx = aoidxx_f(2:2:end);
            %get the latencies of the events
            train_lat = zeros(size(aoidxx,2),1);
            for lt = 1:length(aoidxx)
                train_lat(lt,:) = EEG.event(aoidxx(lt)).latency + TriggerSoundDelay*EEG.srate;
                %generate the feature vector
                stim_oddirr(stl,round(train_lat(lt,:))) = 1;
            end
            
        end
        stim = stim_oddirr';
        
    case 'oddalarm'
        stim_oddalarm = zeros(2,EEG.pnts);
        
        stim_label = {'odd_tone','alarm_tone'};
        for stl = 1:length(stim_label)
            %find the indicies of events
            aoidx = strncmp(stim_label{stl},{EEG.event.type} ,strlength(stim_label{stl}));
            aoidxx_f = find(double(aoidx)~=0);
            
            aoidxx = aoidxx_f(2:2:end);
            %get the latencies of the events
            train_lat = zeros(size(aoidxx,2),1);
            for lt = 1:length(aoidxx)
                train_lat(lt,:) = EEG.event(aoidxx(lt)).latency + TriggerSoundDelay*EEG.srate;
                %generate the feature vector
                stim_oddalarm(stl,round(train_lat(lt,:))) = 1;
            end

        end
        
        stim = stim_oddalarm';
        
    case 'irralarm'
        stim_irralarm = zeros(2,EEG.pnts);
        
        stim_label = {'irr_sound','alarm_tone'};
        for stl = 1:length(stim_label)
            %find the indicies of events
            aoidx = strncmp(stim_label{stl},{EEG.event.type} ,strlength(stim_label{stl}));
            aoidxx_f = find(double(aoidx)~=0);
            
            aoidxx = aoidxx_f(2:2:end);
            %get the latencies of the events
            train_lat = zeros(size(aoidxx,2),1);
            for lt = 1:length(aoidxx)
                train_lat(lt,:) = EEG.event(aoidxx(lt)).latency + TriggerSoundDelay*EEG.srate;
                %generate the feature vector
                stim_irralarm(stl,round(train_lat(lt,:))) = 1;
            end
        end
        stim = stim_irralarm';
        
    case 'oddirrons'
        stim_ons = stims(1,:)';
        ons_idx = find(stim_ons==1);

        stim_oddirr = zeros(2,EEG.pnts);
        
        stim_label = {'odd_tone','irr_sound'};
        for stl = 1:length(stim_label)
            %find the indicies of events
            aoidx = strncmp(stim_label{stl},{EEG.event.type} ,strlength(stim_label{stl}));
            aoidxx_f = find(double(aoidx)~=0);
            
            aoidxx = aoidxx_f(2:2:end);
            %get the latencies of the events
            train_lat = zeros(size(aoidxx,2),1);
            for lt = 1:length(aoidxx)
                train_lat(lt,:) = EEG.event(aoidxx(lt)).latency + TriggerSoundDelay*EEG.srate;
                %generate the feature vector
                stim_oddirr(stl,round(train_lat(lt,:))) = 1;
            end
            
            
           
        end
        %concat the stimuls vectors
        stim = cat(2,stim_oddirr',stim_ons);
        
        
    case 'oddalarmons'
        stim_ons = stims(1,:)';
        ons_idx = find(stim_ons==1);
        
        stim_oddalarm = zeros(2,EEG.pnts);
        
        stim_label = {'odd_tone','alarm_tone'};
        for stl = 1:length(stim_label)
            %find the indicies of events
            aoidx = strncmp(stim_label{stl},{EEG.event.type} ,strlength(stim_label{stl}));
            aoidxx_f = find(double(aoidx)~=0);
            
            aoidxx = aoidxx_f(2:2:end);
            %get the latencies of the events
            train_lat = zeros(size(aoidxx,2),1);
            for lt = 1:length(aoidxx)
                train_lat(lt,:) = EEG.event(aoidxx(lt)).latency + TriggerSoundDelay*EEG.srate;
                %generate the feature vector
                stim_oddalarm(stl,round(train_lat(lt,:))) = 1;
            end
            
            
            
        end
        %concat the stimuls vectors
        stim = cat(2,stim_oddalarm',stim_ons);
        
        
    case 'irralarmons'
        stim_ons = stims(1,:)';
        ons_idx = find(stim_ons==1);
        
        stim_irralarm = zeros(2,EEG.pnts);
        
        stim_label = {'irr_sound','alarm_tone'};
        for stl = 1:length(stim_label)
            %find the indicies of events
            aoidx = strncmp(stim_label{stl},{EEG.event.type} ,strlength(stim_label{stl}));
            aoidxx_f = find(double(aoidx)~=0);
            
            aoidxx = aoidxx_f(2:2:end);
            %get the latencies of the events
            train_lat = zeros(size(aoidxx,2),1);
            for lt = 1:length(aoidxx)
                train_lat(lt,:) = EEG.event(aoidxx(lt)).latency + TriggerSoundDelay*EEG.srate;
                %generate the feature vector
                stim_irralarm(stl,round(train_lat(lt,:))) = 1;
            end
            
            
            
        end
        %concat the stimuls vectors
        stim = cat(2,stim_irralarm',stim_ons);
        
    case 'oddirralarm'
                
        stim_irralarmodd = zeros(3,EEG.pnts);
        
        stim_label = {'odd_tone','irr_sound','alarm_tone'};
        for stl = 1:length(stim_label)
            %find the indicies of events
            aoidx = strncmp(stim_label{stl},{EEG.event.type} ,strlength(stim_label{stl}));
            aoidxx_f = find(double(aoidx)~=0);
            
            aoidxx = aoidxx_f(2:2:end);
            %get the latencies of the events
            train_lat = zeros(size(aoidxx,2),1);
            for lt = 1:length(aoidxx)
                train_lat(lt,:) = EEG.event(aoidxx(lt)).latency + TriggerSoundDelay*EEG.srate;
                %generate the feature vector
                stim_irralarmodd(stl,round(train_lat(lt,:))) = 1;
            end
           
        end
        %concat the stimuls vectors
        stim = stim_irralarmodd';
        
        
    case 'oddirralarmons'
        stim_ons = stims(1,:)';
        ons_idx = find(stim_ons==1);
        
        stim_irralarmodd = zeros(3,EEG.pnts);
        
        stim_label = {'odd_tone','irr_sound','alarm_tone'};
        for stl = 1:length(stim_label)
            %find the indicies of events
            aoidx = strncmp(stim_label{stl},{EEG.event.type} ,strlength(stim_label{stl}));
            aoidxx_f = find(double(aoidx)~=0);
            
            aoidxx = aoidxx_f(2:2:end);
            %get the latencies of the events
            train_lat = zeros(size(aoidxx,2),1);
            for lt = 1:length(aoidxx)
                train_lat(lt,:) = EEG.event(aoidxx(lt)).latency + TriggerSoundDelay*EEG.srate;
                %generate the feature vector
                stim_irralarmodd(stl,round(train_lat(lt,:))) = 1;
            end
            
            
            
        end
        %concat the stimuls vectors
        stim = cat(2,stim_irralarmodd',stim_ons);
        
%% Mel area
    case 'mel'
        
        mfilt = mirfilterbank(maudio,'Mel')
        mirverbose(0);
        mirwaitbar(0);
        menv = mirenvelope(maudio,'Spectro','Mel','Sampling',100);
        mel_spec = squeeze(mirgetdata(menv));
        
        
        save(sprintf('mel_spec_%s_%s.mat',sbj,task{k}),'mel_spec');
        stim = mel_spec;
        
        load(sprintf('mel_spec_%s_%s.mat',sbj,task{k}),'mel_spec');
        temp_dat=mel_spec;
        stim = (temp_dat - min(temp_dat,[],'all')) / ( max(temp_dat,[],'all') - min(temp_dat,[],'all') );
       
    case 'melons'
        load(sprintf('mel_spec_%s_%s.mat',sbj,task{k}),'mel_spec');
        stim_ons= stims(1,:)';
        temp_dat = mel_spec;
        stim_mel = (temp_dat - min(temp_dat,[],'all')) / ( max(temp_dat,[],'all') - min(temp_dat,[],'all') );
        
        stim = [stim_mel stim_ons(1:size(mel_spec,1),:)];
        
    case 'melenv'
        load(sprintf('mel_spec_%s_%s.mat',sbj,task{k}),'mel_spec');
        temp_dat = mel_spec;
        mel_spec = (temp_dat - min(temp_dat,[],'all')) / ( max(temp_dat,[],'all') - min(temp_dat,[],'all') );
        temp_dat = stims(2,:);
        stim_env = (temp_dat - min(temp_dat,[],'all')) / ( max(temp_dat,[],'all') - min(temp_dat,[],'all') );
        stim = cat(2,mel_spec,stim_env(1,1:size(mel_spec,1))');
        
    case 'melmenv'
        load(sprintf('mel_spec_%s_%s.mat',sbj,task{k}),'mel_spec');
        temp_dat = mel_spec;
        mel_spec = (temp_dat - min(temp_dat,[],'all')) / ( max(temp_dat,[],'all') - min(temp_dat,[],'all') );
        temp_dat = stims(4,:);
        stim_env = (temp_dat - min(temp_dat,[],'all')) / ( max(temp_dat,[],'all') - min(temp_dat,[],'all') );
        stim = cat(2,mel_spec,stim_env(1,1:size(mel_spec,1))');
        
    case 'melonsenv'
        load(sprintf('mel_spec_%s_%s.mat',sbj,task{k}),'mel_spec');
        stim_ons = stims(1,:)';
        temp_dat = stims(2,:)';
        stim_env = (temp_dat - min(temp_dat,[],'all')) / ( max(temp_dat,[],'all') - min(temp_dat,[],'all') );
        temp_dat = mel_spec;
        stim_mel = (temp_dat - min(temp_dat,[],'all')) / ( max(temp_dat,[],'all') - min(temp_dat,[],'all') );
        
        stim = [stim_mel stim_ons(1:size(mel_spec,1),:) stim_env(1:size(mel_spec,1),:)];
        
    case 'melonsmenv'
        load(sprintf('mel_spec_%s_%s.mat',sbj,task{k}),'mel_spec');
        stim_ons = stims(1,:)';
        temp_dat = stims(4,:)';
        stim_env = (temp_dat - min(temp_dat,[],'all')) / ( max(temp_dat,[],'all') - min(temp_dat,[],'all') );
        temp_dat = mel_spec;
        stim_mel = (temp_dat - min(temp_dat,[],'all')) / ( max(temp_dat,[],'all') - min(temp_dat,[],'all') );
        
        stim = [stim_mel stim_ons(1:size(mel_spec,1),:) stim_env(1:size(mel_spec,1),:)];
        
    case 'alarmmel'
        load(sprintf('mel_spec_%s_%s.mat',sbj,task{k}),'mel_spec');
        temp_dat = mel_spec;
        stim_mel = (temp_dat - min(temp_dat,[],'all')) / ( max(temp_dat,[],'all') - min(temp_dat,[],'all') );
        stim = zeros(1,EEG.pnts);
        
        stim_label = 'alarm_tone_post';
        %find the indicies of events
        aoidx = strncmp(stim_label,{EEG.event.type} ,strlength(stim_label));
        aoidxx = find(double(aoidx)~=0);
        
        %get the latencies of the events
        train_lat = zeros(size(aoidxx,2),1);
        for lt = 1:length(aoidxx)
            train_lat(lt,:) = EEG.event(aoidxx(lt)).latency + TriggerSoundDelay*EEG.srate;
            %generate the feature vector
            stim(1,round(train_lat(lt,:))) = 1;
        end
        stim = stim';
        stim = cat(2,stim_mel,stim(1:size(mel_spec,1),:));
        
    case 'irrmel'
        load(sprintf('mel_spec_%s_%s.mat',sbj,task{k}),'mel_spec');
        temp_dat = mel_spec;
        stim_mel =  (temp_dat - min(temp_dat,[],'all')) / ( max(temp_dat,[],'all') - min(temp_dat,[],'all') );
        stim = zeros(1,EEG.pnts);
        
        stim_label = 'irr_sound';
        
        %find the indicies of events
        aoidx = strncmp(stim_label,{EEG.event.type} ,strlength(stim_label));
        aoidxx_f = find(double(aoidx)~=0);
        
        %remove double values
        aoidxx = aoidxx_f(2:2:end);

        %get the latencies of the events
        train_lat = zeros(size(aoidxx,2),1);
        for lt = 1:length(aoidxx)
            train_lat(lt,:) = EEG.event(aoidxx(lt)).latency + TriggerSoundDelay*EEG.srate;
            %generate the feature vector
            stim(1,round(train_lat(lt,:))) = 1;
        end
        stim = stim';
        stim = cat(2,stim_mel,stim(1:size(mel_spec,1),:));
        
    case 'oddmel'
        load(sprintf('mel_spec_%s_%s.mat',sbj,task{k}),'mel_spec');
        temp_dat = mel_spec;
        stim_mel = (temp_dat - min(temp_dat,[],'all')) / ( max(temp_dat,[],'all') - min(temp_dat,[],'all') );
        stim = zeros(1,EEG.pnts);

        stim_label = 'odd_tone';
        
        %find the indicies of events
        aoidx = strncmp(stim_label,{EEG.event.type} ,strlength(stim_label));
        aoidxx_f = find(double(aoidx)~=0);
        
        %remove double values
        aoidxx = aoidxx_f(2:2:end);
        
        %get the latencies of the events
        train_lat = zeros(size(aoidxx,2),1);
        for lt = 1:length(aoidxx)
            train_lat(lt,:) = EEG.event(aoidxx(lt)).latency + TriggerSoundDelay*EEG.srate;
            %generate the feature vector
            stim(1,round(train_lat(lt,:))) = 1;
        end
        stim = stim';
        stim = cat(2,stim_mel,stim(1:size(mel_spec,1),:));
        
    case 'oddalarmmel'
        load(sprintf('mel_spec_%s_%s.mat',sbj,task{k}),'mel_spec');
        temp_dat = mel_spec;
        stim_mel = (temp_dat - min(temp_dat,[],'all')) / ( max(temp_dat,[],'all') - min(temp_dat,[],'all') );
        
        stim_irralarm = zeros(2,EEG.pnts);
        stim_label = {'irr_sound','alarm_tone'};
        for stl = 1:length(stim_label)
            %find the indicies of events
            aoidx = strncmp(stim_label{stl},{EEG.event.type} ,strlength(stim_label{stl}));
            aoidxx_f = find(double(aoidx)~=0);
            
            aoidxx = aoidxx_f(2:2:end);
            %get the latencies of the events
            train_lat = zeros(size(aoidxx,2),1);
            for lt = 1:length(aoidxx)
                train_lat(lt,:) = EEG.event(aoidxx(lt)).latency + TriggerSoundDelay*EEG.srate;
                %generate the feature vector
                stim_irralarm(stl,round(train_lat(lt,:))) = 1;
            end
        end 
        stim = stim_irralarm';
        stim = cat(2,stim_mel,stim(1:size(mel_spec,1),:));
        
    case 'oddirrmel'
        load(sprintf('mel_spec_%s_%s.mat',sbj,task{k}),'mel_spec');
        temp_dat = mel_spec;
        stim_mel = (temp_dat - min(temp_dat,[],'all')) / ( max(temp_dat,[],'all') - min(temp_dat,[],'all') );
        
        stim_irralarm = zeros(2,EEG.pnts);
        stim_label = {'odd_tone','irr_sound'};
        for stl = 1:length(stim_label)
            %find the indicies of events
            aoidx = strncmp(stim_label{stl},{EEG.event.type} ,strlength(stim_label{stl}));
            aoidxx_f = find(double(aoidx)~=0);
            
            aoidxx = aoidxx_f(2:2:end);
            %get the latencies of the events
            train_lat = zeros(size(aoidxx,2),1);
            for lt = 1:length(aoidxx)
                train_lat(lt,:) = EEG.event(aoidxx(lt)).latency + TriggerSoundDelay*EEG.srate;
                %generate the feature vector
                stim_irralarm(stl,round(train_lat(lt,:))) = 1;
            end
        end
        stim = stim_irralarm';
        stim = cat(2,stim_mel,stim(1:size(mel_spec,1),:));
        
    case 'irralarmmel'
        load(sprintf('mel_spec_%s_%s.mat',sbj,task{k}),'mel_spec');
        temp_dat = mel_spec;
        stim_mel = (temp_dat - min(temp_dat,[],'all')) / ( max(temp_dat,[],'all') - min(temp_dat,[],'all') );
        
        stim_irralarm = zeros(2,EEG.pnts);
        stim_label = {'irr_sound','alarm_tone'};
        for stl = 1:length(stim_label)
            %find the indicies of events
            aoidx = strncmp(stim_label{stl},{EEG.event.type} ,strlength(stim_label{stl}));
            aoidxx_f = find(double(aoidx)~=0);
            
            aoidxx = aoidxx_f(2:2:end);
            %get the latencies of the events
            train_lat = zeros(size(aoidxx,2),1);
            for lt = 1:length(aoidxx)
                train_lat(lt,:) = EEG.event(aoidxx(lt)).latency + TriggerSoundDelay*EEG.srate;
                %generate the feature vector
                stim_irralarm(stl,round(train_lat(lt,:))) = 1;
            end
        end
        stim = stim_irralarm';
        stim = cat(2,stim_mel,stim(1:size(mel_spec,1),:));
        
        
        
    case 'oddirralarmmel'
        load(sprintf('mel_spec_%s_%s.mat',sbj,task{k}),'mel_spec');
        temp_dat = mel_spec;
        stim_mel = (temp_dat - min(temp_dat,[],'all')) / ( max(temp_dat,[],'all') - min(temp_dat,[],'all') );
        stim_oddirralarm = zeros(3,EEG.pnts);
        
        stim_label = {'odd_tone','irr_sound','alarm_tone'};
        for stl = 1:length(stim_label)
            %find the indicies of events
            aoidx = strncmp(stim_label{stl},{EEG.event.type} ,strlength(stim_label{stl}));
            aoidxx_f = find(double(aoidx)~=0);
            
            aoidxx = aoidxx_f(2:2:end);
            %get the latencies of the events
            train_lat = zeros(size(aoidxx,2),1);
            for lt = 1:length(aoidxx)
                train_lat(lt,:) = EEG.event(aoidxx(lt)).latency + TriggerSoundDelay*EEG.srate;
                %generate the feature vector
                stim_oddirralarm(stl,round(train_lat(lt,:))) = 1;
            end
        end
        stim = stim_oddirralarm';
        stim = cat(2,stim_mel,stim(1:size(mel_spec,1),:));
        
        
%% mTRF envelope area 
    case 'alarmmenv'
        stim_env = stims(4,:)';
        stim_env = (stim_env - min(stim_env,[],'all')) / ( max(stim_env,[],'all') - min(stim_env,[],'all') );
        stim = zeros(1,EEG.pnts);
        
        stim_label = 'alarm_tone_post';
        %find the indicies of events
        aoidx = strncmp(stim_label,{EEG.event.type} ,strlength(stim_label));
        aoidxx = find(double(aoidx)~=0);
        
        %get the latencies of the events
        train_lat = zeros(size(aoidxx,2),1);
        for lt = 1:length(aoidxx)
            train_lat(lt,:) = EEG.event(aoidxx(lt)).latency + TriggerSoundDelay*EEG.srate;
            %generate the feature vector
            stim(1,round(train_lat(lt,:))) = 1;
        end
        stim = stim';
        stim = cat(2,stim,stim_env);
        
    case 'oddmenv'
        stim_env = stims(4,:)';
        stim_env = (stim_env - min(stim_env,[],'all')) / ( max(stim_env,[],'all') - min(stim_env,[],'all') );
        stim = zeros(1,EEG.pnts);
        
        stim_label = 'odd_tone';
        
        %find the indicies of events
        aoidx = strncmp(stim_label,{EEG.event.type} ,strlength(stim_label));
        aoidxx_f = find(double(aoidx)~=0);
        
        %remove double values
        aoidxx = aoidxx_f(2:2:end);

        %get the latencies of the events
        train_lat = zeros(size(aoidxx,2),1);
        for lt = 1:length(aoidxx)
            train_lat(lt,:) = EEG.event(aoidxx(lt)).latency + TriggerSoundDelay*EEG.srate;
            %generate the feature vector
            stim(1,round(train_lat(lt,:))) = 1;
        end
        stim = stim';
        stim = cat(2,stim,stim_env);
        
    case 'irrmenv'
        stim_env = stims(4,:)';
        stim_env = (stim_env - min(stim_env,[],'all')) / ( max(stim_env,[],'all') - min(stim_env,[],'all') );
        stim = zeros(1,EEG.pnts);
        
        stim_label = 'irr_sound';
        
        %find the indicies of events
        aoidx = strncmp(stim_label,{EEG.event.type} ,strlength(stim_label));
        aoidxx_f = find(double(aoidx)~=0);
        
        %remove double values
        aoidxx = aoidxx_f(2:2:end);

        %get the latencies of the events
        train_lat = zeros(size(aoidxx,2),1);
        for lt = 1:length(aoidxx)
            train_lat(lt,:) = EEG.event(aoidxx(lt)).latency + TriggerSoundDelay*EEG.srate;
            %generate the feature vector
            stim(1,round(train_lat(lt,:))) = 1;
        end
        stim = stim';
        stim = cat(2,stim,stim_env);
        
    case 'irralarmmenv'
        stim_env = stims(4,:)';
        stim_env = (stim_env - min(stim_env,[],'all')) / ( max(stim_env,[],'all') - min(stim_env,[],'all') );
        stim_irralarm = zeros(2,EEG.pnts);
        
        stim_label = {'irr_sound','alarm_tone'};
        for stl = 1:length(stim_label)
            %find the indicies of events
            aoidx = strncmp(stim_label{stl},{EEG.event.type} ,strlength(stim_label{stl}));
            aoidxx_f = find(double(aoidx)~=0);
            
            aoidxx = aoidxx_f(2:2:end);
            %get the latencies of the events
            train_lat = zeros(size(aoidxx,2),1);
            for lt = 1:length(aoidxx)
                train_lat(lt,:) = EEG.event(aoidxx(lt)).latency + TriggerSoundDelay*EEG.srate;
                %generate the feature vector
                stim_irralarm(stl,round(train_lat(lt,:))) = 1;
            end
        end 
        stim = stim_irralarm';
        stim = cat(2,stim,stim_env);
        
    case 'oddirrmenv'
        stim_env = stims(4,:)';
        stim_env = (stim_env - min(stim_env,[],'all')) / ( max(stim_env,[],'all') - min(stim_env,[],'all') );
        stim_oddirr = zeros(2,EEG.pnts);
        
        stim_label = {'odd_tone','irr_sound'};
        for stl = 1:length(stim_label)
            %find the indicies of events
            aoidx = strncmp(stim_label{stl},{EEG.event.type} ,strlength(stim_label{stl}));
            aoidxx_f = find(double(aoidx)~=0);
            
            aoidxx = aoidxx_f(2:2:end);
            %get the latencies of the events
            train_lat = zeros(size(aoidxx,2),1);
            for lt = 1:length(aoidxx)
                train_lat(lt,:) = EEG.event(aoidxx(lt)).latency + TriggerSoundDelay*EEG.srate;
                %generate the feature vector
                stim_oddirr(stl,round(train_lat(lt,:))) = 1;
            end
        end
        stim = stim_oddirr';
        stim = cat(2,stim,stim_env);
        
    case 'oddalarmmenv'
        stim_env = stims(4,:)';
        stim_env = (stim_env - min(stim_env,[],'all')) / ( max(stim_env,[],'all') - min(stim_env,[],'all') );
        stim_oddalarm = zeros(2,EEG.pnts);
        
        stim_label = {'odd_tone','alarm_tone'};
        for stl = 1:length(stim_label)
            %find the indicies of events
            aoidx = strncmp(stim_label{stl},{EEG.event.type} ,strlength(stim_label{stl}));
            aoidxx_f = find(double(aoidx)~=0);
            
            aoidxx = aoidxx_f(2:2:end);
            %get the latencies of the events
            train_lat = zeros(size(aoidxx,2),1);
            for lt = 1:length(aoidxx)
                train_lat(lt,:) = EEG.event(aoidxx(lt)).latency + TriggerSoundDelay*EEG.srate;
                %generate the feature vector
                stim_oddalarm(stl,round(train_lat(lt,:))) = 1;
            end
        end
        stim = stim_oddalarm';
        stim = cat(2,stim,stim_env);
        
    case 'oddirralarmmenv'
        stim_env = stims(4,:)';
        stim_env = (stim_env - min(stim_env,[],'all')) / ( max(stim_env,[],'all') - min(stim_env,[],'all') );
        stim_oddirralarm = zeros(3,EEG.pnts);
        
        stim_label = {'odd_tone','irr_sound','alarm_tone'};
        for stl = 1:length(stim_label)
            %find the indicies of events
            aoidx = strncmp(stim_label{stl},{EEG.event.type} ,strlength(stim_label{stl}));
            aoidxx_f = find(double(aoidx)~=0);
            
            aoidxx = aoidxx_f(2:2:end);
            %get the latencies of the events
            train_lat = zeros(size(aoidxx,2),1);
            for lt = 1:length(aoidxx)
                train_lat(lt,:) = EEG.event(aoidxx(lt)).latency + TriggerSoundDelay*EEG.srate;
                %generate the feature vector
                stim_oddirralarm(stl,round(train_lat(lt,:))) = 1;
            end
        end
        stim = stim_oddirralarm';
        stim = cat(2,stim,stim_env);
        

end
        
        
        

                