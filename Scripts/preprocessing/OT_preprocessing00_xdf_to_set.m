%transform that data in the EEG.set files

OT_setup
for s = 1:length(sbj)
     %move to the path where the data is to be saved
     cd(DATAPATH)
     %check if a folder structure exists, and if not, make one
     if ~exist(sbj{s},'dir')
         mkdir(sbj{s})
         cd(sbj{s})
     else
         cd(sbj{s});
     end
     
    for k = 1:length(task)

        %load the EEG data
        EEG = pop_loadxdf([raw_data_path,'sub-',sbj{s},'\ses-S001\eeg\sub-',sbj{s},'_ses-S001_task-',task{k},'_run-001_eeg.xdf'], 'streamtype', 'EEG', 'exclude_markerstreams', {});

        %load the channel layout here
        if size(EEG.data,1)>24
            EEG = pop_chanedit(EEG, 'load',{[chan_path,'mobile24_gyro.elp'],'filetype','autodetect'});
        else
            EEG = pop_chanedit(EEG, 'load',{[chan_path,'mobile24.elp'],'filetype','autodetect'});
        end
        
        EEG = pop_select( EEG, 'channel',{'Fp1','Fp2','Fz','F7','F8','FC1','FC2','Cz','C3','C4','T7','T8','CPz','CP1','CP2','CP5','CP6','TP9','TP10','Pz','P3','P4','O1','O2'});
        
        
        
        %load the audio_stream
        audio_strct = pop_loadxdf([raw_data_path, 'sub-',sbj{s},'\ses-S001\eeg\sub-',sbj{s},'_ses-S001_task-',task{k},'_run-001_eeg.xdf'],...
            'streamname', 'TetrisAudio', 'exclude_markerstreams', {});
        
        %process the audio 
        audio1 = double(audio_strct.data);     % Datapoints are stored in int16
        audio2 = (audio1(1,:)+audio1(2,:))/2;  % Stereo to mono
        audio3 = audio2/3276.7;
        
        audio_strct.data = audio3;
        
        %epoch the audio 
        EP_off = (audio_strct.event(find(strcmpi({audio_strct.event.type}, 'game_end'))).latency - audio_strct.event(find(strcmpi({audio_strct.event.type}, [task{k},'_game_start']))).latency)/audio_strct.srate;
        audio_strct = pop_epoch(audio_strct, {[task{k},'_game_start']}, [0  EP_off]);
        
       
        
        %save the audio structure data
        save([sbj{s},'_',task{k},'_audio.mat'],'audio3');
        save([sbj{s},'_',task{k},'_audio_strct.mat'],'audio_strct');
        
        
        
        for e = 1:size(EEG.event,2)
            % Find events with a 'odd_tone' in their name
            if ismember('odd_tone', EEG.event(e).type)
                % Sometimes there are two numbers in a trigger. One related to
                % the stimulus and the other is the time.
                num_in_trigger = regexp(EEG.event(e).type,'\d+\.?\d*','match');
                end_of_num = regexp(EEG.event(e).type,'\d+\.?\d*\w','end');
                num_events = length(EEG.event);
                for i = 1:length(num_in_trigger)
                    % The time is always before 's_post' in the trigger name
                    if strcmpi(EEG.event(e).type(end_of_num(i):end), 's_post')
                        odd_onset = str2double(num_in_trigger{i})*EEG.srate;
                        if str2double(num_in_trigger{i})>.5 || strcmpi(EEG.event(e).type, 'odd_tone_background_0s_post')
                            EEG.event(num_events+1).type = 'odd_tone';
                            EEG.event(num_events+1).latency = EEG.event(e).latency+odd_onset;
                            EEG.event(num_events+1).duration = 0.5;
                        elseif ismember('odd_tone_alarm_tone', EEG.event(e).type)   % Add specific name for odd tones and alarm sound that occured together
                            EEG.event(e).type = 'alarm_with_odd';
                            EEG.event(num_events+1).type = 'odd_in_alarm';
                            EEG.event(num_events+1).latency = EEG.event(e).latency+odd_onset;
                            EEG.event(num_events+1).duration = 0.5;
                        elseif ismember('odd_tone_irr_sound', EEG.event(e).type) % Add specific name for odd tones that occured in irrelevant sounds
                            EEG.event(num_events+1).type = 'odd_in_irr_sound';
                            EEG.event(num_events+1).latency = EEG.event(e).latency+odd_onset;
                            EEG.event(num_events+1).duration = 0.5;
                        end
                    end
                end
                % Irrelevant sounds are analysed together, so its easier if they
                % have the same name
            elseif strcmpi('irr_sound_1_post', EEG.event(e).type) || strcmpi('irr_sound_2_post', EEG.event(e).type)
                EEG.event(e).type = 'irr_sound';
            end
        end
        
        %separate the data set into eyes open, eyes closed and actual game
        %start
        %closed 
        if sum(ismember({EEG.event.type}, 'eyes_closed_end')) > 0 && sum(ismember({EEG.event.type}, 'eyes_closed_start')) > 0
            
            eyes_c = (EEG.event(find(strcmpi({EEG.event.type}, 'eyes_closed_end'))).latency - EEG.event(find(strcmpi({EEG.event.type}, 'eyes_closed_start'))).latency)/EEG.srate;
            EEG_eyesc = pop_epoch(EEG, {'eyes_closed_start'}, [0  eyes_c], 'epochinfo', 'yes');
            %save as new data sets
            %eyes closed
            filename = [task{k},'_eyes_closed'];
            EEG_eyesc.setname = [sbj{i},'_',filename];
            pop_saveset(EEG_eyesc, 'filename',filename);
            
        end
        
        if sum(ismember({EEG.event.type}, 'eyes_open_end')) > 0 && sum(ismember({EEG.event.type}, 'eyes_open_start')) > 0
            %open
            eyes_o = (EEG.event(find(strcmpi({EEG.event.type}, 'eyes_open_end'))).latency - EEG.event(find(strcmpi({EEG.event.type}, 'eyes_open_start'))).latency)/EEG.srate;
            EEG_eyeso = pop_epoch(EEG, {'eyes_open_start'}, [0  eyes_o], 'epochinfo', 'yes');
            %eyes open
            filename = [task{k},'_eyes_open'];
            EEG_eyeso.setname = [sbj{i},'_',filename];
            pop_saveset(EEG_eyeso, 'filename',filename);
            
        end
        
        if sum(ismember({EEG.event.type}, 'game_end')) > 0 && sum(ismember({EEG.event.type}, [task{k},'_game_start'])) > 0
            
            %full data set
            EP_game = (EEG.event(find(strcmpi({EEG.event.type}, 'game_end'))).latency - EEG.event(find(strcmpi({EEG.event.type}, [task{k},'_game_start']))).latency)/EEG.srate;
            EEG_game = pop_epoch(EEG, {[task{k},'_game_start']}, [0  EP_game], 'epochinfo', 'yes');
            
            %game
            filename = [task{k},'_game_added_trigger'];
            EEG_game.setname = [sbj{i},'_',filename];
            pop_saveset(EEG_game, 'filename',filename);
        end
        
        
        
        

        
        %full 
        filename = [task{k},'_eeg'];
        EEG.setname = [sbj{i},'_',filename];
        pop_saveset(EEG, 'filename',filename);
        
    end
    cd('..')
end