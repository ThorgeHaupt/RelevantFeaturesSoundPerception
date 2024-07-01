%% save the stim vectors to safe processing time
%set up the path
MAINPATH = 'O:\projects\thh_ont\auditory-attention-in-complex-work-related-auditory-envrionments\data files'
addpath(genpath(MAINPATH));



sbj =    {'P001', 'P002','P003','P005', 'P006','P007','P008','P009',...
    'P010','P012','P013', 'P014','P015','P016','P017', 'P018',...
    'P019','P020','P021','P022'};

 for s=1:length(sbj)

    for k=1:2

        %% compute the envelopes
        [EEG,PATH] = OT_preprocessing(s,k,sbj,15);

        cd(PATH)
        if k == 1
            
            %% this is the legit version
            %do it using the MIR toolbox
            %resample the audio
            wavnarrow = load(sprintf('%s_narrow_audio_strct.mat',sbj{s}));

            fs_new = 16000;
            audio = resample(double(wavnarrow.audio_strct.data),fs_new,wavnarrow.audio_strct.srate);
            audio_resiz = audio(1:EEG.xmax*fs_new,1);
            
            %my own function
            [novelty, fs_new] = energy_novelty(double(wavnarrow.audio_strct.data)',wavnarrow.audio_strct.srate,'H',441);
            pks = simp_peak(novelty,0.18);
            if length(pks) > EEG.pnts
                pks = pks(1,1:EEG.pnts)';
            elseif length(pks) < EEG.pnts
                x(s,1) = EEG.pnts - length(pks);
                pks = [pks zeros(1,x(s,1))]';
            end

            nw_aeno = zeros(EEG.pnts,1);
            times = linspace(EEG.xmin, EEG.xmax, EEG.pnts);
            
            
            %derive the mel spectrogram 
            S = melSpectrogram(double(wavnarrow.audio_strct.data),wavnarrow.fs,'NumBands',64,'FrequencyRange',[1 8e3])
            
            %compute the app onset based enevelope
            for i=1:length(myOnsets)
                [~,idx11(i)] = min( abs(times - myOnsets(i)) );
                if idx11(i) < EEG.pnts(end)
                    nw_aeno(idx11(i),1) = 1;
                end
            end
            
            
            clear myOnsets
            
            %compute the trf envelope
            nw_entrf = mTRFenvelope(double(wavnarrow.audio_strct.data)',44100,100);
             if length(nw_entrf) > EEG.pnts
                nw_entrf = nw_entrf(1:EEG.pnts,1);
            elseif length(nw_entrf) < EEG.pnts
                x(s,2) = EEG.pnts - length(nw_entrf);
                nw_entrf = [nw_entrf' zeros(1,x(s,2))]';
            end
            
            %compute the standard envelope
            wav_h = abs( hilbert(double(wavnarrow.audio_strct.data)) );
            [b,a] = butter(3,30/(wavnarrow.audio_strct.srate/2),'low');
            wav_hf = filtfilt( b,a,wav_h );
            wav_hfd = downsample(wav_hf, 441)';
             if length(wav_hfd) > EEG.pnts
                nw_en = wav_hfd(1:EEG.pnts,1);
            elseif length(wav_hfd) < EEG.pnts
                x(s,3) = EEG.pnts - length(wav_hfd);
                nw_en = [wav_hfd' zeros(1,x(s,3))]';
            end
            %
            %compute onset based
            nw_eno = diff( nw_en );
            nw_eno( nw_eno < 0 ) = 0;
            nw_eno(end+1) = 0;
%             
            stims = [pks nw_en nw_eno nw_entrf]';
%             %for better comparison just flip it
%             stim_fa = flip(stims,2);
            save(sprintf('stims_narrow_%s',sbj{s}),'stims')

            
        else
            
            %load the audio files from each path
            wavwide = load('wide_audio.mat');

            times = linspace(EEG.xmin, EEG.xmax, EEG.pnts);
            
            wavwide = load(sprintf('%s_wide_audio_strct.mat',sbj{s}));

            fs_new = 16000;
            audio = resample(double(wavwide.audio_strct.data),fs_new,wavwide.audio_strct.srate);
%             audio_resiz = audio(1:EEG.xmax*fs_new,1);


            

            %my own function
            [novelty, fs_new] = energy_novelty(double(wavwide.audio_strct.data)',wavwide.audio_strct.srate,'H',441);
            pks = simp_peak(novelty,0.18);
            pks = pks(1,1:EEG.pnts)';
            
            %compute the trf envelope
            we_entrf = mTRFenvelope(double(wavwide.audio_strct.data),wavwide.audio_strct.srate,100);
            we_entrf = we_entrf(1:EEG.pnts,1);      %adjust the length of the envelope
            
            %compute the standard envelope
            wav_h = abs( hilbert(double(wavwide.audio_strct.data)) );
            [b,a] = butter(3,30/(wavnarrow.audio_strct.srate/2),'low');
            wav_hf = filtfilt( b,a,wav_h );
            wav_hfd = downsample(wav_hf, 441)';
            we_en = wav_hfd( 1:EEG.pnts,1 );
            
            %compute onset based
            we_eno = diff( we_en );
            we_eno( we_eno < 0 ) = 0;
            we_eno(end+1) = 0;
            stims = [pks we_en we_eno we_entrf]';
            save(sprintf('stims_wide_%s',sbj{s}),'stims')

        end
    end
    
end