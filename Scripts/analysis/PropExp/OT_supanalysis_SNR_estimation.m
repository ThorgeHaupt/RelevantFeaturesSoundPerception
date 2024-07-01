%SNR estimation for the different ERPS
OT_setup
auditory = {'onset','alarm','irregular','odd'};

for s = 1:length(sbj)
    
    for k = 1:2
        [EEG,PATH] = OT_preprocessing(s,k,sbj,20);
        fs = EEG.srate;
        t = [-0.1 0.5]
        for ao = 1:length(auditory)
            
            %extract the stimulus
            label = string(auditory(ao));
            stim = extract_stimulus2(EEG, PATH, label, k,sbj{s},task);
            
            %epoch the data
            [epo_dat,epo_stim,epo_idx,stim_z] = OT_epochize(EEG,stim,t,1);
            
            %clean ERP
            ERP_clean = squeeze(mean(mean(epo_dat,3),1));
            
            [val, val_idx] = min(ERP_clean);
            clean_dat = mean(abs(ERP_clean(1,10:end)));
            %Noise estimation
            epo_temp = epo_dat;
            epo_temp(:,:,1:2:end) = -1*epo_temp(:,:,1:2:end);
            ERP_noise = squeeze(mean(mean(epo_temp,3),1));
            noise_dat = mean(abs(ERP_noise(1,10:end)));
            std_estimate = mean(std(squeeze(mean(epo_dat(:,1:10,:),1)),[],2));
%             std_estimate = mean(abs(epo_dat(:,1:10,:)),'all');
            %SNR estimation
            erp_time = -0.1:1/EEG.srate:0.5;
            widx = dsearchn(erp_time',[0.08 0.12]');
            snr2(s,k,ao) = 10*log10(clean_dat /std_estimate);
            
            %     figure
            %     plot(erp_time,abs(ERP_clean),'k')
            %     hold on
            %     plot(erp_time,abs(ERP_noise),'r')
            %     title(auditory{ao});
        end
    end
end
    

dat = squeeze(mean(snr2,2));
figure
boxplot(dat,auditory);
ylabel('snr ratios')