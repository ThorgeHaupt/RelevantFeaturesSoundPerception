%This script computes the ERPs for differnt soundscape contexts: Density

OT_setup


direc = 1; %specifies the forward modeling
lambda = 0.05;
t = [-0.5 1];
base = [-0.2 -0.01]
erp_time = t(1):0.01:t(2) %0.1 is the sample rate here
base_idx = dsearchn(erp_time',base');

auditory = {'onset','alarm','irregular','odd'};
na_idx = zeros(length(sbj),length(task),length(auditory));


hr_dat = cell(length(task),length(auditory));
lr_dat=  cell(length(task),length(auditory));

win_l = 20;
win_h = 5;

for s=1:length(sbj)
    
    
    for k=1:2
        
        EEG = [];
        [EEG,PATH] = OT_preprocessing(s,k,sbj,40);
        
%         EEG.data = detrend(EEG.data)
        
        %get the event density
        fs_new = 16000;

        %compute density
        label = string(auditory);
        stim = extract_stimulus2(EEG, PATH,'onset', k, sbj{s});
        win_lo = win_l*EEG.srate;
        win_ho = win_h*EEG.srate;

        % Compute the number of output samples
        num_output_samples = floor((length(stim) - win_lo)/win_ho) + 1;

        % Initialize the output signal
        ons_dist = zeros(num_output_samples,1);
        ons_dns = zeros(1,num_output_samples);
        

        % Apply the moving average filter to the input signal
        for i = 1:num_output_samples
            start_index = (i-1)*win_ho + 1;
            end_index = start_index + win_lo - 1;
            ons_dns(1,i) = mean(stim(start_index:end_index,1));
            ons_dist(i,1) = mean(diff(find(stim(start_index:end_index,1))));
        end
        audio_dns = ons_dns;
        %select the 5% highest and lowest peaks
        [audio_dns_sort, audio_sort_idx] = sort(audio_dns,'descend'); % <-- crucial change happens here
        audio_rms_sort = audio_dns_sort;
        sel_thrsh = 10;%ceil(length(audio_dns)*0.05);
        %select the persons rms values for the high and low peaks
%         rms_sbj(s,k,:) = [mean(audio_rms_sort(1,1:sel_thrsh)), mean(audio_rms_sort(1,end+1-sel_thrsh:end))];
        ons_dns = ons_dns(1,audio_sort_idx);

        %find the conversion to the real data ....
        time_win = [linspace(0,EEG.xmax-win_l,size(audio_dns,2))' linspace(win_l,EEG.xmax,size(audio_dns,2))'];
        time_eeg = linspace(EEG.xmin,EEG.xmax,EEG.pnts)';
        
        %adjust to the actual EEG data by finding the minimal offset
        time_avec = [];
        for i=1:size(time_win,1)
            [~, time_idx1] = min(abs(time_eeg - time_win(i,1)));
            [~, time_idx2] = min(abs(time_eeg - time_win(i,2)));
            
            time_avec(i,:) = [time_idx1 time_idx2];
        end
        time_avec_sort = time_avec(audio_sort_idx,:);
        
        
        %select the low peaks
        pks_ldns = [time_avec_sort(end-sel_thrsh+1:end,1) time_avec_sort(end-sel_thrsh+1:end,2)];
        pks_ldns = pks_ldns(pks_ldns(:,1)>=0,:);
        
        %select the high peaks
        pks_dns = [time_avec_sort(1:sel_thrsh,1) time_avec_sort(1:sel_thrsh,2)];
        pks_dns = pks_dns(pks_dns(:,1)>=0,:);
        
         for ao=1:length(auditory)
            
            
            %extract the stimulus
            label = string(auditory(ao));
            stim = extract_stimulus2(EEG, PATH, label, k, sbj{s},task);
            [epo_dat,epo_stim,~,stim_z] = OT_epochize(EEG,stim,t,0);
            stim_idx = find(stim_z==1);

            %find the epochs that fall into the range
            h_idx = [];
            l_idx = [];
            h_rank=[];
            l_rank = [];
            for id = 1:length(stim_idx)
                
                %check if the epoch is in the highest segments
                for sl = 1:sel_thrsh
                    if stim_idx(id) > time_avec_sort(sl,1) && stim_idx(id) < time_avec_sort(sl,2)
                        h_idx = [h_idx; find(stim_idx == stim_idx(id))]; 
                        h_rank = [h_rank;sl];
                    elseif stim_idx(id) > time_avec_sort(end-sl+1,1) && stim_idx(id) < time_avec_sort(end-sl+1,2)
                        l_idx = [l_idx; find(stim_idx == stim_idx(id))];
                        l_rank = [l_rank;sl];

                    end
                end
            end
            %delete double values 
            [hu_idx,ih] = unique(h_idx);
            [lu_idx,il] = unique(l_idx);
            
            %adjust rank
            hr_un = h_rank(ih);
            lr_un = l_rank(il);
            
            %sort the ranks
            [hr_sort,sh_idx] = sort(hr_un,'ascend');
            [lr_sort,sl_idx] = sort(lr_un,'ascend');
            
            %sort the epochs idx according to rank 
            h_sort = hu_idx(sh_idx);
            l_sort = lu_idx(sl_idx);

            %correct for different number of onsets -> throw them out
            if size(l_sort,1) > size(h_sort,1)
                l_dif = size(l_sort,1) - ((size(l_sort,1) - size(h_sort,1)));
                l_sort = l_sort(1:l_dif);
            elseif size(l_sort,1) < size(h_sort,1)
                h_dif = size(h_sort,1) - ((size(h_sort,1) - size(l_sort,1)));
                h_sort = h_sort(1:h_dif);
            end
            
            %should i baseline correct either of the two ERPs?
            h_epo = epo_dat(:,:,h_sort);
            l_epo = epo_dat(:,:,l_sort);

            for tr = 1:length(h_sort)
                for ch = 1:EEG.nbchan
                    %remove the base line
                    h_epo(ch,:,tr) = squeeze(h_epo(ch,:,tr)) - mean(h_epo(ch,base_idx(1):base_idx(2),tr),2);
                    l_epo(ch,:,tr) = squeeze(l_epo(ch,:,tr)) - mean(l_epo(ch,base_idx(1):base_idx(2),tr),2);
                end
            end
                

            high_epo = mean(h_epo,3);
            low_epo =  mean(l_epo,3);
            
            if any(isnan(high_epo),'all') || any(isnan(low_epo),'all')
                na_idx (s,k,ao) = 1;
                
            end
            
            h_dat(s,k,ao,:,:) = high_epo;
            l_dat(s,k,ao,:,:) =low_epo;
            
            %get the difference curve
            h_l_dif(s,k,ao,:,:) =high_epo - low_epo;
            
            nr_ons(s,k,ao,:) = length(h_sort);
            

            
         end
    end
end

%% Save your stuff
cd(saving_path)
DNS_epoch = struct()
DNS_epoch.h_dat = h_dat;
DNS_epoch.l_dat = l_dat;
DNS_epoch.auditory = auditory;
DNS_epoch.nr_ons = nr_ons;
DNS_epoch.h_l_dif = h_l_dif;
DNS_epoch.na_idx = na_idx;
DNS_epoch.erp_time = erp_time;
DNS_epoch.t = 'includes overlapping ERPs, sorted according to event density';
save('DNS20_ERP_518_onsdns.mat','-struct','DNS_epoch')

%% fieldtrip conversion 
ft_dat = cell(length(task),length(auditory),2,length(sbj)); %1= low; 2=high
for s=1:length(sbj)
    
    
    for k=1:2
        
        EEG = [];
        [EEG,PATH] = OT_preprocessing(s,k,sbj,40);
        EEG.pnts = length(erp_time);
        EEG.times = erp_time/1000;
        EEG.xmin = min(erp_time);
        EEG.xmax = max(erp_time);
        EEGl = EEG;
        EEGh = EEG;
        
         for ao=1:length(auditory)
             EEGl.data = squeeze(l_dat(s,k,ao,:,:))
             EEGh.data = squeeze(h_dat(s,k,ao,:,:))
             
             ft_dat{k,ao,1,s} = eeglab2fieldtrip(EEGl,'timelock');
             ft_dat{k,ao,2,s} = eeglab2fieldtrip(EEGh,'timelock');
             chan_layout= eeglab2fieldtrip(EEGh,'chanloc');

             

         end
    end
end

DNS_ft = struct()
DNS_ft.ft_dat = ft_dat
DNS_ft.auditory = auditory;
save('DNS20_TRF_ft_518_onsdns.mat','-struct','DNS_ft')


%% plotting
erp_time = linspace(-500,1000,size(h_dat,5));
base_win = [-400 -200];
base_idx = dsearchn(erp_time',base_win');

%window of interest
win_int = [-100 500];
win_idx = dsearchn(erp_time',win_int');

% c_idx = 8; %this is Cz
% EEGl = EEG;
% EEGh = EEG;

t_erp_time = 0:10:500;
%% plotting the average ERPS and standard deviations
le = length(erp_time);
fig_pos = [57   418   749   560]
%figure,clf

for ao = 1:length(auditory)
    figure
    temp_dat = squeeze(mean(mean(h_dat(:,:,ao,:,:),2,'omitnan'),4,'omitnan'));
    dnsMean = squeeze(mean(temp_dat,1,'omitnan'));
    h = plot(dnsMean,'b','linew',2);
    hold on
    N = size(h_dat,1);
    ySEM = std(temp_dat,1)/sqrt(N);
    CI95 = tinv([0.025 0.975],N-1);
    yCI95 = bsxfun(@times,ySEM,CI95(:));
    conv = yCI95 + dnsMean ;
    x2 = [linspace(1,le,le) fliplr(linspace(1,le,le))];
    inbe = [conv(1,:) fliplr(conv(2,:))];
    f = fill(x2,inbe,'b');
    f.FaceAlpha = 0.2;
    f.EdgeAlpha = 0.4;
    f.LineWidth = 0.5;
    hold on
    
    temp_datl = squeeze(mean(mean(l_dat(:,:,ao,:,:),2,'omitnan'),4,'omitnan'));
    dnslMean = squeeze(mean(temp_datl,1,'omitnan'));
    l = plot(dnslMean,'r','linew',3);
    hold on
    ylSEM = std(temp_datl,1)/sqrt(N);
    CI95 = tinv([0.025 0.975],N-1);
    yCI95 = bsxfun(@times,ylSEM,CI95(:));
    conv = yCI95 + dnslMean;
    x2 = [linspace(1,le,le) fliplr(linspace(1,le,le))];
    inbe = [conv(1,:) fliplr(conv(2,:))];
    f = fill(x2,inbe,'r');
    f.FaceAlpha = 0.1;
    f.EdgeAlpha = 0.2;
    f.LineWidth = 0.5;

    l = legend([h,l],'high dns','low dns','Fontsize',24,'Location','southeast')
    set(l,'Box','off')

    set(gca,'XTick', linspace(1,le,16),'XTickLabel',linspace(-500,1000,16),'Xlim',[40 102],'Fontsize',24)
    title(sprintf('ERP %s',auditory{ao}),'FontSize',30)
    
    xlabel('Time (ms)','Fontsize',24)
    ylabel('[\muV]', 'Interpreter', 'tex','Fontsize',24);
    
    set(gcf,'Position',fig_pos)
    box off
    
%     save_fig(gcf,['\\smb.uni-oldenburg.de\home\lorf0331\Documents\MATLAB\Project\OTtracking\Results\rebuttal\DNS\'],sprintf('%s_DNS_onsdns_ERP',auditory{ao}))
    
end



%% Permutation testing
Map_aMINUSi = 1 - gray; % Invert the grayscale colormap and store it in Map_aMINUSi
Temp = Map_aMINUSi;    % Create a temporary variable to store the modified colormap

Map_aMINUSi(1:60, :) = repmat([1 1 1], [60, 1]);
% Set the first 60 rows to pure white (RGB value [1 1 1])

% The rest of the code remains the same, setting specific rows to specific color levels
Map_aMINUSi(61:65, :) = repmat(Temp(3, :), [5, 1]);
Map_aMINUSi(66:70, :) = repmat(Temp(6, :), [5, 1]);
Map_aMINUSi(71:75, :) = repmat(Temp(10, :), [5, 1]);
Map_aMINUSi(76:80, :) = repmat(Temp(20, :), [5, 1]);
Map_aMINUSi(81:85, :) = repmat(Temp(35, :), [5, 1]);
Map_aMINUSi(86:90, :) = repmat(Temp(45, :), [5, 1]);
Map_aMINUSi(91:94, :) = repmat(Temp(55, :), [4, 1]);

%% get chan layout
s = 1;
k=1;
OT_setup


EEG = [];
[EEG,PATH] = OT_preprocessing(s,k,sbj,40);

chan_layout= eeglab2fieldtrip(EEG,'chanloc');


%%

cfg        = [];
cfg.layout = 'EEG1010.lay';
cfg.channel = chan_layout.elec.label
layout_m = ft_prepare_layout(cfg)
latency = [0 0.3]
eeg_chan = {EEG.chanlocs.labels}'

for ao = 1:length(auditory)
%     for k=1:length(task)
        ft_l = []
        ft_h = []
%         for s = 1:length(sbj)
%             ft_l{1,s} = ft_dat{k,ao,1,s};
%             ft_h{1,s} = ft_dat{k,ao,2,s};
%         end
        
        for s=1:length(sbj)
            
            for k = 1:2
                temp_datl = ft_dat{k,ao,1,s};
                ft_tdatl(k,:,:) = temp_datl.avg;
                
                
                temp_dath = ft_dat{k,ao,2,s};
                ft_tdath(k,:,:) = temp_dath.avg;
            end
            temp_datl.avg = squeeze(mean(ft_tdatl,1,'omitnan'));
            
            
            temp_dath.avg = squeeze(mean(ft_tdath,1,'omitnan'));
           
%             %% baseline correction 
%             temp_datl_dat = temp_datl.avg;
%             temp_dath_dat = temp_dath.avg;
%             base = [-0.2 -0.01];
%             
%             
%             base_idx = dsearchn(temp_dath.time',base');
% 
%             %baseline correct
%             for ch = 1:EEG.nbchan
%                 %remove the base line
%                 temp_datl_dat(ch,:) = temp_datl_dat(ch,:) - mean(temp_datl_dat(ch,base_idx(1):base_idx(2)),2);
%                 temp_dath_dat(ch,:) = temp_dath_dat(ch,:) - mean(temp_dath_dat(ch,base_idx(1):base_idx(2)),2);
%             end
%              temp_datl.avg = temp_datl_dat;
%              temp_dath.avg = temp_dath_dat;
%             
            
            %%
            
            ft_l{s,1} = temp_datl;
            
            
            %baseline correct
            ft_h{s,1} = temp_dath;
        end
        if ao == 3 
            ft_l(10,:) = [];
            ft_h(10,:) = [];
        end
        
        
        %%
  
%         cfg = []
%         cfg.layout = layout_m
%         ft_layoutplot(cfg)

        cfg_neighb        = [];
        cfg_neighb.method = 'distance';
        cfg_neighb.neighbourdist = 80;
        neighbours = ft_prepare_neighbours(cfg_neighb, ft_l{1,1});
        
        
        cfg         = [];
        cfg.channel = {'EEG'};
        cfg.latency = latency;
        cfg.method           = 'montecarlo';
        cfg.statistic        = 'depsamplesT';
        cfg.correctm         = 'cluster';
        cfg.clusteralpha     = 0.05;
        cfg.clusterstatistic = 'maxsum';
        cfg.minnbchan        = 2;
        cfg.neighbours       = neighbours;  % same as defined for the between-trials experiment
        cfg.tail             = 0;
        cfg.clustertail      = 0;
        cfg.alpha            = 0.05;
        cfg.numrandomization = 10000;
        
        if ao == 3
            Nsubj  = length(sbj)-1;
        else
            Nsubj  = length(sbj);
        end
        design = zeros(2, Nsubj*2);
        design(1,:) = [1:Nsubj 1:Nsubj];
        design(2,:) = [ones(1,Nsubj) ones(1,Nsubj)*2];
        
        cfg.design = design;
        cfg.uvar   = 1;
        cfg.ivar   = 2;
        [stat] = ft_timelockstatistics(cfg, ft_l{:}, ft_h{:});
        
        %% get the data to plot ready
        cfg = [];
        cfg.channel   = 'all';
        cfg.latency   = latency;
        cfg.parameter = 'avg';
        GA_l        = ft_timelockgrandaverage(cfg, ft_l{:});
        GA_h         = ft_timelockgrandaverage(cfg, ft_h{:});
        
%         figure
        cfg = []
        cfg.layout = layout_m
        cfg.parameter = 'avg'
        cfg.comment = sprintf('%s %s',auditory{ao}, task{k})
        ft_multiplotER(cfg,GA_h,GA_l)
        %save_fig(gcf,fig_path,sprintf('DNS40_multi_%s_%s',auditory{ao}, task{k}))

        
        cfg = [];
        cfg.operation = 'subtract';
        cfg.parameter = 'avg';
        GA_lvsh    = ft_math(cfg, GA_l, GA_h);
        
        %% bojana code for plotting
        
        comps={'aMINUSi'};
        Cluster=[];
        for idxComps = 1:size(comps,2)
            CustMap=eval(['Map_' comps{idxComps}]);
            
            if isfield(stat,'posclusters') && isempty(stat.posclusters) ==0
                pos_cluster_pvals = [stat.posclusters(:).prob];
                pos_signif_clust = find(pos_cluster_pvals < stat.cfg.clusteralpha);%stat.cfg.clusteralpha
                
                if ~isempty(pos_signif_clust)
                    for idxPos = 1:length(pos_signif_clust)
                        pos = stat.posclusterslabelmat == pos_signif_clust(idxPos);
                        sigposCLM = (pos == 1);
                        probpos(idxPos) = stat.posclusters(idxPos).prob;
                        possum_perclus = sum(sigposCLM,1); %sum over chans for each time- or freq-point
                        
                        Cluster{idxComps,idxPos}=sum(pos,2)>0;
                        
                        ind_min = min(find(possum_perclus~=0));
                        ind_max = max(find(possum_perclus~=0));
                        time_perclus = [stat.time(ind_min) stat.time(ind_max)];%
                        ClusterTime{idxComps,idxPos}=time_perclus;
                        
                         % provides the OG chanlocs indicies
                        pos_int = find(sum(pos(:,ind_min:ind_max),2)> 1);
                        pos_chan = eeg_chan(pos_int);

                        
                        figure('Name',['Probability:' num2str(probpos(idxPos))],'NumberTitle','off')
                        cfgPlot=[];
                        cfgPlot.xlim=[time_perclus(1) time_perclus(2)];
                        cfgPlot.highlight = 'labels';
                        cfgPlot.highlightchannel = pos_chan;
                        cfgPlot.highlightsymbol    = 'x'
                        cfgPlot.highlightcolor     = [220 20 60]./255
                        cfgPlot.highlightsize      = 8
                        cfgPlot.parameter='avg';
                        cfgPlot.layout= layout_m;
%                         cfgPlot.colormap = CustMap;
%                         cfgPlot.zlim = [0 4];
                        cfgPlot.colorbar ='EastOutside';
                        cfgPlot.comment = sprintf('%s [%.3f %.3f]s.',auditory{ao},cfgPlot.xlim(1),cfgPlot.xlim(2))
                        cfgPlot.commentpos = 'title';
                        cfgPlot.colorbar ='yes';
                        cfgPlot.style = 'straight';
                        stat.stat = abs(stat.stat);
                        ft_topoplotER(cfgPlot, GA_lvsh);
%                         save_fig(gcf,[fig_path '\OT14eventdns\TRF\'],sprintf('TRF_posc_%s_%.3f_%.3f',auditory{ao},cfgPlot.xlim(1),cfgPlot.xlim(2)))
%                         save_fig(gcf,[fig_path '\OT14eventdns\ERP\'],sprintf('DNS40epo_ovlp_posc_%s_%.3f_%.3f',auditory{ao},cfgPlot.xlim(1),cfgPlot.xlim(2)))
%                         save_fig(gcf,fig_path,sprintf('TRF_%s_pos',auditory{ao}))
                    end
                end
            end
            
            if isfield(stat,'negclusters') && isempty(stat.negclusters) ==0
                neg_cluster_pvals = [stat.negclusters(:).prob];
                neg_signif_clust = find(neg_cluster_pvals < stat.cfg.clusteralpha);%stat.cfg.clusteralpha
                
                if isempty(neg_signif_clust) == 0
                    for idxNeg = 1:length(neg_signif_clust)
                        neg = stat.negclusterslabelmat == neg_signif_clust(idxNeg);
                        signegCLM = (neg == 1);
                        probneg(idxNeg) = stat.negclusters(idxNeg).prob;
                        negsum_perclus = sum(signegCLM,1); %sum over chans for each time- or freq-point
                        
                        Cluster{idxComps,length(neg_signif_clust)+idxNeg}=sum(neg,2)>0;
                        
                        ind_min = min(find(negsum_perclus~=0));
                        ind_max = max(find(negsum_perclus~=0));
                        time_perclus = [stat.time(ind_min) stat.time(ind_max)];%
                        ClusterTime{idxComps,length(neg_signif_clust)+idxNeg}=time_perclus;
                        
                        % provides the OG chanlocs indicies
                        neg_int = find(sum(neg(:,ind_min:ind_max),2)> 1);
                        neg_chan = eeg_chan(neg_int);
                        
                        figure('Name',['Probability:' num2str(probneg(idxNeg))],'NumberTitle','off')
                        cfgPlot=[];
                        cfgPlot.xlim=[time_perclus(1) time_perclus(2)];
                        cfgPlot.highlight = 'labels';
                        cfgPlot.highlightchannel = neg_chan;
                        cfgPlot.highlightsymbol    = 'x'
                        cfgPlot.highlightcolor     = [220 20 60]./255
                        cfgPlot.highlightsize      = 8
                        cfgPlot.parameter='avg';
                        cfgPlot.layout=layout_m;
%                         cfgPlot.colormap = CustMap;
%                         cfgPlot.zlim = [0 4];
                        cfgPlot.comment = sprintf('%s [%.3f %.3f]s.',auditory{ao},cfgPlot.xlim(1),cfgPlot.xlim(2))
                        cfgPlot.colorbar ='EastOutside';
                        cfgPlot.commentpos = 'title';
                        cfgPlot.colorbar ='yes';
                        cfgPlot.style = 'straight';
                        stat.stat = abs(stat.stat);                        
                        ft_topoplotER(cfgPlot, GA_lvsh);
%                         save_fig(gcf,[fig_path '\OT14eventdns\TRF\'],sprintf('TRF_negc_%s_%.3f_%.3f',auditory{ao},cfgPlot.xlim(1),cfgPlot.xlim(2)))

%                         save_fig(gcf,[fig_path '\OT14eventdns\ERP\'],sprintf('DNS40epo_ovlp_negc_%s_%.3f_%.3f',auditory{ao},cfgPlot.xlim(1),cfgPlot.xlim(2)))
%                         save_fig(gcf,fig_path,sprintf('TRF_%s_neg',auditory{ao}))

                    end
                end
            end
            
            
            
        end
%     end
end
           





            
