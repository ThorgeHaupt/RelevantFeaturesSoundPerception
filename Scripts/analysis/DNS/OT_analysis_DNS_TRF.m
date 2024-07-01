% script to assess the overall event density function 
% data driven i.e. periods with greater event density should differ with
% respect to the neural response compared to periods with less events
OT_setup


direc = 1; %specifies the forward modeling
tmin = -100;
tmax = 500;
lambda = 0.05;
N = 100;
perm = 100;
auditory = {'mTRF envelope','onset','alarm','irregular','odd'}%,'mel'};
% auditory = {'onset'}


win_l = 20;
win_h = 5;
for s=1:length(sbj)
    
    
    for k=1:2
        
        EEG = [];
        [EEG,PATH] = OT_preprocessing(s,k,sbj,20);
        
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
        
        
        
        
        [audio_dns_sort, audio_sort_idx] = sort(audio_dns,'descend'); % <-- crucial change happens here
        audio_rms_sort = audio_dns_sort;
        
        
        
        
        sel_thrsh = 10;%ceil(length(audio_dns)*0.08);
        %select the persons rms values for the high and low peaks
        rms_sbj(s,k,:) = [mean(audio_rms_sort(1,1:sel_thrsh)), mean(audio_rms_sort(1,end+1-sel_thrsh:end))];
        
        
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
            
            %set up the saving container
            resp = cell(size(pks_dns,1),1);
            stims = cell(size(pks_dns,1),1);
            respl = cell(size(pks_dns,1),1);
            stimsl = cell(size(pks_dns,1),1);
            
            %% low pks
            %should i integrate a function that select the throughs instead
            %of just a random subsample
            
            for i = 1:size(pks_ldns,1)
                if pks_ldns(i,2) <= EEG.pnts
                    respl{i,1} = zscore(EEG.data(:,pks_ldns(i,1):pks_ldns(i,2))',[],'all');
                    
                    if ~isbinary(stim)
                        stimsl{i,1} = normalize(stim(pks_ldns(i,1):pks_ldns(i,2),:),'range');%/max(stim(pks_ldns(i,1):pks_ldns(i,2),:));      %important as TRF is not invariant to scaling
                    else
                        stimsl{i,1} = stim(pks_ldns(i,1):pks_ldns(i,2),1);
                        nr_onsl(s,k,ao,i) = sum(stimsl{i,1});
                        
                        
                    end
                else
                end
                
            end
            %remove the empty cells
            null_l(s,k,ao,:) = size(stimsl,1) - sum(cellfun(@sum,stimsl)>0);
            respl = respl(cellfun(@sum,stimsl)>0,:);
            stimsl = stimsl(cellfun(@sum,stimsl)>0,1);
            
            tot_sum = sum(cellfun(@sum,stimsl))
            
            %% high peaks
            
            rep_trf = zeros(EEG.nbchan,length(tmin:10:tmax),N);
            for rep = 1:N
                for i = 1:size(pks_dns,1)
                    if pks_dns(i,2) <= EEG.pnts
                        resp{i,1} = zscore(EEG.data(:,pks_dns(i,1):pks_dns(i,2))',[],'all');
                        if ~isbinary(stim)
                            stims{i,1} = normalize(stim(pks_dns(i,1):pks_dns(i,2),:),'range');%/max(stim(pks_dns(i,1):pks_dns(i,2),1));     %important as TRF is not invariant to scaling
                            
                        elseif contains(auditory{ao},'onset')
                            stims{i,1} = stim(pks_dns(i,1):pks_dns(i,2),1);
                            nr_onsh(s,k,ao,i) = sum(stims{i,1});
                            %adjust the number of trials to that of the low
                            %onset stuff
                            %random segment selection
                            ix = randi(10);
                            if nr_onsh(s,k,ao,i) > nr_onsl(s,k,ao,ix)
                                
                                %set up the new stims container
                                stims_ons = zeros(size(stims{i,1}));
                                
                                %get the differences between high and low
                                %event density events
                                dif_ons = nr_onsh(s,k,ao,i) - nr_onsl(s,k,ao,ix);
                                
                                %select random ones to delete
                                del_ons = randperm(nr_onsh(s,k,ao,i),dif_ons);
                                
                                %find the onset indices
                                ons_idx = find(stims{i,1} ==1);
                                
                                %select the ones to delete
                                ons_del_idx = ons_idx(del_ons);
                                
                                %get the onsets that should not be deleted
                                stims_ons(setdiff(ons_idx,ons_del_idx),1) = 1;
                                stims{i,1} = stims_ons;
                            else
                            end
                        else
                            stims{i,1} = stim(pks_dns(i,1):pks_dns(i,2),1);
                            
                        end
                    else
                    end
                    
                end
                
                model_dns = mTRFtrain(stims,resp,EEG.srate,1,tmin,tmax,lambda,'verbose',0);
                rep_trf(:,:,rep) =  squeeze(model_dns.w)';
                if ~contains(auditory{ao},'onset')
                    disp('break')
                    break;
                end
            end
            %remove the empty cells, fail switch in case something is all
            %zeros
            if contains(auditory{ao},'onset')
                trf_hdns(s,k,ao,:,:) = squeeze(mean(rep_trf,3));
            else
                if sum(cellfun(@sum,stims))<=0 || sum(cellfun(@sum,stimsl))<=0
                    trf_hdns(s,k,ao,:) = NaN;
                    trf_ldns(s,k,ao,:) = NaN;
                    dif_trf(s,k,ao,:) = NaN;
                    break;
                else
                    null_h(s,k,ao,:) = size(stims,1) - sum(cellfun(@sum,stims)>0);
                    resp = resp(cellfun(@sum,stims)>0,:);
                    stims = stims(cellfun(@sum,stims)>0,:);
                    
                    model_dns = mTRFtrain(stims,resp,EEG.srate,1,tmin,tmax,lambda,'verbose',0);
                    trf_hdns(s,k,ao,:,:) = squeeze(model_dns.w)';
                end
            end
            nr_onsh(s,k,ao,:) = sum(cellfun(@sum,stims));
            
            
            %% train trf
            model_ldns = mTRFtrain(stimsl,respl,EEG.srate,1,tmin,tmax,lambda,'verbose',0);
            trf_ldns(s,k,ao,:,:) = squeeze(model_ldns.w)';
            
            dif_trf(s,k,ao,:,:) =  squeeze(model_ldns.w)'-squeeze(model_dns.w)';
            
            
        end
    end
   
end

%save as structure 
cd(saving_path)
trf_dns= struct;
trf_dns.trf_hdns = trf_hdns
trf_dns.trf_ldns = trf_ldns
trf_dns.dif_trf = dif_trf
%trf_hdns.dns = dns_sbj;
trf_dns.nr_onsl = nr_onsl;
trf_dns.nr_onsh = nr_onsh;
trf_dns.auditory = auditory;
trf_dns.tmin = tmin;
trf_dns.tmax = tmax;
trf_dns.txt = '10 most extreme segments were used, TRF weights'
save('trf_dns_518_20.mat','-struct','trf_dns')

%% Make it fieldtrip useable
erp_time = linspace(tmin/1000,tmax/1000,size(trf_hdns,5));
ft_dat = cell(length(task),length(auditory),2,length(sbj)); %1= low; 2=high
for s=1:length(sbj)
    
    
    for k=1:2
        
        EEG = [];
        [EEG,PATH] = OT_preprocessing(s,k,sbj,40);
        EEG.pnts = length(erp_time);
        EEG.times = erp_time;
        EEG.xmin = min(erp_time);
        EEG.xmax = max(erp_time);
        EEGl = EEG;
        EEGh = EEG;
        
         for ao=1:length(auditory)
             EEGl.data = squeeze(trf_ldns(s,k,ao,:,:))
             EEGh.data = squeeze(trf_hdns(s,k,ao,:,:))
             
             ft_dat{k,ao,1,s} = eeglab2fieldtrip(EEGl,'timelock');
             ft_dat{k,ao,2,s} = eeglab2fieldtrip(EEGh,'timelock');
             chan_layout= eeglab2fieldtrip(EEGh,'chanloc');

             

         end
    end
end
trf_dns_ft = struct()
trf_dns_ft.ft_dat = ft_dat
trf_dns_ft.auditory = auditory;
save('DNS_ft_20_518.mat','-struct','trf_dns_ft')


%% plot the average high and low dns TRFs
le = size(trf_hdns(:,1,1,:,:),5);
erp_time = linspace(tmin/1000,tmax/1000,size(trf_hdns,5));
fig_pos = [57   418   749   560]
for ao = 1:length(auditory)
    %% baseline correct 
    base = [-0.2 -0.01];
    base_idx = dsearchn(erp_time',base');
    temp_dath = squeeze(mean(trf_hdns(:,:,ao,:,:),2,'omitnan'));
    temp_datl = squeeze(mean(trf_ldns(:,:,ao,:,:),2,'omitnan'));
    
    %% baseline correct
    for s = 1:length(sbj)
        for ch = 1:EEG.nbchan
            %remove the base line
            temp_datl(s,ch,:) = temp_datl(s,ch,:) - squeeze(mean(temp_datl(s,ch,base_idx(1):base_idx(2)),3));
            temp_dath(s,ch,:) = temp_dath(s,ch,:) - squeeze(mean(temp_dath(s,ch,base_idx(1):base_idx(2)),3));
        end
    end
%% 
    temp_datl = squeeze(mean(temp_datl,2));
    temp_dath = squeeze(mean(temp_dath,2));
    figure,clf
   
    dnsMean = squeeze(mean(temp_dath,1,'omitnan'));
    h = plot(dnsMean,'b','linew',2);
    hold on
    N = size(trf_hdns,1);
    ySEM = std(temp_dath,1)/sqrt(N);
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
    
    dnslMean = squeeze(mean(temp_datl,1,'omitnan'));
    l = plot(dnslMean,'r','linew',2);
    hold on
    ylSEM = std(temp_datl,1)/sqrt(N);
    CI95 = tinv([0.025 0.975],N-1);
    yCI95 = bsxfun(@times,ylSEM,CI95(:));
    conv = yCI95 + dnslMean;
    x2 = [linspace(1,le,le) fliplr(linspace(1,le,le))];
    inbe = [conv(1,:) fliplr(conv(2,:))];
    f = fill(x2,inbe,'r');
    f.FaceAlpha = 0.2;
    f.EdgeAlpha = 0.4;
    f.LineWidth = 0.5;
    if ~strcmp(auditory{ao},'mTRF envelope')
%         leg = legend([h,l],sprintf('high avg. %0.2f',mean(squeeze(nr_onsh(:,:,ao,:)),'all')) ,sprintf('low avg. %0.2f',mean(squeeze(nr_onsl(:,:,ao,:)),'all')),'Fontsize',24,'Location','southeast')
    else
        leg = legend([h,l],'high dns','low dns','Fontsize',24)
        set(leg,'Box','off')
    end
    set(gca,'XTick', linspace(1,le,9),'XTickLabel',linspace(tmin,tmax,9),'Xlim',[9 le-10],'Fontsize',24)%,'ylim',[-59 35])
    
    title(sprintf('%s',auditory{ao}),'FontSize',30)
    box off
    xlabel('Time (ms)','Fontsize',24)
    ylabel('a.u.','Fontsize',24)
    
    set(gcf,'Position',fig_pos)
    
%     save_fig(gcf,fig_path,sprintf('rebutta_ownons_%s',auditory{ao}))
    
end

%% Cluster based permutation testing
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
           
            %% baseline correction 
            temp_datl_dat = temp_datl.avg;
            temp_dath_dat = temp_dath.avg;
            base = [-0.2 -0.01];
            
            
            base_idx = dsearchn(temp_dath.time',base');

            %baseline correct
            for ch = 1:EEG.nbchan
                %remove the base line
                temp_datl_dat(ch,:) = temp_datl_dat(ch,:) - mean(temp_datl_dat(ch,base_idx(1):base_idx(2)),2);
                temp_dath_dat(ch,:) = temp_dath_dat(ch,:) - mean(temp_dath_dat(ch,base_idx(1):base_idx(2)),2);
            end
             temp_datl.avg = temp_datl_dat;
             temp_dath.avg = temp_dath_dat;
            
            
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