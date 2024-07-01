OT_setup
%CV values
nfold = 6;
testfold = 1;

%TRF parameters
Dir = 1; %specifies the forward modeling
tmin = -200;
tmax = 800;
lambdas = linspace(10e-4,10e4,10);
%expmrk 
auditory= {'alarm','odd'}
fig_path = [fig_path '\mininf\prop_exp\']
t = [-0.2 0.8];

base = [-0.2 -0.01]
erp_time = t(1):0.01:t(2) %0.1 is the sample rate here
base_idx = dsearchn(erp_time',base');
%all models
dat_epo = [];
for s=1:length(sbj)
    
    for k=1:2
        
        EEG = [];
        [EEG,PATH] = OT_preprocessing(s,k,sbj,20);
        
        reg = zeros(length(auditory),3,EEG.nbchan);
        for ao = 1:length(auditory)
            
            
            %extract the stimulus
            label = string(auditory(ao));
            
            %extract stimulus
            stim = extract_stimulus2(EEG, PATH, label, k,sbj{s},task);
            
            [STRAIN,RTRAIN] = mTRFpartition(stim,EEG.data',6);
            for i = 1:length(RTRAIN)
                model_train = mTRFtrain(STRAIN, RTRAIN,EEG.srate,1,tmin,tmax,0.05,'verbose',0);
                %baseline correction 
                temp_w = permute(squeeze(model_train.w),[2,1]);
                temp_wc = [];

                for ch = 1:EEG.nbchan
                    temp_wc(ch,:) = temp_w(ch,:) - mean(temp_w(ch,base_idx(1):base_idx(2)));
                end
                dat_epo(s,k,ao,:,:,i) = temp_wc;
            end
            %            [epo_dat,epo_stim,epo_idx,stim_z] = OT_epochize(EEG,stim,t,0);
            %            if k == 1
            %                cl = 'b'
            %            else
            %                cl= 'r'
            %            end
            %baseline correction
            %            for tr = 1:size(epo_dat,3)
            %                 for ch = 1:EEG.nbchan
            %                     %remove the base line
            %                     temp_epo(ch,:,tr) = squeeze(epo_dat(ch,:,tr)) - mean(epo_dat(ch,base_idx(1):base_idx(2),tr),2);
            %                 end
            %             end
            %
            %
            %            dat_epo(s,k,ao,:,:,:) = mean(temp_epo,3);
            %            dat_epo(s,k,ao,:,:) = squeeze(model_train.w);
            
            EEG2 = EEG;

            EEG2.pnts = length(erp_time);
            EEG2.times = model_train.t;
            EEG2.xmin = min(erp_time);
            EEG2.xmax = max(erp_time);
            EEG2.trials = 1;

            EEG2.data = squeeze(mean(dat_epo(s,k,ao,:,:,:),6));
            epo_dat(s,k,ao,:,:,:) = dat_epo(s,k,ao,:,:,:);    
            ft_dat{ao,k,s} = eeglab2fieldtrip(EEG2,'timelock');

        end
    end
end
cd(saving_path)         
chan_layout= eeglab2fieldtrip(EEG,'chanloc');

DNS_ft = struct()
DNS_ft.ft_dat = ft_dat
DNS_ft.auditory = auditory;
DNS_ft.epo_dat = epo_dat;
save('alarm_odd_ft_basecorr.mat','-struct','DNS_ft')
           
time= model_train.t;
dat_epo = epo_dat;
figure
plot(time,squeeze(mean(mean(dat_epo(:,1,1,:,:,:),6)))','b')
hold on
plot(time,squeeze(mean(mean(dat_epo(:,2,1,:,:,:),6)))','r')

figure
plot(time,squeeze(mean(mean(dat_epo(:,1,2,:,:,:),6)))','b')
hold on
plot(time,squeeze(mean(mean(dat_epo(:,2,2,:,:,:),6)))','r')





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



%%

cfg        = [];
cfg.layout = 'EEG1010.lay';
cfg.channel = chan_layout.elec.label
layout_m = ft_prepare_layout(cfg)
latency = [0 0.8]
eeg_chan = {EEG.chanlocs.labels}'
fig_path ='\\smb.uni-oldenburg.de\home\lorf0331\Documents\MATLAB\Project\OTtracking\Results\rebuttal\MARC_ERP\'
for ao = 1:length(auditory)
%     for k=1:length(task)
        ft_nw = []
        ft_wi = []
%         for s = 1:length(sbj)
%             ft_l{1,s} = ft_dat{k,ao,1,s};
%             ft_h{1,s} = ft_dat{k,ao,2,s};
%         end
        
        for s=1:length(sbj)
            
            
            temp_datl = ft_dat{ao,1,s};
            ft_tdatl = temp_datl.avg;
            
            
            temp_dath = ft_dat{ao,2,s};
            ft_tdath = temp_dath.avg;
            
            ft_nw{s,1} = temp_datl;
            
            ft_wi{s,1} = temp_dath;
        end

        
        %%
  
%         cfg = []
%         cfg.layout = layout_m
%         ft_layoutplot(cfg)

        cfg_neighb        = [];
        cfg_neighb.method = 'distance';
        cfg_neighb.neighbourdist = 80;
        neighbours = ft_prepare_neighbours(cfg_neighb, ft_nw{1,1});
        
        
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
        cfg.alpha            = 0.025;
        cfg.numrandomization = 10000;
        
        
        Nsubj  = length(sbj);
        
        design = zeros(2, Nsubj*2);
        design(1,:) = [1:Nsubj 1:Nsubj];
        design(2,:) = [ones(1,Nsubj) ones(1,Nsubj)*2];
        
        cfg.design = design;
        cfg.uvar   = 1;
        cfg.ivar   = 2;
        [stat] = ft_timelockstatistics(cfg, ft_nw{:}, ft_wi{:})
        
        %% get the data to plot ready
        cfg = [];
        cfg.channel   = 'all';
        cfg.latency   = latency;
        cfg.parameter = 'avg';
        GA_l        = ft_timelockgrandaverage(cfg, ft_nw{:});
        GA_h         = ft_timelockgrandaverage(cfg, ft_wi{:});
        
%         figure
        cfg = []
        cfg.layout = layout_m
        cfg.parameter = 'avg'
        cfg.comment = sprintf('%s %s',auditory{ao}, task{k})
        ft_multiplotER(cfg,GA_h,GA_l)
%         save_fig(gcf,fig_path,sprintf('DNS40_multi_%s_%s',auditory{ao}, task{k}))

        
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
                        save_fig(gcf,fig_path,sprintf('DNS40_posc_%s_%.3f_%.3f',auditory{ao},cfgPlot.xlim(1),cfgPlot.xlim(2)))
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
                        save_fig(gcf,fig_path,sprintf('DNS40_negc_%s_%.3f_%.3f',auditory{ao},cfgPlot.xlim(1),cfgPlot.xlim(2)))

                    end
                end
            end
            
            
            
        end
%     end
end


%% RR prep
%get the peak amplitude values
trf_t = tmin:10:tmax; %is in ms
%time windows 
time_w = [300 400;
    330 430];
chan{1,1} = {'Pz','P3','P4','CPz','CP1','CP2'};
chan{2,1}= {'Pz','P3','P4','CPz','CP1','CP2'};
res_dat = [];

for ao = 1:length(auditory)
   
    
    %find the window of interst
    w_idx = dsearchn(time',time_w(ao,:)')
    
    
    for ss = 1:20
        temp_nw(ss,ao,:,:,:) = zscore( squeeze(epo_dat(ss,1,ao,:,:,:)),0,'all');
        temp_wd(ss,ao,:,:,:) = zscore( squeeze(epo_dat(ss,2,ao,:,:,:)),0,'all');
    end
    
%     plot(squeeze(mean(mean(temp_nw(:,ao,:,:,:),5),1))');
%     hold on
%     plot(squeeze(mean(mean(temp_wd(:,ao,:,:,:),5),1))');
% %     get marcs data
%     get the channels
    chan_sub = chan{ao,1};
    chan_idx = ismember({EEG.chanlocs.labels},chan_sub)
    
    
    temp_res_nw(:,ao) = squeeze(mean(mean(mean(temp_nw(:,ao,chan_idx,w_idx(1):w_idx(2),:),3),4),5));
    temp_res_wd(:,ao) = squeeze(mean(mean(mean(temp_wd(:,ao,chan_idx,w_idx(1):w_idx(2),:),3),4),5));
    
    res_dat = [res_dat temp_res_nw(:,ao) temp_res_wd(:,ao)]
    
    
    nw_dat = squeeze(mean(mean(temp_nw(:,ao,chan_idx,:,:),3),5));
    wd_dat = squeeze(mean(mean(temp_wd(:,ao,chan_idx,:,:),3),5));
    
    nw_mean = mean(nw_dat,1)';
    nw_std  = (std(nw_dat,1)/sqrt(length(sbj)))'.*1.96;
    
    wd_mean = mean(wd_dat,1)';
    wd_std  = (std(wd_dat,1)/sqrt(length(sbj)))'.*1.96;
    
    figure
    set(gcf,'pos',fig_pos);
    h1=plot(trf_t,nw_mean,'r');
    hold on 
    h2=plot(trf_t,wd_mean,'b');
    
    % Plot shaded areas (error regions)
    fill([trf_t fliplr(trf_t)], [nw_mean + nw_std; flipud(nw_mean - nw_std)], 'r', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    fill([trf_t fliplr(trf_t)], [wd_mean + wd_std; flipud(wd_mean - wd_std)], 'b', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    
    % Highlight the specific region (300 to 400)
    highlight_x = [time_w(ao,:) flip(time_w(ao,:))];
    highlight_y = [-4 -4 4 4];
    fill(highlight_x, highlight_y, [0.5 0.5 0.5], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    
    title(auditory{ao},'FontSize',16)
        
    % Labels and title
    xlabel('Time');
    ylabel('Amplitude [uV]');
    set(gca,'FontSize',14)
    
    if ao == 1
       line([180 310], [-3.5 -3.5], 'Color', 'k', 'LineWidth', 2);
       hold on
       line([360 800], [-3.5 -3.5], 'Color', 'k', 'LineWidth', 2);
    else
        line([400 800], [-3.5 -3.5], 'Color', 'k', 'LineWidth', 2);
    end
    legend([h1,h2],{'narrow','wide'})


    
    save_fig(gcf,'\\smb.uni-oldenburg.de\home\lorf0331\Documents\MATLAB\Project\OTtracking\Results\rebuttal\MARC_ERP\',auditory{ao});


   
end

csvwrite('peakavg_basecorr.csv',res_dat)



           