%the trifecta model 
%test whether there is an effect of acoustic, sound identitiy, or mental
%state knowledge

OT_setup

Dir = 1; %specifies the forward modeling
tmin = -100;
tmax = 500;
lambdas = linspace(10e-4,10e4,10);

%partition the data set
nfold = 10;
testfold = 1;
base_audi = {'mTRF envelope'} %'onset','mel'; %needs to be done manually for all features
auditory = {'alarm','irregular','odd'};
audi_re = {'alarm','irrelevant','beep'};
audi_cmb{1,1} = {'alarmmenv','irrmenv','oddmenv'}
audi_cmb{1,2} = {'alarmons','irrons','oddons'}
audi_cmb{1,3} = {'alarmmel','irrmel','oddmel'}




for s=1:length(sbj)
    EEG_nw = [];
    [EEG1,PATH] = OT_preprocessing(s,1,sbj,20);
    
    
    EEG_wi = [];
    [EEG2,PATH] = OT_preprocessing(s,2,sbj,20);
    
    for ao = 1:length(auditory)
        EEG_nw = EEG1;
        EEG_wi = EEG2;
        EEG_nw.data = EEG_nw.data';
        EEG_wi.data = EEG_wi.data';

%         %extract the stimulus
%         label = string(auditory(ao));
%         stim_nw = extract_stimulus2(EEG_nw, PATH, label, 1,sbj{s},task);
%         
%         %extract the stimulus
%         label = string(auditory(ao));
%         stim_wi = extract_stimulus2(EEG_wi, PATH, label, 2,sbj{s},task);
%         
%         base_nw = extract_stimulus2(EEG_nw, PATH, base_audi{1,1}, 1,sbj{s},task);
%         base_wi = extract_stimulus2(EEG_wi, PATH, base_audi{1,1}, 2,sbj{s},task);


        %extract the stimulus
        label = string(auditory(ao));
        base_nw = extract_stimulus2(EEG_nw, PATH, label, 1,sbj{s},task);
        
        %extract the stimulus
        label = string(auditory(ao));
        base_wi = extract_stimulus2(EEG_wi, PATH, label, 2,sbj{s},task);
        
        stim_nw  = extract_stimulus2(EEG_nw, PATH, base_audi{1,1}, 1,sbj{s},task);
        stim_wi = extract_stimulus2(EEG_wi, PATH, base_audi{1,1}, 2,sbj{s},task);
        
        if ~isbinary(base_nw)
            %norm stim
            base_nw = (base_nw - min(base_nw,[],'all')) / ( max(base_nw,[],'all') - min(base_nw,[],'all') );
            base_wi = (base_wi - min(base_nw,[],'all')) / ( max(base_wi,[],'all') - min(base_wi,[],'all') );
        end
        
        
        if size(EEG_nw.data,1)>size(base_nw,1)
            EEG_nw.data = EEG_nw.data(1:size(base_nw,1),:);
            stim_nw = stim_nw(1:size(base_nw,1),:);
            
        end
        if size(EEG_wi.data,1)>size(base_wi,1)
            EEG_wi.data = EEG_wi.data(1:size(base_wi,1),:);
            stim_wi = stim_wi(1:size(base_wi,1),1);
        end
        EEG_nw.data = EEG_nw.data';
        EEG_wi.data = EEG_wi.data';
       
        EEG_whole = pop_mergeset(EEG_nw,EEG_wi,0);
        stim_whole = cat(1,stim_nw,stim_wi);
        stim_cond = [];
        base_whole = [];
        stim_cond(:,:) = [stim_nw' zeros(size(stim_wi))']';
        stim_cond(:,:) = [zeros(size(stim_nw))' stim_wi']';
        base_whole = cat(1,base_nw,base_wi);
        
        %start the nested model computations
        for rnd = 1:3
            if rnd == 1 %purely acoustic
                [spart,rpart] = mTRFpartition(base_whole,EEG_whole.data',nfold);
            elseif rnd == 2 %acoustic plus soundidentity
                stim_rns = [base_whole stim_whole];
                [spart,rpart] = mTRFpartition(stim_rns,EEG_whole.data',nfold);
            elseif rnd == 3 %acoustic plus condition informaiton
                stim_rns = [base_whole stim_cond];
                [spart,rpart] = mTRFpartition(stim_rns,EEG_whole.data',nfold);
            end
            model_weights = [];
            
            for tr = 1:length(spart)
                stest = spart(tr,:);
                rtest = rpart(tr,:);
                
                train_idx = setdiff(linspace(1,nfold,nfold),tr,'stable');
                
                strain = spart(train_idx,:);
                rtrain = rpart(train_idx,:);

                

                
                %% use cross-validation
                fs = EEG_nw.srate;
                
                cv = mTRFcrossval(strain,rtrain,fs,Dir,tmin,tmax,lambdas,'Verbose',0);
                
                %get the optimal regression parameter
                l = mean(cv.r,3); %over channels
                [l_val,l_idx] = max(mean(l,1));
                l_opt = lambdas(l_idx);
                
                
                %train the neural model on the optimal regularization parameter
                model_train = mTRFtrain(strain,rtrain,fs,Dir,tmin,tmax,l_opt,'verbose',0);
                
                %predict the neural data
                [PRED,STATS] = mTRFpredict(stest,rtest,model_train,'verbose',0);
                
                reg = STATS.r;
                %% compute the raw scores
                model_train2 = mTRFtrain(strainz,rtrainz,fs,Dir,tmin,tmax,0.05,'verbose',0);
                
                [PRED3,STATS3] = mTRFpredict(stestz,rtestz,model_train2,'verbose',0);
                raw(ao,:) = STATS3.r;
                result_reg(s,ao,rnd,tr,:) = reg;
                %save the weights
                model_weights(:,:,:,tr) = model_train.w;
            end
            
            mod_w{s,ao,rnd} = squeeze(mean(model_weights,4));
        end

    end
end
cd(saving_path);
level_analysis = struct();
level_analysis.result_reg = result_reg;
level_analysis.auditory = auditory
level_analysis.base_audi = base_audi;
level_analysis.mod_w = mod_w;
level_analysis.txt = 'see if adding more levels of information yields better model estimates, added looped training stuff'
save('level_analysis_envmaxrep.mat','-struct','level_analysis')





%% subject the models to permutation testing
erp_time = model_train.t;
ft_dat = cell(length(auditory),length(task),length(sbj)); %1= low; 2=high
load('level_analysis_envmaxrep.mat');
for s=1:length(sbj)
    
    
    for k=1:2
        
        EEG = [];
        [EEG,PATH] = OT_preprocessing(s,k,sbj,20);
        EEG.pnts = length(erp_time);
        EEG.times = erp_time;
        EEG.xmin = min(erp_time);
        EEG.xmax = max(erp_time);
        
        
         for ao=1:length(auditory)
             EEG.data = squeeze(mod_w{s,ao,3}(k+1,:,:));
             
             ft_dat{ao,k,s} = eeglab2fieldtrip(EEG,'timelock');
             chan_layout= eeglab2fieldtrip(EEG,'chanloc');

             

         end
    end
end

DNS_ft = struct()
DNS_ft.ft_dat = ft_dat
DNS_ft.auditory = auditory;
save('env_ft.mat','-struct','DNS_ft')
%% Grey scale
load('env_ft.mat')

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
[EEG,PATH] = OT_preprocessing(1,1,sbj,20);
chan_layout= eeglab2fieldtrip(EEG,'chanloc');

cfg        = [];
cfg.layout = 'EEG1010.lay';
cfg.channel = chan_layout.elec.label
layout_m = ft_prepare_layout(cfg)
latency = [0 0.5]
eeg_chan = {EEG.chanlocs.labels}'
base_idx = dsearchn((tmin:1000/EEG.srate:tmax)',[-100 -10]')
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
            temp_dat = temp_datl.avg';
            for i =1:size(temp_dat,1)
                temp_dat(i,:) = temp_dat(i,:) - mean(temp_dat(i,base_idx(1):base_idx(2)));
            end
             
            temp_datl.avg = temp_dat
            temp_datl.var = temp_datl.var';
            temp_datl.time = temp_datl.time./1000;
            

            
            temp_dath = ft_dat{ao,2,s};
            
            temp_dat = temp_dath.avg';
            for i =1:size(temp_dat,1)
                temp_dat(i,:) = temp_dat(i,:) - mean(temp_dat(i,base_idx(1):base_idx(2)));
            end
            temp_dath.avg = temp_dat;
            temp_dath.var = temp_dath.var';
            temp_dath.time = temp_dath.time./1000;

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
%                         save_fig(gcf,[fig_path '\OT14eventdns\ERP\'],sprintf('DNS40_posc_%s_%.3f_%.3f',auditory{ao},cfgPlot.xlim(1),cfgPlot.xlim(2)))
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
%                         save_fig(gcf,[fig_path '\OT14eventdns\ERP\'],sprintf('DNS40_negc_%s_%.3f_%.3f',auditory{ao},cfgPlot.xlim(1),cfgPlot.xlim(2)))

                    end
                end
            end
            
            
            
        end
%     end
end
        
        

%% get the model weights
%TRF time
trf_t = tmin:10:tmax; %is in ms
%time windows 
time_w = [300 400;
    90 140;
    330 430];
chan{1,1} = {'Pz','P3','P4','CPz','CP1','CP2'};
chan{2,1} = {'Fz','FC1','FC2','Cz','C3','C4'};
chan{3,1}= {'Pz','P3','P4','CPz','CP1','CP2'};

for ao = 1:3
   
    
    %find the window of interst
    w_idx = dsearchn(trf_t',time_w(ao,:)')
    w = squeeze(mod_w(:,ao,3))
    
    for ss = 1:20
        temp_nw = zscore( squeeze(w{ss,1}(2,:,:)),0,'all');
        temp_wd = zscore( squeeze(w{ss,1}(3,:,:)),0,'all');
        temp_w(ss,:,:,:) = permute(cat(3,temp_nw,temp_wd),[3 1 2])%squeeze(w{ss,1}(2:3,:,:));%squeeze(mean(temp_dat{s,rnd}(:,11));%
    end
    
    figure
    plot(squeeze(mean(mean(temp_w(:,1,:,:),4),1)));
    hold on
    plot(squeeze(mean(mean(temp_w(:,2,:,:),4),1)));
%     get marcs data
%     get the channels
    chan_sub = chan{ao,1};
    chan_idx = ismember({EEG.chanlocs.labels},chan_sub)
    
    
    temp_res(:,1,ao) = squeeze(mean(mean(temp_w(:,1,w_idx(1):w_idx(2),chan_idx),3),4));
    temp_res(:,2,ao) = squeeze(mean(mean(temp_w(:,2,w_idx(1):w_idx(2),chan_idx),3),4));
    
    
end


%get the results
for_marc = reshape(temp_res,20,[]);
csvwrite('peakavg.csv',for_marc)

