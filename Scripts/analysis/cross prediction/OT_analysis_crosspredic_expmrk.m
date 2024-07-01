direc = 1; %specifies the forward modeling
tmin = -100;
tmax = 500;
lambdas = linspace(10e-4,10e4,10);
result= [];
fig_pos = [240.200000000000,313,1484.40000000000,692.800000000000];
OT_setup

    
auditory = {'oddirr','oddalarm','irralarm'};
n_fold = 6;
acc_ons_x = zeros(length(sbj),2,length(auditory),n_fold);
acc_x_ons = zeros(length(sbj),2,length(auditory),n_fold);
acc_ons_ons = zeros(length(sbj),2,length(auditory),n_fold);
acc_x_x = zeros(length(sbj),2,length(auditory),n_fold);

for s=1:length(sbj)
   
   
    for k=1:2
        EEG = [];
        [EEG,PATH] = OT_preprocessing(s,k,sbj,20);
       
        fs = EEG.srate; 
        
        for ao = 1:length(auditory)%size(auditory,1)%
            
            resp = EEG.data';
            
            if size(auditory,1)>1
                
                %get the two stimulus information
                label = string(auditory{ao,1});
                stim_x = extract_stimulus2(EEG, PATH, label, k,sbj{s},task);
                
                label = string(auditory{ao,2});
                stim_ons = extract_stimulus2(EEG, PATH, label, k,sbj{s},task);
                
            else
                label = string(auditory{ao});
                stim = extract_stimulus2(EEG, PATH, label, k,sbj{s},task);
                
                %split the stimulus up
                stim_x = stim(:,1);
                stim_ons = stim(:,2);
 
            end
      
            %% Partition the two stimulus data sets
            format longG
            %specifiy training and testing data sets
            [spart_x,rpart] = mTRFpartition(stim_x, resp,n_fold);
            [spart_ons,rpart] = mTRFpartition(stim_ons, resp,n_fold);

            %z-score the neural data
            rpart_z = cellfun(@(x) zscore(x,[],'all'),rpart,'UniformOutput',false);
            train_trials = linspace(1,size(rpart,1),size(rpart,1));
            for tr=1:size(rpart,1)
                
                %select the testing segments
                x_ts = spart_x(tr,1);
                ons_ts = spart_ons(tr,1);
                rep_ts = rpart_z(tr,1);
                
                %select the training segments
                train_idx = setdiff(train_trials,tr,'stable');
                
                x_tr = spart_x(train_idx,1);
                ons_tr = spart_ons(train_idx,1);
                rep_tr = rpart_z(train_idx,1);
                
                cv_ons = mTRFcrossval(ons_tr,rep_tr,fs,direc,tmin,tmax,lambdas,'Verbose',0);
                
                %get the optimal regression parameter over all channels
                l = mean(cv_ons.r,3); %over channels
                [l_val,l_idx] = max(mean(l,1));
                lons_opt = lambdas(l_idx);
                
                cv_x = mTRFcrossval(x_tr,rep_tr,fs,direc,tmin,tmax,lons_opt,'Verbose',0);
                
                %get the optimal regression parameter over all channels
                l = mean(cv_x.r,3); %over channels
                [l_val,l_idx] = max(mean(l,1));
                lx_opt = lambdas(l_idx);
                
                %start the training
                model_x = mTRFtrain(x_tr,rep_tr, fs,direc,tmin,tmax,lx_opt,'verbose',0);
                model_ons = mTRFtrain(ons_tr,rep_tr, fs,direc,tmin,tmax,lons_opt,'verbose',0);

                %do the cross prediction test
                [~,stats_xo] = mTRFpredict(ons_ts,rep_ts,model_x,'verbose',0); %test the x model on onset data
                acc_x_ons(s,k,ao,tr) = mean(stats_xo.r);
                
                [~,stats_onsx] = mTRFpredict(x_ts,rep_ts,model_ons,'verbose',0); %test the onset model on x data data
                acc_ons_x(s,k,ao,tr) = mean(stats_onsx.r);
                
                %do the own prediction test
                [~,stats_xx] = mTRFpredict(x_ts,rep_ts,model_x,'verbose',0); %test the x model on x data
                acc_x_x(s,k,ao,tr) = mean(stats_xx.r);
                
                [~,stats_onon] = mTRFpredict(ons_ts,rep_ts,model_ons,'verbose',0); %test the onset model on onset data data
                acc_ons_ons(s,k,ao,tr) = mean(stats_onon.r);
                
            end
        end
    end
end
%save as structure
cd(saving_path)

OT16_results = struct;
OT16_results.auditory = auditory;
OT16_results.acc_ons_x = acc_ons_x;
OT16_results.acc_x_ons = acc_x_ons;
OT16_results.acc_ons_ons = acc_ons_ons;
OT16_results.acc_x_x = acc_x_x;
OT16_results.txt = 'sound identity marker model'
save('OT16_allmodlsg_518_20.mat','-struct','OT16_results')

fig_path= [fig_path '\cross_pred\combined\']


%% Plottting cross prediction values descriptives
x = {'alarm','irregular','odd','odd','odd','irregular'}
ons = {'mTRF envelope','onset','onset','irregular','alarm','alarm'}
audi = {'odd','alarm','irregular','onset','mTRF envelope'}

audi_full = {'alarm','odd','irregular','onset'}
ylim= [-0.05 0.15]
xlim = [-0.05 0.15]
sz = 20;


figure,clf
tiledlayout(size(audi,2),size(audi,2)-1)

for a = 1:length(audi)
    %get the respective target data
%     target_idx = contains(auditory,audi{a})
    xt_idx = contains(x,audi{a})
    onst_idx = contains(ons,audi{a})
    
    %get the target data
    xt_dat = squeeze(acc_x_x(:,:,xt_idx,:));
    onst_dat = squeeze(acc_ons_ons(:,:,onst_idx,:));
    target_dat = cat(3, reshape( permute(xt_dat,[1 2 4 3]),length(sbj),n_fold*length(task),[]), reshape( permute(onst_dat,[1 2 4 3]),length(sbj),n_fold*length(task),[]));
    
    %     target_dat = cat(3, squeeze(permute(xt_dat,[1 3 2])),squeeze(permute(onst_dat,[1 3 2])) );
    target_dat = reshape(target_dat,[],size(target_dat,3));
    %get the partner data
    xp_dat = squeeze(acc_ons_x(:,:,xt_idx,:));
    onsp_dat = squeeze(acc_x_ons(:,:,onst_idx,:));
    targep_dat = cat(3, reshape( permute(xp_dat,[1 2 4 3]),length(sbj),n_fold*length(task),[]), reshape( permute(onsp_dat,[1 2 4 3]),length(sbj),n_fold*length(task),[]));
    
%     targep_dat = cat(3, squeeze(permute(xp_dat,[1 3 2])),squeeze(permute(onsp_dat,[1 3 2])));
    targep_dat = reshape(targep_dat,[],size(targep_dat,3));
    %get the labels
    xv_la = audi{a};
    yv_la = [ons(xt_idx), x(onst_idx)];
    
    for i = 1:size(target_dat,2)
        nexttile

        
        %get the data
        t_dat = targep_dat(:,i);
        p_dat = target_dat(:,i);
        [cor,p] = corrcoef(t_dat,p_dat)
        p_s(a,i) = p(1,2);
        
        xv = xlim;
        yv = ylim;
        d = (p_dat - xv(1)).*(yv(2) - yv(1)) - (t_dat - yv(1)).*(xv(2) - xv(1) );
        dis = abs(t_dat - p_dat)/sqrt(2);
        
        %plot
        s = scatter(p_dat,t_dat,sz,audi_colorsrgb(yv_la{i}),'filled');
        
        %             s.AlphaData = dis;
        %             s.MarkerFaceAlpha = 'flat';
        hold on
        ls = lsline(gca)
        ls.LineWidth = 2;
        ls.Color = 'k';
        
        s1 = scatter(p_dat(d > 0),t_dat(d>0),sz,audi_colorsrgb(xv_la),'filled')
        %             s1.AlphaData = dis(d<0);
        %             s1.MarkerFaceAlpha = 'flat';
        
        hold on
        plot(xlim,ylim,'--','linew',2)
        ylabel(yv_la{i})%,'FontSize',28)
        xlabel(xv_la)%,'FontSize',28)
        legend(ls,sprintf('r=%0.2f || p=%0.2f',cor(1,2),p(1,2)))
        title(sprintf('%s on %s',yv_la{1,i},xv_la))%,'FontSize',30)
        set(gca,'Ylim',ylim,'Xlim',xlim)%,'FontSize',26)
        
        axis square
    end
    
end

% save_fig(gcf,fig_path,'onsirrcross')
