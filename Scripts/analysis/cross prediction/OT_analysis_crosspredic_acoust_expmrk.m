direc = 1; %specifies the forward modeling
tmin = -100;
tmax = 500;
lambda1 =  0.05;
resukt= [];
n_fold = 6
OT_setup
acc_env = zeros(length(sbj),2,10);
acc_ons = zeros(length(sbj),2,10);
%manual specification of features to be compared
auditory = {'alarm','mTRF envelope'};


for s=1:length(sbj)
   
   
    for k=1:2
        EEG = [];
        [EEG,PATH] = OT_preprocessing(s,k,sbj,20);
       
        fs = EEG.srate;

       
        %% Just do it yourself

        %get the two stimulus information
        label_ons = string(auditory(2));
        stim_ons = extract_stimulus2(EEG, PATH, label_ons, k,sbj{s});
        label_env = string(auditory(1));
        stim_env = extract_stimulus2(EEG, PATH, label_env, k,sbj{s});
        
        %% determine the peak shift
        ons_t = mTRFtrain(stim_ons,EEG.data',EEG.srate,1,-150,550,0.05,'verbose',0)
        env_t = mTRFtrain(stim_env,EEG.data',EEG.srate,1,-150,550,0.05,'verbose',0)
        
        [r,lags] = xcorr(squeeze(ons_t.w(:,:,1)),squeeze(env_t.w(:,:,1)),20,'coef');
        
        %save the weights
        weights_ons(s,k,1,:,:) = squeeze(ons_t.w(:,6:66,:));
        weights_env(s,k,1,:,:) = squeeze(env_t.w(:,6:66,:));
        
        %find the shift for the TRFs
        [~,ind] = max(r)
        opt_lag = lags(ind)
        
       
       
        %% Partition the two stimulus data sets
        format longG

        %specifiy training and testing data sets
        [spart_ons,rpart] = mTRFpartition(stim_ons, EEG.data',n_fold);
        [spart_env,rpart] = mTRFpartition(stim_env, EEG.data',n_fold);
        
        rpart_z = cellfun(@(x) zscore(x,[],'all'),rpart,'UniformOutput',false);
        train_trials = linspace(1,size(rpart,1),size(rpart,1));
        for tr=1:size(rpart,1)
           
            %select the testing segments
            ons_ts = spart_ons(tr,1);
            env_ts = spart_env(tr,1);
            rep_ts = rpart_z(tr,1);
           
            %select the training segments
            train_idx = setdiff(train_trials,tr,'stable');
           
            ons_tr = spart_ons(train_idx,1);
            env_tr = spart_env(train_idx,1);
            rep_tr = rpart_z(train_idx,1);

            %start the training
            model_ons = mTRFtrain(ons_tr,rep_tr, fs,direc,tmin,tmax,lambda1,'verbose',0);
            model_env = mTRFtrain(env_tr,rep_tr, fs,direc,tmin,tmax,lambda1,'verbose',0);
            
            model_env_cross = model_env;
            model_ons_cross = model_ons;
            %shift the little shits, but only for cross prediction
            x_w = circshift(model_env_cross.w,opt_lag,2);
            model_env_cross.w = x_w(:,abs(opt_lag)+1:end-abs(opt_lag),:);
            model_env_cross.t = model_env_cross.t(1,abs(opt_lag)+1:end-abs(opt_lag));
            
            ons_w = circshift(model_ons_cross.w,-opt_lag,2);
            model_ons_cross.w = ons_w(:,abs(opt_lag)+1:end-abs(opt_lag),:);
            model_ons_cross.t = model_ons_cross.t(1,abs(opt_lag)+1:end-abs(opt_lag));
            
            
            %do the cross prediction test
            [~,stats_xo] = mTRFpredict(ons_ts,rep_ts,model_env_cross,'verbose',0); %test the x model on onset data
            acc_x_ons(s,k,tr) = mean(stats_xo.r);
            
            [~,stats_onsx] = mTRFpredict(env_ts,rep_ts,model_ons_cross,'verbose',0); %test the onset model on x data data
            acc_ons_x(s,k,tr) = mean(stats_onsx.r);
            
            %do the own prediction test
            [~,stats_xx] = mTRFpredict(env_ts,rep_ts,model_env,'verbose',0); %test the x model on x data
            acc_x_x(s,k,tr) = mean(stats_xx.r);
            
            [~,stats_onon] = mTRFpredict(ons_ts,rep_ts,model_ons,'verbose',0); %test the onset model on onset data data
            acc_ons_ons(s,k,tr) = mean(stats_onon.r);

           
        end

       
    end
end

cd(saving_path)
% %save the data 
OT_13crosspred_onsenvshift = struct();
OT_13crosspred_onsenvshift.acc_x_ons = acc_x_ons
OT_13crosspred_onsenvshift.acc_ons_ons = acc_ons_ons
OT_13crosspred_onsenvshift.acc_ons_x = acc_ons_x
OT_13crosspred_onsenvshift.acc_x_x = acc_x_x
OT_13crosspred_onsenvshift.auditory = auditory

save('OT_13crosspred_oddalarmshift.mat','-struct','OT_13crosspred_onsenvshift')

%reshape data

fig_path = [fig_path '\cross_pred\ons_env_shift\']
%% Plotting
% scatter
fig_pos =[240.200000000000,313,1484.40000000000,692.800000000000];
x = auditory{1,1};
ons = auditory{1,2};
sz = 80;


ylim= [-0.05 0.15]
xlim = [-0.05 0.15]

figure,clf
set(gcf,'pos',fig_pos)
tiledlayout(1,2)
a = 1


%get the first data set
xox_dat = reshape(acc_x_x,[],1);
onsox_dat = reshape(acc_ons_x,[],1);
x_la = x;
y_la = ons;

%get the data
t_dat = onsox_dat;
p_dat = xox_dat;
[cor,p] = corrcoef(t_dat,p_dat)

xv = xlim;
yv = ylim;
d = (p_dat - xv(1)).*(yv(2) - yv(1)) - (t_dat - yv(1)).*(xv(2) - xv(1) );
dis = abs(t_dat - p_dat)/sqrt(2);

%plot
nexttile,

s = scatter(p_dat,t_dat,sz,audi_colorsrgb(y_la),'filled');
%             s.AlphaData = dis;
%             s.MarkerFaceAlpha = 'flat';
hold on
ls = lsline(gca)
ls.LineWidth = 2;
ls.Color = 'k';

s1 = scatter(p_dat(d > 0),t_dat(d>0),sz,audi_colorsrgb(x_la),'filled')
%             s1.AlphaData = dis(d<0);
%             s1.MarkerFaceAlpha = 'flat';
hold on
plot(xlim,ylim,'--','linew',2)
ylabel(y_la,'Fontsize',28)
xlabel(x_la,'Fontsize',28)
set(gca,'Ylim',ylim,'Xlim',xlim,'Fontsize',26)
legend(ls,sprintf('r=%0.2f || p=%0.2f',cor(1,2),p(1,2)))

title(sprintf('%s on %s',y_la,x_la),'FontSize',30)

axis square

%get the second data set
onsoons_dat = reshape(acc_ons_ons,[],1);
xoons_dat = reshape(acc_x_ons,[],1);

x_la = ons;
y_la = x;

%get the data
t_dat = xoons_dat;
p_dat = onsoons_dat;
[cor,p] = corrcoef(t_dat,p_dat)

xv = xlim;
yv = ylim;
d = (p_dat - xv(1)).*(yv(2) - yv(1)) - (t_dat - yv(1)).*(xv(2) - xv(1) );
dis = abs(t_dat - p_dat)/sqrt(2);

%plot
nexttile,

s = scatter(p_dat,t_dat,sz,audi_colorsrgb(y_la),'filled');
%             s.AlphaData = dis;
%             s.MarkerFaceAlpha = 'flat';
hold on
ls = lsline(gca)
ls.LineWidth = 2;
ls.Color = 'k';

s1 = scatter(p_dat(d > 0),t_dat(d>0),sz,audi_colorsrgb(x_la),'filled')
%             s1.AlphaData = dis(d<0);
%             s1.MarkerFaceAlpha = 'flat';
hold on
plot(xlim,ylim,'--','linew',2)
ylabel(y_la,'Fontsize',28)
xlabel(x_la,'Fontsize',28)
set(gca,'Ylim',ylim,'Xlim',xlim,'Fontsize',26)
legend(ls,sprintf('r=%0.2f || p=%0.2f',cor(1,2),p(1,2)))

title(sprintf('%s on %s',y_la,x_la),'FontSize',30)
axis square

            
     
% save_fig(gcf,fig_path,'cross_onsmenv')
