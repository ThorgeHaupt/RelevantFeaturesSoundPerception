%% 3D modelle
OT_setup

%CV values
nfold = 6;
testfold = 1;
auditory = {'oddirrons','oddalarmons','irralarmons','oddirralarm','oddalarmmenv','oddirrmenv','irralarmmenv'};%'melonsmenv'}%,

%TRF parameters
Dir = 1; %specifies the forward modeling
tmin = -100;
tmax = 500;
lambdas = linspace(10e-4,10e4,10);

for s=1:length(sbj)
    
    
    for k=1:2
        
        EEG = [];
        [EEG,PATH] = OT_preprocessing(s,k,sbj,20);
        
        
        for ao = 1:length(auditory)

            %extract the stimulus
            label = string(auditory(ao));
            
            %extract stimulus
            stim = extract_stimulus2(EEG, PATH, label, k,sbj{s},task);
            
            
            
            %get the neural data
            resp = EEG.data';
            
            if size(resp,1)>size(stim,1)
                resp = resp(1:size(stim,1),:);
            end
            
            %stimuli follow the name i.e. oddirrons 1. odd, 2. irr, 3.ons
            mod{1,:} = stim(:,1); %i.e. odd
            mod{2,:} = stim(:,2); %i.e irr
            mod{3,:} = stim(:,3); %i.e. ons
            mod{4,:} = stim(:,1:2); %i.e. oddirr
            mod{5,:} = stim(:,[1 3]); %i.e. oddons
            mod{6,:} = stim(:,[2 3]); %i.e. irrons
            mod{7,:} = stim; %i.e. oddirrons
            reg_3D = zeros(size(mod,2),EEG.nbchan);
            raw_3D = zeros(size(mod,2),EEG.nbchan);
            %split the model training into the different model parts
            for st = 1:size(mod,1)
                stim = mod{st,:};
                

                %partition the data set
                [strain,rtrain,stest,rtest] = mTRFpartition(stim,resp,nfold,testfold);
                
                %% z-score the input and output data
                strainz = strain;
                stestz = stest;
                
                
                rtrainz = cellfun(@(x) zscore(x,[],'all'),rtrain,'UniformOutput',false);
                rtestz = zscore(rtest,[],'all');
                
                %% use cross-validation
                fs = EEG.srate;
                
                
                cv = mTRFcrossval(strainz,rtrainz,fs,Dir,tmin,tmax,lambdas,'Verbose',0);
                
                %get the optimal regression parameter
                l = mean(cv.r,3); %over channels
                [l_val,l_idx] = max(mean(l,1));
                l_opt = lambdas(l_idx);
                
                %train the neural model on the optimal regularization parameter
                model_train = mTRFtrain(strainz,rtrainz,fs,Dir,tmin,tmax,l_opt,'verbose',0);
                
                %predict the neural data
                [PRED,STATS] = mTRFpredict(stestz,rtestz,model_train,'verbose',0);
                
                reg_3D(st,:) = STATS.r;
                %% compare to no regularization
                model_train2 = mTRFtrain(strain,rtrain,fs,Dir,tmin,tmax,0.05,'verbose',0);
                
                [PRED3,STATS3] = mTRFpredict(stest,rtest,model_train2,'verbose',0);
                raw_3D(st,:) = STATS3.r;
            end
            result_reg_3D(s,k,ao,:,:) = reg_3D;
            result_raw_3D(s,k,ao,:,:) = raw_3D;
        end
        
    end
    
end
cd('\\smb.uni-oldenburg.de\home\lorf0331\Documents\MATLAB\Project\OTtracking\Results\variance_part\3D\')
OT15_3D = struct()
OT15_3D.auditory = auditory;
OT15_3D.result_reg_3D = result_reg_3D;
OT15_3D.result_raw_3D = result_raw_3D;
OT15_3D.txt = 'only model not included the melonsmenv model'

save('3D_vp_allmodl_518_20','-struct','OT15_3D')


%% average over channels
reg_sum = squeeze(mean(result_reg_3D,5));
raw_sum = squeeze(mean(result_raw_3D,5));
dat_sum = cat(5,raw_sum,reg_sum);

%average over conditions
reg_avg = squeeze(mean(reg_sum,2));
raw_avg = squeeze(mean(raw_sum,2));
dat_avg = cat(4,raw_avg,reg_avg);


%%
audi_3D{1,:} = {'odd','irregular','onset','oddirr','oddons','irrons','oddirrons'};
audi_3D{2,:} = {'odd','alarm','onset','oddalarm','oddons','alarmons','oddalarmons'};
audi_3D{3,:} = {'irregular','alarm','onset','irralarm','irrons','alarmons','irralarmons'};
audi_3D{4,:} = {'odd','irregular','alarm','oddirr','oddalarm','irralarm','oddirralarm'};
audi_3D{5,:} = {'odd','alarm','mTRF envelope','oddalarm','oddmenv','alarmmenv','oddalarmmenv'};
audi_3D{6,:} = {'odd','irregular','mTRF envelope','oddirr','oddmenv','irrmenv','oddirrmenv'};
audi_3D{7,:} = {'irregular','alarm','mTRF envelope','irralarm','irrmenv','alarmmenv','irralarmmenv'};
audi_3D{1,:} = {'mel','onset','mTRF envelope','melons','melmenv','onsmenv','melonsmenv'};


ylim = [-0.025 0.2];

for a= 1:size(audi_3D,1)
    %get z data
    for i = 1:size(audi_3D{a,:},2)
        data_cl(i,:) = audi_colorsrgb(audi_3D{a,:}{:,i})
    end
    figure
    h = violinplot(squeeze(reg_avg(:,a,:)),audi_3D{a,:},...
        'ViolinColor',data_cl,...
        'ViolinAlpha',0.45,...
        'ShowMean', true)
    hold on 
    xline(3.5,'--','linew',2)
    ylabel('Accuracy')
    set(gca,'Ylim',ylim)
    title(sprintf('model: %s',auditory{a}))
%     save_fig(gcf,fig_path,sprintf('vp_3D_%s',auditory{a}))
end

%% Variance partitioning averaged conditions
auditory_3D = {'odd','irregular','alarm','onset','mTRF envelope','oddalarm','oddons','oddirr','irralarm','alarmons','irrons','oddmenv','irrmenv','alarmmenv','','oddalarmons','irralarmons','oddirrons','oddirralarm','oddirrenv','oddalarmmenv','irralarmmenv'};
typ = {'raw','regularized'};
fig_pos = [276.200000000000,306.200000000000,980,653.200000000000];
for a = 1:size(audi_3D,1)
    
    
    %search idx of features
    for i=1:size(audi_3D{a,1},2)
        oe(i,:) = strcmp(audi_3D{a,1}{1,i},auditory_3D);
    end
    
    %compute the relevant data
    figure,clf
    set(gcf,'pos',fig_pos)
    t = tiledlayout(1,1)
    
    %compute the unique explained variance
    r2 = squeeze(reg_avg(:,a,:)).^2;   %compute the three overlap models
    ovlp1 = squeeze(r2(:,1)) + squeeze(r2(:,2)) - squeeze(r2(:,4));
    ovlp2 = squeeze(r2(:,1)) + squeeze(r2(:,3)) - squeeze(r2(:,5));
    ovlp3 = squeeze(r2(:,2)) + squeeze(r2(:,3)) - squeeze(r2(:,6));
    ovlp4 = squeeze(r2(:,7)) + squeeze(r2(:,1)) + squeeze(r2(:,2)) + squeeze(r2(:,3)) - squeeze(r2(:,4)) -squeeze(r2(:,5))-squeeze(r2(:,6));
    nexttile
    axis equal
    [H,S] = venn([mean(r2(:,1)),mean(r2(:,2)), mean(r2(:,3))],[mean(ovlp1),mean(ovlp2),mean(ovlp3),mean(ovlp4)],'FaceColor',{audi_colors(audi_3D{a,1}{1,1}),audi_colors(audi_3D{a,1}{1,2}),audi_colors(audi_3D{a,1}{1,3})},'ErrMinMode','TotalError')
    legend(auditory_3D{logical(sum(oe))},'Fontsize',24)
    title('Variance Partitioning','Fontsize',30)
    axis off
    

    sgtitle(auditory{a})
%     save_fig(gcf,[fig_path '\\variance_part\acoustic\'],sprintf('3d_vp_venn_avgc_%s',auditory{a}))
end

