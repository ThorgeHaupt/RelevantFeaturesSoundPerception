OT_setup
fig_path = [fig_path '\variance_part\exp_mrk\'];
%CV values
nfold = 6;
testfold = 1;

%TRF parameters
Dir = 1; %specifies the forward modeling
tmin = -100;
tmax = 500;
lambdas = linspace(10e-4,10e4,10);
%onset
auditory_ons = {'random','onset','alarm','irregular','odd','alarmons','irrons','oddons','oddirrons','oddalarmons','irralarmons','oddirralarm','oddirralarmons'};
%envelope
auditory_env = {'random','envelope','mTRF envelope','alarm','irregular','odd','alarmenv','alarmmenv','oddmenv','irrmenv','irralarmmenv','oddirrmenv','oddalarmmenv','oddirralarmmenv'}
%mel
auditory_mel = {'random','mel','alarm','irregular','odd','alarmmel','irrmel','oddmel','oddirrmel','oddalarmmel','irralarmmel','oddirralarmmel'};
%acoustic comb
auditory_acc = {'onset','envelope','mTRF envelope','mel','onsenv','onsmenv','onsenvmenv','melons','melenv','melmenv','melonsmenv'};
%expmrk 
auditory_exp = {'random','odd','irregular','alarm','oddirr','oddalarm','irralarm','oddirralarm'};
%all models
auditory = unique([auditory_ons, auditory_env auditory_mel auditory_acc auditory_exp],'stable');

for s=1:length(sbj)
    
    
    for k=1:2
        
        EEG = [];
        [EEG,PATH] = OT_preprocessing(s,k,sbj,20);
        
        reg = zeros(length(auditory),EEG.nbchan);
        raw = zeros(length(auditory),EEG.nbchan);
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
            
            reg(ao,:) = STATS.r;
            %% compare to no regularization 
            model_train2 = mTRFtrain(strainz,rtrainz,fs,Dir,tmin,tmax,0.05,'verbose',0);
            %adjust model weights 
            
            [PRED3,STATS3] = mTRFpredict(stestz,rtestz,model_train2,'verbose',0);
            raw(ao,:) = STATS3.r;

        end
        result_reg(s,k,:,:) = reg;
        result_raw(s,k,:,:) = raw;
    end
    
end

cd(saving_path)

OT15_results = struct;
% OT15_results.ons_raw = ons_raw;
% OT15_results.ons_reg = ons_reg;
OT15_results.result_raw = result_raw;
OT15_results.result_reg = result_reg;
OT15_results.auditory = auditory;
OT15_results.txt = 'all of the features are normalized'
% 
save('OT15_result_518_20_withnorm.mat','-struct','OT15_results')

%average over channels
reg_sum = squeeze(mean(result_reg,4));
raw_sum = squeeze(mean(result_raw,4));

%average over conditions
reg_avg = squeeze(mean(reg_sum,2));
raw_avg = squeeze(mean(raw_sum,2));

%load setup
OT_setup

fig_path = [fig_path '\variance_part\exp_mrk\'];
%% Compare model prediction accuracy
%What do i expect
%1 based on preliminary testing unreg>reg data
%2 specific sounds > whole data/ automatic data generation 
%plot the outcome
audi = auditory_acc;
ylim = [-0.03 0.2];
fig_pos = [71,87,1203,907]
figure,clf
set(gcf,'pos',fig_pos)
clear audi_cmp cmp_idx data data_cl
for i = 1:length(audi)
    audi_cmp(i,:) = strcmp(audi{i},auditory);
    cmp_idx(i,:) = find(audi_cmp(i,:)==1)
    data(:,i) = reg_avg(:,audi_cmp(i,:));
    data_cl(i,:) = audi_colorsrgb(audi{i})
end


b = boxplot(data,audi)
box_vars = findall(b,'Tag','Box');
hold on 
xline(4.5,'--','linew',2)

% Fill boxes
for j=1:length(box_vars)
    patch(get(box_vars(j),'XData'),get(box_vars(j),'YData'),audi_colorsrgb(audi{j}),'FaceAlpha',.5);
end
ylabel('Correlation')
set(gca,'Ylim',ylim,'Fontsize',18)
title('Model Performance Accoustic','Fontsize',30)

% save_fig(gcf,'\\smb.uni-oldenburg.de\home\lorf0331\Documents\MATLAB\Project\OTtracking\Results\figs\supplementary\','var_part_allmodl_accoustic_518_20')


%% statistical testing
%nr of unique tests
n_tests = (length(auditory)-1)*length(auditory)/2;
for i = 1:length(auditory)
    for j = 1:length(auditory)
        %use t-test or signrank
        if i ~= j
            [p(i,j),h(i,j),stats] = signrank(reg_avg(:,i),reg_avg(:,j));
            if isnan(h(i,j)); h(i,j) = 0; end
            if isfield(stats,'zval')
                z(i,j) = stats.zval;
                efsz(i,j) = z(i,j)/sqrt(length(sbj));
                rank(i,j) = stats.signedrank;
            end
%             if h(i,j)
%                 % Compute Cohen's d based on the Wilcoxon effect size
%                 dif_dat = reg_avg(:,i)-reg_avg(:,j);
%                 wilcox_es(i,j) = stats.zval / sqrt(n_tests);
%             end
        end
    end
end
%remove double tests, as these make the FDR too conservative
p_real = reshape(triu(p,1),[],1);
p_real(p_real == 0) = [];
%adjust for multiple comparison
[h_c, crit_p, adj_ci_cvrg, adj_p] = fdr_bh(p_real,0.05,'dep','yes');
adj_p(adj_p>1) = 1;

%reassign these values to a matrix
h_c_m = zeros(length(auditory));
adj_p_m = ones(length(auditory));

upper_indices = triu(true(length(auditory)), 1); % Indices of upper triangular part
h_c_m(upper_indices) =h_c;
adj_p_m(upper_indices) = adj_p;

h_c_m = triu(h_c_m)+triu(h_c_m,1)';
adj_p_m = triu(adj_p_m)+triu(adj_p_m,1)';


figure
imagesc(adj_p_m)
set(gca,'XTick',1:1:length(auditory),'XTickLabel',auditory,'YTick',1:1:length(auditory),'YTickLabel',auditory,'FontSize',16)
hold on
sc_idx = [];
sc_idx_p = [];
for i =1:length(auditory)
    for j =1:length(auditory)
        if h_c_m(i,j) ==1
            text(j, i, num2str(adj_p_m(i,j), '%.3f'), 'Color', 'g', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
        end
    end
end
title("Model Comparison (p\_val)")
stats_pv = struct();
stats_pv.adj_p_m = adj_p_m;
stats_pv.labels = auditory;
% save('OT15_stats_pv_withnorm.mat','-struct','stats_pv')
        
% save_fig(gcf,'\\smb.uni-oldenburg.de\home\lorf0331\Documents\MATLAB\Project\OTtracking\Results\figs\supplementary\','vp_st_modcmp_pv_518_20')


%% plot the specific graph 
fig_pos = [139,176.600000000000,1578.80000000000,768];

%feature selection occurs here
audi_s = auditory_exp%{'random','onset','envelope','mTRF envelope','mel','onsmenv','melons','melmenv','melonsmenv'}%;
clear audi_cmp cmp_idx data data_cl audi_idx
%find the indicies
for i=1:length(audi_s)
    audi_cmp(i,:) = strcmp(audi_s{i},auditory);
    cmp_idx(i,:) = find(audi_cmp(i,:)==1)
    data(:,i) = reg_avg(:,audi_cmp(i,:));
    data_cl(i,:) = audi_colorsrgb(audi_s{i})
    audi_idx(i) = find(strcmp(audi_s{i},auditory));
end

%get the p_values
p_int = adj_p_m(audi_idx, audi_idx);
h_int = h_c_m(audi_idx, audi_idx);

z_int = z(audi_idx, audi_idx);
rank_int = rank(audi_idx, audi_idx);
ef_int = efsz(audi_idx, audi_idx);

figure,
tiledlayout('flow')
nexttile
vo = violinplot(data ,audi_s,...
    'ViolinColor',data_cl,...
    'ViolinAlpha',0.45,...
    'ShowMean', true)
hold on 
xline(4.5,'--k','linew',2)
hold on

%manually determine the pairs to be compared
sig_pair = [1,2;...
    2,3;...
    3,4% ;...
    4,5;...
    4,6;...
    5,7;...
    5,8;...
    8,9];
    %6,10;...
    %7,9];    

% Initialize a cell array to store the reshaped data
sig_idx = cell(size(sig_pair, 1), 1);

% Loop to reshape the data
for i = 1:size(sig_pair, 1)
    sig_idx{i} = sig_pair(i, :);
end


sigstar(sig_idx,diag(p_int(sig_pair(:,1),sig_pair(:,2))))

ylabel('Correlation','Fontsize',24)
set(gca,'Xlim',[0 length(audi_s)+1],'Ylim',[-0.025 0.23],'FontSize',22)
set(gcf,'pos',fig_pos)
title('Acoustic Features','FontSize',30)
box off

nexttile
imagesc(p_int)
set(gca,'XTick',1:1:length(audi_s),'XTickLabel',audi_s,'YTick',1:1:length(audi_s),'YTickLabel',audi_s,'FontSize',22)
hold on
sc_idx = [];
sc_idx_p = [];
for i =1:length(audi_s)
    for j =1:length(audi_s)
        %if h_int(i,j) ==1
            text(j, i, num2str(p_int(i,j), '%.3f'), 'Color', 'g', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
        %end
    end
end
title('p\_values','FontSize',30)

nexttile
imagesc(z_int)
set(gca,'XTick',1:1:length(audi_s),'XTickLabel',audi_s,'YTick',1:1:length(audi_s),'YTickLabel',audi_s,'FontSize',22)
hold on
sc_idx = [];
sc_idx_p = [];
for i =1:length(audi_s)
    for j =1:length(audi_s)
        %if h_int(i,j) ==1
            text(j, i, num2str(z_int(i,j), '%.3f'), 'Color', 'g', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
        %end
    end
end
title('z\_values','FontSize',30)

nexttile
imagesc(rank_int)
set(gca,'XTick',1:1:length(audi_s),'XTickLabel',audi_s,'YTick',1:1:length(audi_s),'YTickLabel',audi_s,'FontSize',22)
hold on
sc_idx = [];
sc_idx_p = [];
for i =1:length(audi_s)
    for j =1:length(audi_s)
        %if h_int(i,j) ==1
            text(j, i, num2str(rank_int(i,j), '%.3f'), 'Color', 'g', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
        %end
    end
end
title('rank\_values','FontSize',30)

nexttile
imagesc(ef_int)
set(gca,'XTick',1:1:length(audi_s),'XTickLabel',audi_s,'YTick',1:1:length(audi_s),'YTickLabel',audi_s,'FontSize',22)
hold on
sc_idx = [];
sc_idx_p = [];
for i =1:length(audi_s)
    for j =1:length(audi_s)
        %if h_int(i,j) ==1
            text(j, i, num2str(ef_int(i,j), '%.3f'), 'Color', 'g', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
        %end
    end
end
title('effect\_values','FontSize',30)

% save_fig(gcf,fig_path,char('vp_modlsubcmp_expmrk_518_20'))


      
%% Variance partitioning 
% step1: get the variance explained by all the channels
%is it done for all participants separately, or averaged over participants?
%should be averaged over participants -> reduce individual noise
for k = 1:2
    r2_raw(:,:,k) = squeeze(raw_sum(:,k,:)).^2;
    r2_reg(:,:,k) = squeeze(reg_sum(:,k,:)).^2;
end
r2 = cat(4,r2_raw,r2_reg);
%mean over conditions
r2_m = squeeze(mean(r2,3));

%%  2 Comp model comaprison averaged over conditions
audi_2D{1,:} = {'onset','mTRF envelope','onsmenv'};
audi_2D{2,:} = {'onset','irregular','irrons'};
audi_2D{3,:} = {'onset','odd','oddons'};
audi_2D{4,:} = {'alarm','odd','oddalarm'};
audi_2D{5,:} = {'irregular','odd','oddirr'};
audi_2D{6,:} = {'irregular','alarm','irralarm'};
audi_2D{7,:} = {'mTRF envelope','alarm','alarmmenv'};
audi_2D{8,:} = {'mTRF envelope','irregular','irrmenv'};
audi_2D{9,:} = {'mTRF envelope','odd','oddmenv'};


audi_2D{1,:} = {'mTRF envelope','onset','onsmenv'};
audi_2D{2,:} = {'mTRF envelope','mel','melmenv'};
audi_2D{3,:} = {'onset','mel','melons'};

for a = 1:size(audi_2D,1)
    
    %search idx of features
    for i=1:length(audi_2D{a,1})
        oe(i,:) = strcmp(audi_2D{a,1}{1,i},auditory);
    end
    figure
    set(gcf,'pos',fig_pos)
    t = tiledlayout(1,1)
    for s=2%1:2 %-> raw vs reg
        
        ovlp = squeeze(r2_m(:,oe(1,:),s)) + squeeze(r2_m(:,oe(2,:),s)) - squeeze(r2_m(:,oe(3,:),s));
        nexttile
        axis equal
        [H,S] = venn([sum(r2_m(:,oe(1,:),s)),sum(r2_m(:,oe(2,:),s))],sum(ovlp),'FaceColor',{audi_colors(audi_2D{a,1}{1,1}),audi_colors(audi_2D{a,1}{1,2})})
        
        
        legend(audi_2D{a,1}{1,1},audi_2D{a,1}{1,2},'FontSize',24)
        axis off
        title('Variance Partitioning','FontSize',30)
        
    end
    
%     save_fig(gcf,[fig_path '\variance_part\acoustic\'],sprintf('var_part_avgcond_518_20_%s_vs_%s',audi_2D{a,1}{1,1},audi_2D{a,1}{1,2}))
end