%prediction measured by amount of data explainable
OT_setup
%CV values
nfold = 6;
testfold = 1;

%TRF parameters
Dir = 1; %specifies the forward modeling
tmin = -100;
tmax = 500;
lambdas = linspace(10e-4,10e4,10);
%expmrk 
auditory= {'odd','irregular','alarm','onset','random','oddirralarm','oddirralarmons','mTRF envelope','envelope','mel','onsmenv','melons','melmenv','oddirralarmmel','oddirralarmmenv'};

fig_path = [fig_path '\mininf\prop_exp\']

%all models
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
            
            reg(ao,1,:) = STATS.r;
            %normalize to find zero values easier
            pred = PRED./max(abs(PRED));
            pred_m = pred- mean(pred);
            for ch = 1:EEG.nbchan
                 c = arrayfun(@(x) length(find(pred_m(:,ch) == x)), unique(pred_m(:,ch)), 'Uniform', false);
                 prop_exp(ch,:) = 1-max(cell2mat(c),[],'all')/length(pred_m(:,ch)); %should be the same for every channel                 
            end
            if unique(prop_exp) > 2
                disp('you fucked it up')
            end
            reg(ao,2,:) = prop_exp;
            
            %and now for the data reduced part
            %find onsets
            if sum(strcmp(label,{'alarm','irregular','odd','random','oddirralarm','oddirralarmons'}))>0
                ons_idx = find(sum(stestz,2) == 1);
            elseif sum(strcmp(label,{'onset'}))>0
                ons_idx = find(sum(stestz,2) == 1);
                onsac_idx = ons_idx;
            else
                ons_idx= onsac_idx;
            end
            ons_epo = [];   
            ons_epo(:,1) = ons_idx - abs(((tmin/1000)*EEG.srate));
            ons_epo(:,2) = ons_idx + abs(((tmax/1000)*EEG.srate));
            
            for ep = 1:length(ons_idx)
                %epoch them
                if ons_epo(ep,2) < length(stestz) && ons_epo(ep,1)>0
                    stim_epo{ep,1} = stestz(ons_epo(ep,1):ons_epo(ep,2),:);
                    resp_epo{ep,1} = rtestz(ons_epo(ep,1):ons_epo(ep,2),:);
                end
            end
            stim_epo(cellfun('isempty',stim_epo)) = [];
            resp_epo(cellfun('isempty',resp_epo)) = [];
            
            %predict the neural data
            [PRED,STATS] = mTRFpredict(stim_epo,resp_epo,model_train,'verbose',0);
            clear stim_epo resp_epo ons_epo
            reg(ao,3,:) = mean(STATS.r,1);
            

            

        end
        result_reg(s,k,:,:,:) = reg;
    end
    
end
cd(saving_path)
OT21_propexp = struct();
OT21_propexp.auditory = auditory;
OT21_propexp.result_reg = result_reg;
OT21_propexp.txt = 'envelope is computed on segments where we have onsets and the tmax was adjusted to be the fucking same as in all the other analyses';
save('OT21_propexp_epo_tmax500.mat','-struct','OT21_propexp')


%% plot the results
%get the data
temp_dat = squeeze(mean(mean(result_reg,2),5));
%select subsets of data
audi_part = {'alarm','irregular','odd','onset','oddirralarmons','oddirralarm','random'}
audi_env = {'mTRF envelope','oddirralarmmenv'}
audi_mel = {'mel','oddirralarmmel'}
audi_acc = {'onsmenv','melons','melmenv'}
audi_norm = unique([audi_part,audi_env,audi_mel, audi_acc],'stable');
audi_re = {'odd','irregular','alarm','onset','mTRF envelope','mel','onsmenv','oddirralarm','oddirralarmons','oddirralarmmel','oddirralarmmenv','random'}
audi = audi_norm;%
clear audi_cmp cmp_idx data data_cl sc1 sc2 sc3 sc4
for i = 1:length(audi)
    audi_cmp(i,:) = strcmp(audi{i},auditory);
    cmp_idx(i,:) = find(audi_cmp(i,:)==1)
    data(:,i,:) = temp_dat(:,audi_cmp(i,:),:);
    data_cl(i,:) = audi_colorsrgb(audi{i})
end

xlim_val = [];
xlim_val3 = [];
xlim_val4 = [];

mrk = ['v','o','>','d','p','h','_','<','s','*','+','x','|','^']

%start figure
fig = figure
set(fig,'color','w')
ylim = [-0.15 0.73];

for ao = 1:length(audi)
    temp_dat = permute(data,[1 3 2]);
    x_m = squeeze(mean(temp_dat(:,2,ao)));
    x_std = std(temp_dat(:,2,ao));
    if ismember(audi{ao},audi_part)
        y_m = squeeze(mean(temp_dat(:,3,ao)));
        y_std = std(squeeze(temp_dat(:,3,ao)));
    else
        y_m = squeeze(mean(temp_dat(:,1,ao)));
        y_std = std(squeeze(temp_dat(:,1,ao)));
    end
    
    fprintf('%s explains %0.2f of data; accuracy %0.2f, std:%0.2f \n',audi{ao},x_m,y_m,y_std);
    if ismember(audi{ao},audi_part) 
        sb2=subplot(1,4,2);
        sc1(ao) = scatter(temp_dat(:,2,ao),temp_dat(:,3,ao),40,mrk(ao));
        
        sc1(ao).MarkerFaceColor = data_cl(ao,:);
        sc1(ao).MarkerEdgeColor = data_cl(ao,:);
        sc1(ao).MarkerFaceAlpha = 0.5;
        
        hold on
        l = plot(sb2,x_m,y_m,'Marker',mrk(ao),'MarkerFaceColor',data_cl(ao,:),'MarkerEdgeColor','k','MarkerSize',20);

    elseif ismember(audi{ao},audi_env)
        sb3=subplot(1,4,3);
%         %adjust prop exp of the envelope to the segments of the onset
%         ons_exp = squeeze(temp_dat(:,:,strcmp(audi,'onset')))
%         temp_dat(:,2,ao) = ons_exp(:,2); 
        sc2(ao) = scatter(temp_dat(:,2,ao),temp_dat(:,1,ao),40,mrk(ao));
        sc2(ao).MarkerFaceColor = data_cl(ao,:);
        sc2(ao).MarkerEdgeColor = data_cl(ao,:);
        sc2(ao).MarkerFaceAlpha = 0.5;

        hold on
        
        l = plot(sb3,x_m,y_m,'Marker',mrk(ao),'MarkerFaceColor',data_cl(ao,:),'MarkerEdgeColor','k','MarkerSize',20);
        xlim_val = [xlim_val; x_m x_std]; 
    elseif ismember(audi{ao},audi_mel)
        sb4=subplot(1,4,4);
%         %adjust prop exp of the envelope to the segments of the onset
%         ons_exp = squeeze(temp_dat(:,:,strcmp(audi,'onset')))
%         temp_dat(:,2,ao) = ons_exp(:,2); 
        sc3(ao) = scatter(temp_dat(:,2,ao),temp_dat(:,1,ao),40,mrk(ao));
        sc3(ao).MarkerFaceColor = data_cl(ao,:);
        sc3(ao).MarkerEdgeColor = data_cl(ao,:);
        sc3(ao).MarkerFaceAlpha = 0.5;

        hold on
        
        l = plot(sb4,x_m,y_m,'Marker',mrk(ao),'MarkerFaceColor',data_cl(ao,:),'MarkerEdgeColor','k','MarkerSize',20);
        xlim_val3 = [xlim_val3; x_m x_std];
%     elseif ismember(audi{ao},audi_acc)
%         sb5=subplot(1,5,5);
%         %         %adjust prop exp of the envelope to the segments of the onset
%         %         ons_exp = squeeze(temp_dat(:,:,strcmp(audi,'onset')))
%         %         temp_dat(:,2,ao) = ons_exp(:,2);
%         sc4(ao) = scatter(temp_dat(:,2,ao),temp_dat(:,1,ao),40,mrk(ao));
%         sc4(ao).MarkerFaceColor = data_cl(ao,:);
%         sc4(ao).MarkerEdgeColor = data_cl(ao,:);
%         sc4(ao).MarkerFaceAlpha = 0.5;
%         
%         hold on
%         
%         l = plot(sb5,x_m,y_m,'Marker',mrk(ao),'MarkerFaceColor',data_cl(ao,:),'MarkerEdgeColor','k','MarkerSize',20);
%         xlim_val4 = [xlim_val4; x_m x_std];
%         
%     
    end
    
end


%remove empty cells
arefigures = arrayfun(@(x) ~isempty(fieldnames(x)), sc2);
sc2 = sc2(arefigures);

%remove empty cells
arefigures = arrayfun(@(x) ~isempty(fieldnames(x)), sc3);
sc3 = sc3(arefigures);

% %remove empty cells
% arefigures = arrayfun(@(x) ~isempty(fieldnames(x)), sc4);
% sc4 = sc4(arefigures);
%%

%get first position
pos1 = get(sb2,'Position')
pos1(3) = pos1(3)*2
pos1(1) = pos1(1)*1.3
set(sb2,'Position', pos1);

%set position for the #3 subplot
pos2 = get(sb3,'Position')
pos2(3) = pos2(3)*0.5
pos2(1) = pos1(1)+0.28
pos(2) = pos1(2)
set(sb3,'Position', pos2);

%set position for the #4 subplot
pos3 = get(sb4,'Position')
pos3(3) = pos3(3)*0.5
pos3(1) = pos2(1)+0.1
pos(2) = pos1(2)
set(sb4,'Position', pos3);

% %set position for the #5 subplot
% pos4 = get(sb5,'Position')
% pos4(3) = pos4(3)*0.5
% pos4(1) = pos3(1)+0.1
% pos4(2) = pos1(2)
% set(sb5,'Position', pos4);

%first plot
set(sb2,'XTick',linspace(0,1,11),'XTickLabel',linspace(0,100,11),'Xlim',[0 0.45],'Ylim',ylim,'YColor','none','FontSize',12)
%second plot
set(sb3,...
    'Xlim',[mean(xlim_val(:,1))-mean(xlim_val(:,2))*4 mean(xlim_val(:,1))+mean(xlim_val(:,2))*4],...
    'Ylim',ylim,...
    'YTickLabel',[],...
    'XTick',mean(xlim_val(:,1)),...
    'XTickLabel',100,...
    'FontSize',12,...
    'YColor','none')

set(sb4,...
    'Xlim',[mean(xlim_val3(:,1))-mean(xlim_val3(:,2))*4 mean(xlim_val3(:,1))+mean(xlim_val3(:,2))*4],...
    'Ylim',ylim,...
    'YTickLabel',[],...
    'XTick',mean(xlim_val3(:,1)),...
    'XTickLabel',100,...
    'FontSize',12,...
    'YColor','none')

% set(sb5,...
%     'Xlim',[mean(xlim_val4(:,1))-mean(xlim_val4(:,2))*4 mean(xlim_val4(:,1))+mean(xlim_val4(:,2))*4],...
%     'Ylim',ylim,...
%     'YTickLabel',[],...
%     'XTick',mean(xlim_val4(:,1)),...
%     'XTickLabel',100,...
%     'FontSize',12,...
%     'YColor','none')


% Create a common legend outside of subplots
lgd = legend(sb4,[sc1,sc2,sc3],audi,'Location','northeast','box','off','FontSize',16)
lgd.ItemTokenSize = [30,80]


% % Set the position of the legend
% legendPosition = get(lgd, 'Position');
% legendWidth = 0.15; % Adjust the width of the legend
% legendPosition(1) = 1 - legendWidth; % Set x position of legend
% legendPosition(3) = legendWidth; % Set width of legend
% set(lgd, 'Position', legendPosition);

%reshuffle the labels for the violin plot
% temp_dat = squeeze(mean(reshape(result_reg,[],length(auditory),3,EEG.nbchan),4));
temp_dat = squeeze(mean(mean(result_reg,2),5));
audi_re = {'odd','irregular','alarm','oddirralarm','onset','mTRF envelope','mel','oddirralarmons','oddirralarmmenv','oddirralarmmel','random'}
audi = audi_re;%
clear audi_cmp cmp_idx data data_cl 
for i = 1:length(audi)
    audi_cmp(i,:) = strcmp(audi{i},auditory);
    cmp_idx(i,:) = find(audi_cmp(i,:)==1)
    if ismember(audi{i},audi_part)
        data(:,i) = squeeze(temp_dat(:,audi_cmp(i,:),3));
    else
        data(:,i) = squeeze(temp_dat(:,audi_cmp(i,:),1));
    end
    data_cl(i,:) = audi_colorsrgb(audi{i})
    audi_idx(i) = find(strcmp(audi{i},auditory));

end


%violinplot of the distributions
sb1 = subplot(1,4,1)
vo = violinplot(data ,audi,...
    'ViolinColor',data_cl,...
    'ViolinAlpha',0.45,...
    'ShowMean', true)
set(sb1,'Ylim',ylim,'FontSize',16)
ylabel(sb1,'Prediction Accuracy','FontSize',22)

% statistical testing
n_tests = (length(audi)-1)*length(audi)/2;
for i = 1:length(audi)
    for j = 1:length(audi)
        %use t-test or signrank
        if i ~= j
            [p(i,j),h(i,j),stats] = signrank(data(:,i),data(:,j));
            if isnan(h(i,j)); h(i,j) = 0; end
            if isfield(stats,'zval')
                z(i,j) = stats.zval;
                efsz(i,j) = z(i,j)/(sqrt(length(sbj)*2));
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

%reassign these values to a matrix
h_c_m = zeros(length(audi));
adj_p_m = ones(length(audi));

upper_indices = triu(true(length(audi)), 1); % Indices of upper triangular part
h_c_m(upper_indices) =h_c;
adj_p_m(upper_indices) = adj_p;

h_c_m = triu(h_c_m)+triu(h_c_m,1)';
adj_p_m = triu(adj_p_m)+triu(adj_p_m,1)';
hold on 


%get the p_values
sig_pair = [1,2;...
    2,3;...
    3,4;...
    5,6;...
    6,7;...
    5,7;...
    7,8;...
    9,10]%;...
    %7,9];    

% Initialize a cell array to store the reshaped data
sig_idx = cell(size(sig_pair, 1), 1);

% Loop to reshape the data
for i = 1:size(sig_pair, 1)
    sig_idx{i} = sig_pair(i, :);
end


sigstar(sig_idx,diag(adj_p_m(sig_pair(:,1),sig_pair(:,2))))

box off

%set position for the #5 subplot
pos1 = get(sb1,'Position')
pos1(3) = pos1(3)*1.7
pos1(1) = pos1(1)-0.07
set(sb1,'Position', pos1);

%save_fig(gcf,[fig_path '\mininf\prop_exp\'],'pexp_onsenv_3')


%% stats figure
figure
tiledlayout('flow')
nexttile
imagesc(adj_p_m)
set(gca,'XTick',1:1:length(audi),'XTickLabel',audi,'YTick',1:1:length(audi),'YTickLabel',audi,'FontSize',16)
hold on
sc_idx = [];
sc_idx_p = [];
for i =1:length(audi)
    for j =1:length(audi)
%         if h_c_m(i,j) ==1
            text(j, i, num2str(adj_p_m(i,j), '%.3f'), 'Color', 'g', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
%         end
    end
end
title("Model Comparison (p\_val)")
% save_fig(gcf,[fig_path '\mininf\prop_exp\'],'pexp_onsfair_statspv')

nexttile
imagesc(z)
set(gca,'XTick',1:1:length(audi),'XTickLabel',audi,'YTick',1:1:length(audi),'YTickLabel',audi,'FontSize',22)
hold on
sc_idx = [];
sc_idx_p = [];
for i =1:length(audi)
    for j =1:length(audi)
        %if h_int(i,j) ==1
            text(j, i, num2str(z(i,j), '%.3f'), 'Color', 'g', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
        %end
    end
end
title('z\_values','FontSize',30)

nexttile
imagesc(rank)
set(gca,'XTick',1:1:length(audi),'XTickLabel',audi,'YTick',1:1:length(audi),'YTickLabel',audi,'FontSize',22)
hold on
sc_idx = [];
sc_idx_p = [];
for i =1:length(audi)
    for j =1:length(audi)
        %if h_int(i,j) ==1
            text(j, i, num2str(rank(i,j), '%.3f'), 'Color', 'g', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
        %end
    end
end
title('rank\_values','FontSize',30)

nexttile
imagesc(efsz)
set(gca,'XTick',1:1:length(audi),'XTickLabel',audi,'YTick',1:1:length(audi),'YTickLabel',audi,'FontSize',22)
hold on
sc_idx = [];
sc_idx_p = [];
for i =1:length(audi)
    for j =1:length(audi)
        %if h_int(i,j) ==1
            text(j, i, num2str(efsz(i,j), '%.3f'), 'Color', 'g', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
        %end
    end
end
title('effect\_values','FontSize',30)




