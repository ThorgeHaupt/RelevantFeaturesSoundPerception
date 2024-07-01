%This is the main script

%% Pre_processing
%Step1: extract data streams from XDF files and transform into set files
OT_preprocessing00_xdf_to_set.m

%Step2: estimate ICA weights
OT_preprocessing0_ICA.m 

%Step3: extract features
OT_preprocessing1_savStim.m

%Step4: transform audio into wav
OT_preprocessing2_onsetwav.m

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
edit OT_setup.m %set your paths for this script to work 
%% Nested Model Analysis
saving_path = '' %%indicate where the data is supposed to be saved
OT_analysis_nestedmod.m

%plot the data 
%folder of save results
cd(saving_path);
%get the data of all the stuff
env = load('level_analysis_envmaxrep.mat');
ons = load('level_analysis_onsmaxrep.mat');
mel = load('level_analysis_melmaxrep.mat');
base_audi = {'mTRF envelope','onset','mel'};
auditory = base_audi;
new_dat = cat(4,squeeze(mean(mean(env.result_reg,5),4)),...
    squeeze(mean(mean(ons.result_reg,5),4)),...
    squeeze(mean(mean(mel.result_reg,5),4)));
    
fig_pos = [236 110 1004 868]
%do sig testing
for ba = 1:3
    for ao =1:3
        data = squeeze(new_dat(:,ao,:,ba));
        
        [p(ao,1,ba),~,stats] = signrank(data(:,1),data(:,2))
        z(ao,1,ba) = stats.zval
        w(ao,1,ba) = stats.signedrank
        [p(ao,2,ba),~,stats] = signrank(data(:,2),data(:,3))
        z(ao,2,ba) = stats.zval
        w(ao,2,ba) = stats.signedrank
        [p(ao,3,ba),~,stats] = signrank(data(:,1),data(:,3))
        z(ao,3,ba) = stats.zval
        w(ao,3,ba) = stats.signedrank
    end
end
[h, crit_p, adj_ci_cvrg, adj_p] = fdr_bh(reshape(p,[],1),0.05,'dep')
p_re = reshape(adj_p,3,3,[])
z_re = reshape(z,3,3,[])
w_re = reshape(w,3,3,[])
% dat_temp = reshape(permute(dat_temp,[1 4 3 2]),[],3,3);
for ba = 1:3
    for ao = 1:length(auditory)
        data = squeeze(new_dat(:,ao,:,ba));
        figure
        set(gcf,'pos',fig_pos)
        axh{ao} = gca;
        h = violinplot(data ,{'acoustic','sound identity','condition know'},...
            'ViolinAlpha',0.45,...
            'ShowMean', true)
        h(1,1).BoxColor = audi_colorsrgb(base_audi{1,ba})
        h(1,2).BoxColor = audi_colorsrgb(audi_cmb{1,ba}{1,ao})
        h(1,3).BoxColor = audi_colorsrgb(audi_cmb{1,ba}{1,ao})
        h(1,1).ViolinColor{1,1} = audi_colorsrgb(base_audi{1,ba})
        h(1,2).ViolinColor{1,1} = audi_colorsrgb(audi_cmb{1,ba}{1,ao})
        h(1,3).ViolinColor{1,1} = audi_colorsrgb(audi_cmb{1,ba}{1,ao})
        
        
        title(sprintf('%s and %s',base_audi{1,ba},audi_re{ao}),'FontSize',30)
        hold on
        sigstar({[1,2],[2,3],[1,3]},squeeze(p_re(ao,:,ba)))
        ylabel('Correlation')
        set(gca,'ylim',[-0.01 0.20],'FontSize',22)

%         save_fig(gcf,[fig_path '\nestedDesign\'],sprintf('%s and %s',base_audi{1,ba},audi_re{ao}))
    end
    sav_p{ba,1} = squeeze(p_re(:,:,ba));

end

% plot model weights
rnd_tit = {'acoustic','sound identity','cond. knowledge'}
clim = [-650 300]
fig_pos = [357 121 1169 857]
auditory = {'onset'} % needs to be manually changed
load('level_analysis_onsmaxrep.mat');
for ao =1:length(auditory)
    temp_dat = squeeze(mod_w(:,ao,:));
    figure
    set(gcf,'pos',fig_pos)
    t = tiledlayout(3,2)
    for rnd = 1:3
        temp_dat2 = []

        if rnd < 2
            for s = 1:20 
                temp_dat2(s,:) = squeeze(mean(temp_dat{s,rnd},2));%squeeze(mean(temp_dat{s,rnd}(:,11));
            end
            plot_dat=mean(temp_dat2,1);
            nexttile
            imagesc(plot_dat);
            set(gca,'XTick',linspace(1,61,7),'XTickLabel',linspace(-100,500,7),'clim',clim,'FontSize',16)
            title(rnd_tit{rnd})
            nexttile
            p1 = plot(plot_dat, 'Color',audi_colorsrgb(base_audi{1,1}),'Linew',2);
            set(gca,'XTick',linspace(1,61,7),'XTickLabel',linspace(-100,500,7))
            
        else
            for s = 1:20 
                temp_dat2(s,:,:) = squeeze(mean(temp_dat{s,rnd},3));%squeeze(temp_dat{s,rnd}(:,:,11));
            end
            plot_dat=squeeze(mean(temp_dat2,1));
            nexttile
            imagesc(plot_dat);
            set(gca,'XTick',linspace(1,61,7),'XTickLabel',linspace(-100,500,7),'clim',clim,'FontSize',16)
            title(rnd_tit{rnd})
            nexttile
            plot(plot_dat(1,:)','Color',audi_colorsrgb(base_audi{1,1}),'Linew',2);
            hold on
            plot(plot_dat(2:end,:)','Linew',2)
            set(gca,'XTick',linspace(1,61,7),'XTickLabel',linspace(-100,500,7))
            if rnd == 2
                legend({'acoustic','sound identity'},'Location','southeast')
            else
                legend({'acoustic','narrow','wide'},'Location','southeast')
            end

            
        end
        set(gca,'FontSize',16,'Ylim',clim)

    end
    xlabel(t,'Time (ms.)','FontSize',22)
    ylabel(t,'a.u.','FontSize',22)
    title(t,sprintf('%s and %s (model weights)',base_audi{1,1},audi_re{ao}),'FontSize',30)
%     save_fig(gcf,[fig_path '\nestedDesign\'],sprintf('modw_%s_%s',base_audi{1,1},audi_re{ao}))
end

% comparison to Rosenkranz et. al.,(2023)
OT_analysis_Rosenkranz_Cmp.m

%% Envelope Comparison 
saving_path = '' %indicate where the data is supposed to be saved
OT_analysis_envelope_cmp.m

%% Variance Partitioning
saving_path = '';

%requires manual selection of features of interest for plotting
OT_analysis_variancepart.m

%3D variance partitioning
OT_analysis_variancepart_3D.m

%% Proportion Explained
saving_path = '';
OT_analysis_PropExp.m

%simulation of the proportion explained
%estimate actual data SNR
saving_path = '';
OT_supanalysis_SNR_estimation.m 

%simulate the proportion explained
saving_path = '';
OT_supanalysis_PropExp_simulation.m

%% DNS analysis

%TRF
saving_path = '';
OT_analysis_DNS_TRF.m

%ERP
saving_path = '';
OT_analysis_DNS_ERP.m

%% Cross Prediction 

%Sound identity CrossPred
saving_path = '';
OT_analysis_crosspredic_expmrk.m

%Cross prediction acounting for temporal shift
saving_path = '';
OT_analysis_crosspredic_acoust_expmrk.m






