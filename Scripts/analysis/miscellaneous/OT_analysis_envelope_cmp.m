%% Model Comparison using variance paritioning
OT_setup

Dir = 1; %specifies the forward modeling
tmin = -100;
tmax = 500;
lambdas = linspace(10e-4,10e4,10);

%partition the data set
nfold = 6;
testfold = 1;
auditory = {'envelope','mTRF envelope'};
% auditory= {'melmenv','melonsmenv'}
for s=1:length(sbj)
    
    
    for k=1:2
        
        EEG = [];
        [EEG,PATH] = OT_preprocessing(s,k,sbj,20);
        
        reg = zeros(length(auditory),EEG.nbchan);
        raw = zeros(length(auditory),EEG.nbchan);
        for ao = 1:length(auditory)
            
            %extract the stimulus
            label = string(auditory(ao));
            stim = extract_stimulus2(EEG, PATH, label, k,sbj{s},task);

            %get the neural data
            resp = double(EEG.data');
            
            if size(resp,1)>size(stim,1)
                resp = resp(1:size(stim,1),:);
            end
            
            [strain,rtrain,stest,rtest] = mTRFpartition(stim,resp,nfold,testfold);
            strainz = strain; %normalization occurs at the stim extraction fun
            stestz = stest;   %normalization occurs at the stim extraction fun


            rtrainz = cellfun(@(x) zscore(x,0,'all'),rtrain,'UniformOutput',false);
            rtestz = zscore(rtest,[],'all');
            
            %% use cross-validation
            fs = EEG.srate;
            
            
            cv = mTRFcrossval(strainz,rtrainz,fs,Dir,tmin,tmax,lambdas,'Verbose',0);
            
            %get the optimal regression parameter
            l = mean(cv.r,3); %over channels
            [l_val,l_idx] = max(mean(l,1));
            l_opt = lambdas(l_idx);
            
            %save the lambda
%             l_save(s,k,ao,:) = lambda_opt;

            %train the neural model on the optimal regularization parameter
            model_train = mTRFtrain(strainz,rtrainz,fs,Dir,tmin,tmax,l_opt,'verbose',0);
            weights(s,k,ao,:) = squeeze(mean(model_train.w,3));

%             w = [];
%             b = [];
%             rtrain_chan = [];
%             for chan=1:EEG.nbchan
%                 for i=1:size(rtrainz,1)
%                     rtrain_chan{i,1} = rtrainz{i,1}(:,chan);
%                 end
%                 w = cat(3,model_train.w,w);
%                 b = cat(2,model_train.b,b);
%             end
%             model_train.w = w;
%             model_train.b = b;
%             model_train.w = model_train.w(1,13:66,:);
%             model_train.t = model_train.t(1,13:66);
            %predict the neural data
            [PRED,STATS] = mTRFpredict(stestz,rtestz,model_train,'verbose',0);
            
            reg(ao,:) = STATS.r;            
            %% compute the raw scores          
%             model_train2 = mTRFtrain(strainz,rtrainz,fs,Dir,tmin,tmax,0.05,'verbose',0);
%             
%             [PRED3,STATS3] = mTRFpredict(stestz,rtestz,model_train2,'verbose',0);
%             raw(ao,:) = STATS3.r;
        end
        result_reg(s,k,:,:) = reg;
%         result_raw(s,k,:,:) = raw;
    end
    
end

cd(saving_path)

%save it as a structure with the corresponding labels
analysis08_results = struct;
analysis08_results.result_reg = result_reg;
analysis08_results.auditory = auditory;
analysis08_results.descp = 'cmp envelope vs. mTRFenvelope';
save('envcmp.mat','-struct','analysis08_results')




%average over channels
reg_sum = squeeze(mean(result_reg,4));

%average over conditions
reg_avg = squeeze(mean(reg_sum,2));


OT_setup
fig_path = [fig_path '\variance_part\acoustic\']
%% averaged over conditions
ylim = [-0.025 0.16];
figure,clf
subplot(1,3,1)
h = boxplot(reg_avg,auditory)
% Find all boxes
box_vars = findall(h,'Tag','Box');

% Fill boxes
for j=1:length(box_vars)
    patch(get(box_vars(j),'XData'),get(box_vars(j),'YData'),audi_colorsrgb(auditory{j}),'FaceAlpha',.5);
end
ylabel('Correlation')
set(gca,'Ylim',ylim)
[p,~,stats] = signrank(reg_avg(:,1),reg_avg(:,2));
sigstar( {[1,2]},[p])
title('averaged over cond reg (norm)')
set(gca,'FontSize',14)
box off
legend(auditory,'Box','off')
title('Prediction Accuracy')

subplot(2,1,2)
h = boxplot(raw_avg,auditory)
% Find all boxes
box_vars = findall(h,'Tag','Box');

% Fill boxes
for j=1:length(box_vars)
    patch(get(box_vars(j),'XData'),get(box_vars(j),'YData'),audi_colorsrgb(auditory{j}),'FaceAlpha',.5);
end
ylabel('Correlation')
set(gca,'Ylim',ylim)
title('averaged over cond raw (norm)')

%% envelope comparison

subplot(1,3,2)
plot(model_train.t,squeeze(mean(mean(weights(:,:,1,:),2),1))','Color',audi_colorsrgb('envelope'),'linew',2)
hold on
plot(model_train.t,squeeze(mean(mean(weights(:,:,2,:),2),1))','Color',audi_colorsrgb('mTRF envelope'),'linew',2)
xlabel('Time in ms.')
ylabel('a.u.')
title('Model Weights')
set(gca,'FontSize',14)
box off

%% Frequency content 
%comparison of the two different envelopes
OT_setup

Dir = 1; %specifies the forward modeling
tmin = -100;
tmax = 500;
lambdas = linspace(10e-4,10e4,10);



s=1;
k=2;

EEG = [];
[EEG,PATH] = OT_preprocessing(s,k,sbj,20);


%extract stimulus
stim_env = extract_stimulus2(EEG, PATH, 'envelope', k,sbj{s},task);
stim_menv = extract_stimulus2(EEG, PATH, 'mTRF envelope', k,sbj{s},task);

%normalize this to range between 0 and 1
stim_env_norm = (stim_env - min(stim_env,[],'all')) / ( max(stim_env,[],'all') - min(stim_env,[],'all') );

stim_menv_norm = (stim_menv - min(stim_menv,[],'all')) / ( max(stim_menv,[],'all') - min(stim_menv,[],'all') );

model_env = mTRFtrain(stim_env_norm,EEG.data',100,1,tmin,tmax,0.05,'verbose',0);
model_menv = mTRFtrain(stim_menv_norm,EEG.data',100,1,tmin,tmax,0.05,'verbose',0);
fs = EEG.srate;
%compute the frequency representations
x = [stim_env_norm stim_menv_norm];
% Compute FFT
N = length(x); % Length of the time series
X = fft(x); % Compute FFT
X_mag = abs(X(1:N/2+1,:)); % Extract the magnitude spectrum (discard negative frequencies)
frequencies = (0:N/2) * fs / N; % Frequency axis (Hz)

% % Plot magnitude spectrum
% figure;
% plot(frequencies, X_mag);
% xlabel('Frequency (Hz)');
% ylabel('Magnitude');
% title('Magnitude Spectrum (FFT) of Time Series');
% set(gca,'Ylim',[0 500]);
% 
% %plot the difference
% subplot(1,3,3)
% plot(frequencies,X_mag(:,2)-X_mag(:,1),'linew',1.5)
% set(gca,'Xlim',[-0.1 1],'FontSize',14)
% xlabel('Frequency in Hz.')
% ylabel('\Delta mtrf env - env')
% title('Frequency Difference')
% box off
% set(gca,'FontSize',14)


% wavelet parameters
num_frex = 10;
min_freq =  0.1;
max_freq = 8;


% set range for variable number of wavelet cycles
range_cycles = [ 1 10 ];

% parameters (notice using logarithmically spaced frequencies!)
frex  = linspace(min_freq,max_freq,num_frex);
nCycs = linspace(range_cycles(1),range_cycles(end),num_frex);
time  = -0.2:1/EEG.srate:0.2;
half_wave = (length(time)-1)/2;

% FFT parameters
nWave = length(time);
nData = EEG.pnts;
nConv = nWave+nData-1;


% FFT of data (doesn't change on frequency iteration)
dataX = fft( stim_menv_norm - stim_env_norm,nConv);

% initialize output time-frequency data
tf = zeros(num_frex,EEG.pnts);
% loop over frequencies
for fi=1:num_frex
    
    % create wavelet and get its FFT
    s = nCycs(fi)/(2*pi*frex(fi));
    wavelet  = exp(2*1i*pi*frex(fi).*time) .* exp(-time.^2./(2*s^2));

    waveletX = fft(wavelet,nConv);
    
    % question: is this next line necessary?
    waveletX = waveletX./max(waveletX);
    
    % run convolution
    as = ifft(waveletX.*dataX',nConv);
    as = as(half_wave+1:end-half_wave);
    
    % compute ITPC
    tf(fi,:) = abs(as).^2;
end

% plot results
subplot(1,3,3)
contourf(EEG.times(1,1:10000),frex,tf(:,1:10000),40,'linecolor','none')
colormap parula
xlabel('Sample Points')
ylabel('Frequencies (Hz)')
title('\Delta mtrf env - env')
set(gca,'FontSize',14)