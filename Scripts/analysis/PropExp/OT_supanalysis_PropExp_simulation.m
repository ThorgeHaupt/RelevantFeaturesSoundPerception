%% neural response simulation
samplingFrequency = 1000; % Sampling frequency (Hz)
time_aep = 0:1/samplingFrequency:0.3; % Time vector (1 second duration)
numSamples = length(time_aep);
% Generate a simulated AEP waveform with Gaussian peaks
baseline = 0; % Baseline amplitude
peak_amplitudes = [0.8, -1.5, 1]; % Amplitudes of the peaks
peak_latencies = [0.05, 0.1, 0.19]; % Latencies of the peaks
peak_widths = [0.01, 0.02, 0.02]; % Widths of the Gaussian peaks
% Generate Gaussian peaks
aep = baseline * ones(1, numSamples);
for i = 1:length(peak_amplitudes)
    peak = peak_amplitudes(i) * exp(-(time_aep - peak_latencies(i)).^2 / (2 * peak_widths(i)^2));
    aep = aep + peak;
end
aep = aep/max(abs(aep));

% Plot simulated AEP waveform
figure;
plot(time_aep, aep, 'b', 'LineWidth', 2);
xlabel('Time (s)');
ylabel('Amplitude');
title('Simulated Auditory Evoked Potential (AEP) with Gaussian Peaks');
grid on;
%% generate noise
srate = 1000; % sampling rate in Hz
time  = 0:1/srate:601;
pnts  = length(time);
hz    = linspace(0,srate/2,floor(length(time)/2)+1);
% generate 1/f amplitude spectrum
ed = 5000; % exponential decay parameter
as = rand(1,floor(pnts/2)-1) .* exp(-(1:floor(pnts/2)-1)/ed);
as = [as(1) as 0 0 as(:,end:-1:1)];
% Fourier coefficients
fc = as .* exp(1i*2*pi*rand(size(as)));
% inverse Fourier transform to create the noise
noise = real(ifft(fc)) * pnts;
noise_norm = noise/max(abs(noise));
figure(5), clf
subplot(211)
plot(time,noise,'k')
set(gca,'fontsize',15)
title('Pink noise: Time domain')
xlabel('Time (s)'), ylabel('Amplitude')
subplot(223)
[y,x] = hist(noise,100);
plot(x,y,'k','linew',2)
xlabel('Values'), ylabel('N per bin')
title('Signal histogram (distribution)')
set(gca,'fontsize',15)
subplot(224)
amp = abs(fft(noise)/pnts);
amp(2:end) = 2*amp(2:end);
plot(hz,amp(1:length(hz)),'k')
title('Frequency domain')
set(gca,'fontsize',15)
grid on;



%%
srate =1000;
dur = 10; %in seconds
nh = 30;
nl = 4;
time = 0:1/srate:dur;
numSamples = length(time);
burst_l = dur*srate;
burst_buffer = 0.01*srate;  %spacing of onset marker
low_l = dur*srate;
low_buffer = 1*srate;
stim = [];
time_durl = [];
time_durb= [];
%randomly place ones that are at least spaced 50ms
% Generate random markers at least 50 ms apart
for i = 1:600/dur
    if mod(i, 10) == 0%iseven(i)
        %create the segments
        onset = zeros(size(time));
%         n = 1;
%         while true
%             % Randomly select a position for the marker
%             position = randi([1, numSamples]);
%             % Check if the selected position is at least 50 ms away from the previous marker
%             if all(diff(find(onset))) && (isempty(find(onset(position:min(position+round(burst_buffer),numSamples)))) || ~any(onset(position:min(position+round(burst_buffer),numSamples))))
%                 onset(position) = 1/n; % Set the marker at the selected position
%                 n = n+1;
%             end
%             % Check if the desired number of markers is reached
%             if n == 20 % 50 markers
%                 break;
%             end
%         end
        
        onset(randperm(numSamples,nh)) = 1;
        while min(diff(find(onset))) < burst_buffer
            onset = zeros(size(time));
            ons_idx = randperm(numSamples,nh);
            onset(ons_idx) = 1;
        end
        
%         %simulate neural adaptation
%         ons_idx = find(onset);
%         for i =1:sum(onset)
%             onset(ons_idx(i)) = onset(ons_idx(i))/log(i+2);
%         end
            
        time_durb = [time_durb;length(stim) length(stim)+burst_l];
        stim = [stim onset];
    else
        %create the segments
        onset = zeros(size(time));
%         while true
%             % Randomly select a position for the marker
%             position = randi([1, numSamples]);
%             % Check if the selected position is at least 50 ms away from the previous marker
%             if all(diff(find(onset))) && (isempty(find(onset(position:min(position+round(low_buffer),numSamples)))) || ~any(onset(position:min(position+round(low_buffer),numSamples))))
%                 onset(position) = 1; % Set the marker at the selected position
%             end
%             % Check if the desired number of markers is reached
%             if sum(onset)== 2 % 50 markers
%                 break;
%             end
%         endonset(randperm(numSamples,20)) = 1;
        onset(randperm(numSamples,nl)) = 1;
        while min(diff(find(onset))) < burst_buffer
            onset = zeros(size(time));
            ons_idx = randperm(numSamples,nl);
            onset(ons_idx) = 1;
        end
        
        time_durl = [time_durl;length(stim) length(stim)+burst_l];
        stim = [stim onset];
    end
    disp(i)
    
end


stims = stim;

%generate the signal 
snr = 10;
%instead of addition, convolve the signal
y = conv([zeros(1,length(aep)) stims zeros(1,length(aep))],aep,'valid');

%% prop exp simulation 
nfold = 5;
testfold = 1;
snrs = [-30:2.5:10]%round(logspace(log(1),2,10),0);
%test different SNRs 

%filter the signal
fs = 1000; % Sampling frequency (Hz)
fpass = [0.3, 20]; % Passband frequencies (Hz)

% Define the filter order
filter_order = 100;

% Calculate normalized frequencies
fpass_norm = fpass / (fs / 2);

% Design the Hanning FIR filter
b = fir1(filter_order, fpass_norm, 'bandpass', hann(filter_order + 1));
filtered_data = filter(b, 1, signal);
for i=1:length(snrs)
    %     p_signal = mean(abs(y).^2);
    %
    %     p_noise = p_signal / (10^(snrs(i) / 10));
    %     snr = snrs(i);
    
    desired_snr_db = snrs(i); % Desired SNR level in dB
    desired_snr_linear = 10^(desired_snr_db / 10); % Convert dB to linear scale
    desired_std = mean(abs(aep(1,100:end)))/desired_snr_linear;
    scaling = desired_std/std(noise_norm);
    
    signal = y + (scaling .* noise_norm(1,1:length(y))); % Adjusted noisy signal
    
    
    
    
%     signal = y+ (noise_norm .*sqrt(p_noise));
%     std(noise_norm .*sqrt(p_noise))

    

    
    %% model training
    onset = [stim zeros(1,length(y)-length(stim))];
    [strain,rtrain,stest,rtest] = mTRFpartition(onset',signal',nfold,testfold);
    strainz = strain;
    stestz = stest;
    
    
    rtrainz = cellfun(@(x) zscore(x,[],'all'),rtrain,'UniformOutput',false);
    rtestz = zscore(rtest,[],'all');
    
    fs = srate;
    
    %TRF parameters
    Dir = 1; %specifies the forward modeling
    tmin = -100;
    tmax = 500;
    lambdas = linspace(10e-4,10e4,10);
    cv = mTRFcrossval(strainz,rtrainz,fs,Dir,tmin,tmax,lambdas,'Verbose',0);
    
    %get the optimal regression parameter
    l = mean(cv.r,3); %over channels
    [l_val,l_idx] = max(mean(l,1));
    l_opt = lambdas(l_idx);
    
    %train the neural model on the optimal regularization parameter
    model_train = mTRFtrain(strainz,rtrainz,fs,Dir,tmin,tmax,l_opt,'verbose',0);
    
    %predict the neural data
    [PRED,STATS] = mTRFpredict(stestz,rtestz,model_train,'verbose',0);
    
    %that is the normal prediction accuracy
    reg(i,1,:) = STATS.r;
    
    %normalize to find zero values easier
    pred = PRED./max(abs(PRED));
    pred_m = pred- mean(pred);
    ch = 1;
    c = arrayfun(@(x) length(find(pred_m(:,ch) == x)), unique(pred_m(:,ch)), 'Uniform', false);
    prop_exp(ch,:) = 1-max(cell2mat(c),[],'all')/length(pred_m(:,ch)); %should be the same for every channel
    
    if unique(prop_exp) > 2
        disp('you fucked it up')
    end
    reg(i,2,:) = prop_exp; %percentage of explained data
    
    
    
    %and now for the data reduced part
    %find onsets
    ons_idx = find(sum(stestz,2) == 1);
    onsac_idx = ons_idx;
    
    ons_epo = [];
    ons_epo(:,1) = ons_idx - abs(((tmin/1000)*srate));
    ons_epo(:,2) = ons_idx + abs(((tmax/1000)*srate));
    
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
    reg(i,3,:) = mean(STATS.r,1);
    
end

figure
subplot(1,2,1)

p(1) = plot(snrs,reg(:,1),'-o') 
hold on 
p(2) = plot(snrs,reg(:,3),'-x','Color',"#D95319")
xlabel('snr level','FontSize', 16)
ylabel('correlation','FontSize', 16)
set(gca,'FontSize', 12)
legend({'whole segment','explainable data'},'FontSize', 12,'box','off')
title('Simulation: Correlation as a function of SNR','FontSize', 18)
box off
auditory = {'onset','alarm','irregular','odd'}
for t = 1:4%length(snr_val)
    x_val = mean(dat(:,t))
    y_query = interp1(snrs,reg(:,3),mean(dat(:,t)),'linear')
    x(t) = line([x_val x_val],[-0.1 y_query],'Color',audi_colorsrgb(auditory{t}),'linew',2)
    hold on
    line([-30 x_val],[y_query y_query],'Color',audi_colorsrgb(auditory{t}),'LineStyle','--','linew',2)
    

end
% legend(x(t),auditory,'Location','northwest')
legend(p,{'whole segment','explainable data'},'Location','northwest','FontSize',16)
set(gca,'Ylim',[-0.1 1],'FontSize',16)

%get the actual data
cd('') %where you save the propexp data 
load('')%load the data 
reg_avg = squeeze(mean(mean(result_reg(:,:,:,3,:),2),5));
audi_int = {'onset','odd','irregular','alarm'}
for i=1:length(audi_int)
    audi_cmp(i,:) = strcmp(audi_int{i},auditory);
    cmp_idx(i,:) = find(audi_cmp(i,:)==1)
    data(:,i) = reg_avg(:,audi_cmp(i,:));
    data_cl(i,:) = audi_colorsrgb(audi_int{i})
    audi_idx(i) = find(strcmp(audi_int{i},auditory));
end

subplot(1,2,2)
vo = violinplot(data ,audi_int,...
    'ViolinColor',data_cl,...
    'ViolinAlpha',0.45,...
    'ShowMean', true)
set(gca,'Ylim',[-0.1 1],'FontSize',16)
box off
axis off
% save_fig(gcf,'\\daten.w2kroot.uni-oldenburg.de\home\lorf0331\Documents\MATLAB\Project\OTtracking\Results\rebuttal\simulate propexp\','snrstim_snrreal_cmp')

