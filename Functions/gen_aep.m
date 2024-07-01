function aep_w = gen_aep(t_l,srate)%create the onset model 
% Parameters
samplingFrequency = srate; % Sampling frequency (Hz)
time_aep = 0:1/samplingFrequency:t_l; % Time vector (1 second duration)
numSamples = length(time_aep);

% Generate a simulated AEP waveform with Gaussian peaks
baseline = 0; % Baseline amplitude
peak_amplitudes = [0.8, -1.5, 1]; % Amplitudes of the peaks
peak_latencies = [0.05, 0.1, 0.2]; % Latencies of the peaks
peak_widths = [0.01, 0.02, 0.02]; % Widths of the Gaussian peaks

% Generate Gaussian peaks
aep = baseline * ones(1, numSamples);
for i = 1:length(peak_amplitudes)
    peak = peak_amplitudes(i) * exp(-(time_aep - peak_latencies(i)).^2 / (2 * peak_widths(i)^2));
    aep = aep + peak;
end

aep_w = aep/max(abs(aep));

% Plot simulated AEP waveform
figure;
plot(time_aep, aep, 'b', 'LineWidth', 2);
xlabel('Time (s)');
ylabel('Amplitude');
title('Simulated Auditory Evoked Potential (AEP) with Gaussian Peaks');
grid on;