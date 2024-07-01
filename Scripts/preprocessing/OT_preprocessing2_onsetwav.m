%% script to prepare the audio files for the onset detection master
%MAINPATH = '\\daten.uni-oldenburg.de\psychprojects$\Neuro\Thorge Haupt\data\Thorge\';

MAINPATH = 'O:\projects\thh_ont\auditory-attention-in-complex-work-related-auditory-envrionments\data files'
addpath(genpath(MAINPATH));
sbj =    {'P001', 'P002','P003','P005', 'P006','P007','P008','P009',...
    'P010','P012','P013', 'P014','P015','P016','P017', 'P018',...
    'P019','P020','P021','P022'};

file = {'narrow_audio.mat','wide_audio.mat'};

%the loop was only done to derive the ICA weights
for s=1:length(sbj)
    
    sb = sbj{s};
    %setup path file
    PATH = fullfile(MAINPATH, sb);
    cd(PATH)   
    
    %loads the audio strct with the game audio sample
    wavnarrow = load(sprintf('%s_narrow_audio_strct.mat',sbj{s}));
    wavwide = load(sprintf('%s_wide_audio_strct.mat',sbj{s}));
    
    %write these files into wav files
    audiowrite(sprintf('narrow_audio_game_%s.wav',sb),wavnarrow.audio_strct.data,wavnarrow.audio_strct.srate,'BitsPerSample',64)
    audiowrite(sprintf('wide_audio_game_%s.wav',sb),wavwide.audio_strct.data,wavwide.audio_strct.srate,'BitsPerSample',64)
    
    
end
