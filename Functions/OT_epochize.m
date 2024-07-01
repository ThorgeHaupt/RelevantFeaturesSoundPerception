function [epo_dat,epo_stim,epo_idx,stim_z] = OT_epochize(EEG,stim,t,ovlp)
%a neat function to extract non-overlapping onsets
%Input:
%EEG        just your normal everyday EEG structure
%stim       your binary stimulus vector
%t          the time window over which the overlap will be computed
%ovlp       binary: 1=non overlapping epochs, 0= all epochs
%
%Output:
%epo_dat    your data structure in non overlapping epoch format
%epo_stim   you binary epochized stimulus vector




%get time
t = t*EEG.srate;
t_min = t(1,1);
t_max = t(1,2);

ons_idx = find(stim==1);

%create zero ons vec
stim_z = zeros(size(stim));

%epoch the data such that 100ms and 500ms after the onset the
%data is extracted
erp_r(:,1) = ons_idx - abs(t_min);
erp_r(:,2) = ons_idx + t_max;

%determine whether epochs are overlapping or not
if ovlp
    %select non-overlapping ERPs
    erp_idx = ones(size(erp_r,1),1);
    for i=1:size(erp_r,1)-1
        if erp_r(i,2) > erp_r(i+1,1)
            erp_idx(i+1,1) = 0;
        else
            erp_idx(i+1,1)= 1;
        end
    end
    epo_idx = erp_r(logical(erp_idx),:);
else
    epo_idx = erp_r;
    
end

%fail safe to not go above and beyond
epo_idx = epo_idx(epo_idx(:,2)<EEG.pnts,:);
epo_idx = epo_idx(epo_idx(:,1)>=0,:);
ons_idx = epo_idx(:,1) + abs(t_min);
%select the non-overlapping stimulus
stim_z(ons_idx,1) = 1;

%select the ERP intervals
for ep = 1:length(epo_idx)
    epo_dat(:,:,ep) = EEG.data(:,epo_idx(ep,1):epo_idx(ep,2));
    epo_stim(:,ep) = stim_z(epo_idx(ep,1):epo_idx(ep,2),:);
end
