load('br177yw112_g_day88_ch6_premotor_200ms_spiketimes_acoustics_2023-11-04.mat', 'neuralcase')
spikes = neuralcase.spiketrains; 
onsets = neuralcase.syl_ons_offs(:,1);
for i = 1:length(spikes)
    trialspikes = spikes{i}; 
    n_spikes(i) = length(trialspikes); 
end
mean(n_spikes)
%%

for i = 1:length(spikes)
    trialspikes = spikes{i}; 
    trial_onset = onsets(i); 
    pre_40 = trial_onset - 0.04; % time, 40 ms prior to syl onset 
    trialspikes = trialspikes(trialspikes < trial_onset); 
    trialspikes = trialspikes(trialspikes > pre_40); 
    n_spikes(i) = length(trialspikes); 
end

mean(n_spikes)
