function [all_binned_FRs,mean_binned_FRs, bin_centers] = bin_firing_rates(adjusted_spikes,adj_syl_times,window_to_analyze_s,dt,pre_onset_time,postmotor_cushion)
    
    
    all_binned_FRs = []; % populate later
    
    % the minimum value we will estimate firing rate for is the largest sequence (trial) onset (e.g. the largest amount of time for which we have data for in *every* trial)
    min_syl_onsets = cellfun(@(x) min(x,[],'all'), adj_syl_times, 'UniformOutput', true);
    min_syl_onset = round(min(min_syl_onsets),2);
    
    
    if window_to_analyze_s(1) < min_syl_onset - pre_onset_time %if our analysis window extends into the premotor window, we want to bin our spikes starting there
        minT = window_to_analyze_s(1);
    else
        minT = min_syl_onset - pre_onset_time;
    end
    
    
    % the maximum value we will estimate firing rate for is the smallest sequence (trial) offset (e.g. the largest amount of time for which we have data for in *every* trial)
    max_syl_offset = cellfun(@(x) max(x,[],'all'), adj_syl_times, 'UniformOutput', true); % the largest syllable offset is the end of the trial
    maxT = round(max(max_syl_offset),2); % rounded to the nearest tenth
    
    ts = minT:dt:maxT; % bin edges
    L = length(ts)-1; % number of bins
    
    
    % calculating our bin edges is a little tricky because we always want
    % one edge to = 0, but our motif length is variable. so we need to
    % adjust our minT accordingly.
    
    if ~rem(minT/dt,1)==0
        bins_before_zero = length(find(ts < 0));
        new_minT =  dt*bins_before_zero*-1;
        %disp(['minT adjusted by ',num2str((minT-new_minT)*1000),'ms'])
        ts = new_minT:dt:maxT; % bin edges
        L = length(ts)-1; % number of bins
    end
    
    %new ts
    
    
    
    if window_to_analyze_s(1) < ts(2) % first and last bin numbers can be inaccurate - if analysis window falls into either you need a larger data segment or a smaller analysis window
           warning('Analysis window is in first bin. First and last bin spike counts are not reliable. Adjust analysis window or use a larger data segment.')
    elseif window_to_analyze_s(2) > ts(end-1)
           warning('Analysis window is in last bin. First and last bin spike counts are not reliable. Adjust analysis window or use a larger data segment.')        
    end
    
    bin_centers = ts(1:end-1) + dt*.5; % for plotting
    
    for k = 1:length(adjusted_spikes)
        [binned_spikes, ~] = histcounts(adjusted_spikes{k},ts);
        all_binned_FRs = vertcat(all_binned_FRs,binned_spikes);
    end
    
    mean_binned_FRs = mean(all_binned_FRs)/dt; % spikes per bin and convert to hz
end