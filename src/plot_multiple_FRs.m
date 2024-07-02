% This script analyzes and plots data generated with
% get_motif_spikes_2020 and plots a raster.

% Other functions called: evsoundin (if using cbins), align_spikes,
% timewarping_lyndie_2019 (if timewarp = 1 in parameters), spect_from_waveform,
% distinguishable_colors (https://www.mathworks.com/matlabcentral/fileexchange/29702-generate-maximally-perceptually-distinct-colors)


%clearvars

%% Data files
function [bursts,data_trials] = plot_multiple_FRs(song_file, motif_sumfile,align_syl_num_all,which_trials,trial_order, burstDetector)
% song_file = 'br177yw112_221223_130500_songbout1.mat'; % rhd or cbin for plotting song. leave empty ([]) if you don't want to plot.
trial_no = '';
varargin = {[motif_sumfile]};
%% Set parameters

% align_syl_num_all = [1]; % syllable NUMBER in each sequence to align to (e.g. if you are comparing
% the 'c' in 'abcef' to the first 'c' in 'abcabc' enter [3,3])
align_position = 'onset'; % 'onset' or 'offset',
dt = .005; % window size for binning spike times in seconds
timewarp =0; %1 for yes, 0 for no


%% Main script - loop through files (varargin) and align spike times, bin firing rates, etc


for i = 1:length(varargin)
    clean_trials = {};   

    %choose which sorting order and/or subsets of trials to be plotted
    if strcmp(trial_order,'chron')
        load(varargin{i},'clean_trials_table','pre_onset_time_s','syl_or_seq');
        sorted = clean_trials_table;
        clean_trials = table2cell(sorted);
    elseif strcmp(trial_order,'pitch')
        load(varargin{i},'clean_trials_table','pre_onset_time_s','syl_or_seq')
        sorted = sortrows(clean_trials_table,'pitch_wgt_avg','ascend');
        clean_trials = table2cell(sorted);
    elseif strcmp(trial_order,'amplitude')
        load(varargin{i},'clean_trials_table','pre_onset_time_s','syl_or_seq')
        sorted = sortrows(clean_trials_table,'amp_at_pitchquant','ascend');
        clean_trials = table2cell(sorted);
    elseif strcmp(trial_order,'spectral entropy')        
        load(varargin{i},'clean_trials_table','pre_onset_time_s','syl_or_seq')
        sorted = sortrows(clean_trials_table,'spect_entropy','ascend');
        clean_trials = table2cell(sorted);
    elseif strcmp(trial_order,'wiener entropy')        
        load(varargin{i},'clean_trials_table','pre_onset_time_s','syl_or_seq')
        sorted = sortrows(clean_trials_table,'wiener_entropy','ascend');
        clean_trials = table2cell(sorted);
    end

    if strcmp(which_trials,'all')
        which_trials = 1:size(clean_trials,1);
    end
    
    data_trials = sorted(which_trials,:);
    clean_trials = clean_trials(which_trials,:);
    
    syl_or_seq_all{i} = syl_or_seq;
    pre_onset_time_all(i) = pre_onset_time_s;

    % Extract relevant info for the plot title
    pattern = '^(.*?)_(.*?)_(.*?)_';
    matches = regexp(motif_sumfile, pattern, 'tokens');
    
    % Extract the bird name, syllable, and neuron name
    if ~isempty(matches)
        bird_name = matches{1}{1};
        syllable = matches{1}{2};
        neuron_name = matches{1}{3};
    else
        disp('No bird info matches found');
    end

    
    
    % align spike times. this includes timewarping if timewarp = 1
    if ~exist('postmotor_cushion')
    postmotor_cushion = 0.05;
    end
    
    [adj_spiketimes, adj_syl_times] = align_spikes(clean_trials,align_syl_num_all(i),align_position,pre_onset_time_all(i),syl_or_seq_all{i},timewarp,postmotor_cushion); %
    adj_spiketimes_all{i} = adj_spiketimes;
    adj_syl_times_all{i} = adj_syl_times;
    
    % Get binned firing rates with error bars
    [all_binned_FRs,mean_binned_FRs, bin_centers] = bin_firing_rates(adj_spiketimes,adj_syl_times,[0 0],dt,pre_onset_time_all(1));
  
    all_sequences_binned_FRs{i} = all_binned_FRs;
    mean_binned_FRs_all{i} = mean_binned_FRs;
    bin_centers_all{i} = bin_centers;

    
    clear  syl_or_seq pre_onset_time_s mean_binned_FRs errorbars_2sd adj_spiketimes adj_syl_times align_syl_num bin_centers
end



% % Check that premotor times match
% if length(unique(pre_onset_time_all)) > 1
%     error('Premotor window mismatch!!! Re-run get_motif_spikes with correct window.')
% end

% Check inputs
if length(align_syl_num_all) ~= length(varargin)
    error('Mismatch between number of files and alignment syllable numbers')
end



%% Plot
if ~isempty(song_file)
    n_fig_panels = length(varargin) + 2; %one raster panel, song panels for each input file, plus two panels (will be merged) for mean FRs
    spect_panels = 1; %plot spectrograms in odd-numbered panels (exlcuding the last one)
    raster_panels = 2:1:n_fig_panels-1; %plot rasters in even-numbered panels
    colors = distinguishable_colors(length(varargin));
    
    %check time warping settings
    if timewarp == 1
        disp('Plotting data with time warping...')
    elseif timewarp == 0
        disp('Plotting data without time warping...')
    else
        warning('Timewarp does not equal 0 or 1. Figure is being plotted without time warping. Set timewarp to 1 if this is not correct.')
    end
    
    
    % load song file
    
    if contains(song_file,'rhd')
        
        load([song_file(1:end-3),'mat'],'board_adc_data','frequency_parameters') % assumes file has been converted to .mat after spike sorting etc
        fs = frequency_parameters.amplifier_sample_rate;
        
    elseif contains(song_file,'mat')
        load([song_file(1:end-3),'mat'],'board_adc_data','frequency_parameters') % assumes file has been converted to .mat after spike sorting etc
        fs = frequency_parameters.amplifier_sample_rate;

    elseif contains(song_file,'cbin')
        
        [board_adc_data, fs] = evsoundin('cd',song_file,'obs0');
        
    end
    
    load([song_file,'.not.mat']) %labels etc
    
    figure;
    
        
        % get song data for spectrogram
        time_vector = (0:numel(board_adc_data)-1)/fs; % in seconds
        
        try
            seq_start_idx = strfind(labels,syl_or_seq_all{1});
            seq_end_idx = seq_start_idx(1) + length(syl_or_seq_all{1}) -1;
            plotseq_onset = onsets(seq_start_idx(1))/1000 - pre_onset_time_all(1); % note that syllable times need to be converted to seconds
            plotseq_offset = (offsets(seq_end_idx)/1000)+postmotor_cushion;
            
            time_minus_onset = time_vector - plotseq_onset; % find closest sample to onset time
            plotseq_onset_sample = find(time_minus_onset > (min(abs(time_minus_onset))-5e-06) & time_minus_onset < (min(abs(time_minus_onset)+5e-06)));
            time_minus_offset = time_vector - plotseq_offset; % find closest sample to offset time
            plotseq_offset_sample = find(time_minus_offset > (min(abs(time_minus_offset))-5e-06) & time_minus_offset < (min(abs(time_minus_offset)+5e-06)));
            
            song_in_seq = board_adc_data(plotseq_onset_sample:plotseq_offset_sample+5000);
            align_syl_time_plotting_trial_idx = seq_start_idx(1) + align_syl_num_all(1) -1;
            
            if strcmp(align_position,'onset')
                align_syl_adj = onsets(align_syl_time_plotting_trial_idx)/1000 - plotseq_onset;
            else
                align_syl_adj = offsets(align_syl_time_plotting_trial_idx)/1000 - plotseq_onset;
            end
            
            
            ax(spect_panels(1)) = subplot(n_fig_panels,1,spect_panels(1)); % spectrogram
            [S1,F1_whole,T1,P1_whole] =spect_from_waveform(song_in_seq,fs,0,[0.9 16]);
       
        
        T1 = T1 - align_syl_adj;
        
        imagesc(T1,F1_whole,log(P1_whole)); set(gca,'YD','n');
        set(gca,'XTick',[],'YTick',[3000 6000 9000],'YTickLabel',{'    3','    6','    9'})
        ylabel('Frequency (kHz)');
        clim([-19 -10])
        colormap("turbo")
%         title(['Trial #' num2str(trial_no)])
        title([bird_name ', ' neuron_name])
        subtitle(['syllable ' syllable])
        end


        %Plot Raster
        for b = 1:length(raster_panels)
        ax(raster_panels(b)) = subplot(n_fig_panels,1,raster_panels(b)); %raster
        try
            LineFormat.LineWidth = 2;
%             plotSpikeRaster(adj_spiketimes_all{b},'PlotType','vertline','VertSpikeHeight', 0.8,'XLimForCell',[-0.1 0.1],'LineFormat',LineFormat);
            plotSpikeRaster(adj_spiketimes_all{b},'PlotType','vertline','VertSpikeHeight',1.0,'XLimForCell',[-0.2 0.6],'LineFormat',LineFormat);
            set(gcf, 'color', 'none');    
            set(gca, 'color', 'none');
            xline(0.0, '--r')
            rast = gca; copygraphics(rast,'ContentType','vector','BackgroundColor','none')
        catch ME
            warning(['Cannot plot RASTER'])
            continue
        end
        set(gca,'XTick',[])
        ylabel('Trials')
        rast = gca; copygraphics(rast);
        end
        set(gcf, 'color', 'none');    
        set(gca, 'color', 'none');
        title(['Trial Order: ' trial_order])
        % xline(-0.00, '--r')
%         xline(-0.02, '-.b')
%         xline(0.04, '-.b')

        if burstDetector ==1 
            % Burst Detector
            bursts = cell(length(clean_trials),1);
            for c = 1:length(clean_trials)
                allSpikes = adj_spiketimes_all{1};
                spikeTrain = allSpikes{c}; 
                [allBurstData,tooShort] = maxInterval(spikeTrain);
                burst_spks = cell2mat(allBurstData);
                if isempty(burst_spks)
                    bursts{c} = 0;
                else 
                    bursts{c} = burst_spks;
                end
            end
            bursts = {bursts};
    
            %Plot Detected bursts on Raster plot
            LineFormat = struct();
            LineFormat.color = [1 0 0];
            LineFormat.LineWidth =0.35;
            for b = 1:length(raster_panels)
            ax(raster_panels(b)) = subplot(n_fig_panels,1,raster_panels(b)); %raster
            plotSpikeRaster(bursts{1},'PlotType','vertline','VertSpikeHeight',0.9,LineFormat',LineFormat);
            end
        else bursts = 'NaN';
        end

    % plot firing rates
    ax(n_fig_panels) = subplot(n_fig_panels,1,n_fig_panels); 
    hold on
    for v = 1:length(varargin)
        
        %plot
        x1 = bin_centers_all{v};
        y1 = mean_binned_FRs_all{v};
        plot(x1,y1,'LineStyle','-','LineWidth',1.5,'Color',colors(v,:))
        
        ylabel('Firing rate (Hz)')
        xlabel('Time (s)')
       
    end


    
    linkaxes(ax,'x');
    ylim([0 500])
%     xline(-0.04, '--r')
    xline(-0.0, '--r')
    xlim([-0.04 0.15])
%     xlim([bin_centers_all{1}(1) bin_centers_all{1}(end)])
    set(findall(gcf,'-property','FontSize'),'FontSize',11)
    %legend(h,'Location','southeast','FontSize',10);
    set(gcf,'color','w');
    scrsz = get(0,'ScreenSize'); %below: x, y, width, height, relative to screen size
    set(gcf,'Position',[scrsz(1)/100 scrsz(2)/4 scrsz(3)*(30/100),scrsz(4)/1]);
%     set(gcf,'Position',[scrsz(3)/100 scrsz(4)/4 scrsz(3)*(98/100),scrsz(2)/2]);
    % set(gcf, 'color', 'none');    
    % set(gca, 'color', 'none');
end