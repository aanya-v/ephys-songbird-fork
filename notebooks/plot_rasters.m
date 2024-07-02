birdID= 'br173gr56';
neuronID = 'day69';
syls = 'bcdghjkrm';
neural_channel = 9;
pre_onset_time_s = 0.04;
filelist = 'batch_passed_ch9_clust3.txt';

for i = 1:length(syls)
    currsyl = syls(i);
    % try
    %     quantify_acoustics(currsyl,1,'mat')
    % catch ME
    %     warning(['quantify_acoustics failed for syl ', currsyl])
    %     % rethrow(ME)
    % end
    compile_case(birdID,neuronID,currsyl,neural_channel,pre_onset_time_s,filelist)
end
%% 
song_file = 'br173gr56_220627_125200_songbout2.mat';
align_syl_num_all = 2;
which_trials = 1:100;
trial_order = 'chron';
burstDetector = 0;

for ii = 1:length(syls)
    currsyl = syls(ii);
    % motif_sumfile = ['br177yw112_' currsyl '_day62_ch2_premotor_40ms_spiketimes_acoustics_2024-06-05.mat'];
    motif_sumfile = 'br173gr56_bc_day78_ch11_premotor_40ms_spiketimes_acoustics_2024-07-01.mat'
    try
        plot_multiple_FRs(song_file, motif_sumfile,align_syl_num_all,which_trials,trial_order, burstDetector)
    catch ME
        warning(['Cannot plot current syllable/motif. Motif file may not exist.'])
    end
end
%% 
cd compile_cases

x = dir('*2024-07-01.mat'); 
casenames = {x.name};

cd ..

song_file = 'pk180pu31_220511_095817_songbout2.mat';
align_syl_num_all = 1;
which_trials = 'all';
trial_order = 'chron';
burstDetector = 0;
for j = 1:length(casenames)
    motif_sumfile = casenames{j}
    try
        plot_multiple_FRs(song_file, motif_sumfile,align_syl_num_all,which_trials,trial_order, burstDetector)
    catch ME
        warning(['Cannot plot current syllable/motif. Motif file may not exist.'])
    end
end