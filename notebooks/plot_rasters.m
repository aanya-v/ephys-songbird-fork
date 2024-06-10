birdID= 'br177yw112';
neuronID = 'day62';
syls = 'bkcr';
neural_channel = 2;
pre_onset_time_s = 0.04;
filelist = 'batch_passed_files.txt';

for i = 1:length(syls)
    currsyl = syls(i);
    try
        quantify_acoustics(currsyl,1,'mat')
    catch ME
        warning(['quantify_acoustics failed for syl ', currsyl])
        % rethrow(ME)
    end
    compile_case(birdID,neuronID,currsyl,neural_channel,pre_onset_time_s,filelist)
end
%% 
song_file = 'br177yw112_221127_194802_songbout2.mat';
align_syl_num_all = 1;
which_trials = 1:75;
trial_order = 'chron';
burstDetector = 0;

for ii = 1:length(syls)
    currsyl = syls(ii);
    motif_sumfile = ['br177yw112_' currsyl '_day62_ch2_premotor_40ms_spiketimes_acoustics_2024-06-05.mat'];
    try
        plot_multiple_FRs(song_file, motif_sumfile,align_syl_num_all,which_trials,trial_order, burstDetector)
    catch ME
        warning(['Cannot plot current syllable/motif. Motif file may not exist.'])
    end
end