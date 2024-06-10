function plot_spec_rasters(birdname,days,syls,premotor_win)
% loops through chosen days and syllables to create case files via
% get_motif_spikes_from_list.m
% raster_params.m is used to get channel, filelist, and directory

warning('off','all')
failedSyls ={};


for d= 1:length(days)
    day = days(d)
    neuron_name = ['day' num2str(day)];
    
    % Query the intan channel with isolated neuron, list of passed files
    % change to folder containing files
    [ch, songfile,filelist,directory] = raster_params(birdname,day);
    cd(directory)

    for i = 1:length(syls)
        % If looking for individual syllable
        if length(syls) ==1
            syl_seq = char(syls)
        
        else % If looking for a multi-syllable motif
            syl_seq = syls(i);
        end
          
        if file_exist == 0
            
            try   
                plot_multiple_FRs(song_file, motif_sumfile,align_syl_num_all,which_trials,trial_order, burstDetector)
            catch exception
                    message = ['failed for: ' num2str(day) ' -syl ' syl_seq];
                    failedSyls = [failedSyls; message];
                    disp(message)
            end

        else
            cd compile_cases\
            casefile = strcat(birdname,'_',syl_seq,'_',neuron_name,'_ch',num2str(ch),...
                '_premotor_',num2str(premotor_win*1000),'ms_spiketimes_acoustics_',file_exist,'.mat');

            if exist(casefile,"file")
                 load(casefile,'neuralcase')
                eval(sprintf('allCases.%s_%s_day%d = neuralcase;',birdname,syl_seq,day))            
            end
            
            cd ..
        end
    end
end
end