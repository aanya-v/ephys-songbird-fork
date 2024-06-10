function [f_cutoff,t_assay,spect_params,syl_name_in_notmat]=syllable_params_by_bird_leila(bname,syl)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% f_cutoff: the frequency range with which to quantify the pitch
% t_assay: the time at which to quantify the acoustics


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
syl_params = struct;

if strcmp(bname,'br173gr56')
    switch(syl)
        case 'a'
            f_cutoff=[500 1000];
            t_assay=50;
            spect_params=[0.9 16];
        case 'b'
            f_cutoff=[900 2100];
            t_assay=50;
            spect_params=[.8 10];
        case 'c'
            f_cutoff=[900 2100];
            t_assay=50;
            spect_params=[.8 10];
        case 'd' 
            f_cutoff=[1100 1400];
            t_assay=85;
            spect_params=[0.9 16];
        case 'e'
            f_cutoff=[500 1000];
            t_assay=50;
            spect_params=[0.9 16];
        case 'f'
            disp('f')
            f_cutoff=[1200 2000];
            t_assay=.05;
            spect_params=[.5 8];
        case 'g'
            f_cutoff=[1100 2400];
            t_assay=.04;
            spect_params=[0.8 8];
        case 'h'
            f_cutoff=[800 1800];
            t_assay=80;
            spect_params=[0.8 8];
        case 'j'
            f_cutoff=[500 1900];
            t_assay=80;
            spect_params=[0.8 10];
        case 'k'
            f_cutoff=[1000 3000];
            t_assay=40;
            spect_params=[0.8 10];
        case 'o'
            f_cutoff=[1100 1400];
            t_assay=85;
            spect_params=[0.9 16];
        case 'p'
            f_cutoff=[1000 2000];
            t_assay=50;
            spect_params=[0.8 10];
        case 'q'
            f_cutoff=[1000 2700];
            t_assay=50;
            spect_params=[0.8 10];
        case 'r'
            f_cutoff=[2100 3500];
            t_assay=50;
            spect_params=[0.8 10];
        case 's'
            f_cutoff=[2000 3000];
            t_assay=50;
            spect_params=[0.5 16];
        case 't'
            f_cutoff=[3300 4300];
            t_assay=50;
            spect_params=[0.8 16];
        case 'x'
            f_cutoff=[4000 6000];
            t_assay= 20;
            spect_params=[0.8 8];
        case 'y'
            f_cutoff=[4000 6000];
            t_assay= 20;
            spect_params=[0.8 8];
    end


end


if strcmp(bname,'br177yw112')
    switch(syl)
        case 'a'
            f_cutoff=[1800 3200];
            t_assay=.025;
            spect_params=[.8 10];
        case 'b'
            f_cutoff=[1000 2500];
            t_assay=50;
            spect_params=[.8 16];
        case 'c'
            f_cutoff=[2500 4500];
            t_assay=0.04;
            spect_params=[.8 16];
        case 'd'
            % SECOND NOTE
            f_cutoff=[2000 3500];
            t_assay=.08;
            % FIRST NOTE
%             f_cutoff=[1000 2800];
%             t_assay=.016;
            spect_params=[.8 16];
        case 'f'
            f_cutoff=[1800 3200];
            t_assay=.025;
            spect_params=[.8 10];
        case 'g'
            f_cutoff=[3500 5000];
            t_assay=.025;
            spect_params=[.8 10];
        case 'h'
            f_cutoff=[2000 5000];
            t_assay=70;
            spect_params=[.8 16];
        case 'j'
            f_cutoff=[2000 4000];
            t_assay=30;
            spect_params=[.8 16];
        case 'k'
%             f_cutoff=[1300 3200];
            f_cutoff=[1000 2500];
            t_assay=50;
            spect_params=[.8 16];
        case 'r'
            f_cutoff=[2500 4500];
            t_assay=0.04;
            spect_params=[.8 16];
        case 'm'
            f_cutoff=[000 800];
            t_assay=50;
            spect_params=[.8 16];
    end
end

if strcmp(bname,'gy124wh97')
    switch(syl)
        case 'a'
            f_cutoff=[1000 2500];
            t_assay=50;
            spect_params=[0.9 16];
        case 'b'
            f_cutoff=[1000 2500];
            t_assay=50;
            spect_params=[.8 10];
        case 'c'
            f_cutoff=[1000 2500];
            t_assay=50;
            spect_params=[.8 10];
        case 'j'
            f_cutoff=[1700 2500];
            t_assay=.03;
            spect_params=[.8 10];
        case 'k'
            f_cutoff=[1000 1800];
            t_assay=.02;
            spect_params=[.8 10];
        case 'i'
            f_cutoff=[1000 2500];
            t_assay=50;
            spect_params=[.8 10];
        case 'u'
            f_cutoff=[2500 5000];
            t_assay=50;
            spect_params=[.8 10];
        case 'm'
            f_cutoff=[1000 2500];
            t_assay=50;
            spect_params=[.8 10];
    end
end

if ~exist('f_cutoff')    % if undefined
    warning(['function SYLLABLE_PARAMS_BY_BIRD did not find an if statement with bird ' bname ', this could cause calling functions to crash or not work correctly.']);
    f_cutoff='undefined';
    t_assay='undefined';
    spect_params='undefined';
end