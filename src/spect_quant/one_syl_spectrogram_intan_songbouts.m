% Modified by Lyndie 9/2020 to use Intan files converted to mat
%Purpose:
%A. Fine-tune SYLLABLE_PARAMS_BY_BIRD.m settings to help you quantify 
%syllable pitch correctly.  Run this function before calling the
%latest version of the HEADPHONES_QUANTIFY_PITCH function.
%
%B. A good function to let you look at many iterations of the same
%syllable.

%Arguments:
%SYL_TO_QUANT:  which syllable to quantify.  You can optionally decide to
%  show spectrograms of a syllable that is within a larger motif, such as
%  the second syllable d in the sequence ddbc. Examples: 'a', 'aaaa', 'ccaadd'
%N_SYL_IN_SEQUENCE: optional argument.  Which syllable in SYL_TO_QUANT
%  to show spectrogram of. Example: if SYL_TO_QUANT='abcd', and
%  N_SYL_IN_SEQUENCE=3, it will show only spectrograms of syllable c and 
%  only when c occurs in the motif 'abcd'.
%PITCH_RANGE: optional argument.  minimum and maximum frequency in Hertz.  
%  Use this to look at only those syllable instances where pitch was quantified 
%  in this range. This only works if HEADPHONES_QUANTIFY_PITCH function was 
%  run previously and saved the PEAK_PINTERP_LABELVEC vector.  Otherwise it's 
%  meaningless since no pitch was quantified yet.  Example: [1000 2000]
%
%HOW TO RUN
%1. obs0 should be the bird's own song (if cbin).  See note below when CH_WITH_SONG is defined
%2. In Matlab, navigate to directory containing songfiles that you already labeled
%   using function EVSONGANALY.  The directory must contain both .cbin and
%   .not.mat files.
%3. Open function SYLLABLE_PARAMS_BY_BIRD and edit it to add your bird
%   and the syllables that you labeled.  This function will crash if the
%   bird or syllable is not found in SYLLABLE_PARAMS_BY_BIRD.
%4. Run this function.
%       Hit the [right arrow] key to go to the next set of syllables.
%       Function will exit when you have gone through all files in folder.
%       Hit the [q] key to quit.  You can only zoom in on figures once this
%       function exits, otherwise it just loops waiting for input.
%5. Tweak the SYLLABLE_PARAMS_BY_BIRD parameters
%6. Go to step 4 and see if your tweaks lead to improvement, until you feel
%   confident that the SYLLABLE_PARAMS_BY_BIRD parameters are good.

%OUTPUT FIGURE GUIDE
%Top plots show spectrograms of the desired syllable. Title is the time and 
%number of the songfile from which the syllable originated.  Horizontal
%lines are the minimum and maximum F_CUTOFF.  Function HEADPHONES_QUANTIFY_PITCH 
%will quantify all pitches between these two values, so make sure that the
%harmonic you choose is comfortably within this range.   Make sure there are
%no other harmonics in this range.  This could could lead to erroneous pitch 
%quantification if the second harmonic happens to be at a higher power than
%the one you originally had wanted to use.  Star shows the bin with maximum
%power at time T_ASSAY.
%
%Bottom plots show power at each frequency for the spectrogram column at time
%T_ASSAY in the frequency range F_CUTOFF.  Y axis units are not super important.
%X axis units are frequency in Hertz.  Title is the quantified pitch. If
%there is no quantified pitch (no PEAK_PINTERP_LABELVEC) it will say "0 Hz".  
%Vertical dotted line is the quantified pitch.  Notice how the top plot
%colors (at time T_ASSAY and frequencies between F_CUTOFF) correspond to 
%the bottom plot values.
%
%In the top plot, your goal is that T_ASSAY is always during the same
%part of the syllable.  "Same part" means that the star doesn't
%sometimes fall in the noisy intro to the syllable, and sometimes in the part
%that is a harmonic upsweep, and sometimes in the flat part.  You 
%always want the star to be in a place within the syllable that is
%sung with similar frequencies.
%
%If at all possible, pick a flat, harmonic part of the syllable with high power.
%
%For frequency-sweepy syllables, pick the beginning of the sweep.
%
%In the bottom plot, your goal is to get a sharp, well-defined, single peak 
%in the middle of the F_CUTOFF range.  You don't want to see multiple
%peaks, or peaks that are at the edge of F_CUTOFF range, or peaks that are
%shallow.  If you pick a strong harmonic region that always occurs at the same 
%time, and set your F_CUTOFF to neatly bracket that harmonic, you should be
%fine.
%
%WARNING
%Some syllables will be difficult to find good settings for.  If so, make a
%comment in function SYLLABLE_PARAMS_BY_BIRD and your lab notebook that the
%syllable was difficult and pitch quantification may not always be accurate.
%Remember this and take this into account when you do your analysis of that
%syllable later on.

function one_syl_spectrogram(syl_to_quant,n_syl_in_sequence,pitch_range)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Section 1: Initialization steps before plotting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Defaults: use first syllable in sequence, any pitch
if nargin==1;n_syl_in_sequence=1;pitch_range=[];end
if nargin==2;pitch_range=[];end

fid = getfilelist();fn = fgetl(fid);fclose(fid);
fid = getfilelist(); %not redundant since we want to load the first file again later on

% Discovering "birdname" - assumed to be the beginning of the filename
% string up to the first underscore - e.g. all files should be named
% "birdname_xxx.cbin"
underscore_id = strfind(fn,'_');birdname = fn(1:underscore_id-1);
if isempty(birdname)
   error('BIRDNAME not found.  One possible cause: you are not in a directory containing songfiles');
end

%CH_WITH_SONG FOR CBINS should contain evtaf_amp recording of bird's song
%channel or minimike channel.  Evtaf_amp setting "AI Channel" first entry
%is obs0, 2nd entry is obs1, ... and wires connect to these AI channels
%from the various audio signals. Not used for RHDs.
ch_with_song = 'obs0';

disp(['This is for bird ' birdname])
disp(['Plotting syllable ' syl_to_quant(n_syl_in_sequence) ' in sequence ' syl_to_quant]);
disp(['Channel with bird''s actual song is assumed to be ' ch_with_song]);

%Get the syllable parameters for bird BIRDNAME and syllable SYL_TO_QUANT.
%F_CUTOFF = minimum and maximum frequency between which to do pitch quantification
%T_ASSAY = see below
%SPECT_PARAMS = [(percent spectrogram bin overlap) (spectrogram bin size in milliseconds)]
[f_cutoff,t_assay,spect_params] = syllable_params_by_bird_leila(birdname,syl_to_quant(n_syl_in_sequence));
if strcmp(f_cutoff,'undefined') || strcmp(t_assay,'undefined') || strcmp(spect_params,'undefined')
   error(['Syllable ' syl_to_quant(n_syl_in_sequence) ' does not exist or does not have all ',...
       'three parameters F_CUTOFF, T_ASSAY, SPECT_PARAMS defined in function SYLLABLE_PARAMS_BY_BIRD']);
end

if t_assay>1 %if T_ASSAY>1, then it is the percentage (1-100) of the time through the syllable at which to quantify
    use_pct = 1;
    t_pct = t_assay*.01;
    disp(['Assaying at t = ' num2str(t_assay) '% through syllable with frequency  range [' num2str(f_cutoff) '] Hz'])
else %if T_ASSAY<=1, then it is the number of seconds since syllable start at which to quantify.
    use_pct = 0;
    disp(['Assaying at t = ' num2str(t_assay) 'ms with frequency range [' num2str(f_cutoff) '] Hz'])
end

n_syls = 8; %how many syllables to plot
ct_plotted = 0; %number of syllables plotted on figure

%make fullscreen figure and set up two textbox messages to the user
f = figure; screen_size = get(0, 'ScreenSize'); set(f, 'Position', [0 0 screen_size(3) screen_size(4)]);
textbox1 = annotation('textbox',[0 0 1 1],'String',...
['Right arrow key = next ' num2str(n_syls) ' syllables, q = quit function.'],...
'fontsize',15,'HorizontalAlignment','center','VerticalAlignment','middle');
set(textbox1,'Visible','on'); %show now
textbox2 = annotation('textbox',[0 0 1 1],'String',... %later, lets user know function is done
'one\_syl\_spectrogram exited','fontsize',20,'fontweight','bold',...
'HorizontalAlignment','center','VerticalAlignment','middle');
set(textbox2,'Visible','off'); %only show when this function exits
ax = []; ax2 = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Section 2: Make plots and let user page through syllables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ct = 1;
while 1 %open each file in turn.  Loop will exit using "return" statements in two places below
    %load raw song
    fn = fgetl(fid);
    cbin_fn = fn(1:end-8);
    if (~ischar(fn)) %end of batchfile, no more song files, exit function
        cleanup(fid, textbox1, textbox2, f, ax, ax2)
        return;
    end
    clear peak_pinterp_labelvec
    disp(fn)
    eval(sprintf('load %s',fn))
    if strcmp(cbin_fn(end-3:end),'cbin')
        [rawsong, Fs] = evsoundin('.', cbin_fn,ch_with_song); %load song from current .cbin file
    elseif strcmp(cbin_fn(end-2:end),'rhd')
        load([cbin_fn(1:end-3),'mat'])
        rawsong = board_adc_data;
        Fs = frequency_parameters.amplifier_sample_rate;
    elseif strcmp(cbin_fn(end-2:end),'mat') && length(char(cbin_fn)) < 39
        load([cbin_fn(1:end-3),'mat'])
        rawsong = board_adc_data;
        Fs = frequency_parameters.amplifier_sample_rate;
    elseif strcmp(cbin_fn(end-2:end),'wav') %if file is .wav
        [rawsong] = audioread(cbin_fn);
        Fs = 30000;
        rawsong = resample(rawsong,30000,Fs);
    else ct = ct + 1;
    end

    %get onsets, offsets, sample times, and ids of matching syllables or syllable sequences
    onsets = onsets/1000; %get syllable onsets and offsets into sec, not msec
    offsets = offsets/1000;
    t_samples = 1:length(rawsong); 
    t = t_samples/Fs; %get time in seconds for each value in rawsong
    id = strfind(labels,syl_to_quant); %find either single syllable 'a' or motif 'aabbc'
    id = id + n_syl_in_sequence - 1; %change ID to match # syllable in sequence 'aabbcc', if single syllable 'a' ID stays the same.
    disp(['Found ' num2str(length(id)) ' examples']);

    %vector of pitches of the target syllable, if pitch was previously
    %quantified using headphones_quantify_pitch.
    if ~exist('peak_pinterp_labelvec','var')
        peak_pinterp_labelvec = zeros(1,length(labels));
    end
    if numel(pitch_range)>0 %find ids of all pitches of target syllable that are within the set frequency range
        id = intersect(find(peak_pinterp_labelvec>pitch_range(1) & peak_pinterp_labelvec<pitch_range(2)),id);
    end

    for x=1:length(id)    %Make two plots for each syllable  
        %Keypress code block
        if ct_plotted==n_syls %if all syllables are plotted already on current plot
            set(f,'Color', [0.8 0.8 0.8]);
            %wait for user to hit a key after current plots are finished
            %right arrow key = go to next syllables. q = quit
            key = '';
            while ~strcmp(key,'rightarrow') && ~strcmp(key,'q')
                i = 0;
                while i==0, i = waitforbuttonpress;  end
                key = lower(get(gcf, 'CurrentKey'));
            end
            if strcmp(key,'q') %exit function
                cleanup(fid, textbox1, textbox2, f, ax, ax2);
                return;
            end
            set(f,'Color',[0.7 0.7 0.7]); %make figure background darker so user knows it is processing next syllables
            drawnow %make darker figure background visible immediately
            for i=1:n_syls*2 %clear everything on all subplots before plotting next syllables
                subplot(2,n_syls,i);
                cla; title(''); xlabel('');
            end
            ax = []; ax2 = [];
            ct_plotted=0; %reset CT_PLOTTED to plot the next set of syllables
        end
        
        ct_plotted = ct_plotted+1; %plot the next syllable

        %Make top plot
        ax(ct_plotted) = subplot(2,n_syls,ct_plotted); 
        hold on
        ylim([100 12000]);
        on = onsets(id(x));
        off = offsets(id(x));
        if use_pct
            t_assay = t_pct*(off-on);
        end
        if (off-on)<spect_params(2)/1000    %if syl is too short
            off = on+spect_params(2)/1000;    % extend end of syl
            disp('Extending offset')
        end
        on_id = find(abs((t-on))==min(abs(t-on))); %which time in RAWSONG is closest to syllable onset?
        off_id = find(abs((t-off))==min(abs(t-off)));%...to syllable offset?
        syl_wav = rawsong(on_id:off_id); %raw waveform for current syllable
        [~,F1,T1,P1] = spect_from_waveform(syl_wav,Fs,0,spect_params); %get syllable spectrogram
        imagesc(T1,F1,log(P1)); set(gca,'YD','n');
        if ct_plotted==1, %only make axis labels in one subplot for cleaner-looking graphs
            xlabel('Time (Seconds)');
            ylabel('Frequency (Hz)'); 
        end
        f_cut_id = find(F1>f_cutoff(1) & F1<f_cutoff(2));
        F1 = F1(f_cut_id); %Frequencies at which power was computed.
        P1 = P1(f_cut_id,:); %P1 rows are power at various frequencies F1.
                           %P1 columns are different times in syllable
                           %at millisecond intervals specified in SPECT_PARAMS(2)
        %the frequency column whose time (since syllable start) is closest to T_ASSAY
        t_id = find(abs((T1-t_assay))==min(abs(T1-t_assay)));
        if length(t_id)>1
            disp('Two equidistant time windows - choosing earliest one')
            t_id = min(t_id);
        end
        spect_slice = P1(:,t_id); %select one frequency column at a time slice containing T_ASSAY
        max_p_id = spect_slice==max(spect_slice); %frequency bin with the highest power
        plot(t_assay,F1(max_p_id),'k*') %plot starred point at the highest power for this syllable at time T_ASSAY
        xl = get(gca,'xlim');
        plot(xl,[f_cutoff(1) f_cutoff(1)],'k'); %f_cutoff horizontal line
        plot(xl,[f_cutoff(2) f_cutoff(2)],'k'); %f_cutoff horizontal line
        ida = find(cbin_fn=='_', 1, 'last' )-6;    %last character '_' +1
        idb = find(cbin_fn=='.', 1, 'last' )-1;    %last character '.' -1
        title(cbin_fn(ida:idb)); %shows only the time and file ID # to make a short readable title

        %Make bottom plot
        ax2(ct_plotted) = subplot(2,n_syls,ct_plotted+n_syls);
        hold on
        xlim(f_cutoff);
        semilogy(F1,spect_slice) %plot frequency power in the selected column as line
        semilogy(F1,spect_slice,'.') %plot frequency power as dots
        yl = get(gca,'ylim');
        val = peak_pinterp_labelvec(id(x)); %pitch for current syllable, if pitch was previously quantified by a version of HEADPHONES_QUANTIFY_PITCH function
        plot([val val],yl,'k:') %dotted line
        title([num2str(val) ' Hz']);
        xlab{1} = ['t=' num2str(round(onsets(id(x))*100)/100) 'sec']; %time in song that syllable occurred
        xlab{2} = ['Context:' labels(max([1 id(x)-5]):id(x)-1) ' ' labels(id(x)) ' '  labels(id(x)+1:min([length(labels) id(x)+5]))]; %what were surrounding syllables
        xlabel(xlab,'fontweight','bold')
    end
end
end

%Utility function: make list of all .not.mat files in current directory
function fid = getfilelist
if isunix, 
    !ls *.not.mat > batchfile
else
    !dir /B *.not.mat > batchfile
end
fid = fopen('batchfile','r');
end

%Utility function: clean up before returning
function cleanup(fid, textbox1, textbox2, f, ax, ax2)
    fclose(fid);
    set(textbox1,'Visible','off'); 
    set(textbox2,'Visible','on'); 
    set(f,'Color', [0.8 0.8 0.8]); 
    try %try-catch so no error message if nothing was plotted
        linkaxes(ax,'xy'); %link axes so zoom-in will zoom in on all plots
        linkaxes(ax2,'x');
    catch e
    end
end