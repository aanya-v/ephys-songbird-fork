%OneSylSpecWithPitchContours
%Show spectrograms of renditions of SYL_TO_QUANT with the estimated pitch
%contour plotted over the spectrogram, to get a sense of how well the
%algorithm for pitch contours is estimating the pitch
%
%
%Arguments:
%SYL_TO_QUANT:  which syllable to quantify.  You can optionally decide to
%  show spectrograms of a syllable that is within a larger motif, such as
%  the second syllable d in the sequence ddbc. Examples: 'a', 'aaaa', 'ccaadd'
%N_SYL_IN_SEQUENCE: optional argument.  Which syllable in SYL_TO_QUANT
%  to show spectrogram of. Example: if SYL_TO_QUANT='abcd', and
%  N_SYL_IN_SEQUENCE=3, it will show only spectrograms of syllable c and 
%  only when c occurs in the motif 'abcd'.


function OneSylSpecWithPitchContours(syl_to_quant,n_syl_in_sequence)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Section 1: Initialization steps before plotting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Defaults: use first syllable in sequence, any pitch
if nargin==1;n_syl_in_sequence=1;end

fid = getfilelist();fn = fgetl(fid);fclose(fid);
fid = getfilelist(); %not redundant since we want to load the first file again later on

% Discovering "birdname" - assumed to be the beginning of the filename
% string up to the first underscore - e.g. all files should be named
% "birdname_xxx.cbin"
underscore_id = strfind(fn,'_');birdname = fn(1:underscore_id-1);
if isempty(birdname)
   error('BIRDNAME not found.  One possible cause: you are not in a directory containing songfiles');
end

%CH_WITH_SONG should contain evtaf_amp recording of bird's song
%channel or minimike channel.  Evtaf_amp setting "AI Channel" first entry
%is obs0, 2nd entry is obs1, ... and wires connect to these AI channels
%from the various audio signals.
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

n_syls = 4; %how many syllables to plot
ct_plotted = 0; %number of syllables plotted on figure

%make fullscreen figure and set up two textbox messages to the user
f = figure; screen_size = get(0, 'ScreenSize'); set(f, 'Position', [0 0 screen_size(3) screen_size(4)]);
textbox1 = annotation('textbox',[0 0 1 1],'String',...
['Right arrow key = next ' num2str(n_syls) ' syllables, q = quit function.'],...
'fontsize',15,'HorizontalAlignment','center','VerticalAlignment','cap');
set(textbox1,'Visible','on'); %show now
textbox2 = annotation('textbox',[0 0 1 1],'String',... %later, lets user know function is done
'one\_syl\_spectrogram exited','fontsize',20,'fontweight','bold',...
'HorizontalAlignment','center','VerticalAlignment','middle');
set(textbox2,'Visible','off'); %only show when this function exits
ax = []; ax2 = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Section 2: Make plots and let user page through syllables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
    else %if file is .wav
        [rawsong,Fs,~] = wavread(cbin_fn);
        rawsong = resample(rawsong,32000,Fs);
        Fs = 32000;
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

        %Makeplot
        ax(ct_plotted) = subplot(1,n_syls,ct_plotted); 
        hold on
        ylim([100 10000]);
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
        if exist('pitchContours','var')
            if ~isempty(pitchContours(id(x)).contour)
                plot(pitchContours(id(x)).timeVector,pitchContours(id(x)).contour,'k.','linewidth',2)
            end
        end
        ida = find(cbin_fn=='_', 1, 'last' )+1;    %last character '_' +1
        idb = find(cbin_fn=='.', 1, 'last' )-1;    %last character '.' -1
        title(cbin_fn(ida:idb)); %shows only the time and file ID # to make a short readable title
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