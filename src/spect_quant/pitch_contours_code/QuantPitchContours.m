function QuantPitchContours(syl,NFFT,overlap,sigma)
% QuantPitchContours
% estimates contour of pitch for syllables, then adds a struct variable
% containing the pitch contours and other information (see below) to each
% .not.mat file in a directory. This function uses the frequency cutoffs
% defined in the syllable_params_by_bird function as the "f_cutoff" 
% variable, so that variable must be defined for the bird ID in the filenames
% of the directory you're in when you run this function.
%
% syntax:
% QuantPitchContours(syl,NFFT,overlap,sigma)
%
% input arguments:
%   syl -- syllable for which function should find pitch contours, e.g., 'h'
%   NFFT -- number of samples from audio file to use for Fast Fourier
%       Transform(FFT)
%   overlap -- number of samples to overlap for successive FFTs on song
%       file
%   sigma -- parameter used in exponential function to generate Gaussian
%       window used for spectrogram
%
% output: struct "pitchContours" with the following fields
%     contour -- n x 1 vector containing 
%     timeVector -- n x 1 vector of time bins returned from FFT
%     NFFT, overlap, sigma -- same as above, saved in struct in case
%       they're needed for analysis
%
% usage example:
%   QuantPitchContours('a',1024,1020,1.0)
%
% last edited: Tomek Jorke, 02.23.15

% first get all .not.mat files into a variable, to loop through
% loading each .not.mat file and associated data from .cbin
notmats = ls('*.not.mat');

spbb_fname='syllable_params_by_bird_leila'; % Specify the name of your syllable_params_by_bird*.m file

% Discovering "birdname" - assumed to be the beginning of the filename
% string up to the first underscore - e.g. all files should be named
% "birdname_xxx.cbin"
fn=notmats(1,:);
birdname=fn(1:strfind(fn,'_')-1);

% get parameters required out of syllable_params_by_bird
eval(sprintf('[f_cutoff,~,spect_params]=%s(birdname,syl);',spbb_fname))
lowFreqCutoff = f_cutoff(1); % f_cutoff is variable loaded from
highFreqCutoff = f_cutoff(2); % syllable_params_by_bird

% loadin' loop
for i = 1:size(notmats,1)

    % open .not.mat
    load(notmats(i,:),'-mat','labels','onsets','offsets','pitchContours')
        % '-mat' argument forces open as .mat file
        % for some reason it kept trying to open as ASCII file
    warning off % suppress warnings if vars are not found in file because
                % this is expected the first time QuantPitchContours is run on a .not.mat file.
    
    % have to clear pitchContours if number of labels has changed, e.g.,
    % because you changed threshold and thus number of syllables to label
    lblsize = numel(labels);
    if exist('pitchContours','var')
        if isstruct('pitchContours') && numel(pitchContours) ~= lblsize
            clear pitchContours
        elseif ~isstruct('pitchContours')
            clear pitchContours
        end
    end
    
        % find ids of targeted syl
    ids = strfind(labels,syl);
    disp(['Found ' num2str(length(ids)) ' examples of ''' syl ''' in ' notmats(i,:)])
    
    % open .cbin, load raw song data
    %   first get filename from .not.mat filename
    notmat_strid = strfind(notmats(i,:),'.not.mat');
    cbin = notmats(i, 1:notmat_strid - 1);
	[rawsong,fs] = evsoundin('',cbin,'obs0'); % fs = sampling frequency
    
    % Dialog box that halts the program if it finds that Fs ~= 32k, and asks
    % you whether you want to resample to that frequency (by clicking 'ok')
    % or quit. Happens if Evtaf_amp sampling frequency changes (Actual AI Rate)
    if fs ~= 32000 && ~exist('resample_song','var')
        old_fs = fs;
        new_fs = 32000;
        [p,q] = rat(new_fs/old_fs,0.0001);
        questntxt = {['Sampling frequency (Fs) for ' cbin_fn 'equals ' num2str(fs)...
            ' instead of ' num2str(new_fs)],...
            ['click ''OK'' to resample with ratio ' num2str(p) ' / '  num2str(q)],...
            ['otherwise click ''Quit'' and copy files with this Fs'],...
            ['to a different folder to prevent crashes']};
        btn = questdlg(questntxt,'SAMPLING ERROR','OK','Quit','OK');
        switch btn
            case 'OK'
                resample_song = 1;
                rawsong = resample(rawsong,p,q);
                if quantify_shifted_song;shiftedsong = resample(shiftedsong,p,q);end
                if quantify_mini_mic;mini_mic_song = resample(mini_mic_song,p,q);end
                fs = new_fs;
            case 'Quit'
                error('non-matching sampling frequencies')
        end
    elseif fs ~= 32000 && resample_song == 1
        rawsong = resample(rawsong,p,q);
        fs = new_fs;
    end   
    
    %now that we have fs from current file, create window
    t = -NFFT/2+1 : NFFT/2;
    sigma = (sigma/1000) * fs;
    % Gaussian and first derivative as windows.
    window = exp(-(t/sigma).^2);

    % get onsets and offsets
    sylOnsets = onsets(ids);
    sylOffsets = offsets(ids);
    % convert on- and offsets to time in units of samples
    % (have to first convert from s to ms)
    sylOnsets = floor((sylOnsets / 1000) * fs);  
    sylOffsets = ceil((sylOffsets / 1000) * fs);
    
    % initialize variable that holds output from GetPitchContours function
    if ~exist('pitchContours','var')
        pitchContours = struct('contour',[],...
            'timeVector',[],...
            'NFFT',[],...
            'overlap',[],...
            'sigma',[]);
    end
    
    for j = 1:length(ids)
        sylFromRawsong = rawsong(sylOnsets(j):sylOffsets(j));
       
        try 
        [~,freqsVector,timeVector,powerVector] = spectrogram(sylFromRawsong,window,overlap,NFFT,fs);
        catch err
            if strcmp(err.identifier,'signal:welchparse:invalidSegmentLength')
                disp(['Rendition ' num2str(currRend) ' could not be converted to spectrogram'])
                pitchContours(currRend,1) = {Nan};
            else % if "catch" caught something besides "welchparse" error
                rethrow(err) % throw that error instead
            end
        end

        %This loop determines the pitch (fundamental frequency) at each time bin by
        %taking the peak of the autocorrelation function within a range specified
        %by the user
        freqBinIds = find(freqsVector < highFreqCutoff & freqsVector > lowFreqCutoff);
        freqsForPitchquant = freqsVector(freqBinIds);
        powerForPitchquant = powerVector(freqBinIds,:);

        pitch = zeros(length(timeVector),1);

        for currentTimeBin = 1:length(timeVector)
            spectSlice=powerForPitchquant(:,currentTimeBin);
            pitch(currentTimeBin,1) = ComputeScalarPitch(spectSlice,freqsForPitchquant);
        end
        pitchContours(ids(j)).contour = pitch;
        pitchContours(ids(j)).timeVector = timeVector(:); % (:) <- makes timeVector a column vector
        pitchContours(ids(j)).NFFT = NFFT;
        pitchContours(ids(j)).overlap = overlap;
        pitchContours(ids(j)).sigma = sigma;
    end % of j loop to go through ids vector
    save(notmats(i,:),'pitchContours','-append')
    
    clear labels onsets offsets pitchContours
    
end % of i loop to go through .not.mat filename 

function [peak_pinterp] = ComputeScalarPitch(spectSlice,freqsForPitchquant)
% Compute a single scalar value of pitch based on spect_slice

% find maximum of power spectrum
maxPowerId = find(spectSlice==max(spectSlice));
% if maximum is at first or last index of spectrum
if maxPowerId==1 || maxPowerId==length(spectSlice)
    % save that maximum for now
    oldMaxPowerId = maxPowerId;
    % and then look for a convex peak in between limits
    % by taking the 2nd derivative of the power spectrum
    spect2ndDerivative = sign(diff(sign(diff(spectSlice))));
    % then find the index of negative values. By the second derivative
    % test, these will be a local maxima
    neg2ndDerivId = find(spect2ndDerivative<0) + 1; % have to add 1
    % because the vector shrinks by one every time you take the derivative
    % (think of it as losing the first and last values of the spectrum --
    % you can ignore that the last value isn't there because it doesn't
    % change the value of the indices in the same way that not having the
    % first value does)
    
    % now take whichever one of the local maxima has the highest power
    maxPowerId = neg2ndDerivId(find(spectSlice(neg2ndDerivId)==max(spectSlice(neg2ndDerivId))));
end

% if NO convex peak found above, use old edge peak
if isempty(maxPowerId)
    peak=freqsForPitchquant(oldMaxPowerId);
    peak_pinterp = peak;
else
    % use convex peak
    peak = freqsForPitchquant(maxPowerId);
    x_pinterp = freqsForPitchquant(maxPowerId-1:maxPowerId+1);
    y_pinterp = spectSlice(maxPowerId-1:maxPowerId+1);
    [peak_pinterp] = pinterp_internal(x_pinterp, y_pinterp); % parabolic interpolation
end

function xmax = pinterp_internal(xs, ys)
% parabolically interpolates the peak given three points
x = [xs.^2, xs, [1;1;1]];
xInv = inv(x);
abc = xInv * ys;
xmax = -abc(2)/(2*abc(1));