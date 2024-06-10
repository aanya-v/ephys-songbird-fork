%syllable_specgram should be your p (ex. p1)
%syllable_specgrams_tmaxsize is the size of your new data (ex. 64 will
%produce a 64 by 64 spectogram
%plot can be 1 if you want to plot the spectograms and 0 if you just want
%the new data

%example usage
%newerdata = intrapolate(p1, 64, 1)


function spec_interpolated = spec_interpolate (syllable_specgram, syllable_specgrams_tmaxsize, plot)
  % state original dimensions based on syllable_specgram
  %disp(['Before: ' num2str(size(syllable_specgram))]);

  % Extract time and frequency dimensions
  specgram_tlen = size(syllable_specgram, 2);
  syllable_specgrams_fsize = size(syllable_specgram, 1);

  % Create zero-filled interpolated spectrogram with desired size
  syllable_specgram_interp = zeros(syllable_specgrams_tmaxsize, syllable_specgrams_tmaxsize);

  % Define arbitrary duration
  arbitrary_duration = 1.0;

  % Create old time and frequency grids
  Told = linspace(0, arbitrary_duration, specgram_tlen);
  Fold = linspace(0, arbitrary_duration, syllable_specgrams_fsize);

  % Create new time and frequency grids for interpolation
  xx = linspace(0, arbitrary_duration, syllable_specgrams_tmaxsize);
  yy = linspace(0, arbitrary_duration, syllable_specgrams_tmaxsize);
  
 [meshX,meshY] = meshgrid(xx,yy);
 gridPoints = {Fold, Told};
 interp = griddedInterpolant(gridPoints, syllable_specgram, 'nearest');
 newdata = interp(meshX,meshY);
 spec_interpolated = rot90(newdata);

  if plot
      clims = [0.00000 0.000001];
      colormap parula;
      figure; imagesc(spec_interpolated, clims); %plots new data
      title('New Spectogram');

      colorbar;
      figure; imagesc(syllable_specgram,clims) %plots original data
      axis xy
      colorbar;
      title('Old Spectogram');

  end
end

