function [reshapedMatrix] = createSpectMatrix(spect_uncut,selected_bins,bin_length_setting)
% Open compiled spectrograms and pad columns with zero columns to the
% desired window size 'maxColumns'

n_trials = length(spect_uncut);

if bin_length_setting == "fixed"
    % will transform spec into a fixed-size n x n spectrogram using interpolation
    n_time_bins = selected_bins;
    freq_binlength = selected_bins;
else
    n_time_bins = diff(selected_bins)+1;
    freq_binlength = length(spect_uncut{1});

    % Pad each cell with zeros to have equal number of columns
    for i = 1:n_trials
        currentColumns = size(spect_uncut{i}, 2);
        if currentColumns < n_time_bins
            % Padding with zeros
            spect_uncut{i} = [spect_uncut{i}, zeros(freq_binlength, n_time_bins - currentColumns)];
        end
    
    end
end

% Initialize an empty matrix to hold the reshaped data
reshapedMatrix = zeros(n_trials, freq_binlength * n_time_bins);

% Counter for the row index in reshapedMatrix
rowIndex = 1;

% Loop through each set and each cell to generate each trial's spectrogram
for i = 1:n_trials
    curr_spect = spect_uncut{i};

    if bin_length_setting == "fixed"
        plot = 0;
        spect_select = spec_interpolate(curr_spect,selected_bins,plot);
    else
        % Select time bins for current spectrogram
        spect_select = curr_spect(:,selected_bins(1):selected_bins(2));
    end
        
    
    % Reshape the current cell into a row vector
    currentRow = reshape(spect_select, 1, []);

    % Assign it to the reshapedMatrix
    reshapedMatrix(rowIndex, :) = currentRow;
    
    % Update the row index
    rowIndex = rowIndex + 1;
end
end