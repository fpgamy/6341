% Project 2. 6341
% TD-DFT
% The purpose of the window in time-dependent Fourier transform is to
% limit the extent of the sequence to be transformed so that the spectral
% characteristics are approximately constant over the duration of the
% window. Based on the waterfall from gqrx the spectrum varies
% significantly only between 0.5 second windows. 

% N is the number of DFT samples (samples in the frequency dimension)
% L is the length of the TD-DFT window
% R is the amount of jumped samples between the signals (samples in the
% time dimension)

% Requirements:
% Need L <= N to guarantee reconstruction of the windowed segments from
% X[k]
% Need R < L so that no samples are missing from the windows

fid = fopen('gqrx_20181124_194330_940000000_2048000_fc.raw','rb');
hw_freq = 940e6;
IQData = fread(fid,'float32');
IQData = IQData(1:2:end) + 1i*IQData(2:2:end);
fclose(fid);

fs = 2.048e6;
L  = 0.125*fs;
R  = L/2;
N  = 2^ceil(log2(L));
pspectrum(abs(IQData), fs, 'spectrogram', 'FrequencyResolution', 4e3, ...
    'Overlap', 0, 'Leakage', 1);%, 'FrequencyLimits', [0.9e6 1e6]);%, FrequencyResolution', 125);

% pspectrum(abs(IQData), fs);