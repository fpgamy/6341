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

fid = fopen('gqrx_20181123_033424_939788000_2048000_fc.raw','rb');
hw_freq = 939788000;
IQData = fread(fid,'float32');
IQData = IQData(1:2:end) + 1i*IQData(2:2:end);
fclose(fid);

fs = 2.048e6;
% signal is bandlimited from 0 to 1.024e6
% freqeuncy resolution ~= 1.024e6/L * 6db mainlobe width = 4k
% for blackman window mainlobe is 12*pi/M
% 4kHz/2.048MHz * 2*pi = 12*pi/(L-1)
% L = 6 * 2.048Mhz / 4kHz = 3073
% for sake of efficiency we shall choose N to be a power of 2 = 2048
% for L = 4096, time resolution = 2048 / fs = 1ms
N  = 4096;
L  = N;
R  = L/2;

window = blackman(L);
TDDFT = td_dft(transpose(IQData), R, L, N, window);
periodogram_db = 20.*log10(abs(TDDFT)); 
X = (hw_freq/1e6) + (((-fs/2):1:(fs/2))/(1e6));
Y_labels = datenum(seconds((1:(ceil(length(IQData)/R)))*R/fs) + datetime(2018, 11, 23, 03, 34, 24));

% Create figure
figure1 = figure;

% Create axes
axes1 = axes('Parent',figure1);
hold(axes1,'on');

imagesc(X, Y_labels, periodogram_db, 'Parent',axes1);
colorbar

xlim(axes1,[(hw_freq - fs/2)/1e6  (hw_freq + fs/2)/1e6]);

set(gca, 'YDir', 'normal');
xlabel('Frequency (MHz)');
datetick('y', 'mmm.dd,yyyy HH:MM:SS');
ylim([Y_labels(1), Y_labels(end)])

figure;
plot(20*log10(abs(fftshift(fft(window)))));

figure;
periodograms_single_time_inst = abs(TDDFT).*abs(TDDFT);
unbias_sf = 1/(sum(abs(window).*abs(window)));
periodograms_single_time_inst = periodograms_single_time_inst.*unbias_sf;
unbiased_average_periodogram = mean(periodograms_single_time_inst);
hold on;
plot(10*log10(unbiased_average_periodogram));

flat_periodogram = reshape(abs(TDDFT).*abs(TDDFT), 1, []);
noise_cutoff_index = ceil(size(TDDFT, 1) * size(TDDFT, 2) / 5);
td_dft_sorted = sort(flat_periodogram, 'ascend');
noise_cutoff = unbias_sf * mean(td_dft_sorted(1:noise_cutoff_index));
plot([0 N], [10*log10(noise_cutoff) 10*log10(noise_cutoff)]);
hold off

noise_cutoff_index = ceil(size(TDDFT, 1) * size(TDDFT, 2) / 1.1);
td_dft_sorted = sort(flat_periodogram, 'ascend');
noise_cutoff = unbias_sf * mean(td_dft_sorted(noise_cutoff_index));

p = unbias_sf .* abs(TDDFT) .* abs(TDDFT);
average_noise_periodogram = zeros(1, N);
number_of_samples = zeros(1, N);
for freq_bin = 1:1:N
    freq_col = p(:,freq_bin);
    number_of_samples(freq_bin) = sum(freq_col < noise_cutoff);
    average_noise_periodogram(freq_bin) = mean(freq_col(freq_col < noise_cutoff));
end

figure;
plot(10*log10(average_noise_periodogram));

% Bayesian average
bayesian_average_noise_periodogram = zeros(1, N);
prior_periodogram = average_noise_periodogram;
prior_periodogram(isnan(prior_periodogram)) = 0;

PRIOR_W = N/32;
for i = 1:10
    for freq_bin = 1:1:N
        if (freq_bin < ((PRIOR_W+1)/2) || freq_bin > (N-((PRIOR_W+1)/2)))
            C = number_of_samples(freq_bin);
            prior = prior_periodogram(freq_bin);
        else
            C = mean(number_of_samples(freq_bin-floor(PRIOR_W/2):freq_bin+(ceil((PRIOR_W+1)/2))));
            prior = mean(prior_periodogram(freq_bin-floor(PRIOR_W/2):freq_bin+(ceil((PRIOR_W+1)/2))));
        end

        freq_col = p(:,freq_bin);
        noise_samples = sum(freq_col < noise_cutoff);
        bayesian_average_noise_periodogram(freq_bin) = (sum(freq_col(freq_col < noise_cutoff))+C*prior)/(noise_samples + C);
    end
    prior_periodogram = bayesian_average_noise_periodogram;
end

figure;
hold on;
plot(10*log10(bayesian_average_noise_periodogram));
plot(10*log10(unbiased_average_periodogram));
hold off;

one_zero_periodogram = (abs(TDDFT).*abs(TDDFT));
eta = 6/unbias_sf;
for i = 1:N
    fft_for_bin_i = movmean(one_zero_periodogram(:,i), 4);
    one_locations = (find(fft_for_bin_i >  (eta*(bayesian_average_noise_periodogram(i))) ) );
    ones_p = zeros(1, length(fft_for_bin_i));
    ones_p(one_locations) = 1;
    one_zero_periodogram(:,i) = ones_p;
end

figure;
imagesc(one_zero_periodogram);
