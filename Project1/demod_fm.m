fid = fopen('gqrx_20181016_105047_98500200_2048000_fc.raw','rb');
central_freq = 98500200;
IQData = fread(fid,'float32');
IQData = IQData(1:2:end) + 1i*IQData(2:2:end);
sampleRate = 2.048e6;

plot_spectrum(IQData, sampleRate, central_freq);
title('$$|X(j\Omega)|$$','interpreter','latex');

sample_no = transpose(1:1:length(IQData));

% plot_spectrum(IQData, sampleRate)
% AA Filter
Fpass = 125000;      % Passband Frequency
Fstop = 128000;      % Stopband Frequency
Apass = 0.05;        % Passband Ripple (dB)
Astop = 50;          % Stopband Attenuation (dB)
match = 'passband';  % Band to match exactly

% Construct an FDESIGN object and call its CHEBY1 method.
h  = fdesign.lowpass(Fpass, Fstop, Apass, Astop, sampleRate);
Hd = design(h, 'cheby1', 'MatchExactly', match);
fvtool(Hd);

aa_IQData = filter(Hd, IQData);
plot_spectrum(aa_IQData, sampleRate, central_freq);
title('Anti-Aliased $$|X(j\Omega)|$$','interpreter','latex');
% clear IQData
downsampled_IQData = downsample(aa_IQData, 8);
sampleRate = sampleRate / 8;

plot_spectrum(downsampled_IQData, sampleRate, central_freq);
title('$$|X_d(j\Omega)|$$','interpreter','latex');

% limiter 
downsampled_IQData = downsampled_IQData./abs(downsampled_IQData);

plot_spectrum(downsampled_IQData, sampleRate, central_freq);
title('Amplitude Limited $$|X_{d}(j\Omega)|$$','interpreter','latex');

% differentiator
N = 100;                       % Order
F = [0 0.95 0.98 1];  % Frequency Vector
A = [0 1 0 0];                 % Amplitude Vector
W = [1 1];                     % Weight Vector

% Calculate the coefficients using the FIRPM function.
b  = firls(N, F, A, W, 'differentiator');
Hd = dfilt.dffir(b);
fvtool(Hd);

[gd, ] = grpdelay(Hd);
gd = mean(gd);
% delayed_IQData = transpose([zeros(1, gd) transpose(downsampled_IQData)]);
delayed_IQData = downsampled_IQData(1:end-gd);

plot_spectrum(delayed_IQData, sampleRate, central_freq);
title('$$|Y_{delay}(j\Omega)|$$','interpreter','latex');

% differentiated_IQData = transpose([ transpose(filter(Hd, downsampled_IQData)) zeros(1, gd) ]);
differentiated_IQData = filter(Hd, downsampled_IQData);
differentiated_IQData(1:gd) = [];
plot_spectrum(differentiated_IQData, sampleRate, central_freq);
title('$$|Y_{differentiated}(j\Omega)|$$','interpreter','latex');

mult_IQData = differentiated_IQData.*conj(delayed_IQData);

plot_spectrum(mult_IQData, sampleRate, central_freq);
title('$$|Y_{2}(j\Omega)|$$','interpreter','latex');

% clear downsampled_IQData
% clear delayed_IQData
% clear differentiated_IQData
discrim_out = imag(mult_IQData);
% clear mult_IQData
plot_spectrum(discrim_out, sampleRate, central_freq);
title('$$|M(j\Omega)|$$','interpreter','latex');

% discriminator
tau = 75e-6;
T = 1/sampleRate;
alpha = 1/(tan(T/(2*tau)));
% fvtool(numd, dend);
deemphasis_out = filter([1 1], [1+alpha (1-alpha)], discrim_out);
% clear discrim_out

plot_spectrum(deemphasis_out, sampleRate, central_freq);
title('$$|M_{d}(j\Omega)|$$','interpreter','latex');

%low pass filter around 15kHz

Fpass = 15;          % Passband Frequency
Fstop = 18;          % Stopband Frequency
Apass = 1;           % Passband Ripple (dB)
Astop = 100;         % Stopband Attenuation (dB)
match = 'stopband';  % Band to match exactly

% Construct an FDESIGN object and call its CHEBY2 method.
h  = fdesign.lowpass(Fpass, Fstop, Apass, Astop, sampleRate/1000);
Hd = design(h, 'cheby2', 'MatchExactly', match);

% fvtool(Hd);

output_256 = filter(Hd, deemphasis_out);

plot_spectrum(deemphasis_out, sampleRate, central_freq);
title('$$|A_{m}(j\Omega)|$$','interpreter','latex');

% clear demphasis_out
audio_out = downsample(output_256, 4);

% clear output_256
sampleRate = sampleRate / 4;

plot_spectrum(audio_out, sampleRate, central_freq);
title('$$|A(j\Omega)|$$','interpreter','latex');

soundsc(audio_out, sampleRate);