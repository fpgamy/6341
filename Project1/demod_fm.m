fid = fopen('gqrx_20181016_105150_106700200_2048000_fc.raw','rb');
IQData = fread(fid,'float32');
IQData = IQData(1:2:end) + 1i*IQData(2:2:end);
sampleRate = 2.048*10^6;

% plot_spectrum(IQData, sampleRate)

sample_no = transpose(1:1:length(IQData));

IQData = IQData .* exp(-1i.*sample_no.*2*pi*(106.7/2048));
% plot_spectrum(IQData, sampleRate)
% AA Filter
Fpass = 0.05625;          % Passband Frequency
Fstop = 0.0625;           % Stopband Frequency
Dpass = 0.0057563991496;  % Passband Ripple
Dstop = 0.0001;           % Stopband Attenuation
dens  = 20;               % Density Factor

% Calculate the order from the parameters using FIRPMORD.
[N, Fo, Ao, W] = firpmord([Fpass, Fstop], [1 0], [Dpass, Dstop]);

% Calculate the coefficients using the FIRPM function.
b  = firpm(N, Fo, Ao, W, {dens});
Hd = dfilt.dffir(b);
% fvtool(Hd);

aa_IQData = filter(Hd, IQData);
% clear IQData
downsampled_IQData = downsample(aa_IQData, 8);
sampleRate = sampleRate / 8;

% figure;
% plot(abs(real(downsampled_IQData) + 1i*(imag(downsampled_IQData))));

% limiter 
mask = ones(size(downsampled_IQData));
mask = mask./abs(downsampled_IQData);

% dRL = limiter(-10, 'SampleRate', sampleRate);
downsampled_IQData = (downsampled_IQData .* mask);
% downsampled_IQData = dRL(real(downsampled_IQData)) + 1i*dRL(imag(downsampled_IQData));
% figure;
% plot(abs(real(downsampled_IQData) + 1i*(imag(downsampled_IQData))));

% differentiator
N = 100;                       % Order
F = [0 0.95 0.98 1];  % Frequency Vector
A = [0 1 0 0];                 % Amplitude Vector
W = [1 1];                     % Weight Vector

% Calculate the coefficients using the FIRPM function.
b  = firls(N, F, A, W, 'differentiator');
Hd = dfilt.dffir(b);
% fvtool(Hd);

[gd, ] = grpdelay(Hd);
delayed_IQData = transpose([zeros(1, gd(1)) transpose(downsampled_IQData)]);
differentiated_IQData = transpose([ transpose(filter(Hd, downsampled_IQData)) zeros(1, gd(1)) ]);
mult_IQData = differentiated_IQData.*delayed_IQData;
% clear downsampled_IQData
% clear delayed_IQData
% clear differentiated_IQData
discrim_out = imag(mult_IQData);
% clear mult_IQData
plot_spectrum(discrim_out, sampleRate);
% discriminator
tau = 75e-6;
T = 1/sampleRate;
alpha = 1/(tan(T/(2*tau)));
% fvtool(numd, dend);
deemphasis_out = filter([1 1], [1+alpha (1-alpha)], discrim_out);
% clear discrim_out

plot_spectrum(deemphasis_out, sampleRate);

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
% clear demphasis_out
audio_out = downsample(output_256, 4);
% clear output_256
sampleRate = sampleRate / 4;
soundsc(audio_out, sampleRate);