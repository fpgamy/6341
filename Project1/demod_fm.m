fid = fopen('gqrx_20181029_051700_105409000_2048000_fc_105p7MHz.raw','rb');
central_freq = 105.7e6;
hw_freq = 105409000;
sampleRate = 2.048e6;

% read the data from the file
IQData = fread(fid,'float32');
% convert it to complex values
IQData = IQData(1:2:end) + 1i*IQData(2:2:end);
% plot the spectrum
plot_spectrum(IQData, sampleRate, hw_freq);
title('$$|X(j\Omega)|$$','interpreter','latex');
y1=get(gca,'ylim');

% plot the hw are central frequency
hold on
plot([hw_freq/1e6 hw_freq/1e6],y1, 'Color','red','LineWidth', 3,'LineStyle','--');
plot([central_freq/1e6 central_freq/1e6],y1, 'Color','green','LineWidth', 3,'LineStyle','--');
legend('Spectrum','HW Frequency', 'Station Freqency');
hold off

% frequency shift so the central frequency is in the centre of the
% frequency band
sample_no = transpose(1:1:length(IQData));
IQData = IQData .* exp(1i.*sample_no.*2*pi*((hw_freq-central_freq)/2.048e6));

% plot of the shifted frequency response
plot_spectrum(IQData, sampleRate, central_freq);
title('$$|X(j\Omega)|$$','interpreter','latex');
y1=get(gca,'ylim');

% plot the hw are central frequency
hold on
plot([hw_freq/1e6 hw_freq/1e6],y1, 'Color','red','LineWidth', 3,'LineStyle','--');
plot([central_freq/1e6 central_freq/1e6],y1, 'Color','green','LineWidth', 3,'LineStyle','--');
legend('Spectrum','HW Frequency', 'Station Freqency');
hold off

% anti-aliasing Filter
Fpass = 126000;           % Passband Frequency
Fstop = 128000;           % Stopband Frequency
Dpass = 0.057501127785;   % Passband Ripple
Dstop = 0.0031622776602;  % Stopband Attenuation
dens  = 20;               % Density Factor

% calculate the order from the parameters using FIRPMORD.
[N, Fo, Ao, W] = firpmord([Fpass, Fstop]/(sampleRate/2), [1 0], [Dpass, Dstop]);

% calculate the coefficients using the FIRPM function.
b  = firpm(N, Fo, Ao, W, {dens});
Hd = dfilt.dffir(b);
disp(Hd);
% plot the frequency response
fvtool(Hd);

% apply the filter and plot the resulting spectrum
aa_IQData = filter(Hd, IQData);
plot_spectrum(aa_IQData, sampleRate, central_freq);
title('Anti-Aliased $$|X(j\Omega)|$$','interpreter','latex');
y1=get(gca,'ylim');

% plot the edge frequencies
hold on
plot([(central_freq-128e3)/1e6 (central_freq-128e3)/1e6],y1, 'Color','red','LineWidth', 3,'LineStyle','--')
plot([(central_freq+128e3)/1e6 (central_freq+128e3)/1e6],y1, 'Color','red','LineWidth', 3,'LineStyle','--')
legend('Spectrum', 'Stopband Frequency');
hold off

% decimate by factor of 8
downsampled_IQData = downsample(aa_IQData, 8);
sampleRate = sampleRate / 8;

% plot the resulting spectrum
plot_spectrum(downsampled_IQData, sampleRate, central_freq);
title('$$|X_d(j\Omega)|$$','interpreter','latex');

% plot comparing the input and output of the decimator
figure;
subplot(2, 1, 1);
plot_discrete_spectrum(aa_IQData, 8*sampleRate);
title('Anti-Aliased $$|X(e^{j\omega})|$$','interpreter','latex');

subplot(2, 1, 2);
plot_discrete_spectrum(downsampled_IQData, sampleRate);
title('$$|X_d(e^{j\omega})|$$','interpreter','latex');

% limiter 
downsampled_IQData = downsampled_IQData./abs(downsampled_IQData);

% plot the spectrum of the limiter
plot_spectrum(downsampled_IQData, sampleRate, central_freq);
title('Amplitude Limited $$|X_{d}(j\Omega)|$$','interpreter','latex');

% differentiator
M = 32;
beta = 2.4;

% use kaiser window method
w = kaiser(M, beta);
n = -((M-1)/2):1:((M-1)/2);
h_diff = (cos(pi * (n)) ./ (n)) - (sin(pi * (n)) ./ (pi * (n).^2));
b = w.*transpose(h_diff);

% plot the frequency response
fvtool(b, 1);
disp('Differentiator Type:');
disp(firtype(b));

differentiated_IQData = filter(b, 1, downsampled_IQData);

% plot the spectrum of the differentiated data
plot_spectrum(differentiated_IQData, sampleRate, central_freq);
title('$$|Y_{differentiated}(j\Omega)|$$','interpreter','latex');

% gd = (M-1)/2;
[gd, ] = grpdelay(b, 1);
gd = mean(gd);

% non-integer delay

% bandlimited interpolation filter
bi = intfilt(2, 2, 1);

% plot the response
fvtool(bi, 1);
% upsample and interpolate
upsampled_IQData = upsample(downsampled_IQData, 2);

%group delay = 3
interpolated = filter(bi, 1, upsampled_IQData);

% anti-aliasing filter
% ensure group delay is 28
N     = 56;    % Order
Fpass = 0.48;  % Passband Frequency
Fstop = 0.5;   % Stopband Frequency
Wpass = 1;     % Passband Weight
Wstop = 1;     % Stopband Weight
dens  = 20;    % Density Factor

% calculate the coefficients using the FIRPM function.
b  = firpm(N, [0 Fpass Fstop 1], [1 1 0 0], [Wpass Wstop], {dens});
Hd = dfilt.dffir(b);

% plot the spectrum
fvtool(Hd);

% apply the filter
aa_interpolated = filter(b, 1, interpolated);

% in an ideal system the magnitude of the samples would still be limited
aa_interpolated((1+(gd*2)):end) = aa_interpolated((1+(gd*2)):end)./abs(aa_interpolated((1+(gd*2)):end));

% decimate
delayed_IQData = downsample(aa_interpolated, 2);

% plot the spectrum of the delayed samples. Magnitude should be
% approximately unchanged
plot_spectrum(delayed_IQData, sampleRate, central_freq);
title('$$|Y_{delay}(j\Omega)|$$','interpreter','latex');


% plot the time domain samples to show the delay
figure;
subplot(2, 1, 1);
samples = 1:length(upsampled_IQData);
samples = samples./2;
% must interpolate by factor 2
plot(samples, real(filter(bi, 1, upsample(downsampled_IQData, 2))));
hold on;
plot(samples, real(filter(bi, 1, upsample(delayed_IQData, 2))), 'color', 'red');
hold off;
grid minor;
title('Input and Output of Delay System');
ylabel('Real Part');
xlabel('Sample Number');
xlim([100 400]);
legend('Input', 'Delayed');

subplot(2, 1, 2);
plot(samples, imag(filter(bi, 1, upsample(downsampled_IQData, 2))));
hold on;
plot(samples, imag(filter(bi, 1, upsample(delayed_IQData, 2))), 'color', 'red');
hold off;
grid minor;
ylabel('Imaginary Part');
xlabel('Sample Number');
xlim([100 400]);

% multiply the output of the differentiator and the delay
mult_IQData = differentiated_IQData.*conj(delayed_IQData);

% plot the spectrum
plot_spectrum(mult_IQData, sampleRate, central_freq);
title('$$|Y_{2}(j\Omega)|$$','interpreter','latex');

% take the imaginary part
discrim_out = imag(mult_IQData);
% plot the spectrum
plot_spectrum(discrim_out, sampleRate, central_freq);
title('$$|M(j\Omega)|$$','interpreter','latex');

% discriminator
tau = 75e-6;
T = 1/sampleRate;

% calculated using the bilinear transform
alpha = 1/(tan(T/(2*tau)));

% plot the frequency response
fvtool([1 1], [1+alpha (1-alpha)]);

% apply the filter
deemphasis_out = filter([1 1], [1+alpha (1-alpha)], discrim_out);

% plot the spectrum
plot_spectrum(deemphasis_out, sampleRate, central_freq);
title('$$|M_{d}(j\Omega)|$$','interpreter','latex');

% low pass filter around 15kHz to get mono audio
Fpass = 14.5;    % Passband Frequency
Fstop = 15;      % Stopband Frequency
Apass = 0.5;     % Passband Ripple (dB)
Astop = 50;      % Stopband Attenuation (dB)
match = 'both';  % Band to match exactly

% construct an FDESIGN object and call its ELLIP method.
h  = fdesign.lowpass(Fpass, Fstop, Apass, Astop, sampleRate/1000);
Hd = design(h, 'ellip', 'MatchExactly', match);

% plot the frequency response
fvtool(Hd);

% apply the filter
output_256 = filter(Hd, deemphasis_out);

% plot the spectrum
plot_spectrum(output_256, sampleRate, central_freq);
title('$$|A_{m}(j\Omega)|$$','interpreter','latex');

% decimate by a factor of 4 to get a playable sampling rate
audio_out = downsample(output_256, 4);

sampleRate = sampleRate / 4;

%plot the spectrum
plot_spectrum(audio_out, sampleRate, central_freq);
title('$$|A(j\Omega)|$$','interpreter','latex');

% play the audio with normalised volume
soundsc(audio_out, sampleRate);