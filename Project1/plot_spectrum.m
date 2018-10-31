function plot_spectrum(IQData, sampleRate)
    IQData_spectrum = fftshift(fft(IQData));
    N = size(IQData_spectrum, 1);
    f = -sampleRate:2*sampleRate/N:sampleRate;
    f = f(1:end-1);
    plot(f, abs(IQData_spectrum)/N);
end