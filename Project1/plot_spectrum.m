function plot_spectrum(IQData, sampleRate, central_freq)
    spectrum_seg = IQData(1:(sampleRate/200));
    spectrum_seg = spectrum_seg .* blackman((sampleRate/200));
    IQData_spectrum = fftshift(fft(spectrum_seg));
    N = size(IQData_spectrum, 1);
    f = -(sampleRate/2):(sampleRate/N):(sampleRate/2 - (sampleRate/N));
    f = f + central_freq;
    f = f./1e6;
    figure;
    plot(f, 20*log10(abs(IQData_spectrum)));
    xlim([(central_freq  -(sampleRate/2))/1e6, (central_freq + (sampleRate/2))/1e6]);
    xlabel('Frequency (MHz)');
    ylabel('Magnitude (dB)');
    grid minor
    set(findall(gcf,'-property','FontSize'),'FontSize',20)
end