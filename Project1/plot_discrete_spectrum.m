function plot_discrete_spectrum(IQData, sampleRate)
    spectrum_seg = IQData(1:(sampleRate/1));
    spectrum_seg = spectrum_seg .* blackman((sampleRate/1));
    IQData_spectrum = fftshift(fft(spectrum_seg));
    N = size(IQData_spectrum, 1);
    f = -(sampleRate/2):(sampleRate/N):(sampleRate/2 - (sampleRate/N));
    f = f./(sampleRate/2);
    figure;
    plot(f, 20*log10(abs(IQData_spectrum)));
    xlim([-1 1]);
    xlabel('\omega (\pi radians/sample)');
    ylabel('Magnitude (dB)');
    grid minor
    set(findall(gcf,'-property','FontSize'),'FontSize',20)
end