function X = td_dft(x, R, L, N, w)
    X = zeros(ceil(length(x)/R), N);
    x_padded = [x zeros(1, L-1)];
    for r = 1:R:length(x)
        x_r = transpose(x_padded(r:(L-1+r))).*w(1:L);
        X_rR = fftshift(fft(x_r, N));
        X(floor((r/R) + 1),:) = X_rR;
    end
end
