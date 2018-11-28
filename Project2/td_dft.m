function X = td_dft(x, R, L, N, w)
    X = zeros(ceil(length(x)/R), N);
    x_padded = [x zeros(1, L-1)];
    for r = 1:R:length(x)
        x_r = transpose(x_padded(r:(L-1+r))).*w(1:L);
        X_rR = fftshift(fft(x_r, N));
        X(floor((r/R) + 1),:) = X_rR;
        if ((r == 1) || (r >= -R+length(x)))
            figure; plot(20*log10(abs(X_rR)))
        end
    end
end
% x_padded = x;
%     X = zeros(ceil((length(x)+L-1)/R), N);
%     if (N < L)
%         error("N must be greater than or equal to L");
%     else
%         if (N > L)
%             x_padded = [x zeros(1, N-length(x))];
%         end
%     end
%     for l = 1:N
%         w_minus_n = w(length(w):-1:1);
%         W = transpose(exp((1i*2*pi*l).*(1:length(w))./N));
%         h_k = w_minus_n .* W;
%         X_k = conv(x_padded, h_k);
%         X_k = downsample(X_k, R);
%         X(:,l) = X_k;
%     end
%     X = transpose(X);
% end