function p = td_dft_point(x, n, k, L, N, w)
    acc = complex(0, 0);
    for m = 0:L-1
        acc = acc + x(n+m)*w(m)*exp(-1i*2*pi*k/N);
    end
    p = acc;
end