function s = surrogate(x,K)
    phases = exp(2i*pi*randn(size(x,1),K));
    s = real(ifft(phases.*fft(x).*flipud(conj(phases))));
end