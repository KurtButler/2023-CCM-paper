function y = surrogate(x)
    N = size(x,1);
    if mod(N,2) % Odd N
        m = floor(0.5*N) + 1;
        u = rand(m-1,1);
        phi = exp( 2i*pi* [0; u; -flipud(u)]);
    else
        m = floor(0.5*N);
        u = rand(m-1,1);
        phi = exp( 2i*pi* [0; u; 0; -flipud(u)]);
    end
    Y = fft(x).*phi;
    y = real(ifft(Y));
end