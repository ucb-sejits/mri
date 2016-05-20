function f = fftc(F)

f = ifftshift(fft(fftshift(F)));


