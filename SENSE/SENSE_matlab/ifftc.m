function f = ifftc(F)

f = fftshift(ifft(ifftshift(F)));


