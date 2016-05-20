function f = ifft2c(F)

f = fftshift(ifft2(ifftshift(F)));


