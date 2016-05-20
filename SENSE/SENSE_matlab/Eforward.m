function outp = Eforward(image,csm,k_traj,weights,n,coils,wg,osf,scale)

        outp = repmat(reshape(image,size(csm,1),size(csm,2)),[1 1 size(csm,3)]) .* csm;
        outp = NuFFT(outp,coils,k_traj,weights,n,osf,wg)/scale;
        outp = outp(:);