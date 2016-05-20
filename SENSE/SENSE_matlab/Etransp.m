function outp = Etransp(data,csm,k_traj,weights,n,coils,wg,osf,scale)

        outp = NuFFT_adj(data,coils,k_traj,weights,n,osf,wg)*scale;   
        outp = sum(conj(csm) .* outp,3);
        outp = outp(:);
        
        
        
        