function m = NuFFT(image,coils,k_traj,w,n,osf,wg)



%% pad image with zeros to the oversampling size

image = padarray(image, [n*(osf-1)/2,n*(osf-1)/2,0] );

%complex representation of k-space trajectory
k_traj = k_traj(:,1)+1i*k_traj(:,2);


%% Kernel computation 
% width of the kernel on the original grid
kw = wg/osf;
kosf = floor(0.91/(osf*1e-3));
kwidth = osf*kw/2;
beta = pi*sqrt((kw*(osf-0.5)).^2-0.8);

% compute kernel
om = [0:kosf*kwidth]/(kosf*kwidth);
p = besseli(0,beta*sqrt(1-om.*om));
p = p./p(1);
p(end) = 0;

   
    %% Divide image by Fourier transform of the kernel:
    % compute deapodization function
    x = [-osf*n/2:osf*n/2-1]/n;
    sqa = sqrt(pi^2*kw^2*x.^2-beta^2);
    dax = sin(sqa)./(sqa);
    % normalize by DC value
    dax = dax/dax(osf*n/2);
    % make it a 2D array
    da = dax'*dax;
    da = repmat(da,[1,1,8]);
    image = image./da;
   
[A,B] = size(k_traj);
m = zeros(A,coils);

%Do the NuFFT for each coil
for j=1:coils
    
%Compute FFT of each coil image
kspace(:,:,j) = fft2c(image(:,:,j));

% convert k-space trajectory to matrix indices
nx = (n*osf/2+1) + osf*n*real(k_traj);
ny = (n*osf/2+1) + osf*n*imag(k_traj);
    
% loop over samples in kernel at grid spacing
    for lx = -kwidth:kwidth,
      for ly = -kwidth:kwidth,

        % find nearest samples
        nxt = round(nx+lx);
        nyt = round(ny+ly);

        % seperable kernel value
        kkr = min(round(kosf*sqrt(abs(nx-nxt).^2+abs(ny-nyt).^2)+1), floor(kosf*kwidth)+1);
        kwr = p(kkr);

        % if data falls outside matrix, put it at the edge, zero out below
        nxt = max(nxt,1); nxt = min(nxt,osf*n);
        nyt = max(nyt,1); nyt = min(nyt,osf*n);

        for i = 1:length(nxt)
            vector(i) = kspace(nxt(i),nyt(i),j);
        end

        % accumulate gridded data

        m(:,j) = m(:,j) + vector.'.*kwr';
      end
    end

end

m = m.*repmat(sqrt(w),[1,8]);










