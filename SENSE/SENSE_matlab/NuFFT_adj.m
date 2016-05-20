function image = NuFFT_adj(data,coils,k_traj,w,n,osf,wg)

    
    k = k_traj(:,1)+1i*k_traj(:,2);
    % convert to single column
    w = repmat(w,[coils 1]);
    % width of the kernel on the original grid
    kw = wg/osf;
    % preweight
    dw = data.*sqrt(w);
    %dw = data.*w;
        
    n_k = length(data)/coils;
    dw = reshape(dw,[n_k,coils]);
    
    % compute kernel, assume e1 is 0.001, assuming nearest neighbor
    kosf = floor(0.91/(osf*1e-3));
    % half width in oversampled grid units
    kwidth = osf*kw/2;
    % beta from the Beatty paper
    beta = pi*sqrt((kw*(osf-0.5)).^2-0.8);

    % compute kernel
    om = [0:kosf*kwidth]/(kosf*kwidth);
    p = besseli(0,beta*sqrt(1-om.*om));
    p = p./p(1);
    % last sample is zero so we can use min() below for samples bigger than kwidth
    p(end) = 0;

    % convert k-space samples to matrix indices
    nx = (n*osf/2+1) + osf*n*real(k);
    ny = (n*osf/2+1) + osf*n*imag(k);

    m = zeros(osf*n,osf*n,coils);

for j=1:coils
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

        % accumulate gridded data
        m(:,:,j) = m(:,:,j)+sparse(nxt,nyt,dw(:,j).*kwr',osf*n,osf*n);
      end;
    end;

% zero out data at edges, which is probably due to data outside mtx
m(:,1,j) = 0; m(:,osf*n,j) = 0;
m(1,:,j) = 0; m(osf*n,:,j) = 0;

im(:,:,j) = ifft2c(m(:,:,j));

end

    % compute deapodization function
    x = [-osf*n/2:osf*n/2-1]/n;
    sqa = sqrt(pi^2*kw^2*x.^2-beta^2);
    dax = sin(sqa)./(sqa);
    % normalize by DC value
    dax = dax/dax(osf*n/2);
    % make it a 2D array
    da = dax'*dax;

    da = repmat(da,[1,1,8]);
    % deapodize
    im = im./da;

    %% trim the images
    %correct image size: 
    image = zeros(n,n,coils);
    for j=1:coils
    image(:,:,j) = im((osf-1)*n/2+1:n+(osf-1)*n/2,(osf-1)*n/2+1:n+(osf-1)*n/2,j);
    end
    %return the result
