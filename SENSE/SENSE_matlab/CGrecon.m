clear all 
close all
clc;
% get data from http://hansenms.github.io/sunrise/sunrise2013/
load('hansen_exercises.mat')
%Decorrelate noise
dmtx = ismrm_calculate_noise_decorrelation_mtx(noise_spiral);
data_spiral = ismrm_apply_noise_decorrelation_mtx(data_spiral,dmtx);
smaps_prew = ismrm_apply_noise_decorrelation_mtx(smaps,dmtx);
csm = smaps_prew;

%Prepare NUFFT
N = [size(smaps,1) size(smaps,2)];
osf =2;
wg = 5;
K = N*osf;
n=N(1);
coils = size(csm,3);
scale = sqrt(prod(K));

E_forward = @(x) Eforward(x,csm,k_spiral,w_spiral,n,coils,wg,osf,scale);
E_transp = @(y) Etransp(y,csm,k_spiral,w_spiral,n,coils,wg,osf,scale);

%Generate RHS 
coildata = data_spiral(:).*repmat(sqrt(w_spiral),[coils 1]);

%% Algo
x=zeros(size(csm,1)*size(csm,2),1);
b = E_transp(coildata);

r = b-E_transp(E_forward(x));
    p=r;
    count=0;
for i = 1:13
        rsold=r'*r;
        EHEp = E_transp(E_forward(p));
        alpha=rsold/(p'*EHEp);
        x=x+alpha*p;
        r=r-alpha*EHEp;   
        rsnew=r'*r;
        beta = (rsnew)/(rsold);
        p=r+beta*p;
        rsold=rsnew;
        count = count+1;
%         filename = sprintf('%s_%d.png','iteration',i);
%         img = reshape(x,size(csm,1),size(csm,2));
%         img = img./max(max(abs(img)));
%         imwrite(abs(img),filename)
end
relresidual = norm(b-E_transp(E_forward(x)),2)/norm(b);

img = reshape(x,size(csm,1),size(csm,2));

figure,imshow(abs(img),[]);

