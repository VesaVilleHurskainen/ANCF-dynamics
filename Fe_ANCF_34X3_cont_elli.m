function Fe=Fe_ANCF_34X3_cont_elli(ElemDofs,E,nu,ks2,ks3,L,H,W,ee,nxi)
% Implemented for elliptic cross section by Henrik Ebel on 26.10.2015.
% adapted for element 34X3 by VVH

Fe0=zeros(1,ElemDofs);

% in xi direction use Gaussian quadrature
[xiv,wxi]=gauleg2(-1,1,nxi);
% in eta/zeta directions use circle integration formula from Abrahamowitz and Stegun (1972)
[pcirc,wcirc] = circleInteg(5);

detF0=(L*H*W)/8;
% detF0=L/2;

prefac = pi; % pre-factor from integration formula
for kk1=1:length(wcirc),
        for ii1=1:nxi,
            xxi=xiv(ii1);
            eetaM=pcirc(kk1,1);
            zzetaM=pcirc(kk1,2);
            Fe0=Fe0+prefac*FedV_ANCF_34X3_cont_mex(E,nu,NaN,ks2,ks3,H,W,L,NaN,NaN,NaN,ee,xxi,eetaM,zzetaM)*wxi(ii1)*wcirc(kk1);
        end
end

Fe=detF0*Fe0;