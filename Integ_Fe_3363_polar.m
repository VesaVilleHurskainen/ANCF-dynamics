function Fe = Integ_Fe_3363_polar(ElemDofs,E,nu,G,ks2,ks3,H,W,L,ee,nxi,na,nr)
% Function: Integrate elastic forces using Gaussian quadrature

Fe=zeros(1,ElemDofs);

Wi = H;
Wo = W;

lam = (nu*E)/((1+nu)*(1-2*nu));

% Get quadrature points and weights
% nxi = 4;
% na = 6;
% nr = 2;
 
[xiv,wxi]=gauleg2(-1,1,nxi);
[av,wa]=gauleg2(-1,1,na);
[rv,wr]=gauleg2(-1,1,nr);


% Compute quadratures
for kk1=1:nr       
for jj1=1:na
for ii1=1:nxi
    xxi=xiv(ii1)*L/2;
    aa=av(jj1)*pi+pi;
    rr=rv(kk1)*((Wo-Wi)/4)+(Wi+Wo)/4;
    eeta=rr*cos(aa);
    zzeta=rr*sin(aa);
    Fe=Fe+Fedx_3363c(E,G,lam,L,ks2,ks3,xxi,eeta,zzeta,ee)'*rr*wxi(ii1)*wa(jj1)*wr(kk1);
end
end
end

%for iiL=1:nxiL
%    xxi=xivL(iiL)*L/2;
%    Fev=Fev+FeVdx_3363c(E,nu,L,A,xxi,ee)'*wxiL(iiL);
%end

% Sum up components (with scaling factors for interval change)
detF0=(L/2)*pi*((Wo-Wi)/4);
%detFv=L/2;
Fe=(Fe*detF0);

% for ii1 =1:nxi
%     xxi=xiv(ii1)*L/2;
%     fi = @(y,z) Fedx_3363c(E,G,lam,L,ks2,ks3,xxi,y,z,ee);
%     Fe = Fe+(L/2)*(circleInteg(fi,Wo/2,5,ElemDofs)-circleInteg(fi,Wi/2,5,ElemDofs))'*wxi(ii1);
% end