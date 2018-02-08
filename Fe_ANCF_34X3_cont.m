function Fe=Fe_ANCF_34X3_cont(~,ElemDofs,E,G,lam,ks2,ks3,L,H,W,ee,nxi,neta,nzeta)
% Funktio: Laskee elastiset voimat k‰ytt‰en selektiivist‰ integrointia

%Fe=zeros(1,ElemDofs);
%Fev=zeros(1,ElemDofs);
Fe0=zeros(1,ElemDofs);


% gauleg laskee gaussin integrointipisteet ja painot

 %  detF0=1/8*L*H*W;     % f(x)=f(xi)|J|, only scaling (no rot, no def) between x and xi  
          [xiv,wxi]=gauleg2(-1,1,nxi);
          [etav,weta]=gauleg2(-1,1,neta);
          [zetav,wzeta]=gauleg2(-1,1,nzeta);
         
%          nxiL=3;
%           
%       [xivL,wxiL] = Lobatto(nxiL);
%        [xivL,wxiL] = gauleg2(-1,1,nxiL);
%           

% 3D palkkia varten
detF0=1/8*L*H*W;    % f(x)=f(xi)|J|, only scaling (no rot, no def) between x and xi  

%detF0v=L/2;          % dx=detF0 dxi

%for kk1=1:nzeta,        
%for jj1=1:neta,


% % Gaussin integrointia 
% for ii1=1:nxi,
%     x=xv(ii1);
%     Fe=Fe+detF0*Fedx_ANCF_3333(E,G,ks2,ks3,H,W,L,A,Iy,Iz,ee,x)*wxi(ii1);  
% end


% % Gaussian quadratures for the elstic forces due to Simo and Vu-Quoc approach and Lobatto for the thickness terms
% % nyt k‰ytet‰‰n yht‰ monta integroitnipistett‰ moelmmissa 
for kk1=1:nzeta,        
for jj1=1:neta,
for ii1=1:nxi,
    xxi=L/2*xiv(ii1);
    eeta=H/2*etav(jj1);
    zzeta=W/2*zetav(kk1);
    Fe0=Fe0+FedV_34X3(E,G,lam,L,ks2,ks3,xxi,eeta,zzeta,ee)*wxi(ii1)*weta(jj1)*wzeta(kk1);      
%     Fe0=Fe0+FedV_ANCF_331011101103_cont_Marko(E,nu,G,ks2,ks3,H,W,L,A,Iy,Iz,ee,xxi,eeta,zzeta)*wxi(ii1)*weta(jj1)*wzeta(kk1);
end
end
end

% for iiL=1:nxiL,
%      xxi=xivL(iiL);
%      Fev=Fev+FeVdx_ANCF_3333_ench(E,nu,G,ks2,ks3,H,W,L,A,Iy,Iz,ee,xxi)*wxiL(iiL);
% end

Fe=Fe0*detF0;
%Fe=(Fe0*detF0+Fev*detF0v);