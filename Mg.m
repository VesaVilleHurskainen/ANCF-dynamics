function M = Mg(~,nloc,nl,nx,Element,rho,H,W,L,CrossSec)      
% Creates mass matrix of system using function MassM to get mass matrix of
% single element.
% For the students written by MKM
% adapted by Henrik Ebel on 10.11.2015
% modified by VVH

M=zeros(nx,nx);

if Element==3403
    xloc = xlocAllANCF_3403(nloc);
else
    disp('****** Please choose valid element. Error in function ''Mg''. ******');
end

% loop over all elements
for k = 1:nl
    xlock = xloc(k,:);
    %eek=ee(xlock);        
    Mk=MassM(Element,rho,H,W,L/nl,CrossSec);
%     Mk=Me_34X3_elli(L/nl,H,W,rho);
    M(xlock,xlock) = M(xlock,xlock) + Mk; 
end