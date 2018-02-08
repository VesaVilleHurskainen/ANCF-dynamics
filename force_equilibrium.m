% Force equilibrium equation for computation of static equilibrium position
% Simplified version of eom_Dynam.

function Feq = force_equilibrium(eebc,data)

% --- Extract parameters from struct --- 
xloc=data.xloc;             % Matrix of DOFs per element
%nloc=data.nloc;
n=data.n;                   % Number of elements
nn=data.nn;                 % Number of nodes
nx=data.nx;                 % Total DOFs of system
ndof=data.ndof;             % Free DOFs of system
ElemDofs=data.ElemDofs;     % DOFs per element
E=data.E;                   % Young's modulus
G=data.G;                   % Shear modulus
nu=data.nu;                 % Poisson's constant
%lambda=data.lambda;
ks=data.ks;
rho=data.rho;               % Density
g=data.g;                   % Gravitational constant (used in cases 1-3)
L=data.L;                   % Beam length
bc=data.bc2;                 % Vector indicating linearly constrained DOF
elementID=data.elementID;   % ID number of element being used
CrossSec=data.CrossSec;     % Cross section number
Case=data.Case;             % Case number
%IntegFormu=data.IntegFormu;
if CrossSec == 5
    Wi1 = data.Wi1;
    Wi2 = data.Wi2;
    Wo = data.Wo;
else
    H=data.H;                   % Beam height
    W=data.W;                   % Beam width
end
alpha = data.alpha;         % Baumgarte stabilization constant alpha
beta = data.beta;           % Baumgarte stabilization constant beta
A=data.A;                   % Beam cross-sectional area
Iy=data.Iy;                 % Beam square moments of inertia
Iz=data.Iz;
    
Fe=zeros(nx,1);
% M=zeros(nx,nx);
Fext=zeros(nx,1);

ee = data.ee0;
ee(bc) = eebc;

% --- Loop over elements - calculate mass matrix and external forces ---
for k = 1:n
        xlock = xloc(k,:);
        eek=ee(xlock);

        if elementID==3363
            ks2=ks;
            ks3=ks;
            %%%%%%
            nxi=5; %5
            nyi=5; %5
            nzi=5; %5
            if CrossSec == 1
                % rectangular cross section
                Fek=element3363_mex('Fe_ANCF_3363_cont',elementID,ElemDofs,E,nu,NaN,ks2,ks3,H,W,L/n,NaN,NaN,NaN,eek,nxi,nyi,nzi)';
            elseif CrossSec == 2
                % circular cross section
                % only uses nxi from above, for circular cross-section
                % integration: see function:
                Fek=element3363_mex('Fe_ANCF_3363_cont_circ',elementID,ElemDofs,E,nu,NaN,ks2,ks3,H,W,L/n,NaN,NaN,NaN,eek,nxi,nyi,nzi)';
%                 Fek=Fe_ANCF_3363_cont_circ(Element,ElemDofs,E,nu,G,ks2,ks3,H,W,L,A,Iy,Iz,eek,nxi,nyi,nzi)';
                % look at number of Gauss points in fcn above!
            elseif CrossSec == 3
                % only uses nxi from above, for elliptic cross-section
                % integration: see function:
                Fek=element3363_mex('Fe_ANCF_3363_cont_elli',elementID,ElemDofs,E,nu,NaN,ks2,ks3,H,W,L/n,NaN,NaN,NaN,eek,nxi,nyi,nzi)';
%                 Fek=Fe_ANCF_3363_cont_elli(Element,ElemDofs,E,nu,G,ks2,ks3,H,W,L,A,Iy,Iz,eek,nxi,nyi,nzi)';
                % look at number of Gauss points in fcn above!
            elseif CrossSec == 4
%                 Fek=Integ_Fe_3363_polar_mex(ElemDofs,E,nu,G,ks2,ks3,H,W,L/n,eek,4,8,4)';
                  Fek=element3363_mex('Fe_ANCF_3363_cont_circ',elementID,ElemDofs,E,nu,NaN,ks2,ks3,W,W,L/n,NaN,NaN,NaN,eek,nxi,nyi,nzi)'...
                    -element3363_mex('Fe_ANCF_3363_cont_elli',elementID,ElemDofs,E,nu,NaN,ks2,ks3,H,H+0.05,L/n,NaN,NaN,NaN,eek,nxi,nyi,nzi)';
            elseif CrossSec == 5
			    % hollow circular cross section
				% integrate in polar coordinates (experimental)
%                 Fek=Integ_Fe_3363_polar_mex(ElemDofs,E,nu,G,ks2,ks3,H,W,L/n,eek,4,8,4)';
                 Fek=element3363_mex('Fe_ANCF_3363_cont_elli',elementID,ElemDofs,E,nu,NaN,ks2,ks3,Wo,Wo,L/n,NaN,NaN,NaN,eek,nxi,nyi,nzi)'...
                    -element3363_mex('Fe_ANCF_3363_cont_elli',elementID,ElemDofs,E,nu,NaN,ks2,ks3,Wi1,Wi2,L/n,NaN,NaN,NaN,eek,nxi,nyi,nzi)';
            end
        elseif elementID==3333
            ks2=ks;
            ks3=ks;
            %%%%%%
            % St.-Venant-Kirchhoff approach with splitting of D
            nxi=4;
            nyi=3;
            nzi=3;
            if CrossSec==1
                % rectangular cross section
%                 Fek=Fe_ANCF_3333_ench(Element,ElemDofs,E,nu,G,ks2,ks3,H,W,L,A,Iy,Iz,eek,nxi,nyi,nzi)';
%                 Fek=Fe_ANCF_3333csplit(Element,ElemDofs,E,nu,G,ks2,ks3,H,W,L,A,Iy,Iz,eek,nxi,nyi,nzi)';
                Fek=element3333_mex('Fe_ANCF_3333csplit',elementID,ElemDofs,E,nu,NaN,ks2,ks3,H,W,L/n,NaN,NaN,NaN,eek,nxi,nyi,nzi)';
            elseif CrossSec==2
                % circular cross section
                % only uses nxi from above (for 1d integral) -> for 2d
                % circular integral: see function itself:
                %Fek=Fe_ANCF_3333_ench_circ(Element,ElemDofs,E,nu,G,ks2,ks3,H,W,L,A,Iy,Iz,eek,nxi,nyi,nzi)';
%                 Fek=Fe_ANCF_3333c_split_circ(Element,ElemDofs,E,nu,G,ks2,ks3,H,W,L,A,Iy,Iz,eek,nxi,nyi,nzi)';
                Fek=element3333_mex('Fe_ANCF_3333c_split_circ',elementID,ElemDofs,E,nu,NaN,ks2,ks3,H,W,L/n,NaN,NaN,NaN,eek,nxi,nyi,nzi)';
            elseif CrossSec == 3
                % only uses nxi from above, for elliptic cross-section
%                 integration: see function:
%                 Fek=Fe_ANCF_3333_ench_elli(Element,ElemDofs,E,nu,G,ks2,ks3,H,W,L,A,Iy,Iz,eek,nxi,nyi,nzi)';
%                 Fek=Fe_ANCF_3333c_split_elli(Element,ElemDofs,E,nu,G,ks2,ks3,H,W,L,A,Iy,Iz,eek,nxi,nyi,nzi)';
                Fek=element3333_mex('Fe_ANCF_3333c_split_elli',elementID,ElemDofs,E,nu,NaN,ks2,ks3,H,W,L/n,NaN,NaN,NaN,eek,nxi,nyi,nzi)';
                % look at number of Gauss points in fcn above!
            else
                error('Please choose valid cross section! Error in function ''eom_Dynam''.\n')
            end
        elseif elementID==33330
            ks2=ks;
            ks3=ks;
            %%%%%%
            nxi=2;
            nyi=3;
            nzi=3;
            if CrossSec==1
                % rectangular cross section
                error('Cross section not implemented for 33330')
            elseif CrossSec==2
                % circular cross section
                % only uses nxi from above (for 1d integral) -> for 2d
                Fek=element3333_mex('Fe_ANCF_3333s',elementID,ElemDofs,E,nu,G,ks2,ks3,H,W,L/n,A,Iy,Iz,eek,nxi,nyi,nzi)';
            elseif CrossSec == 3
                % only uses nxi from above, for elliptic cross-section
                error('Cross section not implemented for 33330')
            elseif CrossSec == 4
                % hollow circular cross section
                % only uses nxi from above (for 1d integral) -> for 2d
                Fek=element3333_mex('Fe_ANCF_3333s',elementID,ElemDofs,E,nu,G,ks2,ks3,H,W,L/n,A,Iy,Iz,eek,nxi,nyi,nzi)';
            else
                error('Please choose valid cross section! Error in function ''eom_Dynam''.\n')
            end
        elseif elementID==3403
%             Mk=MassM(elementID,rho,H,W,L/n);
            ks2=ks;
            ks3=ks;
            %%%%%%
            nxi=5; %5
            nyi=5; %5
            nzi=5; %5
            if CrossSec == 1
				error('')
            elseif CrossSec == 2
				Fek = Fe_ANCF_34X3_cont_circ(elementID,ElemDofs,E,G,data.lambda,ks2,ks3,L/n,H,W,eek,nxi,nyi,nzi)';
            elseif CrossSec == 3
				Fek = Fe_ANCF_34X3_cont_elli(elementID,ElemDofs,E,G,data.lambda,ks2,ks3,L/n,H,W,eek,nxi,nyi,nzi)';
            elseif CrossSec == 4
				error('')
            elseif CrossSec == 5
				error('')
            end
        else
            error('Please select a valid element. Error in ''eom_dynam.m''.\n')
        end
        
        if g~=0
            Fextk=FextBF(elementID,nx,rho,g,L/n,H,W,Case,CrossSec);
        else
            Fextk=zeros(length(xlock),1);
        end
        Fe(xlock)=Fe(xlock)+Fek;

%         M(xlock,xlock) = M(xlock,xlock) + Mk;
        Fext(xlock)=Fext(xlock)+Fextk;
end

if Case == 10 %!!!!!!!!!!!!!!!
    %Fext(xlocANCF_3363(data.midnode,5)) = -1e11*Wo;
    Fext(xlocANCF_3403(1,1:120)) = 1e10*[
            0
            0
            0
            0
         H/16
            0
            0
            0
            0
            0
            0
            0
            0
            0
            0
            0
            0
            0
            0
            0
            0
            0
            0
            0
            0
      H^3/384
            0
            0
            0
            0
            0
            0
            0
            0
    -(9*H)/16
            0
            0
            0
            0
            0
            0
            0
            0
            0
            0
            0
            0
            0
            0
            0
            0
            0
            0
            0
            0
 -(3*H^3)/128
            0
            0
            0
            0
            0
            0
            0
            0
    -(9*H)/16
            0
            0
            0
            0
            0
            0
            0
            0
            0
            0
            0
            0
            0
            0
            0
            0
            0
            0
            0
            0
 -(3*H^3)/128
            0
            0
            0
            0
            0
            0
            0
            0
         H/16
            0
            0
            0
            0
            0
            0
            0
            0
            0
            0
            0
            0
            0
            0
            0
            0
            0
            0
            0
            0
      H^3/384
            0
            0
            0
            0];
end

% Add gravitational force from rigid disk or point load in cantilever case     
if elementID == 3333 || elementID == 33330
    if Case == 7
		Fext(xlocANCF_3333(data.midnode,2)) = Fext(xlocANCF_3333(data.midnode,2)) - data.md;
    elseif Case == 9
        Fext(xlocANCF_3333(nn,2)) = data.Fy;
    end
elseif elementID == 3363
	if Case == 7
		Fext(xlocANCF_3363(data.midnode,2)) = Fext(xlocANCF_3363(data.midnode,2)) - data.md;
    elseif Case == 9
        Fext(xlocANCF_3363(nn,2)) = data.Fy;
    end
elseif elementID == 3403
    if Case == 9
        Fext(xlocANCF_3403(nn,2)) = data.Fy;
    end
end

Fextc=Fext(bc);
Fec=Fe(bc);

Feq = Fextc-Fec;