% This function computes the right side of the ODE

function [ydot, flag, new_data] = eom_Dynam(t,y,data)

% --- Extract parameters from struct --- 
xloc=data.xloc;             % Matrix of DOFs per element
n=data.n;                   % Number of elements
nn=data.nn;                 % Number of nodes
nx=data.nx;                 % Total DOFs of system
ndof=data.ndof;             % Free DOFs of system
ElemDofs=data.ElemDofs;     % DOFs per element
E=data.E;                   % Young's modulus
% G=data.G;                   % Shear modulus
nu=data.nu;                 % Poisson's ratio
ks=data.ks;                 % Shear correction factor
L=data.L;                   % Beam length
bc=data.bc;                 % Vector indicating linearly constrained DOF
elementID=data.elementID;   % ID number of element being used
CrossSec=data.CrossSec;     % Cross section number
Wi1 = data.Wi1;
Wi2 = data.Wi2;
Wo = data.Wo;
ee0 = data.ee0;
ee0dot = data.ee0dot;
Case = data.Case;
    
ee=zeros(nx,1);             % nodal positions, for all DOFs
eedot=zeros(nx,1);          % nodal velocities, -"-

% --- Positions ---
ee(~bc)=ee0(~bc);
ee(bc)=y(1:ndof);

% --- Velocities ---
eedot(~bc)=ee0dot(~bc);
eedot(bc)=y(ndof+1:ndof*2);

Fe=zeros(nx,1);
Fext=zeros(nx,1);

% ---- Constraints ---
% Rotation imposed on 1st end of beam

if Case ~= 3
    [a,ad,~] = motion_case10(t);

    q5a=cos(a);
    q6a=sin(a);
    q8a=-q6a;
    q9a=q5a;

    q5ad=-ad*sin(a);
    q6ad=ad*cos(a);
    q8ad=-q6ad;
    q9ad=q5ad;

    ee(5) = q5a;
    ee(6) = q6a;
    ee(8) = q8a;
    ee(9) = q9a;
    eedot(5) = q5ad;
    eedot(6) = q6ad;
    eedot(8) = q8ad;
    eedot(9) = q9ad;

    ee(xlocANCF_3403(nn,5)) = q5a;
    ee(xlocANCF_3403(nn,6)) = q6a;
    ee(xlocANCF_3403(nn,8)) = q8a;
    ee(xlocANCF_3403(nn,9)) = q9a;
    eedot(xlocANCF_3403(nn,5)) = q5ad;
    eedot(xlocANCF_3403(nn,6)) = q6ad;
    eedot(xlocANCF_3403(nn,8)) = q8ad;
    eedot(xlocANCF_3403(nn,9)) = q9ad;
end
	
% --- Loop over elements - calculate mass matrix and external forces ---
for k = 1:n
        xlock = xloc(k,:);
        eek=ee(xlock);

        if elementID==3403
            ks2=ks;
            ks3=ks;
            %%%%%%
            nxi=6; %5
            if CrossSec == 1
				error('')
            elseif CrossSec == 2
				Fek = Fe_ANCF_34X3_cont_elli(ElemDofs,E,nu,ks2,ks3,L/n,Wo,Wo,eek,nxi)';
            elseif CrossSec == 3
				Fek = Fe_ANCF_34X3_cont_elli(ElemDofs,E,nu,ks2,ks3,L/n,Wi1,Wo,eek,nxi)';
            elseif CrossSec == 4
				error('')
            elseif CrossSec == 5
				Fek = Fe_ANCF_34X3_cont_elli(ElemDofs,E,nu,ks2,ks3,L/n,Wo,Wo,eek,nxi)'...
                    - Fe_ANCF_34X3_cont_elli(ElemDofs,E,nu,ks2,ks3,L/n,Wi1,Wi2,eek,nxi)';
            end
        else
            error('Please select a valid element. Error in ''eom_dynam.m''.\n')
        end
               
        Fe(xlock)=Fe(xlock)+Fek; 
end

if Case == 3
    Fext(end-30+2) = -1000;
end

% --- Apply linear constraints ---
Fextc=Fext(bc);
Fec=Fe(bc);

% --- Force equation (external - internal) ---
Fextcall=Fextc-Fec;
y1dot=eedot(bc);            % velocities

% --- No nonlinear constraints, use sparse matrices ---
Fextcsall=sparse(Fextcall);
y2sdot=data.MInvSparse*Fextcsall;
y2dot=full(y2sdot);

% --- Results to vector ---
ydot = zeros(2*ndof,1);
ydot(1:ndof)=y1dot(1:ndof); 
ydot(ndof+1:2*ndof)=y2dot(1:ndof);
ydot=ydot(:);

flag = 0;
new_data = [];
