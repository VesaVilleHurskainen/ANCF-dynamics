% Beam element large displacement dynamic test code
% Coded by MKM, 2011
% Modified by Henrik Ebel, 2015
% Modified further by VVH, 2016-2017

clc;
clear;
close all;
format long;

resolution=10; % resolution of visualization - points per length unit

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                              PARAMETERS                               %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Type of analysis                      %
% 1: dynamic large displacement problem
Analysis=1;

%%%%%%%%%%%%%%%%%%%
% Element options %
% 3333 : three-node beam element with second-order approximation in 
%        x-direction and linear approximation in transversal directions. 
%        Uses  splitting-approach for Poisson-locking-free calculation of 
%        the elastic forces (see Nachgebauer et al. paper, in principal 
%        selective reduced integration). 
%        Uses St. Venant-Kirchhoff material law.
% 3363 : three-node beam element with second-order approximation in all
%        directions, St. Venant-Kirchhoff material law
elementID=3403;
n=1; % number of elements

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Case                      %
% 1: standard dyn. pendulum
% 2: circular dyn. pendulum (currently elements 3333 and 3363 only)
% 3: elliptic dyn. pendulum (currently elements 3333 and 3363 only)
% 7: rotating beam (currently elements 3333 and 3363 only)
Case=4;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Integrator options                                               %
% 1 : sundialsTB - needs proper installation of sundials toolbox
%                  https://computation.llnl.gov/casc/sundials/main.html
% 7 : radau IIa ('radau5' by  Ch. Engstler - see function radau5.m)
% 15: ode15s
% 23: ode23s
% 24: ode23t
% 45: ode45
IntegFormu=15;
t_end=5;                            % end time of simulation
t0=0;                                % w.l.o.g. t0=0
deltatspan=0.001;                     % time step (only relevant for output)
tspan=(0:deltatspan:t_end)';
% Tolerances - for both Matlab, radau5 and sundialsTB solvers
RelTol=1e-3; %1e-3
AbsTol=1e-5; %1e-6
% Options for sundials integrator only, see documentation below:
LinearSolver='BiCGStab';
LMM='BDF';
NonlinearSolver='Newton';
MaxOrder=5; % 5
MaxStep=1e-4;
% MaxStep=1e-3;
MaxNumSteps=100;

% from sundialsTB doc:
%
% LinearSolver - Linear solver type [Dense|Diag|Band|GMRES|BiCGStab|TFQMR]
% Specifies the type of linear solver to be used for the Newton nonlinear
% solver (see NonlinearSolver). Valid choices are: Dense (direct, dense
% Jacobian), Band (direct, banded Jacobian), Diag (direct, diagonal Jacobian),
% GMRES (iterative, scaled preconditioned GMRES), BiCGStab (iterative, scaled
% preconditioned stabilized BiCG), TFQMR (iterative, scaled transposen-free QMR).
% The GMRES, BiCGStab, and TFQMR are matrix-free linear solvers.
%
% LMM - Linear Multistep Method [ ’Adams’ | ’BDF’ ]
% This property specifies whether the Adams method is to be used instead
% of the default Backward Differentiation Formulas (BDF) method.
% The Adams method is recommended for non-stiff problems, while BDF is
% recommended for stiff problems.
%
% NonlinearSolver - Type of nonlinear solver used [ Functional | Newton ]
% The ’Functional’ nonlinear solver is best suited for non-stiff
% problems, in conjunction with the ’Adams’ linear multistep method,
% while ’Newton’ is better suited for stiff problems, using the ’BDF’
% method.
%
% MaxNumSteps - Maximum number of steps [positive integer | 500]
% CVode will return with an error after taking MaxNumSteps internal steps
% in its attempt to reach the next output time.
%
% MaxStep - Maximum stepsize [ positive scalar | inf ]
% Defines an upper bound on the integration step size.
%
% MaxOrder - Maximum method order [ 1-12 for Adams, 1-5 for BDF | 5 ]
% Defines an upper bound on the linear multistep method order.
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                      END OF PARAMETER DEFINITIONS                     %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Element definitions
if elementID==3403
    Element='34X3c';
    elemDim=3;
    ElemNodes=4;
    DofsAtNode=30;
    ElemDofs=ElemNodes*DofsAtNode;      % degrees of freedom per element
else
    error('Please choose a valid element type. Error in script ''run_Dynam''.')
end

% Geometry parameters etc.
% Element assembly
if Case == 1   	% RADIAL EXPANSION ROTOR EXAMPLE
    CrossSec=2;
    L=2;
    rho=7850;
    E=210e9;
    nu=0.3;
    Wo=0.2;
    Wi1=0;
    Wi2=0;
    ks=1;
    G=E/(2*(1+nu));
    lambda=E*nu/((1+nu)*(1-2*nu));
    g=0;
elseif Case == 2
    CrossSec=5;     	% hollow circular cross section
    L=7.9;				% length
    Wo=2188e-3;          	% outer diameter
    Wi1=Wo-2*(30e-3+2e-3);	% inner diameter
    Wi2=Wo-2*(30e-3-2e-3);
%     Wi1=0;
%     Wi2=0;
    g=0;				% gravity
    rho=7000;
    E=150e9;			% Young's modulus
    nu=0.33;				% Poisson's ratio
    ks=1;				% shear coefficient
    G=E/(2*(1+nu));		% shear modulus
    lambda=E*nu/((1+nu)*(1-2*nu));	% Lamé constant
elseif Case == 3   	% STATIC
    CrossSec=3;
    L=2;
    rho=7850;
    E=210e9;
    nu=0.3;
    Wo=0.2;
    Wi1=0.1;
    Wi2=0;
    ks=1;
    G=E/(2*(1+nu));
    lambda=E*nu/((1+nu)*(1-2*nu));
    g=0;
elseif Case == 4   	% solid elliptic
    CrossSec=3;
    L=2;
    rho=7850;
    E=210e9;
    nu=0.3;
    Wo=0.2;
    Wi1=0.15;
    Wi2=0;
    ks=1;
    G=E/(2*(1+nu));
    lambda=E*nu/((1+nu)*(1-2*nu));
    g=0;
else
    error('Please select a valid case. Error in ''run_Dynam.m''.')
end

% Element assembly
if elementID==3403
    [P0,nloc] = linemesh3nodes_v10_ANCF_3403(n,L);
    xloc=xlocAllANCF_3403(nloc);
end

[nn,~] = size(P0);
nx = DofsAtNode*nn+0;  % number of degrees of freedom of all elements,
midnode = ceil(nn/2);
% constraints not yet eliminated

% nodal coordinates are formed over the entire structure, from P0
for jj=1:nn,
    ee0((jj-1)*DofsAtNode+1:(jj-1)*DofsAtNode+DofsAtNode)=P0(jj,:);
end
ee0=ee0(:);
% ee0(end+1:end+3) = [0;0.05;0];

% Boundary conditions initialization
bc = true(1,nx);
if elementID==3403
    bc(xlocANCF_3403(1,1:30))=0;
    if Case ~= 3
        for j=2:nn-1
            bc(xlocANCF_3403(j,2:3))=0;
        end
        bc(xlocANCF_3403(nn,2:30))=0;
%       bc(xlocANCF_3403(nn,2:3))=0;
    end
end
ndof = sum(bc);      % number of degrees of freedom, with eliminated constraints

ee0dot=zeros(nx,1);  % Initial velocity vector, all zeros in our examples,
y0=zeros(2*ndof,1);  % y contains all the positions and velocities
y0(1:ndof)=ee0(bc);  % initial positions
y0(ndof+1:2*ndof)=ee0dot(bc);        % initial velocities

% OPTIONS FOR MATLAB INTEGRATORS:
% MassMat=Mg([],nloc,nl,nx,elementID,rho,H,W,L);
% MassMatc=MassMat(bc,bc);
% MassMatcs=MassMatc;
% MasscS=sparse([eye(ndof),zeros(ndof,ndof);zeros(ndof,ndof),MassMatcs]);
% options = odeset('AbsTol',AbsTol,'RelTol',RelTol,'OutputFcn',@ OutputFcn,'Mass',MasscS);
if IntegFormu == 7 % radau5
    options = odeset('AbsTol',AbsTol,'RelTol',RelTol,'OutputFcn',@ OutputFcnWrapper);
    options.Complex = 'off';
    options.Index=[2*ndof; 0; 0];
else
    options = odeset('AbsTol',AbsTol,'RelTol',RelTol,'OutputFcn',@ OutputFcn,'Stats','on');
end

% Compute mass matrix
if CrossSec == 2
    data.M = Mg([],nloc,n,nx,elementID,rho,Wo,Wo,L,2);
elseif CrossSec == 3
    data.M = Mg([],nloc,n,nx,elementID,rho,Wi1,Wo,L,3);
elseif CrossSec == 5
    data.M = Mg([],nloc,n,nx,elementID,rho,Wo,Wo,L,2) - Mg([],nloc,n,nx,elementID,rho,Wi1,Wi2,L,3);
else
    error('Choose proper cross-section')
end

% Use struct for simple parameter passing
data.xloc=xloc;
data.nloc=nloc;
data.n=n;
data.nn=nn;
data.nx=nx;
data.ndof=ndof;
data.Element=Element;
data.ElemDofs=ElemDofs;
data.E=E;
data.G=G;
data.nu=nu;
data.lambda=lambda;
data.ks=ks;
data.rho=rho;
data.g=g;
data.Wi1 = Wi1;
data.Wi2 = Wi2;
data.Wo = Wo;
data.L=L;
data.bc=bc;
data.elementID=elementID;
data.CrossSec=CrossSec;
data.Case=Case;
data.resolution=resolution;
data.elemDim=elemDim;
data.DofsAtNode=DofsAtNode;
data.ee0=ee0;
data.IntegFormu=IntegFormu;
data.Mc = data.M(bc,bc);
data.MInvSparse=sparse(inv(data.Mc));
data.midnode = midnode;
data.ee0 = ee0;
data.ee0dot = ee0dot;

if Case == 3
    fsopts = optimoptions('fsolve','MaxFunctionEvaluations',1e5);
    fstatic = @(y) eom_Dynam(0,y,data);
    ystatic = fsolve(fstatic,y0,fsopts);
    ystatic(1:ndof)
    return
end

tic;
% profile on

% BEGIN SOLUTION PROCESS
display(['Beginning with time integration, t = 0, t_end = ',num2str(t_end),' s'])
tdisp=0;
if IntegFormu == 45 % ode45
    display('Integrator: ode45 (Matlab), options:')
    options
    [t,y]=ode45(@eom_Dynam,tspan,y0,options,data);
elseif IntegFormu == 15
    display('Integrator: ode15s (Matlab), options:')
    options
    [t,y]=ode15s(@eom_Dynam,tspan,y0,options,data);
elseif IntegFormu == 23.
    display('Integrator: ode23s (Matlab), options:')
    options
    [t,y]=ode23s(@eom_Dynam,tspan,y0,options,data);
elseif IntegFormu == 24
    display('Integrator: ode23t (Matlab), options:')
    options
    [t,y]=ode23t(@eom_Dynam,tspan,y0,options,data);
elseif IntegFormu == 7 % radau II a
    display('Integrator: radau5 (Ch. Engstler, based on E. Hairer and G. Wanner [Fortran]), options:')
    if elementID==3363
        warning('off','MATLAB:nearlySingularMatrix')
        warning('warning ''MATLAB:nearlySingularMatrix'' turned off!')
    end
    options
    [t,y]=radau5('odefile',tspan,y0,options,data);
else
    disp('Please select a valid integrator.  Error in file ''run_Dynam''.');
end
CPUtime=toc;
disp(['CPU-time: ' num2str(CPUtime)]);
% profile report


% POSTPROCESSING
% Recover solution for full nodal coordinate vector
eedat=zeros(length(t),nx);

for ii=1:length(t)
    eedat(ii,bc)=ee0(bc)';
end

eedat(:,bc)=y(:,1:ndof);

% Recover solution for time derivative of the nodal coordinate vector
eedotdat=zeros(length(t),nx);
eedotdat(:,bc)=y(:,ndof+1:2*ndof);

reslocy=zeros(length(tspan),3);
reslocz=zeros(length(tspan),3);
refdata=zeros(length(tspan),nx);
diffvec=zeros(length(tspan),nx);

for i=1:length(tspan)
    [a,ad,~] = motion_case10(tspan(i));

    q5a=cos(a);
    q6a=sin(a);
    q8a=-q6a;
    q9a=q5a;

    q5ad=-ad*sin(a);
    q6ad=ad*cos(a);
    q8ad=-q6ad;
    q9ad=q5ad;

    eedat(i,5) = q5a;
    eedat(i,6) = q6a;
    eedat(i,8) = q8a;
    eedat(i,9) = q9a;
    eedotdat(i,5) = q5ad;
    eedotdat(i,6) = q6ad;
    eedotdat(i,8) = q8ad;
    eedotdat(i,9) = q9ad;

    eedat(i,xlocANCF_3403(nn,5)) = q5a;
    eedat(i,xlocANCF_3403(nn,6)) = q6a;
    eedat(i,xlocANCF_3403(nn,8)) = q8a;
    eedat(i,xlocANCF_3403(nn,9)) = q9a;
    eedotdat(i,xlocANCF_3403(nn,5)) = q5ad;
    eedotdat(i,xlocANCF_3403(nn,6)) = q6ad;
    eedotdat(i,xlocANCF_3403(nn,8)) = q8ad;
    eedotdat(i,xlocANCF_3403(nn,9)) = q9ad;

    refdata(i,:) = ee0;
    for j=1:nn
        refdata(i,xlocANCF_3403(j,5:6)) = [q5a,q6a];
        refdata(i,xlocANCF_3403(j,8:9)) = [q8a,q9a];
    end

    reslocy(i,:) = getLoc3403(eedat(i,xlocANCF_3403(1,1:120))',0,1,0,data);
    reslocz(i,:) = getLoc3403(eedat(i,xlocANCF_3403(1,1:120))',0,0,1,data);
    rynorm(i) = norm(reslocy(i,2:3));
    rznorm(i) = norm(reslocz(i,2:3));

    diffvec(i,:) = eedat(i,:)-refdata(i,:);
end

figure
plot(tspan,rynorm)
figure
plot(tspan,rznorm)

figure
for i=1:n
    drawee = refdata(end,(i-1)*90+1:i*90+30)+diffvec(end,(i-1)*90+1:i*90+30);
    hold on
    draw34X3(drawee',data)
end