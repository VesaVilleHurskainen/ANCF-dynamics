function [tout,yout,stats] = radau5(odefile,tspan,y0,options,varargin)
% function [tout,yout,stats] = radau5(odefile,tspan,y0,options,[varargin])
%
% RADAU5 Solve stiff differential equations, fifth order method.
%   [T,Y] = RADAU5('F',TSPAN,Y0) with TSPAN = [T0 TFINAL] integrates the
%   system of first order differential equations y' = F(t,y) from time T0 to
%   TFINAL with initial conditions Y0.  'F' is a string containing the name
%   of an ODE file.  Function F(T,Y) must return a column vector.  Each row
%   in solution array Y corresponds to a time returned in column vector T.
%   To obtain solutions at specific times T0, T1, ..., TFINAL (all
%   increasing or all decreasing), use TSPAN = [T0 T1 ... TFINAL].
%	
%   [T,Y] = RADAU5('F',TSPAN,Y0,OPTIONS) solves as above with default
%   integration parameters replaced by values in OPTIONS, an argument
%   created with the ODESET function.  See ODESET for details.  Commonly
%   used options are scalar relative error tolerance 'RelTol' (1e-3 by
%   default) and vector of absolute error tolerances 'AbsTol' (all
%   components 1e-6 by default).
%   
%   It is possible to specify TSPAN, Y0 and OPTIONS in the ODE file (see
%   ODEFILE).  If TSPAN or Y0 is empty, then RADAU5 calls the ODE file
%   [TSPAN,Y0,OPTIONS] = F([],[],'init') to obtain any values not supplied
%   in the RADAU5 argument list.  Empty arguments at the end of the call
%   list may be omitted, e.g. RADAU5('F').
%   
%   The Jacobian matrix dF/dy can be critical to reliability and efficiency.
%   Use ODESET to set JConstant 'on' if dF/dy is constant.  Set Vectorized
%   'on' if the ODE file is coded so that F(T,[Y1 Y2 ...]) returns [F(T,Y1)
%   F(T,Y2) ...].  Set JPattern 'on' if dF/dy is a sparse matrix and the ODE
%   file is coded so that F([],[],'jpattern') returns a sparsity pattern
%   matrix of 1's and 0's showing the nonzeros of dF/dy.  Set Jacobian 'on'
%   if the ODE file is coded so that F(T,Y,'jacobian') returns dF/dy.
%
%   As an example, the command
%   
%       radau5('vdpode',[0 3000],[2 0],[],1000);
%   
%   solves the system y' = vdpode(t,y) with mu = 1000, using the default
%   relative error tolerance 1e-3 and the default absolute tolerance of 1e-6
%   for each component.  When called with no output arguments, as in this
%   example, RADAU5 calls the default output function ODEPLOT to plot the
%   solution as it is computed.
%	
%   RADAU5 also solves problems M*y' = F(t,y) with a constant mass matrix M
%   that is nonsingular and (usually) sparse.  Use ODESET to set Mass 'on'
%   if the ODE file is coded so that F([],[],'mass') returns M (see
%   FEM2ODE).  Use ODE15S if M is time-dependent.
%  
%   In addition to the options provided by ODESET one can set the
%   following parameters to optimize the performance of the code.
%   WARNING: We have to set these _after_ the last call to 'odeset'
%   since ODESET deletes unrecognized options. :-(
%
%    options.Complex  [  on  | {off} ]  
%       Solution is complex
%    options.Predict  [ {on} |  off  ]  
%       Starting values of newton iterations
%    options.Gustaf   [ {on} |  off  ]    
%       Switch for step size controller
%    options.Colmmd   [ {on} |  off  ]
%       Use column minimum degree permutation when solving sparse systems.
%    options.Index    [ 3 x 1 vector of integers ] {[neq; 0; 0]}
%       Dimension of the index 1,2,3 variables
%    options.MaxNewt  [ integer {7}  ]  
%       Maximum number of Newton iterations for solving the nonlinear systems.  
%    options.Debug    [ integer {0}  ]  
%        = 0   No print output
%        = 1   Short integration monitor
%        = 2   Note stepsize reductions
%        = 3   Explain stepsize reductions
%        = 4   Error estimates (embedded method)
%        = 5   Short convergence monitor (Newton iteration)
%        = 10  Error estimates (Newton iteration)
%
%   See also ODEFILE and
%       other ODE solvers:   ODE15S, ODE45, ODE23, ODE113
%       options handling:    ODESET, ODEGET
%       output functions:    ODEPLOT, ODEPHAS2, ODEPHAS3, ODEPRINT
%       odefile examples:    VDPODE, BRUSSODE, B5ODE, CHM6ODE, FEM2ODE
%       Jacobian functions:  NUMJAC, COLGROUP
%
%   RADAU5 is an implementation of a fifth order implicit Runge-Kutta method
%   with 3 stages (RADAU IIA), stepsize control and continuous output. By
%   default, Jacobians are generated numerically.  

%   Details are to be found in the book
%
%      E. Hairer and G. Wanner, Solving Ordinary Differential
%      Equations II. Stiff and Differential-Algebraic Problems.
%      Springer Series in Computational Mathematics 14,
%      Springer-Verlag 1991, second edition 1996.
%

%   This is essentially a port of the FORTAN code RADAU5 to MATLAB
%   Authors of the FORTRAN code: 
%
%      E. Hairer and G. Wanner
%      Universite de Geneve, Dept. de Mathematiques
%      CH-1211 Geneve 24, SWITZERLAND 
%      email:  hairer@divsun.unige.ch, wanner@divsun.unige.ch
%

%   MATLAB code by
%
%      Ch. Engstler
%      Universitaet Tuebingen, Mathematisches Institut
%      Auf der Morgenstelle 10
%      D-72074 Tuebingen, GERMANY
%      email: engstler@na.uni-tuebingen.de

% version 1.03, Oct 12 1999

% GK: 07.12.10
% GK: Index 2,3 Problematik korrigiert: Da options.Index nicht existiert
% GK:                                   wird mit global g_nind gearbeitet

% 01.12.2015: Modified by Henrik Ebel in order to make code usable with
% simple dynamic ANCF simulation code - now allows for output fcns that
% need additional parameters passed in varargin (like with OutputFcns used 
% for built-in Matlab integrators). 
% Also: removed nmax for long calculations.

global g_nind;         % GK: Eingefuegt

true = logical(1);
false = ~true;

nsteps = 0;
nfailed = 0;
naccpt = 0;
nfevals = 0;
npds = 0;
ndecomps = 0;
nsolves = 0;

if nargin == 0
  error('Not enough input arguments.  See RADAU5.');
elseif ~isstr(odefile) & ~isa(odefile, 'inline')
  error('First argument must be a single-quoted string.  See RADAU5.');
end

if nargin == 1
  tspan = []; y0 = []; options = [];
elseif nargin == 2
  y0 = []; options = [];
elseif nargin == 3
  options = [];
end

% Get default tspan and y0 from odefile if none are specified.
if isempty(tspan) | isempty(y0)
  if (nargout(odefile) < 3) & (nargout(odefile) ~= -1)
    msg = sprintf('Use radau5(''%s'',tspan,y0,...) instead.',odefile);
    error(['No default parameters in ' upper(odefile) '.  ' msg]);
  end
  [def_tspan,def_y0,def_options] = feval(odefile,[],[],'init',varargin{:});
  if isempty(tspan)
    tspan = def_tspan;
  end
  if isempty(y0)
    y0 = def_y0;
  end
  if isempty(options)
    options = def_options;
  else
    options = odeset(def_options,options);
  end
end


% Test that tspan is internally consistent.
tspan = tspan(:);
ntspan = length(tspan);
if ntspan == 1
  t = 0;
  next = 1;
else
  t = tspan(1);
  next = 2;
end
tfinal = tspan(ntspan);
if t == tfinal
  error('The last entry in tspan must be different from the first entry.');
end
tdir = sign(tfinal - t);
if any(tdir * (tspan(2:ntspan) - tspan(1:ntspan-1)) <= 0)
  error('The entries in tspan must strictly increase or decrease.');
end

y = y0(:);
neq = length(y);

sqrtneq = sqrt(neq);
sqrt3neq = sqrt(3*neq);
ns = 3*neq;

% Get options, and set defaults.
rtol = odeget(options,'RelTol',1e-3);
atol = odeget(options,'AbsTol',1e-6);
atol = atol(:);
if (rtol < 0) | any(atol <= 0)  
  error('RelTol must be nonnegative and all entries of AbsTol must be positive.');
end
if ( rtol == 0 )
  fprintf('WARNING: Curious input for RelTol=0.\n');
end
if length(atol) == 1
  atol = atol + zeros(neq,1);
elseif length(atol) ~= neq
  msg = sprintf('an AbsTol vector of length 1 or %d.',neq);
  error(['Solving ' upper(odefile) ' requires ' msg]);
end
% change tolerances
expm=2/3;
if ( rtol > 0 )
  quot = atol/rtol;
  rtol = 0.1*rtol^expm;
  atol = rtol.*quot;
else
  atol = 0.1*atol.^expm;
end

% By default, hmax is 1/2 of the interval.
hmax = min(abs(tfinal-t), abs(odeget(options,'MaxStep',0.5*(tfinal-t))));
if hmax <= 0
  error('Option ''MaxStep'' must be greater than zero.');
end
htry = abs(odeget(options,'InitialStep',1e-2));
if htry <= 0
  error('Option ''InitialStep'' must be greater than zero.');
end
h = min(htry,hmax);

hmin = 16*eps*(abs(t)+1);
hmin = min(hmin, hmax);

if nargout > 0
  outfun = odeget(options,'OutputFcn');
else
  outfun = odeget(options,'OutputFcn','odeplot');
end
if isempty(outfun)
  haveoutfun = false;
else
  haveoutfun = true;
  outputs = odeget(options,'OutputSel',1:neq);
end
haveeventfun = false;

refine = odeget(options,'Refine',1);
printstats = strcmp(odeget(options,'Stats','off'),'on');

% information about the Jacobian
vectorized = strcmp(odeget(options,'Vectorized','off'),'on');
Janalytic = strcmp(odeget(options,'Jacobian','off'),'on');
JConstant = strcmp(odeget(options,'JConstant','off'),'on');
Jpattern = strcmp(odeget(options,'JPattern','off'),'on');
if Jpattern
  Js = feval(odefile,[],[],'jpattern',varargin{:});
else
  Js = [];
end

mass = strcmp(odeget(options,'Mass','off'),'on');
if mass
  havemass = true;
  M = feval(odefile,t,[],'mass',varargin{:});
  Mconstant = strcmp(odeget(options,'MassConstant','on'),'on');
  if ~Mconstant
    error('For a non-constant mass matrix, M(t)*y'', use ODE15S.');
  end
else
  havemass = false;
  M = sparse((1:neq)',(1:neq)',1,neq,neq); % sparse Identity matrix
  Mconstant = true;
end

% Set the output flag.
if (ntspan > 2) || ((ntspan == 2) && (refine == 0)) % MODIFIED TO MODERN MATLAB (HENRIK EBEL)
  outflag = 1;                          % output only at tspan points
elseif refine <= 1
  outflag = 2;                          % computed points, no refinement
else
  outflag = 3;                          % computed points, with refinement
  S = (1:refine-1)' / refine;
end

% Allocate memory if we're generating output.
if nargout > 0
  if ntspan > 2                         % output only at tspan points
    tout = zeros(ntspan,1);
    yout = zeros(ntspan,neq);
  else                                  % alloc in chunks
    chunk = max(ceil(128 / neq),refine);
    tout = zeros(chunk,1);
    yout = zeros(chunk,neq);
  end
  nout = 1;
  tout(nout) = t;
  yout(nout,:) = y.';
end

% Initialize the output function.
if haveoutfun
  feval(outfun,[t tfinal],y(outputs),'init',varargin{:});
end

%
% Additional options.
%
if isempty(options)
  options.dummy = 0;
end
if isfield(options,'Complex')
  if strcmp('on',options.Complex)
    isreal = false;
    fprintf('------------------------------------\n');
    fprintf('problem requires complex arithmetic.\n');
    fprintf('------------------------------------\n');
  else
    isreal = true;
  end
else
  isreal = true;
  warning('Solution is supposed to be real ... !')
  fprintf('You can inhibit this warning message by setting ')
  fprintf('options.Complex = {on|off}\n');
end
% DEBUG
if isfield(options,'Debug')
  DEBUG = options.Debug;
else
  DEBUG = 0;
end
% Predict: switch for starting values of newton iterations.
if isfield(options,'Predict')
  if strcmp('on',options.Predict)
    start_n = false;
  elseif strcmp('off',options.Predict)
    start_n = true;
  else
    fprintf('WARNING: Curious input for Predict=''%s''.\n',options.Predict);
    start_n = false;
  end
else
  start_n = false; % default
end
% MaxNewt: maximum number of Newton iterations
if isfield(options,'MaxNewt')
  if options.MaxNewt > 0
    max_newt = options.MaxNewt;
  elseif options.MaxNewt == 0
    max_newt = 7;
  else
    fprintf('WARNING: Curious input for MaxNewt=''%d''.\n',options.MaxNewt);
    max_newt = 7;
  end
else
  max_newt = 7; % default
end
% Index: dimension of the index 1,2,3 variables
if isfield(options,'Index')
  if length(options.Index)==3
    nind = options.Index;
    if sum(nind) ~= neq,
      fprintf('curious input for Index: [ %d %d %d ].\n',nind);
      error ('can''t recover from previous errors.' );
    end
  else
    fprintf('WARNING: Curious input for Index.\n');
    nind = [neq; 0; 0];
  end
else
  % nind = [neq; 0; 0]; % default      GK: auskommentiert :GK
  nind = g_nind;        %              GK: Eingefuegt     :GK
end
zero = 0; one = 1; two = 2;
index = [zero(ones(size(1:nind(1)))) one(ones(size(1:nind(2)))) ...
    two(ones(size(1:nind(3))))]';
% Gustaf: switch for step size controller
if isfield(options,'Gustaf')
  if strcmp('on',options.Gustaf)
    Pred = true;
  elseif strcmp('off',options.Gustaf)
    Pred = false;
  else
    fprintf('WARNING: Curious input for Gustaf=''%s''.\n',options.Gustaf);
    Pred = true;
  end
else
  Pred = true; % default
end
% Colmmd 
if isfield(options,'Colmmd')
  if strcmp('on',options.Colmmd)
    Colmmd = true;
  elseif strcmp('off',options.Colmmd)
    Colmmd = false;
  else
    fprintf('WARNING: Curious input for Colmmmd=''%s''.\n',options.Colmmd);
    Colmmd = true;
  end
else
  Colmmd = true; % default
end
%
% Set method parameters.
%
nmax = 100000;

% work(1) = eps, the smallest number satisfying 1+uround > 1 obsolete.
% work(2) = safe the safety factor in step size prediction,
safe = 0.9;
% work(3) = thet: decides whether the jacobian should be recomputed. 
%   increase thet, to 0.1 say, when jacobian evaluations are costly.
%   for small systems thet should be smaller (0.001 say). negativ thet
%   forces the code to compute the jacobian after every accepted step.
%   default 1e-3.
thet = 1e-3;
% work(4) = fnewt: stopping crierion for newton's method, usually chosen <1.
%   smaller values of work(4) make the code slower, but safer.
%   default min(0.03d0,rtol(1)**0.5d0)
fnewt = 0;  % option in a future version.
if (rtol > 0)
  tolst = rtol;
else
  tolst = min(atol);
end
if ( fnewt == 0 ) 
   fnewt = max(10*eps/tolst,min(0.03,tolst^0.5));
else
   if (fnewt <= eps/tolst)
      error([ 'Curious value for fnewt: ' sprintf('%e',fnewt) ]);
   end
end
% work(5:6) = [quot1; quot2], parameters for step size selection
%   if quot1 < hnew/hold < quot2, then the step size is not changed.
%   this saves, together with a large work(3), LU-decompositions and
%   computing time for large systems. for small systems one may have
%   quot1 = 1, quot2=1.2, for large full systems quot1 = 0.99, quot2 =
%   2 might be good. defaults quot1 = 1, quot2 = 1.2 .
quot1 = 1;
quot2 = 1.2;
% work(7:8) = [ fac_lo; fac_hi], parameters for step size selection
%   the new step size is chosen subject to the restriction work(8) <=
%   hnew/hold <= work(9) default values: work(8)=0.2d0, work(9)=8.d0
fac_lo = 5;
fac_hi = 1/8;
%  
% Define the RK coefficients. We have 3 stages.
%
% real eigenvalue of inv(A)
gamma0 = (6+81^(1/3)-9^(1/3))/30;
gamma = 1/gamma0;
% conjugate complex eigenvalue pair of inv(A): alpha +- i beta
alpha = (12-81^(1/3)+9^(1/3))/60;
beta  = (81^(1/3)+9^(1/3))*sqrt(3)/60;
cno = alpha^2+beta^2;
alpha = alpha/cno;
beta = beta/cno;

sqrt6 = sqrt(6);
if isreal
  % inv(T)*inv(A)*T = [ gamma 0 0; 0 alpha -beta; 0 beta alpha]
  % TI = inv(T);

  T11 = 9.1232394870892942792d-02;
  T12 = -0.14125529502095420843d0;
  T13 = -3.0029194105147424492d-02;
  T21 = 0.24171793270710701896d0;
  T22 = 0.20412935229379993199d0;
  T23 = 0.38294211275726193779d0;
  T31 = 0.96604818261509293619d0;
  T32 = 1;
  % 
  TI11 = 4.3255798900631553510d0;
  TI12 = 0.33919925181580986954d0;
  TI13 = 0.54177053993587487119d0;
  TI21 = -4.1787185915519047273d0;
  TI22 = -0.32768282076106238708d0;
  TI23 = 0.47662355450055045196d0;
  TI31 = -0.50287263494578687595d0;
  TI32 = 2.5719269498556054292d0;
  TI33 = -0.59603920482822492497d0;
else
  A=zeros(3,3);
  A(1,1) = (88-7*sqrt6)/360;
  A(2,1) = (296+169*sqrt6)/1800;
  A(3,1) = (16 - sqrt6)/36;
  A(1,2) = (296-169*sqrt6)/1800;
  A(2,2) = (88+7*sqrt6)/360;
  A(3,2) = (16 + sqrt6)/36;
  A(1,3) = (-2+3*sqrt6)/225;
  A(2,3) = (-2-3*sqrt6)/225;
  A(3,3) = 1/9;
  [T,D] = eig(A);
  T = T(:,[3,2,1]);
  TI = T\eye(3);
  temp = alpha;
  alpha = temp+1i*beta; % modified, HE
  beta = temp-1i*beta; % modified, HE
  for i1 = 1:3,
    for i2 = 1:3,
      eval(sprintf('T%d%d = T(%d,%d);',i1,i2,i1,i2));
      eval(sprintf('TI%d%d = TI(%d,%d);',i1,i2,i1,i2));
    end
  end
end

c = zeros(3,1);
c(1) = (4-sqrt6)/10;
c(2) = (4+sqrt6)/10;
c(3) = 1;
c1 = c(1);
c2 = c(2);
c3 = 1;
c1m1 = c1-1;
c2m1 = c2-1;
c1mc2 = c1-c2;

% e is needed for error estimation
% we leave out the multiplication by gamma_0 [HW] since this will simplify
% the computation of the error later.
e = zeros(1,3);
e(1) = -13-7*sqrt6;
e(2) = -13+7*sqrt6;
e(3) = -1;
e = e'/3;


% Initialize
cfac = safe*(2*max_newt+1);
fac_con = 1;

posneg = sign(tfinal-t);
hmaxn = min(abs(hmax),abs(tfinal-t));
if (abs(h) <= 10*eps), h = 1.d-6; end
h = min(abs(h),hmaxn);
h = h*posneg;
h_old = h;
h_fac = h;

reject = false;
first = true;

if ( (t+h*1.0001-tfinal)*posneg >= 0 ),
  h = tfinal-t;
  last = true;
else
  last = false;
end
hhfac = h;

scal_0 = rtol*abs(y)+atol;
scal = scal_0./(h.^index);
cont = zeros(neq,3);
f = zeros(neq,3);

% The input arguments of odefile determine the args to use to evaluate f.
if nargin(odefile) == 2
  args = {};                            % odefile accepts only (t,y)
else
  args = [{''} varargin];               % use (t,y,'',p1,p2,...)
end

f0 = feval(odefile,t,y,args{:});

[m,n] = size(f0);
if n > 1
  error(['Function ' odefile '(t,y) must return a column vector.'])
elseif m ~= neq
  error(['Vector ' odefile '(t,y) must be same length as initial conditions.']);
end
nfevals = nfevals+1;

% compute the Jacobian at (t0,y0)
if Janalytic
  J = feval(odefile,t,y,'jacobian',varargin{:});
else
  [J,facj,gj,nF] = numjac(odefile,t,y,f0,atol,[],vectorized,Js,[],args{:});
  nfevals = nfevals + nF; 		% stats
end

% probably user-supplied in a future version.
perm = [];

SparseMMD = issparse(J) & Colmmd;
if SparseMMD & isempty(perm) 
  perm = colmmd(M/h-J);
end

npds = npds + 1; 			% stats
needNewJ = false;
needNewLU = true;
consth = false;
tooslow = false;
itfail = false;
JCurrent = true;
newt = 0;
err_print = 0;
% --------------------------------------------------------------------
%
%    *argl*
% 
ones1s = ones(3,1);
ones1n = ones(neq,1);
% --------------------------------------------------------------------
radinfo(nind,Pred,SparseMMD,start_n,fnewt);

% while ( nsteps <= nmax & t < tfinal) % MODIFIED - REMOVED nsteps
while ( t < tfinal)  
  % LOOP FOR ADVANCING ONE STEP.
  while true
    
    strategy = '??,??';
    if needNewJ, strategy(1:2)='nJ'; else strategy(1:2)='kJ'; end
    if needNewLU, strategy(4:5)='nD'; else strategy(4:5)='kD'; end
    
    if ( DEBUG >= 1 ),
      fprintf('radau5: nsteps=%d, t=%7.2e, h=%8.2e, newt=%d, err=%9.2e [%s]\n', nsteps, t, h, newt, err_print, strategy);
    end
  
    if abs(h) <= hmin
      msg = sprintf(['Failure at t=%e.  Unable to meet integration ' ...
	    'tolerances without reducing the step size below ' ...
	    'the smallest value allowed (%e) at time t.\n'], ...
	  t,hmin);
      warning(msg);
      if nargout ~= 0
	tout = tout(1:nout);
	yout = yout(1:nout,:);
	stats = [nsteps; nfailed; nfevals; npds; ndecomps; nsolves];
      end
      return;
    end

    newt = 0;
    fac_con = max(fac_con, eps)^0.8;
    theta = abs(thet);
    newnrm = 1e+30;
    
    % LOOP FOR THE COMPUTATION OF YNEW
    tooslow = false;
    itfail = false;

    while true

      nsteps = nsteps+1;

      % compute the Jacobian
      if needNewJ,
	if JConstant
	  error('internal error: needNewJ & JConstant');
	end
	if Janalytic
	  J = feval(odefile,t,y,'jacobian',varargin{:});
	else
	  f0 = feval(odefile,t,y,args{:});
	  [J,facj,gj,nF] = ...
	      numjac(odefile,t,y,f0,atol,facj,vectorized,Js,gj,args{:});
	  nfevals = nfevals + nF + 1; % stats
	end
	npds = npds + 1; 		% stats
	JCurrent = true;
	needNewJ = false;
	needNewLU = true;
      end
      
      % Compute the iteration matrices and their decomposition
      % E1 = gamma/h*M - J
      % E2 = (alpha+i*beta)/h*M - J
      
      if needNewLU,
	a = alpha/h;
	b = beta/h;
	g = gamma/h;
	B = g*M-J;
	if SparseMMD
	  [L1,U1,P1] = lu( B(:,perm) );
	else
	  [L1,U1,P1] = lu( B );
	end
	if (any(abs(diag(U1)))<eps),
	  fprintf('WARNING: g*M-J is singular.\n');
	end
	if isreal
	  B = (a+i*b)*M-J;
	  if SparseMMD
	    [L2,U2,P2] = lu( B(:,perm) );
	  else
	    [L2,U2,P2] = lu( B );
	  end
	  if (any(abs(diag(U2)))<eps),
	    fprintf('WARNING: (a+i*b)*M-J is singular.\n');
	  end
	else
	  if SparseMMD 
	    B = a*M - J;
	    [L2,U2,P2] = lu ( B(:,perm) );
	    B = b*M - J;
	    [L3,U3,P3] = lu ( B(:,perm) );
	  else
	    [L2,U2,P2] = lu ( a*M-J );
	    [L3,U3,P3] = lu ( b*M-J );
	  end
	end
	ndecomps = ndecomps+1;
	needNewLU = false;
      end % needNewLU,
      
      % compute starting values for the Newton iteration
      if ( first | start_n )
	if ( DEBUG >= 10 ),
	  fprintf('    Newton: assuming zero initial values.\n');
	end
	z = zeros(neq,3);
	w = zeros(neq,3);
      else
	if ( DEBUG >= 10 ),
	  fprintf('    Newton: computing initial values.\n');
	end
	cq = c*(h/h_old);
	z(:,1) = cq(1)*( cont(:,1) + (cq(1)-c2m1)*( cont(:,2) + ...
	    (cq(1)-c1m1)*cont(:,3))); 
	z(:,2) = cq(2)*( cont(:,1) + (cq(2)-c2m1)*( cont(:,2) + ...
	    (cq(2)-c1m1)*cont(:,3))); 
	z(:,3) = cq(3)*( cont(:,1) + (cq(3)-c2m1)*( cont(:,2) + ...
	    (cq(3)-c1m1)*cont(:,3))); 
	w(:,1) = TI11*z(:,1)+TI12*z(:,2)+TI13*z(:,3);
	w(:,2) = TI21*z(:,1)+TI22*z(:,2)+TI23*z(:,3);
	w(:,3) = TI31*z(:,1)+TI32*z(:,2)+TI33*z(:,3);
      end % compute starting values

      itfail = true;
      % LOOP FOR THE SIMPLIFIED NEWTON ITERATION
      for newt = 1:max_newt
      
	% function evaluation
	f(:,1) = feval(odefile,t+c1*h,y+z(:,1),args{:});
	f(:,2) = feval(odefile,t+c2*h,y+z(:,2),args{:});
	f(:,3) = feval(odefile,t+h,   y+z(:,3),args{:});
	nfevals = nfevals + 3;

	z(:,1) = TI11*f(:,1)+TI12*f(:,2)+TI13*f(:,3);
	z(:,2) = TI21*f(:,1)+TI22*f(:,2)+TI23*f(:,3);
	z(:,3) = TI31*f(:,1)+TI32*f(:,2)+TI33*f(:,3);

	% solve the linear systems
      
	if havemass,
	  Mw = M*w;
	else
	  Mw = w;
	end
	if isreal
	  z(:,1) = z(:,1)-g*Mw(:,1);
	  z(:,2) = z(:,2)-a*Mw(:,2)+b*Mw(:,3);
	  z(:,3) = z(:,3)-b*Mw(:,2)-a*Mw(:,3);
	  if SparseMMD
	    z(perm,1) = U1\(L1\(P1*z(:,1)));
	    z(perm,2) = U2\(L2\(P2*(z(:,2)+1i*z(:,3)))); % modified, HE
	  else
	    z(:,1) = U1\(L1\(P1*z(:,1)));
	    z(:,2) = U2\(L2\(P2*(z(:,2)+1i*z(:,3)))); % modified, HE
	  end
	  z(:,3) = imag(z(:,2));
	  z(:,2) = real(z(:,2));
	else
	  z(:,1) = z(:,1)-g*Mw(:,1);
	  z(:,2) = z(:,2)-a*Mw(:,2);
	  z(:,3) = z(:,3)-b*Mw(:,3);
	  if SparseMMD
	    z(perm,1) = U1\(L1\(P1*z(:,1)));
	    z(perm,2) = U2\(L2\(P2*z(:,2)));
	    z(perm,3) = U3\(L3\(P3*z(:,3)));
	  else
	    z(:,1) = U1\(L1\(P1*z(:,1)));
	    z(:,2) = U2\(L2\(P2*z(:,2)));
	    z(:,3) = U3\(L3\(P3*z(:,3)));
	  end
	end
	nsolves = nsolves+1;
      
	% estimate the error in the current iteration step
	scal = scal_0./(h.^index);
	newnrm = norm(z./[scal scal scal],'fro')/sqrt3neq;
	if ( DEBUG >= 10 ),
	  fprintf('    newt = %d, newnrm = %e,  theta = %e\n',...
	      newt,newnrm,theta);
	end
      
	% check for slow convergence
	if ( newt > 1 & newt < max_newt ),
	  
	  thq = newnrm/oldnrm;
	  if ( newt == 2 ),
	    theta = thq;
	  else
	    theta = sqrt(thq*thqold);
	  end
	  thqold = thq;
	  
	  if ( theta < 0.99 ), % convergence at least
	    fac_con = theta/(1-theta);
	    dyth = fac_con*newnrm*theta^(max_newt-1-newt)/fnewt;
	    if ( dyth >= 1 ),  
	      % we can not  expect convergence after max_newt steps. 
	      qnewt = max(1.e-4,min(20,dyth));
	      hhfac = .8*qnewt^(-1/(4+max_newt-1-newt));
	      h = hhfac*h;
	      tooslow = true;
	      itfail = true;
	      break;
	    end
	  else 
	    if ( DEBUG >= 5 ),
	      fprintf ('  no convergence. theta >= 0.99\n');
	    end
	    tooslow = true;
	    itfail = true;
	    break;
	  end % if ( theta < 0.99 )
	
	end % check for slow convergence
      
	% convergence of Newton iteration ok. Perform the step.
	oldnrm = max(newnrm,eps);
	w = w+z;
      
	% z = kron(T,eye(n))*w
	z(:,1) = T11*w(:,1)+T12*w(:,2)+T13*w(:,3);
	z(:,2) = T21*w(:,1)+T22*w(:,2)+T23*w(:,3);
	if isreal
	  z(:,3) = T31*w(:,1)+    w(:,2);
	else
	  z(:,3) = T31*w(:,1)+T32*w(:,2)+T33*w(:,3);
	end    
    
	if (fac_con*newnrm <= fnewt),
	  itfail = false;
	  break
	end
    
      end % for newt = 1:max_newt
      
      if ~itfail,
	break
      else
	reject = true;
	last = false;
	if ( theta >= 0.99 )
	  if ( DEBUG >= 3 ),
	    fprintf ('  radau5: unexpected step rejection.\n');
	  end
	  h = h/2;
	  hhfac = 0.5;
	else
	  if ( DEBUG >= 3 ),
	    fprintf('  radau5: convergence too slow.\n');
	  end
	end
	if ( DEBUG >= 2 ),
	  fprintf('  radau5: step aborted. restart with h = %e\n', h);
	end
	if ( JCurrent ),
	  needNewJ = false;
	  needNewLU = true;
	else
	  needNewJ = ~JConstant;
	  needNewLU = true;
	end
      end
    
    end	% loop for the simplified Newton iteration

    % At this point the Newton iteration converged to a solution. 
    % Our next task is to estimate the local error.
    if ( DEBUG >= 5 ),
      fprintf('  success after %d steps, theta = %e\n', newt, theta);
    end
      
    % estimate the local error of the Runge-Kutta method: [HW], p.134
    temp = z*(e/h);
    if havemass,
      temp = M*temp;
    end
    scal = scal_0./(h.^index);
    if SparseMMD
      err_v(perm,1) = U1\(L1\(P1*(f0+temp)));
    else
      err_v = U1\(L1\(P1*(f0+temp)));
    end
    err = norm(err_v./scal,2);
    err = max(err/sqrtneq,1e-10);
    if ( DEBUG >= 4 ),
      fprintf('  radau5: first  error estimate: err = %e\n', err);
    end
    if ( (err >= 1) & (first | reject) ),
      nfevals=nfevals+1;
      err_v = feval(odefile,t,y+err_v,args{:});
      if SparseMMD
	err_v(perm,1) = U1\(L1\(P1*(err_v+temp)));
      else
	err_v = U1\(L1\(P1*(err_v+temp)));
      end
      err = norm(err_v./scal,2);
      err = max(err/sqrtneq,1e-10);
      if ( DEBUG >= 4 ),
	fprintf('  radau5: second error estimate: err = %e\n', err);
      end
    end
    err_print = err;

    % step size selection
    % we require 1/fac_lo <= hnew/h <= 1/fac_hi
    fac = min(safe, cfac/(2*max_newt+newt));
    quot = max(fac_hi,min(fac_lo,(err^0.25)/fac));
    hnew = h/quot;
      
    if ( err <= 1), break, end
      
    reject = true;
    last = false;
    if ( naccpt >= 1 ), nfailed = nfailed+1; end
	
    if first
      h = h/10;
      hhfac = 0.1;
    else
      hhfac = hnew/h;
      h = hnew;
    end
    
    if ( DEBUG >= 2 ),
      fprintf('  radau5: step rejected. restart with h = %e\n', h);
    end

    if JCurrent,
      needNewJ = false;
      needNewLU = true;
    else
      needNewJ = ~JConstant;
      needNewLU = true;
    end
	
    if haveoutfun
      if feval(outfun,t,y(outputs),'failed',varargin{:}) == 1
	if nargout ~= 0
	  tout = tout(1:nout);
	  yout = yout(1:nout,:);
	  stats = [nsteps; nfailed; nfevals; npds; ndecomps; nsolves];
	end
	return;
      end
    end
      
  end  % basic integration step
  
  %
  % step is accepted
  %
  
  first = false;
  naccpt = naccpt+1;
  
  if (Pred)
    if (naccpt > 1)
      facgus = (hacc/h)*(err^2/erracc)^0.25d0/safe;
      facgus = max(fac_hi,min(fac_lo,facgus));
      quot = max(quot,facgus);
      hnew = h/quot;
    end
    hacc = h;
    erracc = max(1.0d-2,err);
  end

  told = t;
  h_old = h;
  t = t+h;
  y = y+z(:,3);

  % collocation polynamial 
  
  cont(:,3) = z(:,1) / c1;
  cont(:,2) = ( z(:,1) - z(:,2) ) / c1mc2;
  cont(:,1) = ( z(:,2) - z(:,3) ) / c2m1;
  cont(:,3) = ( cont(:,2) - cont(:,3) ) / c2;
  cont(:,2) = ( cont(:,2) - cont(:,1) ) / c1m1;
  cont(:,3) = cont(:,2)-cont(:,3);
  
  tnew = t;
  
  if nargout > 0
      oldnout = nout;
      if outflag == 2 			% computed points, no refinement
          nout = nout + 1;
          if nout > length(tout)
              tout = [tout; zeros(chunk,1)];
              yout = [yout; zeros(chunk,neq)];
          end
          tout(nout) = tnew;
          yout(nout,:) = y.';
      elseif outflag == 3 		% computed points, with refinement
          nout = nout + refine;
          if nout > length(tout)
              tout = [tout; zeros(chunk,1)]; 	% requires chunk >= refine
              yout = [yout; zeros(chunk,neq)];
          end
          ii = oldnout+1:nout-1;
          tinterp = told+h*S;
          yinterp = ntrprad(tinterp,[],[],t,y,h,cont);
          tout(ii) = tinterp;
          yout(ii,:) = yinterp.';
          tout(nout) = tnew;
          yout(nout,:) = y.';
      elseif outflag == 1 		% output only at tspan points
          while next <= ntspan
              if tdir * (tnew - tspan(next)) < 0
                  if haveeventfun && done
                      nout = nout + 1;
                      tout(nout) = tnew;
                      yout(nout,:) = y.';
                  end
                  break;
              elseif tnew == tspan(next)
                  nout = nout + 1;
                  tout(nout) = tnew;
                  yout(nout,:) = y.';
                  next = next + 1;
                  break;
              end
              nout = nout + 1;                % tout and yout are already allocated
              tout(nout) = tspan(next);
              s = (tspan(next) - tnew) / h;
              yout(nout,:) = ntrprad(tspan(next),[],[],t,y,h,cont).';
              next = next + 1;
          end
      end
      
      if haveoutfun
          if outflag == 3
              for ii = nout-refine+1:nout-1
                  feval(outfun,tout(ii),yout(ii,:).','interp',varargin{:});
              end
          elseif outflag == 1
              for ii = oldnout+1:nout
                  feval(outfun,tout(ii),yout(ii,:).','interp',varargin{:});
              end
          end
          if feval(outfun,tnew,y,'',varargin{:}) == 1
              tout = tout(1:nout);
              yout = yout(1:nout,:);
              stats = [nsteps; nfailed; nfevals; npds; ndecomps; nsolves];
              return;
          end
      end
      
  elseif haveoutfun
      if outflag == 2
          if feval(outfun,tnew,y(outputs),'',varargin{:}) == 1
              return;
          end
      elseif outflag == 3                 % computed points, with refinement
          for ii = 1:refine-1
              s = S(ii);
              yinterp = y + s*(cont(:,1)+(s-c2m1)*(cont(:,2)+(s-c1m1)*cont(:,3)));
              feval(outfun,tnew + h*S(ii),yinterp,'interp',varargin{:});
          end
      elseif outflag == 1 		% output only at tspan points
          ninterp = 0;
          while next <= ntspan
              if tdir * (tnew - tspan(next)) < 0
                  if haveeventfun && done
                      ninterp = ninterp + 1;
                      tinterp(ninterp,1) = tnew;
                      yinterp(:,ninterp) = y;
                  end
                  break;
              elseif tnew == tspan(next)
                  ninterp = ninterp + 1;
                  tinterp(ninterp,1) = tnew;
                  yinterp(:,ninterp) = y;
                  next = next + 1;
                  break;
              end
              ninterp = ninterp + 1;
              tinterp(ninterp,1) = tspan(next);
              s = (tspan(next) - tnew) / h;
              yinterp(:,ninterp) = y + s*(cont(:,1)+(s-c2m1)*(cont(:,2)+(s-c1m1)*cont(:,3)));
              next = next + 1;
          end
          if ninterp > 0
              if feval(outfun,tinterp(1:ninterp),yinterp(outputs,1:ninterp),'',varargin{:}) == 1
                  return;
              end
          end
      end
  end

  JCurrent = JConstant;
  
  if (last),
    h = hopt;
    idid = 1;
    break;
  end
  
  f0 = feval(odefile,t,y,args{:});
  nfevals = nfevals+1;
  
  hnew = posneg*min(abs(hnew),hmaxn);
  % ??? hopt = hnew;
  hopt = min(h,hnew);
  
  if (reject)
    hnew = posneg*min(abs(hnew),abs(h)); 
  end
  reject = false;
  
  if ( (t+hnew/0.99-tfinal)*posneg >= 0 ),
    h = tfinal-t;
    hhfac = h;
    last = true;
  else
    qt = hnew/h;
    hhfac = h;
  end
  
  consth = (qt >= quot1) & (qt <= quot2);
  consth = consth & (theta <= thet | JConstant);

  if consth,
    
    needNewJ = false;
    needNewLU = false;
    
  else
    
    if ( ~last ), h = hnew; end
    if (theta <= thet),
      needNewJ = false;
      needNewLU = true;
    else
      needNewJ = ~JConstant;
      needNewLU = true;
    end
    
  end
  
  scal_0 = abs(y)*rtol + atol;
  scal = scal_0./(h.^index);

end

if ( DEBUG >=1 ),
  fprintf('radau5: nsteps=%d, t=%7.2e, h=%8.2e, newt=%d, err=%9.2e [%s]\n', nsteps, t, h, newt, err_print, strategy);
end

if nargout~= 0
  tout = tout(1:nout);
  yout = yout(1:nout,:);
  stats=[ naccpt, nfailed, nfevals, npds, ndecomps, nsolves]';
end

if printstats
  fprintf('%g steps\n', nsteps);
  fprintf('%g successful steps\n', naccpt);
  fprintf('%g failed attempts\n', nfailed);
  fprintf('%g aborted attempts\n', nsteps-naccpt-nfailed);
  fprintf('%g calls to odefile\n', nfevals);
  fprintf('%g partial derivatives\n', npds);
  fprintf('%g LU decompositions\n', ndecomps);
  fprintf('%g solutions of linear systems\n', nsolves);
end



function yinterp = ntrprad(tinterp,t,y,tnew,ynew,h,cont)
%NTRPRAD Interpolation helper function for RADAU5.
%   YINTERP = NTRPRAD(TINTERP,T,Y,TNEW,YNEW,H,F) uses data computed in RADAU5
%   to approximate the solution at time TINTERP.
%   
%   See also RADAU5.
%fprintf('yinterp = ntrprad(tinterp,t,y,tnew,ynew,h,cont)\n'); % commented
%by Henrik Ebel on 01.12.2015
c1m1 = (4-sqrt(6))/10-1;
c2m1 = (4+sqrt(6))/10-1;

s = ((tinterp - tnew) / h)';       % may be a row vector
m = length(tinterp);
n = length(ynew);

yinterp = cont(:,ones(m,1)+2);
yinterp = yinterp.*( s(ones(n,1),:)-c1m1 ) + cont(:,ones(m,1)+1);
yinterp = yinterp.*( s(ones(n,1),:)-c2m1 ) + cont(:,ones(m,1));
yinterp = yinterp.*  s(ones(n,1),:)        + ynew(:,ones(m,1));



function radinfo(nind,Pred,Colmmd,start_n,fnewt)
%RADINFO Helper function for RADAU5.
fprintf('-----------------------------------------------------------------------\n')
if nind(2)==0 & nind(3)==0
  fprintf('System type: ordinary differential equation or index 1 system\n')
elseif nind(3)==0
  fprintf('System type: index 2 with %d/%d index 1/2 equations\n',nind(1),nind(2))
else
  fprintf('System type: index 3 with %d/%d/%d index 1/2/3 equations\n',nind)
end

if Pred
  fprintf(' - predictive step size controller.\n');
else
  fprintf(' - classical step size control.\n');
end

if Colmmd
  fprintf(' - mimimum degree ordering.\n');
else
  fprintf(' - standard (non)sparse solver.\n');
end

if ~start_n
  fprintf(' - trying to guess initial values for the Newton iteration.\n');
else
  fprintf(' - starting Newton iteration from zero initial values\n');
end
fprintf(' - stopping criterion for Newton''s method: %9.3e\n',fnewt);

fprintf('-----------------------------------------------------------------------\n')
