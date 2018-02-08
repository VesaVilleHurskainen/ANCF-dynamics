function [r]=kinematic_beams(Element,H,W,L,XI,e,DIM)
% Fiilussa kuvataan elementin kinematiikka. Pyrit‰‰n esitt‰m‰‰n muotofunktiot 
% xi-koordinaatistossa [-1..1]...ainakin k‰ytett‰ess‰ isoparametrisia
% ancf-elementtej‰
% Fiilu on jatkoa shapefunc-koodeihin
% Recoded by MKM 2007

% Kevennetty versio, ainoastaan elementti B12 

% Elements 3343 and 3343s111 added by Henrik Ebel in October 2015

e=e(:);


if Element==3333 || Element==3334  %  dof beam
    xi=XI(1);
    eta=XI(2);
    zeta=XI(3);
    
    Svec=[1/2*xi*(-1+xi);
        1/4*H*xi*eta*(-1+xi);
        1/4*W*zeta*xi*(-1+xi);
        -(-1+xi)*(xi+1);
        -1/2*H*eta*(-1+xi)*(xi+1);
        -1/2*W*zeta*(-1+xi)*(xi+1);
        1/2*xi*(xi+1);
        1/4*H*xi*eta*(xi+1);
        1/4*W*zeta*xi*(xi+1)];
    
    % Muotofunktiovektori matriisiksi
    for ii=1:3,
        for jj=1:9,
            jjj=((jj-1)*3+1)+(ii-1);
            S(ii,jjj)=Svec(jj);
        end
    end
    
    r=S*e;
    
elseif Element==3343
    xi=XI(1);
    eta=XI(2);
    zeta=XI(3);
    % shape function vector
    Svec=[1/2*xi^2-1/2*xi;
          1/4*eta*H*xi^2-1/4*eta*H*xi;
          1/4*zeta*W*xi^2-1/4*zeta*W*xi;
          1/8*eta*H*zeta*W*xi^2-1/8*eta*H*zeta*W*xi;
          -xi^2+1;
          -1/2*eta*H*xi^2+1/2*eta*H;
          -1/2*zeta*W*xi^2+1/2*zeta*W;
          -1/4*eta*H*zeta*W*xi^2+1/4*eta*H*zeta*W;
          1/2*xi^2+1/2*xi;
          1/4*eta*H*xi^2+1/4*eta*H*xi;
          1/4*zeta*W*xi^2+1/4*zeta*W*xi;
          1/8*eta*H*zeta*W*xi^2+1/8*eta*H*zeta*W*xi];
    % determine shape function matrix
    % Using kronecker tensor product - because it is fast AND fancy
    S = kron(Svec',eye(3));
    r=S*e;
elseif Element==3343111
    xi=XI(1);
    eta=XI(2);
    zeta=XI(3);
    % shape function vector
    Svec=[1/2*xi^2-1/2*xi;
          1/4*eta*H*xi^2-1/4*eta*H*xi;
          1/4*zeta*W*xi^2-1/4*zeta*W*xi;
          1/48*L*eta*H*zeta*W*xi^3-1/32*L*eta*H*zeta*W*xi^2+5/96*L*eta*H*zeta*W;
          -xi^2+1;
          -1/2*eta*H*xi^2+1/2*eta*H;
          -1/2*zeta*W*xi^2+1/2*zeta*W;
          -1/24*L*eta*H*zeta*W*xi^3+1/8*L*eta*H*zeta*W*xi+1/12*L*eta*H*zeta*W;
          1/2*xi^2+1/2*xi;
          1/4*eta*H*xi^2+1/4*eta*H*xi;
          1/4*zeta*W*xi^2+1/4*zeta*W*xi;
          1/48*L*eta*H*zeta*W*xi^3+1/32*L*eta*H*zeta*W*xi^2-1/96*L*eta*H*zeta*W];
    % determine shape function matrix
    % Using kronecker tensor product - because it is fast AND fancy
    S = kron(Svec',eye(3));
    r=S*e;
elseif Element==3353
    xi=XI(1);
    eta=XI(2);
    zeta=XI(3);
    % shape function vector
    Svec=[1/2*xi^2-1/2*xi;
          1/4*eta*H*xi^2-1/4*eta*H*xi;
          1/4*zeta*W*xi^2-1/4*zeta*W*xi;
          1/16*eta^2*H^2*xi^2-1/16*eta^2*H^2*xi;
          1/16*zeta^2*W^2*xi^2-1/16*zeta^2*W^2*xi;
          -xi^2+1;
          -1/2*eta*H*xi^2+1/2*eta*H;
          -1/2*zeta*W*xi^2+1/2*zeta*W;
          -1/8*eta^2*H^2*xi^2+1/8*eta^2*H^2;
          -1/8*zeta^2*W^2*xi^2+1/8*zeta^2*W^2;
          1/2*xi^2+1/2*xi;
          1/4*eta*H*xi^2+1/4*eta*H*xi;
          1/4*zeta*W*xi^2+1/4*zeta*W*xi;
          1/16*eta^2*H^2*xi^2+1/16*eta^2*H^2*xi;
          1/16*zeta^2*W^2*xi^2+1/16*zeta^2*W^2*xi];
    % determine shape function matrix
    % Using kronecker tensor product - because it is fast AND fancy
    S = kron(Svec',eye(3));
    r=S*e;
elseif Element==3363
    xi=XI(1);
    eta=XI(2);
    zeta=XI(3);
    % shape function vector
    Svec=[1/2*xi^2-1/2*xi;
          1/4*eta*H*xi^2-1/4*eta*H*xi;
          1/4*zeta*W*xi^2-1/4*zeta*W*xi;
          1/8*eta*H*zeta*W*xi^2-1/8*eta*H*zeta*W*xi;
          1/16*eta^2*H^2*xi^2-1/16*eta^2*H^2*xi;
          1/16*zeta^2*W^2*xi^2-1/16*zeta^2*W^2*xi;
          -xi^2+1;
          -1/2*eta*H*xi^2+1/2*eta*H;
          -1/2*zeta*W*xi^2+1/2*zeta*W;
          -1/4*eta*H*zeta*W*xi^2+1/4*eta*H*zeta*W;
          -1/8*eta^2*H^2*xi^2+1/8*eta^2*H^2;
          -1/8*zeta^2*W^2*xi^2+1/8*zeta^2*W^2;
          1/2*xi^2+1/2*xi;
          1/4*eta*H*xi^2+1/4*eta*H*xi;
          1/4*zeta*W*xi^2+1/4*zeta*W*xi;
          1/8*eta*H*zeta*W*xi^2+1/8*eta*H*zeta*W*xi;
          1/16*eta^2*H^2*xi^2+1/16*eta^2*H^2*xi;
          1/16*zeta^2*W^2*xi^2+1/16*zeta^2*W^2*xi];
    % determine shape function matrix
    % Using kronecker tensor product - because it is fast AND fancy
    S = kron(Svec',eye(3));
    r=S*e;
elseif Element==34103
    xi=XI(1);
    eta=XI(2);
    zeta=XI(3);
    % shape function vector
    Svec=[-9/16*xi^3+9/16*xi^2+1/16*xi-1/16;
          -9/32*eta*H*xi^3+9/32*eta*H*xi^2+1/32*eta*H*xi-1/32*eta*H;
          -9/32*zeta*W*xi^3+9/32*zeta*W*xi^2+1/32*zeta*W*xi-1/32*zeta*W;
          -9/64*eta*H*zeta*W*xi^3+9/64*eta*H*zeta*W*xi^2+1/64*eta*H*zeta*W*xi-1/64*eta*H*zeta*W;
          -9/128*eta^2*H^2*xi^3+9/128*eta^2*H^2*xi^2+1/128*eta^2*H^2*xi-1/128*eta^2*H^2;
          -9/128*zeta^2*W^2*xi^3+9/128*zeta^2*W^2*xi^2+1/128*zeta^2*W^2*xi-1/128*zeta^2*W^2;
          -9/256*eta^2*H^2*zeta*W*xi^3+9/256*eta^2*H^2*zeta*W*xi^2+1/256*eta^2*H^2*zeta*W*xi-1/256*eta^2*H^2*zeta*W;
          -9/256*eta*H*zeta^2*W^2*xi^3+9/256*eta*H*zeta^2*W^2*xi^2+1/256*eta*H*zeta^2*W^2*xi-1/256*eta*H*zeta^2*W^2;
          -3/256*eta^3*H^3*xi^3+3/256*eta^3*H^3*xi^2+1/768*eta^3*H^3*xi-1/768*eta^3*H^3;
          -3/256*zeta^3*W^3*xi^3+3/256*zeta^3*W^3*xi^2+1/768*zeta^3*W^3*xi-1/768*zeta^3*W^3;
          27/16*xi^3-9/16*xi^2-27/16*xi+9/16;
          27/32*eta*H*xi^3-9/32*eta*H*xi^2-27/32*eta*H*xi+9/32*eta*H;
          27/32*zeta*W*xi^3-9/32*zeta*W*xi^2-27/32*zeta*W*xi+9/32*zeta*W;
          27/64*eta*H*zeta*W*xi^3-9/64*eta*H*zeta*W*xi^2-27/64*eta*H*zeta*W*xi+9/64*eta*H*zeta*W;
          27/128*eta^2*H^2*xi^3-9/128*eta^2*H^2*xi^2-27/128*eta^2*H^2*xi+9/128*eta^2*H^2;
          27/128*zeta^2*W^2*xi^3-9/128*zeta^2*W^2*xi^2-27/128*zeta^2*W^2*xi+9/128*zeta^2*W^2;
          27/256*eta^2*H^2*zeta*W*xi^3-9/256*eta^2*H^2*zeta*W*xi^2-27/256*eta^2*H^2*zeta*W*xi+9/256*eta^2*H^2*zeta*W;
          27/256*eta*H*zeta^2*W^2*xi^3-9/256*eta*H*zeta^2*W^2*xi^2-27/256*eta*H*zeta^2*W^2*xi+9/256*eta*H*zeta^2*W^2;
          9/256*eta^3*H^3*xi^3-3/256*eta^3*H^3*xi^2-9/256*eta^3*H^3*xi+3/256*eta^3*H^3;
          9/256*zeta^3*W^3*xi^3-3/256*zeta^3*W^3*xi^2-9/256*zeta^3*W^3*xi+3/256*zeta^3*W^3;
          -27/16*xi^3-9/16*xi^2+27/16*xi+9/16;
          -27/32*eta*H*xi^3-9/32*eta*H*xi^2+27/32*eta*H*xi+9/32*eta*H;
          -27/32*zeta*W*xi^3-9/32*zeta*W*xi^2+27/32*zeta*W*xi+9/32*zeta*W;
          -27/64*eta*H*zeta*W*xi^3-9/64*eta*H*zeta*W*xi^2+27/64*eta*H*zeta*W*xi+9/64*eta*H*zeta*W;
          -27/128*eta^2*H^2*xi^3-9/128*eta^2*H^2*xi^2+27/128*eta^2*H^2*xi+9/128*eta^2*H^2;
          -27/128*zeta^2*W^2*xi^3-9/128*zeta^2*W^2*xi^2+27/128*zeta^2*W^2*xi+9/128*zeta^2*W^2;
          -27/256*eta^2*H^2*zeta*W*xi^3-9/256*eta^2*H^2*zeta*W*xi^2+27/256*eta^2*H^2*zeta*W*xi+9/256*eta^2*H^2*zeta*W;
          -27/256*eta*H*zeta^2*W^2*xi^3-9/256*eta*H*zeta^2*W^2*xi^2+27/256*eta*H*zeta^2*W^2*xi+9/256*eta*H*zeta^2*W^2;
          -9/256*eta^3*H^3*xi^3-3/256*eta^3*H^3*xi^2+9/256*eta^3*H^3*xi+3/256*eta^3*H^3;
          -9/256*zeta^3*W^3*xi^3-3/256*zeta^3*W^3*xi^2+9/256*zeta^3*W^3*xi+3/256*zeta^3*W^3;
          9/16*xi^3+9/16*xi^2-1/16*xi-1/16;
          9/32*eta*H*xi^3+9/32*eta*H*xi^2-1/32*eta*H*xi-1/32*eta*H;
          9/32*zeta*W*xi^3+9/32*zeta*W*xi^2-1/32*zeta*W*xi-1/32*zeta*W;
          9/64*eta*H*zeta*W*xi^3+9/64*eta*H*zeta*W*xi^2-1/64*eta*H*zeta*W*xi-1/64*eta*H*zeta*W;
          9/128*eta^2*H^2*xi^3+9/128*eta^2*H^2*xi^2-1/128*eta^2*H^2*xi-1/128*eta^2*H^2;
          9/128*zeta^2*W^2*xi^3+9/128*zeta^2*W^2*xi^2-1/128*zeta^2*W^2*xi-1/128*zeta^2*W^2;
          9/256*eta^2*H^2*zeta*W*xi^3+9/256*eta^2*H^2*zeta*W*xi^2-1/256*eta^2*H^2*zeta*W*xi-1/256*eta^2*H^2*zeta*W;
          9/256*eta*H*zeta^2*W^2*xi^3+9/256*eta*H*zeta^2*W^2*xi^2-1/256*eta*H*zeta^2*W^2*xi-1/256*eta*H*zeta^2*W^2;
          3/256*eta^3*H^3*xi^3+3/256*eta^3*H^3*xi^2-1/768*eta^3*H^3*xi-1/768*eta^3*H^3;
          3/256*zeta^3*W^3*xi^3+3/256*zeta^3*W^3*xi^2-1/768*zeta^3*W^3*xi-1/768*zeta^3*W^3];
    % determine shape function matrix
    % Using kronecker tensor product - because it is fast AND fancy
    S = kron(Svec',eye(3));
    r=S*e;    
elseif Element==3273
    xi=XI(1);
    eta=XI(2);
    zeta=XI(3);
    % shape function vector
    Svec=[1/4*xi^3-3/4*xi+1/2;
          1/8*L*xi^3-1/8*L*xi^2-1/8*L*xi+1/8*L;
          -1/4*eta*H*xi+1/4*eta*H;
          -1/4*zeta*W*xi+1/4*zeta*W;
          -1/8*eta*H*zeta*W*xi+1/8*eta*H*zeta*W;
          -1/16*eta^2*H^2*xi+1/16*eta^2*H^2;
          -1/16*zeta^2*W^2*xi+1/16*zeta^2*W^2;
          -1/4*xi^3+3/4*xi+1/2;
          1/8*L*xi^3+1/8*L*xi^2-1/8*L*xi-1/8*L;
          1/4*eta*H*xi+1/4*eta*H;
          1/4*zeta*W*xi+1/4*zeta*W;
          1/8*eta*H*zeta*W*xi+1/8*eta*H*zeta*W;
          1/16*eta^2*H^2*xi+1/16*eta^2*H^2;
          1/16*zeta^2*W^2*xi+1/16*zeta^2*W^2];
    % determine shape function matrix
    % Using kronecker tensor product - because it is fast AND fancy
    S = kron(Svec',eye(3));
    r=S*e;
else
    disp('****** Elementill‰ ei indeksi‰ !! (virheilmoitus tiedostosta kinem_v004_Beam.m) ******');
end    

if DIM==2
        r(3)=XI(3);   % Yleisyyden takia, laitetaan Z-koordinaatiksi vaikkapa XI(3) 
end
