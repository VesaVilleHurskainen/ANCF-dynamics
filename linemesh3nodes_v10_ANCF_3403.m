function [P,nloc] = linemesh3nodes_v10_ANCF_3403(n,a)

% make a linemesh with three node ANCF beam element 3363 - adapted by
% Henrik Ebel from code written by Marko Mattikainen for the 
% beam element 3333


% generate nodal coordinates in xy-plane
dx = a/n;
nodes=n*4-(n-1);

% Geospace for spacing
%xk=geospace(0,a,nodes,0.95)';

%xk=geospace(0,a,nodes,1)';
xk=linspace(0,a,nodes)'; % should do the same as the geospace cmd above...
yk=zeros(1,nodes)';
zk=zeros(1,nodes)';

nullmat=zeros(nodes,3);
onesvec=ones(nodes,1);

%dr1dx=nullmat;
%dr1dx(:,1)=onesvec;

dr1dy=nullmat;
dr1dy(:,2)=onesvec;

dr1dz=nullmat;
dr1dz(:,3)=onesvec;

P = [];
Pk = [xk yk zk dr1dy dr1dz nullmat nullmat nullmat nullmat nullmat nullmat nullmat];  
P = [P; Pk]; 

% generate element connectivity
nloc = [];

for k = 1:n
    loc = [(k-1)*3+1, (k-1)*3+2, (k-1)*3+3, (k-1)*3+4];
    nloc = [nloc; loc];
end

%rk by HE: now nloc contains in each row the information, which nodes (i.e.
%with which index, nodes counted from 1) belong to the element
% P contains the initial element degrees of freedom of each node (one line
% per node)
