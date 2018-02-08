function loc = xlocANCF_3403(nodes,comps)
% make location vector in x from nodelist and componentlist
% adapted to beam element 3363 by Henrik Ebel on 20.10.2015

% loc = [];
% for n=1:length(nodes)
%   nn = nodes(n);
%   loc = [loc 18*nn-18+comps];
% end
ln = length(nodes);
lc = length(comps);
loc = zeros(1,ln*lc);
for n=1:ln
    nn = nodes(n);
    loc(1,(n-1)*lc+1:1:n*lc)=30*nn-30+comps;
end