function xlocAll= xlocAllANCF_3403(nloc)
% function xlocAllANCF makes full xloc for all elements 
% row - element, column - dof


[nl,~] = size(nloc);

xlocAll = zeros(nl,30*4); % (now 30 DOFs at each node)

for k = 1:nl
  n1 = [nloc(k,1)*30-29:1:nloc(k,1)*30];
  n2 = [nloc(k,2)*30-29:1:nloc(k,2)*30];
  n3 = [nloc(k,3)*30-29:1:nloc(k,3)*30];
  n4 = [nloc(k,4)*30-29:1:nloc(k,4)*30];
  xlocAll(k,:) = [n1 n2 n3 n4];
end
% rk by HE: now xlocAll contains in each row the information, which indices
% the element DOFs of the respective elements have - first the numbers of
% the 12 element degrees of freedom of the first node of the element, then
% of the second, then of the third