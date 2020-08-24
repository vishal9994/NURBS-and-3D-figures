function [nn, basis] = DeBoor(deg,kvect,strt,endt)

knots = unique(kvect(strt:endt));
div = 30;
nn = zeros((div-1)*(length(knots)-1)+1,1);
for j = 1:length(knots)-1
if j == 1
 nn(1:div) = linspace(knots(j),knots(j+1),div);
else
 nn((div-1)*(j-1)+1:(div-1)*j+1) = linspace(knots(j),knots(j+1),div);
end  
end

nbasis = length(kvect) - deg - 1;
k = nbasis + deg;
m = length(nn);
tbasis = zeros(k,m);

for B = 1:k
for n = 1:m
 t = nn(n);
 if kvect(B) <= t && t < kvect(B+1)
 tbasis(B,n) = 1;
 end
end
end

basis = recurse(deg,kvect,nn,1,tbasis);
%display(Basis);
% -----------------------------------------------------------------------

function ctbasis = recurse(deg,kvect,nn,curr,ptbasis)
k = size(ptbasis,1)-1;
m = size(ptbasis,2);
ctbasis = zeros(k,m);

for B = 1:k
for n = 1:m
   t = nn(n);    
if kvect(B+curr) - kvect(B) ~= 0
  ctbasis(B,n) = ((t-kvect(B))/(kvect(B+curr)-kvect(B))) * ptbasis(B,n);
end
    
if kvect(B+curr+1) - kvect(B+1) ~= 0
  ctbasis(B,n) = ctbasis(B,n) + ((kvect(B+curr+1)-t)/(kvect(B+curr+1)-kvect(B+1))) * ptbasis(B+1,n);
end
 end
end

if curr ~= deg
ctbasis = recurse(deg,kvect,nn,curr+1,ctbasis);
end

% -----------------------------------------------------------------------
