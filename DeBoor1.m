function [basis] = DeBoor1(deg,kvect,strt,ent,get_t)
nn = get_t;
nbasis = length(kvect) - deg - 1;
k = nbasis + deg;
m = length(nn);
tbasis = zeros(k,m);

for ii = 1:k
for n = 1:m
 t = nn(n);
 if kvect(ii) <= t && t < kvect(ii+1)
tbasis(ii,n) = 1;
 end
end
end

basis = recurse(deg,kvect,nn,1,tbasis);

% -----------------------------------------------------------------------

function ctbasis = recurse(deg,kvect,nn,curr,ptbasis)
k = size(ptbasis,1)-1;
m = size(ptbasis,2);
ctbasis = zeros(k,m);

for ii = 1:k
for n = 1:m
 t = nn(n);
if kvect(ii+curr) - kvect(ii) ~= 0
ctbasis(ii,n) = ((t-kvect(ii))/(kvect(ii+curr)-kvect(ii))) * ptbasis(ii,n);
end
if kvect(ii+curr+1) - kvect(ii+1) ~= 0
 ctbasis(ii,n) = ctbasis(ii,n) + ((kvect(ii+curr+1)-t)/(kvect(ii+curr+1)-kvect(ii+1))) * ptbasis(ii+1,n);
end
end
end

if curr ~= deg
ctbasis = recurse(deg,kvect,nn,curr+1,ctbasis);
end