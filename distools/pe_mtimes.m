%PE_MTIMES Pseudo-Euclidean matrix product
%
%  C = PE_MTIMES(A,B)
%
% A and/or B should be datasets. Their PE signature SIG is retrieved,
% C = A*J*B, where J is a diagonal matrix with 1's, followed by -1's.
% J = diag ([ONES(SIG(1),1);  -ONES(SIG(2),1)])

function c = pe_mtimes(a,b)

siga = getsig(a);
sigb = getsig(b); 

if isdataset(a) & isdataset(b)
  if any(siga~=sigb)
    error('PE signatures should be equal')
  end
  J = diag([ones(1,siga(1))  -ones(1,siga(2))]);
  c = a*J*b;
elseif isdataset(a)
  J = diag([ones(1,siga(1))  -ones(1,siga(2))]);
  c = a*J*b;
elseif isdataset(b) % problem! make doubles
  J = diag([ones(1,sigb(1))  -ones(1,sigb(2))]);
  c = a*J*(+b);
else % just doubles, no PE space
  c = a*b;
end
  
  