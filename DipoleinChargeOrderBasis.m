function[P]=DipoleinChargeOrderBasis(L)

A=permn([1 0 -1],L);
s=sum(A,2);
[sorted, row_ids] = sort(s, 'descend');
B=A(row_ids,:);

P=sparse(3^L,3^L);

for k=1:3^L
    for m=1:L
        P(k,k)=B(k,m)*(m-1);
    end
end
