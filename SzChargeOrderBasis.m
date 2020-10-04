function[Sz]=SzChargeOrderBasis(L)

clear i
clear j

A=permn([1 0 -1],L);
s=sum(A,2);
Sz=sparse(3^L,3^L);
blockSize=zeros(2*L+1,1);
bsz=0;

for i=1:2*L+1
   blockSize(i,1)=length(find(s==i-L-1));
   for j=1:blockSize(i,1)
       Sz(j+bsz,j+bsz)=-(i-L-1);       
   end
   bsz=bsz+blockSize(i,1);
end



