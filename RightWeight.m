function[a,rho]=RightWeight(i,L,Ot)

L0=speye(3);
L1=sparse(3,3); L1(1,2)=1; L1(2,1)=1; L1=L1/sqrt(2/3);
L2=sparse(3,3); L2(1,2)=-i; L2(2,1)=i; L2=L2/sqrt(2/3);
L3=sparse(3,3); L3(1,1)=1; L3(3,3)=-1; L3=L3/sqrt(2/3);
L4=sparse(3,3); L4(1,3)=1; L4(3,1)=1; L4=L4/sqrt(2/3);
L5=sparse(3,3); L5(1,3)=-i; L5(3,1)=i; L5=L5/sqrt(2/3);
L6=sparse(3,3); L6(2,3)=1; L6(3,2)=1; L6=L6/sqrt(2/3);
L7=sparse(3,3); L7(2,3)=-i; L7(3,2)=i; L7=L7/sqrt(2/3);
L8=sparse(3,3); L8(1,1)=1/sqrt(3); L8(2,2)=-2/sqrt(3); L8(3,3)=1/sqrt(3); L8=L8/sqrt(2/3);

clear basis
clear indices
clear k

basis=basisStrings2(L);
B=basis(:,i+1:L);

stot=0;
for h=i+1:L
    stot=stot+h;
end

indices=find(sum(B,2)==stot);
for k=1:length(indices)
    if basis(indices(k,1),i)==0
        indices(k,1)=0;
    end
end

indices=indices(indices>0);

a=sparse(9^L,1);
rho=0;

for j=1:length(indices)
    S=basisStringMatrix(basis(indices(j,1),:),L0,L1,L2,L3,L4,L5,L6,L7,L8);
    a(j,1)=trace(S*Ot)/(3^L); 
    rho=rho + conj(a(j,1))*a(j,1);
end

a=full(a);
