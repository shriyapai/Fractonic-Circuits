function[O]=SzSitei(i,L)
%i=str2num(i);
%L=str2num(L);

Sz=speye(3); Sz(2,2)=0; Sz(3,3)=-1/sqrt(2/3); Sz(1,1)=1/sqrt(2/3); 
I=speye(3);

if i==1
    O=Sz;
    for k=2:L
        O=kron(O,I);
    end
else 
    O=I;
    for k=2:i-1
       O= kron(O,I);
    end
    O=kron(O,Sz);
    for k=i+1:L
        O=kron(O,I);
    end
end

