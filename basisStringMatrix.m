function[S]=basisStringMatrix(v,L0,L1,L2,L3,L4,L5,L6,L7,L8)

tic

if v(1)==1
    S=L0;
end
if v(1)==2
    S=L1;
end
if v(1)==3
    S=L2;
end
if v(1)==4
    S=L3;
end
if v(1)==5
    S=L4;
end
if v(1)==6
    S=L5;
end
if v(1)==7
    S=L6;
end
if v(1)==8
    S=L7;
end
if v(1)==9
    S=L8;
end

for j=2:length(v)
    if v(j)==1
        S=kron(S,L0);
    end
    if v(j)==2
        S=kron(S,L1);
    end
    if v(j)==3
        S=kron(S,L2);
    end
    if v(j)==4
        S=kron(S,L3);
    end
    if v(j)==5
        S=kron(S,L4);
    end
    if v(j)==6
        S=kron(S,L5);
    end
    if v(j)==7
        S=kron(S,L6);
    end
    if v(j)==8
        S=kron(S,L7);
    end
    if v(j)==9
        S=kron(S,L8);
    end
end

toc