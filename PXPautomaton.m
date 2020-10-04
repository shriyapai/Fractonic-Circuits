function[U]=PXPautomaton(L,ep)  %L is the length of the spin-1/2 chain
%constructs the unitary under open boundary conditions. PXP for the first
%(last)site is just XP (PX).

X = sparse([0 1 ; 1 0]) ;
Z = sparse([1 0; 0 -1]) ;
P = sparse((speye(2) - Z)/2) ;
th = (pi/2)*(1+ep) ;

U_odd = sparse(zeros(2^L)) ; 
U_even = sparse(zeros(2^L));
U = sparse(zeros(2^L));

odd_exp = kron(X,P) ;
for i=3:L
    odd_exp = kron(odd_exp,speye(2));
end
U_odd = expm(-1i*th*(odd_exp)) ;

for j = 2:L-1
    if (mod(j,2)==0)
        if j==2
            even_exp = kron(P,X);
            even_exp = kron(even_exp,P) ;
            for i=4:L
                even_exp = kron(even_exp,speye(2)) ;
            end
            U_even = expm(-1i*th*(even_exp)) ;
        else
            even_exp = speye(2) ;
            for k=2:j-2
                even_exp = kron(even_exp, speye(2)) ;
            end
            even_exp = kron(even_exp,P);
            even_exp = kron(even_exp,X) ;
            even_exp = kron(even_exp,P) ;
            if j<= (L-2)
                for i=j+2:L
                    even_exp = kron(even_exp,speye(2));
                end 
            end
        end
        U_even = U_even*expm(-1i*th*(even_exp)) ;
    end
    if (mod(j,2)==1)
        odd_exp = speye(2) ;
        if j==3
            odd_exp = kron(odd_exp,P);
            odd_exp = kron(odd_exp,X) ;
            odd_exp = kron(odd_exp,P) ;
            for i=5:L
                odd_exp = kron(odd_exp,speye(2));
            end
        else 
            for k=2:j-2
                odd_exp = kron(odd_exp, speye(2)) ;
            end
            odd_exp = kron(odd_exp,P);
            odd_exp = kron(odd_exp,X) ;
            odd_exp = kron(odd_exp,P) ;
            if j<= (L-2)
                for i=j+2:L
                    odd_exp = kron(odd_exp,speye(2));
                end
            end
        end
        U_odd = U_odd*expm(-1i*th*(odd_exp)) ;
    end
end

if mod(L,2) == 1
    odd_exp = kron(speye(2),speye(2)) ;
    for i=3:L-2
        odd_exp = kron(odd_exp,speye(2));
    end
    odd_exp = kron(odd_exp,P) ;
    odd_exp = kron(odd_exp,X) ;
    U_odd = U_odd*expm(-1i*th*(odd_exp)) ;
else
    even_exp = kron(speye(2),speye(2)) ;
    for i=3:L-2
        even_exp = kron(even_exp,speye(2));
    end
    even_exp = kron(even_exp,P) ;
    even_exp = kron(even_exp,X) ;
    U_even = U_even*expm(-1i*th*(even_exp)) ;
end

U = U_odd*U_even;