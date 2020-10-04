function[P]=projector_Fibonacci(L) 
%this code constructs the projector onto the Fibonacci subspace. This is
%done using two layers of unitaries. The first layer acts on sites (1,2), (3,4),
%(5,6),..., while the second acts on (1), (2,3), (4,5),...Each two-site
%unitary projects into the Fibonacci subspace of two sites

P2 = [0 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1];
P_layer1 = P2;

%layer 1
if mod(L,2) == 0
    for j=3:2:L
        P_layer1 = kron(P_layer1,P2);
    end
else 
    for j=3:2:(L-1) 
        P_layer1 = kron(P_layer1,P2);
    end
    P_layer1 = kron(P_layer1,speye(2));
end

%layer 2
P_layer2 = speye(2);
if mod(L,2) == 1
    for j=2:2:L
        P_layer2 = kron(P_layer2,P2);
    end
else
    for j=2:2:(L-1) 
        P_layer2 = kron(P_layer2,P2);
    end
    P_layer2 = kron(P_layer2,speye(2));
end

P = P_layer1*P_layer2 ;
    