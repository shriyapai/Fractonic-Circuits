function S = Entropy(psi,L)
%computing half-chain entanglement entropy

psi_R = reshape(psi, 2^(L-floor(L/2)), 2^(floor(L/2)));
rho_R = psi_R' * psi_R;
lambda = nonzeros(eig(full(rho_R)));
S = -real(lambda'*log(lambda));