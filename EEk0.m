function[D,QEeval_only,S]=EEk0(L,ep)

if L<=16
    mat=test_fib2(L);
    orbits=sort_mat2(mat);
end
if L==18
    mat = load('mat18.dat');
    orbits = load('orbits18.dat');
end
if L==20
    mat = load('mat20.dat');
    orbits = load('orbits20.dat');
end
if L==22
    mat = load('mat22.dat');
    orbits = load('orbits22.dat');
end
if L==24
    mat = load('mat24.dat');
    orbits = load('orbits24.dat');
end

no_orbits = size(orbits,1);
U=sparse(zeros(no_orbits,no_orbits));
theta = (pi/2)*(1+ep) ;

for j=1:no_orbits
    [coefficients]=test_update(orbits(j,:),theta,orbits);
    U(:,j)=coefficients;
end
D=full(U);

[QEevec,QEeval] = eig(full(U));
QEeval_only=real(i*log(eig(full(U)))); %diagonalize U and calculate the quasienergies
%Sp = (L/2)*log(2)-(1/2) ;
SpFib = (L/2)*(log((1+sqrt(5))/2))-(1/2) ;
states = sparse(zeros(2^L,no_orbits));

for j=1:no_orbits
    states(:,j)=construct_orbit_state(L,orbits(j,:)) ;
end

%entropy calculation of states in the k=0 Fibonacci subspace
EEstates = sparse(zeros(2^L,1));
S = sparse(zeros(length(QEeval_only),1)) ; 
for m=1:no_orbits
    EEstates(:,1) = states*QEevec(:,m);
    EEstates(:,1) = EEstates(:,1)./sqrt(EEstates(:,1)'*EEstates(:,1)) ;
    S(m,1) = Entropy(EEstates,L) ;
end

[~,col_idx] = sort(QEeval_only(:,1),'descend');
QEeval_only = QEeval_only(col_idx,:);
S = S(col_idx,:);

figure(1)
if ep==0
    plot(QEeval_only/pi,S./SpFib,'b.','MarkerSize',12)
end
if ep==0.01
    plot(QEeval_only/pi,S./SpFib,'r.','MarkerSize',12)
end
if ep==0.1
    plot(QEeval_only/pi,S./SpFib,'k.','MarkerSize',12)
end
ylabel('Entanglement entropy, S/S_P')
xlabel('Quasienergy, \epsilon/\pi')
xlim([-1 1])
ylim([0 1])
hold on