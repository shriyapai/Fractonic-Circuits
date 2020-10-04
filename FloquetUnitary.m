function[U]=FloquetUnitary(L)

clearvars -except L 

%timestep 1
I=speye(3);
UA = haarUnitarySpinOne;
w=phaseAction;

if mod(L,3)==0
    for i=1:((L/3)-1)
        U2 = UA;
        UA=kron(UA,U2);
    end
end

if mod(L,3)==1
    for i=1:(((L-1)/3)-1)
        U2 = UA;
        UA=kron(UA,U2);
    end
    UA=kron(UA,w);
end

if mod(L,3)==2
    for i=1:(((L-2)/3)-1)
        U2 = UA;
        UA=kron(UA,U2);
    end
    UA=kron(UA,w);
    UA=kron(UA,w);
end

%timestep 2
UB = phaseAction;
U3 = haarUnitarySpinOne;
w1=phaseAction;

if mod(L,3)==0
    for i=1:((L/3)-1)
        UB=kron(UB,U3);
    end
    UB=kron(UB,w1);
    UB=kron(UB,w1);
end

if mod(L,3)==1
    for i=1:((L-1)/3)
        UB=kron(UB,U3);
    end
end

if mod(L,3)==2
    for i=1:(((L-2)/3))
        UB=kron(UB,U3);
    end
    UB=kron(UB,w1);
end


%timestep 3
w2 = phaseAction;
UC=kron(w2,w2);
U4 = haarUnitarySpinOne;

if mod(L,3)==0
    for i=1:((L/3)-1)
        UC=kron(UC,U4);
    end
    UC=kron(UC,w2);
end

if mod(L,3)==1
    for i=1:(((L-1)/3)-1)
        UC=kron(UC,U4);
    end
    UC=kron(UC,w2);
    UC=kron(UC,w2);
end

if mod(L,3)==2
    for i=1:(((L-2)/3))
        UC=kron(UC,U4);
    end
end

U=UA*UB*UC;

[QEevec,QEeval] = eig(full(U));
QEeval_only1=real(eig(full(U)));
QEeval_only=sort(QEeval_only1);

x=[1:3^L];

figure(1)
plot(x,QEeval_only,'b.','MarkerSize',8)
xlabel('Eigenvalue number')
ylabel('Quasienergy, E')
hold on

%entanglement entropy

A=permn([1 0 -1],L);
s=sum(A,2);
[sorted, row_ids] = sort(s, 'descend');
B=A(row_ids,:);

A2=permn([1 0 -1],floor(L/2));
s2=sum(A2,2);
[sorted2, row_ids2] = sort(s2, 'descend');
B2=A2(row_ids2,:);

compareB=B(:,1:floor(L/2));
S=sparse(3^L,1);
site=[1:L];

for t=1:3^L
    rho=QEevec(:,t)*QEevec(:,t)';
    R=sparse(3^(L-floor(L/2)),3^(L-floor(L/2)));
    [rID,cID] = find(rho);
    rho_other_basis=sparse(3^L,3^L);
    for a=1:length(rID)
        %for b=1:length(cID)
           rho_other_basis(row_ids(rID(a,1)),row_ids(cID(a,1)))=rho(rID(a,1),cID(a,1)); 
        %end
    end 
    SzExp=zeros(L,1);
    RW=sparse(3^L,L);
    for si=1:L
        [aS,RW(t,si)]=RightWeight(si,L,rho_other_basis);
        Sz=SzSitei(si,L);
        SzExp(si,1)=real(trace(rho_other_basis*Sz));
    end 
    figure(3)
    plot(site,SzExp,'k.-','MarkerSize',9)
    
    for z=1:3^(floor(L/2))
        indices = find(ismember(compareB,B2(z,:),'rows'));
        rho_red=sparse(length(indices),length(indices));
        for d=1:length(indices)
            for f=1:length(indices)
                rho_red(d,f)=rho_red(d,f)+rho(indices(d,1),indices(f,1));
            end
        end
        R=R+rho_red;
    end
    R=full(R);
    rho_ev=real(eig(R));
    rho_ev=rho_ev/sum(rho_ev);
    for k=1:length(rho_ev)
        if (rho_ev(k,1)<= 0.01)
            S(t,1)=S(t,1);
        else
            S(t,1)=S(t,1)-rho_ev(k,1)*log(rho_ev(k,1));
        end
    end
end
S=real(S);

figure(2)
plot(QEeval_only1,S,'r.','MarkerSize',8)
ylabel('Entanglement entropy, S')
xlabel('Quasienergy, E')
hold on

%for nb=1:3^L
%    figure(3)
%    plot(site,RW(nb,:)/sum(RW(nb,:)),'k.-','MarkerSize',10)
    %ylim([0 1])
%    hold on
%end










