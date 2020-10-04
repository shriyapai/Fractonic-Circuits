function[coeff,states]=test_update_rule(a,theta,orbits)

L=length(a);
no_orbits = size(orbits,1);
states = [a];
coeff = [1];
%z=L;

for k=2:2:L
    states_int2 = states;
    if (k~=L)
        for j=1:size(states_int2,1)
            states_int = states;
            if (states_int2(j,k-1)==0) & (states_int2(j,k+1)==0)
                c=coeff(j,1);
                coeff(j,1) = coeff(j,1)*(cosh(-1i*theta));
                states(j,k)=1-states(j,k);
                states = [states_int; states(j,:)];
                coeff = [coeff; c*sinh(-1i*theta)];
            end
        end
    end
    
    if (k==L)
        for j=1:size(states_int2,1)
            states_int = states;
            if (states_int2(j,L-1)==0) & (states_int2(j,1)==0)
                c=coeff(j,1);
                coeff(j,1) = coeff(j,1)*(cosh(-1i*theta));
                states(j,L)=1-states(j,L);
                states = [states_int; states(j,:)];
                coeff = [coeff; c*sinh(-1i*theta)];
            end
        end
    end
end

for k=1:2:L
    states_int2 = states;
    if (k~=L) & (k~=1)
        for j=1:size(states_int2,1)
            states_int = states;
            if (states_int2(j,k-1)==0) & (states_int2(j,k+1)==0)
                c=coeff(j,1);
                coeff(j,1) = coeff(j,1)*(cosh(-1i*theta));
                states(j,k)=1-states(j,k);
                states = [states_int; states(j,:)];
                coeff = [coeff; c*sinh(-1i*theta)];
            end
        end
    end
    
    if (k==1)
        for j=1:size(states_int2,1)
            states_int = states;
            if (states_int2(j,L)==0) & (states_int2(j,2)==0)
                c=coeff(j,1);
                coeff(j,1) = coeff(j,1)*(cosh(-1i*theta));
                states(j,1)=1-states(j,1);
                states = [states_int; states(j,:)];
                coeff = [coeff; c*sinh(-1i*theta)];
            end
        end
    end
    
    if (k==L)
        for j=1:size(states_int2,1)
            states_int = states;
            if (states_int2(j,L-1)==0) & (states_int2(j,1)==0)
                c=coeff(j,1);
                coeff(j,1) = coeff(j,1)*(cosh(-1i*theta));
                states(j,L)=1-states(j,L);
                states = [states_int; states(j,:)];
                coeff = [coeff; c*sinh(-1i*theta)];
            end
        end
    end 
end

%idx = find(abs(coeff)<(10^(-z)));
%states(idx,:) = [];
%coeff(idx,:) = [];

bin_array=zeros(size(states,1),1);
for h=1:size(states,1)
    bin_array(h,1) = binary_num_array(states(h,:));
end
[~,col_idx] = sort(bin_array(:,1),'descend');
states = states(col_idx,:);
coeff = coeff(col_idx,:);