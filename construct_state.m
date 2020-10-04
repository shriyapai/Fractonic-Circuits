function[state]=construct_state(L,config)

if config(1,1)==0
    state = sparse([0;1]);
else
    state = sparse([1;0]);
end

for k=2:L
    if config(1,k)==0
        state =kron(state,sparse([0;1]));
    end
    if config(1,k)==1
        state = kron(state,sparse([1;0]));
    end
end