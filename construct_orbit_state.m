function[state]=construct_orbit_state(L,config) 
%config is the representative configuration chosen in constructing U

v=config;
a=circshift(config,2);
config = [config; a];
while ~isequal(a,v)
    a = circshift(a,2);
    config = [config; a];
end
unique_config = unique(config,'rows','stable');

state=sparse(zeros(2^L,1));

for j=1:size(unique_config,1)
    state = state + construct_state(L,unique_config(j,:));
end

state = state/(sqrt(size(unique_config,1)));


