function[orbit_configs,N]=generate_orbits(L,config)

v=config;
a=circshift(config,2);
config = [config; a];
while ~isequal(a,v)
    a = circshift(a,2);
    config = [config; a];
end
orbit_configs = unique(config,'rows','stable');

N=size(orbit_configs,1);