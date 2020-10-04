function[b]=move_cyclic_rep(a)

L=length(a);
[orbit_configs,N]=generate_orbits(L,a);
bin_array=zeros(size(orbit_configs,1),1);

for h=1:size(orbit_configs,1)
    bin_array(h,1) = binary_num_array(orbit_configs(h,:));
end
idx = find(bin_array == max(bin_array));

b = orbit_configs(idx,:);
