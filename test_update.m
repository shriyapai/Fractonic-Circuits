function[coefficients]=test_update(a,theta,orbits)

L=length(a);
no_orbits = size(orbits,1);
coefficients = sparse(zeros(no_orbits,1));

[a_configs,Na]=generate_orbits(L,a);

for p=1:Na
    [coeff,states]=test_update_rule(a_configs(p,:),theta,orbits);
    for q1=1:size(states,1)
        [states_orbits,Ns]=generate_orbits(L,states(q1,:));
        coeff(q1,1) = coeff(q1,1)/(sqrt(Ns)*sqrt(Na));
        v2 = move_cyclic_rep(states(q1,:));
        [q3, index] = ismember(orbits, v2, 'rows');
        index = find(index);
        coefficients(index,1) = coefficients(index,1) + coeff(q1,1);
    end
end

idx = find(abs(coefficients)<(10^(-L)));
coefficients(idx,:) = 0;



