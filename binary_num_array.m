function[n]=binary_num_array(a)

n=0;
L=length(a);
for k=1:length(a)
    n = n + (2^(k-1))*a(L-(k-1));
end
