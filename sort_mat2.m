function[orbits]=sort_mat2(mat)
tosort = mat;
orbits = [];
while ~isempty(tosort)
    circvec = tosort(1, :);
    orbits=[orbits;circvec];
    
    c2=ismember(tosort,circvec,'rows');
    idx2=find(c2);
    tosort(idx2,:)=[];
    
    a=circshift(circvec,2);
    c=ismember(tosort,a,'rows');
    idx=find(c);
    tosort(idx,:)=[];
    
    while ~isequal(a,circvec)
        a=circshift(a,2);
        c=ismember(tosort,a,'rows');
        idx=find(c);
        tosort(idx,:)=[];
    end
end