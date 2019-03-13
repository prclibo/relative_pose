
rays1 = normr(rays1);
rays2 = normr(rays2);
E = E ./ norm(E_4pst0(:));
for i = 1:size(rays1, 1)
    Ex1 = E * rays1(i, :)';
    Etx2 = E' * rays2(i, :)';
    err = (rays2(i, :) * E * rays1(i, :)').^2 ./ ...
        sum([Ex1(1:2).^2; Etx2(1:2).^2])
end



E = E_4pst0; for i = 1:size(rays1, 1); Ex1 = E * rays1(i, :)'; Etx2 = E' * rays2(i, :)'; err(i) = (rays2(i, :) * E * rays1(i, :)').^2 ./  sum([Ex1(1:2).^2; Etx2(1:2).^2]); end