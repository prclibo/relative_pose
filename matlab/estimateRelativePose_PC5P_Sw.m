function E = estimateRelativePose_PC5P_Sw(rays1, rays2)

for i = 1:5
    AE(i, :) = reshape(rays2(i, :)' * rays1(i, :), 1, []);
end
nulle = null(AE);

[x, y, z] = solver_sw5pt(nulle(:, 1), nulle(:, 2), nulle(:, 3), nulle(:, 4));

E = [];
for i = 1:numel(x)
    Evec = nulle * [x(i); y(i); z(i); 1];
    E = [E; reshape(Evec, 3, 3)];
end

end