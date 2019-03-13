rgb_posi = nan(3, numel(rgb_poses));
for i = 1:numel(rgb_poses)
    if ~isempty(rgb_poses{i})
        rgb_posi(:, i) = rgb_poses{i}(1:3, 4);
    end
end

% rgb_posi = rgb_posi(:, 1:100);
plot(rgb_posi(1, :), rgb_posi(3, :), 'b');
hold on; axis equal;

abs_5p = eye(4);
abs_4pst0 = eye(4);

posi_5p = nan(3, numel(rel_5p));
posi_4pst0 = nan(3, numel(rel_4pst0));
for i = 1:numel(rel_5p)
    if isempty(rel_5p{i}) && isempty(rel_4pst0{i}); continue; end
    abs_5p = abs_5p * rel_5p{i};
    abs_4pst0 = abs_4pst0 * rel_4pst0{i};
    posi_5p(1:3, i) = abs_5p(1:3, 4);
    posi_4pst0(1:3, i) = abs_4pst0(1:3, 4);
end

posi_5p = posi_5p(:, isfinite(posi_5p(1, :)));
posi_4pst0 = posi_4pst0(:, isfinite(posi_4pst0(1, :)));

plot(posi_5p(1, :), posi_5p(3, :), 'r');
plot(posi_4pst0(1, :), posi_4pst0(3, :), 'g');



return

abs_5p = eye(4);
abs_4pst0 = eye(4);

pos_5p(1:2, 1) = 0;
pos_4pst0(1:2, 1) = 0;
for i = 2:numel(rel_5p)
    abs_5p = abs_5p * rel_5p{i};
    abs_4pst0 = abs_4pst0 * rel_4pst0{i};
    pos_5p(1:3, i) = abs_5p(1:3, 4);
    pos_4pst0(1:3, i) = abs_4pst0(1:3, 4);
end
