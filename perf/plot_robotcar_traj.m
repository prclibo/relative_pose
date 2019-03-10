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
