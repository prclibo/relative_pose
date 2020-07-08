clear, clc

E1 = sym('E1%d%d', [3 ,3]);
E2 = sym('E2%d%d', [3 ,3]);
E3 = sym('E3%d%d', [3 ,3]);
E4 = sym('E4%d%d', [3 ,3]);
E5 = sym('E5%d%d', [3 ,3]);

syms x y z w
E = x*E1 + y*E2 + z*E3 +  w*E4 + E5;

eqs = sym(zeros(11, 1));
eqs(1) = det(E);

Et = transpose(E);
te = 2*(E*Et)*E - trace(E*Et)*E;
eqs(2:10) = te(:);

eqs(11) = trace(E);

unknown_vars = symvar([x, y, z, w]);
known_vars = setdiff(symvar(eqs), unknown_vars);
known = {};
for var = known_vars
    known = [known {char(var)}];
end
unknown = {};
for var = unknown_vars
    unknown = [unknown {char(var)}];
end

kngroups = ones(9,1)*[1 2 3 4 5];
[res, export] = gbs_CreateCode('pc_4pst0_nulle', eqs, known, unknown, kngroups);


