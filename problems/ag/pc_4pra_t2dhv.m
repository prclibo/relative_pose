clear, clc

syms u1 u2 u3 s
u = [u1; u2; u3];
% all_monos = monomialsUpTo(u, 5);
% null_monos = monomialsUpTo(u(1:2), 5);
% 
% numeq = 36;
% C = sym('c%d%d', [numeq, numel(all_monos)]);
% 
% Cz = sym(zeros(numeq, numel(null_monos)));
% eqs = sym(zeros(numeq, 1));
% for i = 1:numeq
%     eqs(i) = C(i, :) * transpose(all_monos);
%     [Cz(i, :), ~] = coeffs(eqs(i), u(1:2));
% end

R = 2 * (u * transpose(u) - s * skew(u)) + (s * s - transpose(u) * u) * eye(3);

% q = transpose(sym('q%d%d', [4, 3]));
% qq = transpose(sym('qq%d%d', [4, 3]));
q = transpose(sym('q%d%d', [4, 3]));
qq = transpose(sym('qq%d%d', [4, 3]));

% gc_zp_data(4, false, [symvar([q, qq]), sym('s')]);
% q = eval(q);
% qq = eval(qq);
% R = eval(R);

eqs = sym(zeros(5, 1));
eqs(1) = det(Fijk(q, qq, R, 2, 3, 4));
eqs(2) = det(Fijk(q, qq, R, 3, 4, 1));
eqs(3) = det(Fijk(q, qq, R, 4, 1, 2));
eqs(4) = det(Fijk(q, qq, R, 1, 2, 3));
eqs(5) = transpose(u) * u + s^2 - 1;

aeqs = sym([]);
mono_d1 = monomialsUpTo(u, 1);
mono_d3 = monomialsUpTo(u, 3);

for mono = mono_d1
    for i = 1:4
        aeqs(end + 1) = eqs(i) * mono;
    end
end
for mono = mono_d3
    aeqs(end + 1) = eqs(5) * mono;
end

monos = monomialsUpTo(u(1:2), 5);
mono_words = arrayfun(@char, monos, 'UniformOutput', false);
mono_orders = cell2mat(...
    arrayfun(@(x) GetMonomialOrder(x{1}, {'u1', 'u2'}),...
    mono_words, 'UniformOutput', false));

order_map = containers.Map(mono_words, mono_orders);

C = sym([]);
for i = 1:numel(aeqs)
    eq = aeqs(i);
    [coeff, term] = coeffs(eq, u(1:2));
    term_words = arrayfun(@char, term, 'UniformOutput', false);
    term_indices = cell2mat(arrayfun(@(x) order_map(x{1}), term_words,...
        'UniformOutput', false));
    
    C(i, term_indices) = coeff;
end

rem = C;
D = sym(zeros([5, size(C)]));
for i = 5:-1:0
    [D(i + 1, :, :), rem] = quorem(rem, u3^i);
end
    

pc_sample_data(4, false, setdiff([symvar([q, qq, R]), transpose(u), s], []));
u_gt = eval(u);


for r = 1:size(D, 2)
    for c = 1:size(D, 3)
        disp([r, c]);
        eval(D(1, r, c))
    end
end

Dins = zeros(size(D));
for c = 1:size(D, 3)
    Dins(:, :, c) = eval(D(:, :, c));
end

vars = setdiff(symvar(D), u);
Dins = vpa(subs(D, vars, eval(vars)), 32);


