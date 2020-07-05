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

[monos, mono_words, mono_orders] = monomialsUpTo(u, 5);
order_map = containers.Map(mono_words, mono_orders);

C = sym([]);
for i = 1:numel(aeqs)
    eq = aeqs(i);
    [coeff, term] = coeffs(eq, u);
    term_words = arrayfun(@char, term, 'UniformOutput', false);
    term_indices = cell2mat(arrayfun(@(x) order_map(x{1}), term_words,...
        'UniformOutput', false));
    
    C(i, term_indices) = coeff;
end

[v_monos, ~, ~] = monomialsUpTo(u, 4);
x1_monos = v_monos .* u(1);
x1_words = arrayfun(@char, x1_monos, 'UniformOutput', false);
x1_indices = cell2mat(arrayfun(@(x) order_map(x{1}), x1_words,...
        'UniformOutput', false));

x2_monos = setdiff(monos, x1_monos);
x2_words = arrayfun(@char, x2_monos, 'UniformOutput', false);
x2_indices = cell2mat(arrayfun(@(x) order_map(x{1}), x2_words,...
        'UniformOutput', false));

imap = zeros(1, numel(monos));
imap(x1_indices) = 1:numel(x1_indices);
imap(x2_indices) = 1:numel(x2_indices);

pc_sample_data(4, false, setdiff([symvar([q, qq, R])], []));
u_ = [u1; u2; u3];
clear u1 u2 u3;

C = eval(C);

C1 = C(:, x1_indices);
C2 = C(:, x2_indices);
Bbar = -C2 \ C1;

xv_monos = v_monos .* u(2);
B = zeros(numel(v_monos));
for i = 1:numel(xv_monos)
    mono = xv_monos(i);
    [quo, rem] = quorem(mono, u(1), u(1));
    order = order_map(char(mono));
    index = imap(order);
    if quo == 0
        B(i, :) = Bbar(index, :);
    else
        B(i, index) = 1;
    end
end
    
