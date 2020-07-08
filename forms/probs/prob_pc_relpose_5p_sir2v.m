classdef prob_pc_relpose_5p_sir2v < problem
    properties
        N = 5;
    end
    methods
        function [eqs_sym, abbr_subs] = gen_eqs_sym(obj)
            [in, out] = obj.gen_arg_subs();
            v = out.v;
            R = 2 * (v * transpose(v) - skew(v)) + (1 - transpose(v) * v) * eye(3);
            
            eqs_sym = sym([]);
            abbr_subs = struct();
            indices = combinator(obj.N,3,'c');
            for row = 1:size(indices, 1)
                i = indices(row, 1);
                j = indices(row, 2);
                k = indices(row, 3);
                
                [F, F_subs] = obj.Fijk(i, j, k, in.q, in.qq, R);
                eqs_sym(end + 1) = det(F);
                abbr_subs = catstruct(abbr_subs, F_subs);
            end
        end
        function [in_subs, out_subs] = gen_arg_subs(obj)
            in_subs.q = sym('q%d%d', [3, 5]);
            in_subs.qq = sym('qq%d%d', [3, 5]);
            out_subs.v = sym('v%d', [3, 1]);
        end
        function [in_zp, out_zp] = rand_arg_zp(obj, p)
            [in_zp, out_zp] = vzp_rand_relpose_pinhole(5, p, {'q', 'qq'}, {'v'});
        end
        function [in_rl, out_rl] = rand_arg_rl(obj)
            [in_rl, out_rl] = vrl_rand_relpose_pinhole(5, {'q', 'qq'}, {'v'});
        end
        function [F, abbr_subs] = Fijk(obj, i, j, k, q, qq, R)
            pij = sym(sprintf('p%d%d%d_%d_%d_%%d', i, j, k, i, j), [3, 1]);
            pik = sym(sprintf('p%d%d%d_%d_%d_%%d', i, j, k, i, k), [3, 1]);
            ppij = sym(sprintf('pp%d%d%d_%d_%d_%%d', i, j, k, i, j), [3, 1]);
            ppik = sym(sprintf('pp%d%d%d_%d_%d_%%d', i, j, k, i, k), [3, 1]);
            
            F = [transpose(qq(:, j)) * R * pij, transpose(ppij) * R * q(:, j);
                transpose(qq(:, k)) * R * pik, transpose(ppik) * R * q(:, k)];
            pij_val = skew(q(:, i)) * q(:, j);
            pik_val = skew(q(:, i)) * q(:, k);
            ppij_val = skew(qq(:, i)) * qq(:, j);
            ppik_val = skew(qq(:, i)) * qq(:, k);
            
            c = num2cell([pij_val; pik_val; ppij_val; ppik_val]);
            f = cellstr(string([pij, pik, ppij, ppik]));
            abbr_subs = cell2struct(c, f);
        end
    end
end