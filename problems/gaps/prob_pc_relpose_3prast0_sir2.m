classdef prob_pc_relpose_3prast0_sir2__example < problem
    properties
        N = 3;
    end
    methods
        function [in_subs, out_subs] = gen_arg_subs(obj)
            N = obj.N; %#ok<*PROP>
            in_subs.q = sym('q%d%d', [3, N]);
            in_subs.qq = sym('qq%d%d', [3, N]);
            in_subs.s = sym('s');
            out_subs.u = sym('u%d', [3, 1]);
        end
        function [eqs_sym, abbr_subs] = gen_eqs_sym(obj)
            [in, out] = obj.gen_arg_subs();
            u = out.u;
            s = in.s;
            R = 2 * (u * u.' - s * skew(u)) + (s * s - u.' * u) * eye(3);
            
            eqs_sym = sym([]);
            abbr_subs = struct();
            for i = 1:obj.N
                j = mod(i, obj.N) + 1;
                k = mod(i + 1, obj.N) + 1;
                
                [F, F_subs] = obj.Fijk(i, j, k, in.q, in.qq, R);
                F(3, :) = [-transpose(u) * in.q(:, i), transpose(u) * in.qq(:, i)];

                eqs_sym(end + 1) = det(F([1, 2], :));
                eqs_sym(end + 1) = det(F([2, 3], :));
                eqs_sym(end + 1) = det(F([3, 1], :));
                 % CATSTRUCT is a utility function shipped in GAPS
                abbr_subs = catstruct(abbr_subs, F_subs);
            end
            eqs_sym(end + 1) = out.u.' * out.u + in.s .* in.s - 1;
        end
        function [F, abbr_subs] = Fijk(obj, i, j, k, q, qq, R)
            pij = sym(sprintf('p_%d_%d_%%d', i, j), [3, 1]);
            pik = sym(sprintf('p_%d_%d_%%d', i, k), [3, 1]);
            ppij = sym(sprintf('pp_%d_%d_%%d', i, j), [3, 1]);
            ppik = sym(sprintf('pp_%d_%d_%%d', i, k), [3, 1]);
            
            F = [qq(:, j).' * R * pij, ppij.' * R * q(:, j);
                qq(:, k).' * R * pik, ppik.' * R * q(:, k)];
            pij_val = skew(q(:, i)) * q(:, j);
            pik_val = skew(q(:, i)) * q(:, k);
            ppij_val = skew(qq(:, i)) * qq(:, j);
            ppik_val = skew(qq(:, i)) * qq(:, k);
            
            c = num2cell([pij_val; pik_val; ppij_val; ppik_val]);
            f = cellstr(string([pij, pik, ppij, ppik]));
            abbr_subs = cell2struct(c, f);
        end
        function [in_zp, out_zp] = rand_arg_zp(obj, p)
            [in_zp, out_zp] = vzp_rand_relpose_pinhole(3, p, {'q', 'qq', 's'}, {'u'}, 'ZeroScrewTransl', 'inner');
        end
        function [in_rl, out_rl] = rand_arg_rl(obj)
            [in_rl, out_rl] = vrl_rand_relpose_pinhole(3, {'q', 'qq', 's'}, {'u'}, 'ZeroScrewTransl', 'inner');
        end
    end
end