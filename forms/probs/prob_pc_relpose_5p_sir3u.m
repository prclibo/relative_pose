classdef prob_pc_relpose_5p_sir3u < problem
    properties
        N = 5;
    end
    methods
        function [eqs_sym, abbr_subs] = gen_eqs_sym(obj)
            [in, out] = obj.gen_arg_subs();
            R = sym_quat2dcm([out.s; out.u]);
            obj.config.R = R;
            
            eqs_sym = sym([]);
            abbr_subs = struct();
            for i = 1:obj.N
                for j = i + 1:obj.N
                	for k = j + 1:obj.N
                        G = obj.Gijk(i, j, k, in.q, in.qq, R);
                        eqs_sym(end + 1) = det(G);
                        obj.config.G{i} = G;
                    end
                end
            end
            eqs_sym(end + 1) = out.u.' * out.u + out.s .* out.s - 1;
        end
        function [in_subs, out_subs] = gen_arg_subs(obj)
            in_subs.q = sym('q%d%d', [3, 5]);
            in_subs.qq = sym('qq%d%d', [3, 5]);
            out_subs.s = sym('s');
            out_subs.u = sym('u%d', [3, 1]);
        end
        function [in_zp, out_zp] = rand_arg_zp(obj, p)
            [in_zp, out_zp] = vzp_rand_relpose_pinhole(5, p, {'q', 'qq'}, {'s', 'u'});
        end
        function [in_rl, out_rl] = rand_arg_rl(obj)
            [in_rl, out_rl] = vrl_rand_relpose_pinhole(5, {'q', 'qq'}, {'s', 'u'});
        end
        function G = Gijk(obj, i, j, k, q, qq, R)
            G = [-(qq(:, i).') * skew(R * q(:, i));
                -(qq(:, j).') * skew(R * q(:, j));
                -(qq(:, k).') * skew(R * q(:, k))];
        end
    end
end