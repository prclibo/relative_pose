classdef prob_pc_relpose_4pra_p6d < problem
    properties
        N = 4;
        q = sym('q%d%d', [3, 4]);
        qq = sym('qq%d%d', [3, 4]);
        u = sym('u%d', [3, 1]);
        s = sym('s');
        t = sym('t%d', [3, 1]);
    end
    methods
        function [eqs_sym, abbr_subs] = gen_eqs_sym(obj)
            R = sym_quat2dcm([obj.s; obj.u]);            
            E = skew(obj.t) * R;
            eqs_sym = sym([]);
            for i = 1:obj.N
                eqs_sym(i) = obj.qq(:, i).' * E * obj.q(:, i);
            end
            eqs_sym(obj.N + 1) = obj.u.' * obj.u + obj.s .* obj.s - 1;
            eqs_sym(obj.N + 2) = obj.t.' * obj.t - 1;

            abbr_subs = struct([]);
        end
        function [in_subs, out_subs] = gen_arg_subs(obj)
            in_subs = struct();
            in_subs.q = obj.q;
            in_subs.qq = obj.qq;
            in_subs.s = obj.s;
            out_subs.u = obj.u;
            out_subs.t = obj.t;
        end
        function [in_zp, out_zp] = rand_arg_zp(obj, p)
            [in_zp, out_zp] = vzp_rand_relpose_pinhole(4, p, {'q', 'qq', 's'}, {'u', 't'});
        end
        function [in_rl, out_rl] = rand_arg_rl(obj)
            [in_rl, out_rl] = vrl_rand_relpose_pinhole(4, {'q', 'qq', 's'}, {'u', 't'});
        end
    end
end