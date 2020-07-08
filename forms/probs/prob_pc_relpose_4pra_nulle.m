classdef prob_pc_relpose_4pra_nulle < problem
    properties
        N = (4)
        NE = sym('NE%d%d%d', [[9] - (4), 3, 3]);
        w = sym('w%d', [[9] - (4) - 1, 1]);
    end
    methods
        function [eqs_sym, abbr_subs] = gen_eqs_sym(obj)
            NE = permute(obj.NE, [2, 3, 1]);
            NE = reshape(NE, 9, [9] - obj.N);
            tau = sym('tau');
            E = reshape(NE * [obj.w; 1], 3, 3);
            eqs_sym = sym([]);
            eqs_sym(1) = det(E);
            Et = transpose(E);
            te = 2*(E*Et)*E - trace(E*Et)*E;
            eqs_sym(2:10) = te(:);
            eqs_sym(end + 1) = (tau^2 - 1) * trace(E*Et) + 2 * (tau + 1) * trace(E * E) - 2 * tau * trace(E)^2;

            abbr_subs = struct([]);
        end
        function [in_subs, out_subs] = gen_arg_subs(obj)
            in_subs = struct();
            in_subs.NE = obj.NE;
            in_subs.tau = sym('tau');
            out_subs.w = obj.w;
        end
        function [in_zp, out_zp] = rand_arg_zp(obj, p)
            [in_zp, out_zp] = vzp_rand_relpose_pinhole(obj.N, p, {'NE', 'tau'}, {'w'});
        end
        function [in_rl, out_rl] = rand_arg_rl(obj)
            [in_rl, out_rl] = vrl_rand_relpose_pinhole(obj.N, {'NE', 'tau'}, {'w'});
        end
    end
end