classdef prob_pc_relpose_5p_nulle < problem
    properties
        NE = sym('NE%d%d%d', [4, 3, 3]);
        w = sym('w%d', [3, 1]);
%         E1 = sym('NE1%d%d', [3 ,3]);
%         E2 = sym('NE2%d%d', [3 ,3]);
%         E3 = sym('NE3%d%d', [3 ,3]);
%         E4 = sym('NE4%d%d', [3 ,3]);
    end
    methods
        function [eqs_sym, abbr_subs] = gen_eqs_sym(obj)
            NE = permute(obj.NE, [2, 3, 1]);
            NE = reshape(NE, 9, 4);
            E = reshape(NE * [obj.w; 1], 3, 3);

            eqs_sym = sym([]);
            eqs_sym(1) = det(E);
            Et = transpose(E);
            te = 2*(E*Et)*E - trace(E*Et)*E;
            eqs_sym(2:10) = te(:);

            abbr_subs = struct([]);
        end
        function [in_subs, out_subs] = gen_arg_subs(obj)
            in_subs = struct();
            in_subs.NE = obj.NE;
%             in_subs.E1 = obj.E1;
%             in_subs.E2 = obj.E2;
%             in_subs.E3 = obj.E3;
%             in_subs.E4 = obj.E4;
            out_subs.w = obj.w;
        end
        function [in_zp, out_zp] = rand_arg_zp(obj, p)
            [in_zp, out_zp] = vzp_rand_relpose_pinhole(5, p, {'NE'}, {'w'});
%             in_zp.NE = in_zp.NE + 1;
%             in_zp.NE = randi(10000, size(obj.NE));
            out_zp = struct;
        end
        function [in_rl, out_rl] = rand_arg_rl(obj)
            [in_rl, out_rl] = vrl_rand_relpose_pinhole(5, {'NE'}, {'w'});
        end
    end
end