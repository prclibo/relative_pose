classdef prob_pc_relpose_5p_nulle_ne__simple < problem
    % A simple instance to construct a pinhole camera 5-pt relative pose
    % estimation problem solver.
    methods
        function [in, out] = gen_arg_subs(obj)
            % Each field in `in_subs/out_subs` will become an input/output
            % argument in the generated solver.
            %
            % `in_subs.NE` is a 4x3x3 matrix made up of NEijk symbols. The
            % created solver will expect input argument NE to be 4x3x3
            % matrix and fill its element values to NEijk respectively.
            %
            % See your favorite 5-pt paper for details.
            %
            % Base vectors of the null space of the essential matrix.
            in.NE = sym('NE%d%d%d', [4, 3, 3]);
            % The weights
            out.w = sym('w%d', [3, 1]);
        end
        function [eqs_sym, abbr_subs] = gen_eqs_sym(obj)
            [in, out] = gen_arg_subs(obj);
            NE = permute(in.NE, [2, 3, 1]);
            NE = reshape(NE, 9, 4);
            E = reshape(NE * [out.w; 1], 3, 3);
            eqs_sym = sym([]);
            
            % Construct polynomial system as symbolics
            eqs_sym(1) = det(E);
            Et = transpose(E);
            te = 2*(E*Et)*E - trace(E*Et)*E;
            eqs_sym(2:10) = te(:);
            
            % abbr_subs (Abbreviation substitution) is used to declare
            % intermediate variables when computing coefficients and their
            % expansions. This this very simple case we will not use it.
            abbr_subs = struct([]);
        end
        function [in_zp, out_zp] = rand_arg_zp(obj, p)
            % For 5-pt problem, arbitrary value of NE always corresponds to
            % valid polynomial system. Therefore we just instantiate a Zp
            % case by random integer. However, usually this does not hold
            % for a general minimal problem. You need to consider the
            % geometry constraint during instantiation.
            %
            % in_zp/out_zp is expected to have same fields as of in_subs
            % and out_subs returned by `gen_arg_subs`. But the field values
            % are symbolic integers.
            in_zp.NE = sym(randi([1, p - 1], [4, 3, 3]));
            % out_zp can be omitted as it is not used right now.
            out_zp = struct();
        end
        function [in_rl, out_rl] = rand_arg_rl(obj)
            error(['This is similar to `rand_arg_zp` with real field values. ',...
                'We are not using it here as it is for benchmarking']);
        end
    end
end