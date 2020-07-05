classdef prob_pc_relpose_4pra_sir2__example < problem
    properties
        N = 4;
    end
    methods
        function [in_subs, out_subs] = gen_arg_subs(obj)
            in_subs.q = sym('q%d%d', [3, 4]);
            in_subs.qq = sym('qq%d%d', [3, 4]);
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
                eqs_sym(end + 1) = det(F);
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
            quat = zp_rand_unit(4, p);
            s = quat(1);
            u1 = quat(2);
            u2 = quat(3);
            u3 = quat(4);
            u = [u1; u2; u3];
            
            % Make use of some utility functions.
            R = zp_quat2dcm(quat, p);            
            t = zp_rand_unit(3, p);
            
            Q = sym(randi(p, [3, obj.N]));
            % Yes in Matlab 2018 you don't need to REPMAT.
            QQ = mod(R * Q + repmat(t, [1, obj.N]), p);
            
            % We ignore to make Q and QQ homogenous
            in_zp.q = Q; in_zp.qq = QQ;
            in_zp.s = s; out_zp.u = u;
        end
        function [in_rl, out_rl] = rand_arg_rl(obj)
            quat = normc(rand([4, 1]));
            s = quat(1);
            u1 = quat(2);
            u2 = quat(3);
            u3 = quat(4);
            u = [u1; u2; u3];
            
            R = 2 * (u * u.' - s * skew(u)) + (s * s - u.' * u) * eye(3);
            t = normc(rand([3, 1]));
            Q = rand([3, obj.N]);
            QQ = R * Q + repmat(t, [1, obj.N]);
            in_rl.q = Q; in_rl.qq = QQ;
            in_rl.s = s; out_rl.u = u;
        end
    end
end