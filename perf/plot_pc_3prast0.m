errs = quantile(transl_errs_pc4pst0(:, 3, 1, :), 0.25, 1);
errs = reshape(errs, [K2, 1]);
errs = repmat(errs, [1, K3]);

ref_errs = quantile(transl_errs_pc5p(:, 3, 1, 1), 0.25, 1);
ref_errs = reshape(ref_errs, [1, 1]);
ref_errs = repmat(ref_errs, [K2, K3]);

mesh(errs), hold on, mesh(ref_errs);

errs = quantile(transl_errs_pc3prast0(:, 3, :, :), 0.25, 1);
errs = reshape(errs, [K2, K3]);
mesh(errs), hold on


ref_errs = quantile(transl_errs_pc5p(:, 3, 1, 1), 0.25, 1);
ref_errs = reshape(ref_errs, [1, 1]);
ref_errs = repmat(ref_errs, [K2, K3]);
mesh(ref_errs);

figure, 

errs = quantile(rot_errs_pc3prast0(:, 3, :, :), 0.25, 1);
errs = reshape(errs, [K2, K3]);

ref_errs = quantile(rot_errs_pc5p(:, 3, 1, 1), 0.25, 1);
ref_errs = reshape(ref_errs, [1, 1]);
ref_errs = repmat(ref_errs, [K2, K3]);

mesh(errs), hold on, mesh(ref_errs);