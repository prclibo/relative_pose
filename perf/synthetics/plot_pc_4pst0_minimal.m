nstd_px = nstd_ray * focal;

fig_rot_err = figure();
errs = quantile(rot_errs_pc5p(:, :, 1, 1), 0.25, 1);
plot(nstd_px, errs, 'x-'), hold on;
fig_transl_err = figure();
errs = quantile(transl_errs_pc5p(:, :, 1, 1), 0.25, 1);
plot(nstd_px, errs, 'x-'), hold on;
for k3 = 1:K3
    figure(fig_rot_err);
    errs = quantile(rot_errs_pc4pst0(:, :, 1, k3), 0.25, 1);
    plot(nstd_px, errs, 'g-'), hold on;
    figure(fig_transl_err);
    errs = quantile(transl_errs_pc4pst0(:, :, 1, k3), 0.25, 1);
    plot(nstd_px, errs, 'g-'), hold on;
end
figure(fig_transl_err)
legend('PC-5P', 'PC-4PST0');
