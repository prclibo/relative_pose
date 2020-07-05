nstd_px = nstd_ray * focal;

fig_rot_err = figure();
errs = quantile(rot_errs_pc5p(:, :, 1, 1), 0.25, 1);
errs = mean(rot_errs_pc5p(:, :, 1, 1), 1, 'omitnan');
plot(nstd_px, errs, 'x-'), hold on;
fig_transl_err = figure();
errs = quantile(transl_errs_pc5p(:, :, 1, 1), 0.25, 1);
errs = mean(transl_errs_pc5p(:, :, 1, 1), 1, 'omitnan');
plot(nstd_px, errs, 'x-'), hold on;
for k2 = 1:K2
    figure(fig_rot_err);
    errs = quantile(rot_errs_pc4pra(:, :, k2, 1), 0.25, 1);
    errs = mean(rot_errs_pc4pra(:, :, k2, 1), 1, 'omitnan');
    plot(nstd_px, errs, 'g-'), hold on;
    figure(fig_transl_err);
    errs = quantile(transl_errs_pc4pra(:, :, k2, 1), 0.25, 1);
    errs = mean(transl_errs_pc4pra(:, :, k2, 1), 1, 'omitnan');
    plot(nstd_px, errs, 'g-'), hold on;
end
figure(fig_transl_err)
legend('PC-5P', 'PC-4PRA');