nstd_px = nstd_ray * focal;

fig_rot_err = figure();
set(fig_rot_err, 'position', [100,100,200,150]);
errs = mean(rot_errs_pc5p(:, :, 1, 1), 1);
plot(nstd_px, errs, 'xr-'), hold on;

fig_transl_err = figure();
set(fig_transl_err, 'position', [100,100,200,150]);
errs = mean(transl_errs_pc5p(:, :, 1, 1), 1);
plot(nstd_px, errs, 'xr-'), hold on;
for k3 = [1,K3]
    figure(fig_rot_err);
%     errs = quantile(rot_errs_pc4pst0(:, :, 1, k3), 0.25, 1);
    errs = mean(rot_errs_pc4pst0(:, :, 1, k3), 1);
    plot(nstd_px, errs, 'g-'), hold on;
    figure(fig_transl_err);
%     errs = quantile(transl_errs_pc4pst0(:, :, 1, k3), 0.25, 1);
    errs = mean(transl_errs_pc4pst0(:, :, 1, k3), 1);
    plot(nstd_px, errs, 'g-'), hold on;
end
figure(fig_transl_err)
legend('PC-5P', 'PC-4PST0');
xticklabels(arrayfun(@(x) sprintf('%.2f', x), rad2deg(nstd_ray), 'uniformoutput', false));
xlabel('ray turb (px)')
% yticklabels(arrayfun(@(x) sprintf('%.1f', x), nstd_transl * 100, 'uniformoutput', false))
ylabel('transl err (%)')

figure(fig_rot_err)
legend('PC-5P', 'PC-4PST0');
xticklabels(arrayfun(@(x) sprintf('%.2f', x), rad2deg(nstd_ray), 'uniformoutput', false));
xlabel('ray turb (px)')
% yticklabels(arrayfun(@(x) sprintf('%.1f', x), nstd_transl * 100, 'uniformoutput', false))
ylabel('transl err (%)')

set(fig_rot_err,'Units','Inches');
pos = get(fig_rot_err,'Position');
set(fig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(fig, ['./pc_4pst0_ransac', part], '-dpdf', '-r0');
set(fig_transl_err,'Units','Inches');
pos = get(fig_transl_err,'Position');
set(fig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(fig, ['./pc_4pst0_ransac', part], '-dpdf', '-r0');
