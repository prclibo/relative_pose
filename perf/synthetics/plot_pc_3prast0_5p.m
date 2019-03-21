% part = 'rot';
part = 'transl';

errs_pc5p = eval(sprintf('%s_errs_pc5p', part));
errs_pc3prast0 = eval(sprintf('%s_errs_pc3prast0', part));

% errs_pc5p = transl_errs_pc5p;
% errs_pc3prast0 = transl_errs_pc3prast0;

fig = figure;
set(fig, 'position', [100,100,200,200]);

h = {};
for level = 1:numel(nstd_ray)

    ref_errs = mean(errs_pc5p(:, level, 1, 1), 1);
    ref_errs = reshape(ref_errs, [1 ,1]);
%     ref_errs = repmat(ref_errs, [K2, K3]);
    
%     mesh(ref_errs), hold on
    
    errs = mean(errs_pc3prast0(:, level, :, :), 1);
    errs = reshape(errs, [K2, K3]);
    
    [r, c] = find(isnan(errs));
    for i = 1:numel(r)
        bottom = min(r + 1, size(errs, 1));
        top = max(r - 1, 1);
        left = max(c - 1, 1);
        right = min(c + 1, size(errs, 2));
        nhood = errs(top:bottom, left:right);
        errs(r, c) = mean(nhood(:), 'omitnan');
    end
    
    errs = imgaussfilt(errs, 0.8);
%     mesh(errs), hold on

    [C, h{level}] = contour(errs, 'levellist', [ref_errs, ref_errs], 'showtext', 'on',...
        'LineColor', 'b');

    hold on

end

% https://undocumentedmatlab.com/blog/customizing-contour-plots
drawnow();

xticklabels(arrayfun(@(x) sprintf('%.2f', x), rad2deg(nstd_angle(2:2:end)), 'uniformoutput', false));
xlabel('angle turb (deg)')
yticklabels(arrayfun(@(x) sprintf('%.0f', x), nstd_transl * 100, 'uniformoutput', false))
ylabel('transl turb (%)')

for i = 1:numel(h)
    for j = 1:numel(h{i}.TextPrims)
        temp_err = str2double(h{i}.TextPrims(j).String);
        h{i}.TextPrims(j).String = sprintf('%.1f (%.1f)',...
            temp_err, nstd_ray(i) * focal);
    end
end
drawnow();


set(fig,'Units','Inches');
pos = get(fig,'Position');
set(fig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(fig, ['./pc_3prast0_5p_', part], '-dpdf', '-r0');
