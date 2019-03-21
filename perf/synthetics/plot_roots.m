
fig = figure();

boxplot([nroots_pc3prast0(:), nroots_pc4pra(:), nroots_pc4pst0(:), nroots_pc5p(:)], 'OutlierSize',1, 'symbol','g.');

lines = findobj(fig, 'type', 'line', 'Tag', 'Median');
set(lines, 'Color', 'g');

xticklabels({'3P', '4PRA', '4PST0', '5P'});
% xtickangle(45)
ylabel('real root count');
set(fig, 'position', [100,100,200,150]);

set(fig,'Units','Inches');
pos = get(fig,'Position');
set(fig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(fig, ['./pc_roots'], '-dpdf', '-r0');

