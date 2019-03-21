
load num3prast0.mat
load num4pra.mat
load num4pst0.mat
load num4pst0_ess.mat
load num3prast0.mat
load num5p.mat
i = 1:90000;
acc = log10([num3prast0(i)', num4pra(i)', num4pst0_ess(i)', num5p(i)']);

fig = figure();

boxplot(acc, 'OutlierSize',1, 'symbol','g.', 'datalim', [-20, 1]);

lines = findobj(fig, 'type', 'line', 'Tag', 'Median');
set(lines, 'Color', 'g');

xticklabels({'3P', '4PRA', '4PST0', '5P'});
% xtickangle(45)
ylabel('log10 error');
set(fig, 'position', [100,100,200,150]);

set(fig,'Units','Inches');
pos = get(fig,'Position');
set(fig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(fig, ['./pc_numerical'], '-dpdf', '-r0');
