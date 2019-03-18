for seq_i = 5:10
    disp(seq_i);
    perf_pc_kitti;
    disp('median error');
    disp(median([t_err_3prast0(:), t_err_4pra(:), t_err_4pst0(:), t_err_5p(:)], 'omitnan'));
    disp('mean error');
    disp(mean([t_err_3prast0(:), t_err_4pra(:), t_err_4pst0(:), t_err_5p(:)], 'omitnan'));
    disp('2pot rate');
    disp(sum(degen_2pot(:, 1)) / sum(isfinite(t_err_4pst0)));
    clear
end
    