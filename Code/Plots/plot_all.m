function plot_all()

plot_schematic(1);
plot_schematic(8);
plot_schematicAppen(8);
plot_simulation_visual();
plot_simulation_powers();
plot_simulation_heatmaps('1DHeat');
plot_simulation_heatmaps('HDHeat');
% plot_simulation_heatmaps_pvals();
plot_simulation_powerCompare(1);
plot_simulation_powerCompare(0);
%plot_simulation_slopegraph();
%plot_simulation_outliers();
%plot_simulation_performanceProfiles();
plot_realData();

% plot_simulation_permutation();