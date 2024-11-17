%% create gui and visualizations

% create GUI object;
gui = StandardGUI(sim, [25, 75, 1800, 1000], 0, [0, 0.175, 1.0, 0.8], [8, 22], [0.015, 0.0225], [0, 0, 1.0, 0.05], [1, 6], ...
  elementGroups, elementsInGroup);

gui.addVisualization(MultiPlot({'history lookL','history lookLWB', 'history lookC', 'history lookRWB','history lookR'}, {'output', 'output', 'output', 'output', 'output'}, ...
  [1, 1, 1, 1, 1], 'horizontal', ...
  {'XLim', [-historyDuration, 10], 'YTick', [0, 1], 'YLim', [-0.5, 1.5], 'Box', 'on', 'XGrid', 'on', 'YGrid', 'on', 'fontweight', 'bold', 'fontsize', 16}, ...
  {{'r-', 'LineWidth', 2, 'XData', 0:-1:-historyDuration+1}, ...
  {'r-.', 'LineWidth', 2, 'XData', 0:-1:-historyDuration+1}, ...
  {'g-.', 'LineWidth', 2, 'XData', 0:-1:-historyDuration+1}, ...
  {'b-.', 'LineWidth', 2, 'XData', 0:-1:-historyDuration+1}, ...
  {'b-', 'LineWidth', 2, 'XData', 0:-1:-historyDuration+1}}, ...
  'looking time', 'time (8ms/timestep)'), [1, 16], [1.5,3]);
gui.addVisualization(MultiPlot({'word', 'word'}, {'activation', 'output'}, [1, 10], 'horizontal', ...
  {'XTick', [1 20], 'YTick', [-9 9], 'XLim', [1, fieldSize_wd], 'YLim', [-15, 15], 'Box', 'on', 'fontweight', 'bold', 'fontsize', 16}, {{'Color', 'b'}, {'Color', 'r'}}, 'word'), ...
  [2, 4], [1, 3]);


% % gui.addVisualization(MultiPlot({'hword'}, {'output'}, [1], 'horizontal', ...
% %     {'XLim', [1, fieldSize_wd], 'YLim', [0, 1], 'Box', 'on'}, {{'b', 'LineWidth', 1}}, ...
% %   'hebbian word'), [2, 4], [1, 3] );

%commented to be replaced by hwm_f below % % gui.addVisualization(MultiPlot({'hwm_s'}, {'output'}, [1], 'horizontal', ...
% % %   {'XLim', [1, fieldSize_spt], 'YLim', [0, 1], 'Box', 'on'}, {{'b', 'LineWidth', 1}}, 'hwm s'), ...
% % %   [2, 7], [1, 3]);

% % gui.addVisualization(MultiPlot({['hwm_f' num2str(1)]}, {'output'}, [1], 'horizontal', ...
% %   {'XLim', [1, fieldSize_ftr], 'YLim', [0, 1], 'Box', 'on'}, {{'b', 'LineWidth', 1}}, 'hwm f'), ...
% %   [2, 7], [1, 3]);
gui.addVisualization(MultiPlot({'wm_s', 'wm_s'}, {'activation', 'output'}, [1, 10], 'horizontal', ...
  {'XTick', [1 100], 'YTick', [-9 9], 'XLim', [1, fieldSize_spt], 'YLim', [-15, 15], 'Box', 'on', 'fontweight', 'bold', 'fontsize', 16}, {{'Color', 'b'}, {'Color', 'r'}}, 'spatial working memory'), ...
  [2, 10], [1, 3]);


gui.addVisualization(MultiPlot({'atn_sa', 'atn_sa'}, {'activation', 'output'}, [1, 10], 'horizontal', ...
  {'XTick', [1 100], 'YTick', [-9 9], 'XLim', [1, fieldSize_spt], 'YLim', [-15, 15], 'Box', 'on', 'fontweight', 'bold', 'fontsize', 16}, {{'Color', 'b'}, {'Color', 'r'}}, 'spatial attention (scene)'), ...
  [2, 13], [1, 3]);
gui.addVisualization(MultiPlot({'con_s', 'con_s'}, {'activation', 'output'}, [1, 10], 'horizontal', ...
  {'XTick', [1 100], 'YTick', [-9 9], 'XLim', [1, fieldSize_spt], 'YLim', [-15, 15], 'Box', 'on', 'fontweight', 'bold', 'fontsize', 16}, {{'Color', 'b'}, {'Color', 'r'}}, 'spatial contrast field'), ...
  [1, 10], [1, 3]);
% gui.addVisualization(MultiPlot({'hcon_s'}, {'output'}, [1], 'horizontal', ...
%   {'XLim', [1, fieldSize_spt], 'YLim', [0, 1], 'Box', 'on'}, {{'b', 'LineWidth', 1}}, 'hcon s'), ...
%   [1, 7], [1, 3]);
% % gui.addVisualization(MultiPlot({['hcon_f' num2str(1)]}, {'output'}, [1], 'horizontal', ...
% %   {'XLim', [1, fieldSize_ftr], 'YLim', [0, 1], 'Box', 'on'}, {{'b', 'LineWidth', 1}}, 'hcon f'), ...
% %   [1, 7], [1, 3]);

% gui.addVisualization(MultiPlot({'motor', 'motor', 'mrn', 'mrn', 'cos_m', 'cos_m'}, ...
%   {'activation', 'output', 'activation', 'output', 'activation', 'output'}, ...
%   [1, 10, 1, 10, 1, 10], 'horizontal', ...
%   {'XLim', [1, fieldSize_spt], 'YLim', [-15, 15], 'Box', 'on'}, ...
%   {{'Color', 'b'}, {'Color', 'r'}, ...
%   {'bs', 'MarkerSize', 4, 'XDataMode', 'manual', 'XData', 5}, ...
%   {'rs', 'MarkerSize', 4, 'XDataMode', 'manual', 'XData', 5}, ...
%   {'bo', 'MarkerSize', 4, 'XDataMode', 'manual', 'XData', 10}, ...
%   {'ro', 'MarkerSize', 4, 'XDataMode', 'manual', 'XData', 10}}, ...
%   'motor'), [1, 13], [1, 3]);

gui.addVisualization(MultiPlot({'atn_sr', 'atn_sr'}, {'activation', 'output'}, [1, 10], 'horizontal', ...
  {'XTick', [1 100], 'YTick', [-9 9], 'XLim', [1, fieldSize_spt], 'YLim', [-15, 15], 'Box', 'on', 'fontweight', 'bold', 'fontsize', 16}, {{'Color', 'b'}, {'Color', 'r'}},  'spatial attention (retina)'), ...
  [2, 19], [1, 3]);
gui.addVisualization(MultiPlot({'ior_s', 'ior_s', 'pd_c1', 'pd_c1', 'pd_c2', 'pd_c2', 'cos', 'cos'}, ...
  {'activation', 'output', 'activation', 'output', 'activation', 'output', 'activation', 'output'}, ...
  [1, 10, 1, 10, 1, 10, 1, 10], 'horizontal', ...
  {'XTick', [1 100], 'YTick', [-9 9], 'XLim', [1, fieldSize_spt], 'YLim', [-15, 15], 'Box', 'on', 'fontweight', 'bold', 'fontsize', 16}, ...
  {{'Color', 'b'}, {'Color', 'r'}, ...
  {'bs', 'MarkerSize', 4, 'XDataMode', 'manual', 'XData', 1}, ...
  {'rs', 'MarkerSize', 4, 'XDataMode', 'manual', 'XData', 1}, ...
  {'bs', 'MarkerSize', 4, 'XDataMode', 'manual', 'XData', 5}, ...
  {'rs', 'MarkerSize', 4, 'XDataMode', 'manual', 'XData', 5}, ...
  {'bo', 'MarkerSize', 4, 'XDataMode', 'manual', 'XData', 10}, ...
  {'ro', 'MarkerSize', 4, 'XDataMode', 'manual', 'XData', 10}}, ...
  'inhibition of return (IOR)'), [1, 19], [1, 3]);


for i = 1 : 1
  n = num2str(i);
  if i == nFeatures % x-axis labels only for the bottommost axes
    xla = 'space (scene)';
    xlr = 'space (retina)';
    xlw = 'word';
  else
    xla = [];
    xlr = [];
    xlw = [];
  end
  ov = 3.1*i; % vertical offset for placing visualizations of the current feature
  
  gui.addVisualization(ScaledImage(['wf' n], 'activation', [-7.5, 7.5], ...
    {'XTick', [1 20], 'YTick', [10 300], 'CLim', [-7.5, 7.5], 'YAxisLocation', 'right', 'YDir', 'normal', 'fontweight', 'bold', 'fontsize', 16}, {}, ['word-feature field'], xlw), [ov, 4], [3, 3]);
  gui.addVisualization(ScaledImage(['wm_c' n], 'activation', [-7.5, 7.5], ...
    {'XTick', [1 100], 'YTick', [10 300], 'CLim', [-7.5, 7.5], 'YAxisLocation', 'right', 'YDir', 'normal', 'fontweight', 'bold', 'fontsize', 16}, {}, ['working memory (scene)'], xla), [ov, 10], [3, 3]);
  gui.addVisualization(ScaledImage(['hwm_c' n], 'output', [0, 1.0], ...
    {'XTick', [1 100], 'YTick', [10 300], 'CLim', [0,1.0], 'YAxisLocation', 'right', 'YDir', 'normal', 'fontweight', 'bold', 'fontsize', 16}, {}, ['memory trace (scene)'], xla), [ov, 7], [3, 3]);
  gui.addVisualization(ScaledImage(['hwf' n], 'output', [0, 1.0], ...
    {'XTick', [1 20], 'YTick', [10 300], 'CLim', [0,1.0], 'YAxisLocation', 'right', 'YDir', 'normal', 'fontweight', 'bold', 'fontsize', 16}, {}, ['word-feature memory trace'], xlw), [ov, 1], [3, 3]);
  gui.addVisualization(ScaledImage(['atn_c' n], 'activation', [-7.5, 7.5], ...
    {'XTick', [1 100], 'YTick', [10 300], 'CLim', [-7.5, 7.5], 'YAxisLocation', 'right', 'YDir', 'normal', 'fontweight', 'bold', 'fontsize', 16}, {}, ['attention (scene)'], xla), [ov, 13], [3, 3]);
  gui.addVisualization(MultiPlot({['wm_f' n], ['wm_f' n]}, {'activation', 'output'}, [1, 10], 'vertical', ...
    {'XTick', [-9 9], 'YTick', [10 300], 'YLim', [1, fieldSize_ftr], 'XLim', [-15, 15], 'XDir', 'reverse', 'YAxisLocation', 'right', 'Box', 'on', 'fontweight', 'bold', 'fontsize', 16}, ...
    {{'Color', 'b'}, {'Color', 'r'}}, ['VWM']), [ov, 16], [3, 1]);
  gui.addVisualization(MultiPlot({['con_f' n], ['con_f' n]}, {'activation', 'output'}, [1, 10], 'vertical', ...
    {'XTick', [-9 9], 'YTick', [10 300], 'YLim', [1, fieldSize_ftr], 'XLim', [-15, 15], 'XDir', 'reverse', 'YAxisLocation', 'right', 'Box', 'on', 'fontweight', 'bold', 'fontsize', 16}, ...
    {{'Color', 'b'}, {'Color', 'r'}}, ['F con']), [ov, 17], [3, 1]);
  gui.addVisualization(MultiPlot({['atn_f' n], ['atn_f' n]}, {'activation', 'output'}, [1, 10], 'vertical', ...
    {'XTick', [-9 9], 'YTick', [10 300], 'YLim', [1, fieldSize_ftr], 'XLim', [-15, 15], 'XDir', 'reverse', 'YAxisLocation', 'right', 'Box', 'on', 'fontweight', 'bold', 'fontsize', 16}, ...
    {{'Color', 'b'}, {'Color', 'r'}}, ['F atn']), [ov, 18], [3, 1]);
  
  gui.addVisualization(ScaledImage(['vis_f' n], 'activation', [-7.5, 7.5], ...
    {'XTick', [1 100], 'YTick', [10 300], 'CLim', [-7.5, 7.5], 'YAxisLocation', 'right', 'YDir', 'normal', 'fontweight', 'bold', 'fontsize', 16}, {}, ['visual field'], xlr, ['colour']), ...
    [ov, 19], [3, 3]);
end


for i = 2 : 2
  n = num2str(i);
  if i == nFeatures % x-axis labels only for the bottommost axes
    xla = 'space (scene)';
    xlr = 'space (retina)';
    xlw = 'word';
  else
    xla = [];
    xlr = [];
    xlw = [];
  end
  ov = 3*i; % vertical offset for placing visualizations of the current feature
  
  gui.addVisualization(ScaledImage(['wf' n], 'activation', [-7.5, 7.5], ...
    {'XTick', [1 20], 'YTick', [10 300], 'CLim', [-7.5, 7.5], 'YAxisLocation', 'right', 'YDir', 'normal', 'fontweight', 'bold', 'fontsize', 16}, {}, [''], xlw), [ov, 4], [3, 3]);
  gui.addVisualization(ScaledImage(['wm_c' n], 'activation', [-7.5, 7.5], ...
    {'XTick', [1 100], 'YTick', [10 300], 'CLim', [-7.5, 7.5], 'YAxisLocation', 'right', 'YDir', 'normal', 'fontweight', 'bold', 'fontsize', 16}, {}, [''], xla), [ov, 10], [3, 3]);
  gui.addVisualization(ScaledImage(['hwm_c' n], 'output', [0, 1.0], ...
    {'XTick', [1 100], 'YTick', [10 300], 'CLim', [0,1.0], 'YAxisLocation', 'right', 'YDir', 'normal', 'fontweight', 'bold', 'fontsize', 16}, {}, [''], xla), [ov, 7], [3, 3]);
  gui.addVisualization(ScaledImage(['hwf' n], 'output', [0, 1.0], ...
    {'XTick', [1 20], 'YTick', [10 300], 'CLim', [0,1.0], 'YAxisLocation', 'right', 'YDir', 'normal', 'fontweight', 'bold', 'fontsize', 16}, {}, [''], xlw), [ov, 1], [3, 3]);
  gui.addVisualization(ScaledImage(['atn_c' n], 'activation', [-7.5, 7.5], ...
    {'XTick', [1 100], 'YTick', [10 300], 'CLim', [-7.5, 7.5], 'YAxisLocation', 'right', 'YDir', 'normal', 'fontweight', 'bold', 'fontsize', 16}, {}, [''], xla), [ov, 13], [3, 3]);
  gui.addVisualization(MultiPlot({['wm_f' n], ['wm_f' n]}, {'activation', 'output'}, [1, 10], 'vertical', ...
    {'XTick', [-9 9], 'YTick', [10 300], 'YLim', [1, fieldSize_ftr], 'XLim', [-15, 15], 'XDir', 'reverse', 'YAxisLocation', 'right', 'Box', 'on', 'fontweight', 'bold', 'fontsize', 16}, ...
    {{'Color', 'b'}, {'Color', 'r'}}, ['']), [ov, 16], [3, 1]);
  gui.addVisualization(MultiPlot({['con_f' n], ['con_f' n]}, {'activation', 'output'}, [1, 10], 'vertical', ...
    {'XTick', [-9 9], 'YTick', [10 300], 'YLim', [1, fieldSize_ftr], 'XLim', [-15, 15], 'XDir', 'reverse', 'YAxisLocation', 'right', 'Box', 'on', 'fontweight', 'bold', 'fontsize', 16}, ...
    {{'Color', 'b'}, {'Color', 'r'}}, ['']), [ov, 17], [3, 1]);
  gui.addVisualization(MultiPlot({['atn_f' n], ['atn_f' n]}, {'activation', 'output'}, [1, 10], 'vertical', ...
    {'XTick', [-9 9], 'YTick', [10 300], 'YLim', [1, fieldSize_ftr], 'XLim', [-15, 15], 'XDir', 'reverse', 'YAxisLocation', 'right', 'Box', 'on', 'fontweight', 'bold', 'fontsize', 16}, ...
    {{'Color', 'b'}, {'Color', 'r'}}, ['']), [ov, 18], [3, 1]);
  
  gui.addVisualization(ScaledImage(['vis_f' n], 'activation', [-7.5, 7.5], ...
    {'XTick', [1 100], 'YTick', [10 300], 'CLim', [-7.5, 7.5], 'YAxisLocation', 'right', 'YDir', 'normal', 'fontweight', 'bold', 'fontsize', 16}, {}, [''], xlr, ['shape']), ...
    [ov, 19], [3, 3]);
end


% general control buttons
gui.addVisualization(TimeDisplay(), [1, 5], [1, 1], 'control');
gui.addControl(GlobalControlButton('Pause', gui, 'pauseSimulation', true, false, false, 'pause simulation'), [1, 1]);
gui.addControl(GlobalControlButton('Reset', gui, 'resetSimulation', true, false, true, 'reset simulation'), [1, 2]);
gui.addControl(GlobalControlButton('Parameters', gui, 'paramPanelRequest', true, false, false, 'open parameter panel'), [1, 3]);
gui.addControl(GlobalControlButton('Save', gui, 'saveParameters', true, false, true, 'save parameter settings'), [1, 4]);
gui.addControl(GlobalControlButton('Load', gui, 'loadParameters', true, false, true, 'load parameter settings'), [1, 5]);
gui.addControl(GlobalControlButton('Quit', gui, 'quitSimulation', true, false, false, 'quit simulation'), [1, 6]);
