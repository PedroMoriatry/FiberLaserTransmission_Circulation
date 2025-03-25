%% ===================== 主程序 FiberLaserTransmission_Circulation_1.m =====================
function FiberLaserTransmission_Circulation_1()
% 初始化参数结构体
[simParams, pulseParams, plotParams] = initialize_parameters();

% 生成初始脉冲
pulseParams = initialize_pulse(simParams, pulseParams);

% 执行腔循环传播
[pulseParams, simParams] = run_cavity_loop(simParams, pulseParams, plotParams);

% 最终绘图和数据导出
if plotParams.boolean_imagenow
    plot_final_results(simParams, pulseParams, plotParams);
end
if plotParams.boolean_export
    export_data(simParams, pulseParams);
end
end
