%% ===================== 参数初始化函数 initialize_parameters.m =====================
function [simParams, pulseParams, plotParams] = initialize_parameters()
% 模拟参数
simParams = struct(...
    'n_window',       80000 + 1, ...       % 时域采样点数
    't_window_max',  200e-12, ...          % 时域窗口范围（秒）
    'lambda0',       2.78e-6, ...          % 中心波长（米）
    'beta2',         -0.083e-24, ...       % 二阶色散系数（s²/m）
    'beta3',         -0.476e-39, ...       % 三阶色散系数（s³/m）
    'gamma',         1.688e-4 ...          % 非线性系数（W⁻¹·m⁻¹）
);

% 脉冲数据
pulseParams = struct(...
    'At',           [], ...               % 时域电场
    'It_L',         [], ...               % 时域强度矩阵
    'Ilambda_L',    [], ...               % 频域强度矩阵
    'E_L',          [], ...               % 能量数组
    'z_L',          [] ...                % 传播距离数组
);

% 绘图与输出控制
plotParams = struct(...
    'boolean_imagenow',  true, ...        % 实时绘图开关
    'boolean_imagefast', false, ...       % 快速绘图模式
    'boolean_export',    false ...        % 数据导出开关
);
end
