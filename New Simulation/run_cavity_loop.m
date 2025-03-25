%% ===================== 腔循环函数 run_cavity_loop.m =====================
function [pulseParams, simParams] = run_cavity_loop(simParams, pulseParams, plotParams)
N = 95; % 循环次数

% 预分配独立变量（确保每个迭代独立操作）
It_L = pulseParams.It_L;       % n_window × N 矩阵
Ilambda_L = pulseParams.Ilambda_L; % n_lambda × N 矩阵
E_L = pulseParams.E_L;         % N × 1 向量
z_L = pulseParams.z_L;         % N × 1 向量

% 开启并行池
if isempty(gcp('nocreate'))
    parpool('local', 4); % 使用4个CPU核心
end

% 并行循环（每个迭代独立处理）
parfor i = 1:N
    % 复制局部变量（每个worker独立）
    localSim = simParams; 
    localPulse = struct(...
        'At', pulseParams.At, ...      % 初始脉冲数据
        'It', It_L(:,i), ...            % 当前迭代的时域强度切片
        'Ilambda', Ilambda_L(:,i), ...  % 当前迭代的频域强度切片
        'E', E_L(i), ...               % 当前迭代的能量值
        'z', z_L(i) ...                % 当前迭代的位置值
    );
    
    % 光学元件操作（使用局部变量）
    localPulse = apply_OC(localPulse, 0.4);     % OC镜透射率40%
    localPulse = apply_HWP(localPulse, 38);     % 半波片角度38度
    localPulse = apply_ISO(localPulse);         % 隔离器
    localPulse = apply_QWP(localPulse, 180);    % 四分之一波片
    localPulse = propagate_Ge(localPulse, localSim, 0.17); % 锗棒传播
    localPulse = propagate_fiber(localPulse, localSim, 3.0, 1.2e-9, +1); % 光纤传播
    
    % 更新切片变量（仅当前迭代的切片）
    It_L(:,i) = localPulse.It;
    Ilambda_L(:,i) = localPulse.Ilambda;
    E_L(i) = localPulse.E;
    z_L(i) = localPulse.z;
end

% 将结果写回结构体
pulseParams.It_L = It_L;
pulseParams.Ilambda_L = Ilambda_L;
pulseParams.E_L = E_L;
pulseParams.z_L = z_L;
end