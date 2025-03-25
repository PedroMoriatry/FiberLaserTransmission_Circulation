%% ===================== 光纤传播函数 propagate_fiber.m =====================
function pulseParams = propagate_fiber(pulseParams, simParams, L, Esat, Bpump)
% 输入参数:
%   L - 光纤长度 (m)
%   Esat - 饱和能量 (J)
%   Bpump - 泵浦方向 (+1:输入端泵浦, -1:输出端泵浦)

% 分步傅里叶法参数
dL = 0.005; % 初始步长
n_steps = round(L / dL);

% 预分配内存
It_L = zeros(simParams.n_window, n_steps);
E_L = zeros(n_steps, 1);

% 主循环
for i = 1:n_steps
    % 频域处理增益和色散
    Aw = fft(pulseParams.At);
    Aw = Aw .* exp(1i * dL * (simParams.beta2/2 * (simParams.w_window_fft - simParams.w0).^2));
    pulseParams.At = ifft(Aw);
    
    % 非线性效应
    pulseParams.At = pulseParams.At .* exp(1i * dL * simParams.gamma * abs(pulseParams.At).^2);
    
    % 记录数据
    It_L(:,i) = abs(pulseParams.At).^2;
    E_L(i) = sum(It_L(:,i)) * simParams.dt_window;
    
    % 动态调整步长
    if max(It_L(:,i)) > 2e4
        dL = 0.01; % 高光强时减小步长
    else
        dL = 0.1;  % 低光强时增大步长
    end
end

% 更新结构体
pulseParams.It_L = It_L;
pulseParams.E_L = E_L;
pulseParams.z_L = cumsum(ones(n_steps,1)*dL);
end