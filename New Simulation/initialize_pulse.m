%% ===================== 脉冲初始化函数 initialize_pulse.m =====================
function pulseParams = initialize_pulse(simParams, pulseParams)
c = 3e8;
t_window = linspace(-simParams.t_window_max/2, simParams.t_window_max/2, simParams.n_window)';
A0 = sqrt(4e-9 / 300e-15); % 计算初始振幅（能量4nJ，脉宽300fs）

% 生成双曲正割脉冲
pulseParams.At = A0 * sech(t_window / (1/(2*log(sqrt(2)+1))*300e-15)) .* ...
    exp(1i*2*pi*c/simParams.lambda0 * t_window);

% 双偏振初始化
pulseParams.At = [pulseParams.At*sqrt(2)/2, pulseParams.At*sqrt(2)/2]; 
end
