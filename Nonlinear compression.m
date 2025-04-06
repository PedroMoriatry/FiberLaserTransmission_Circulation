%% 参数设置
clear; clc;
% 空芯光纤参数
core_radius = 270e-6;   % 纤芯半径270μm
fiber_length = 2.0;     % 光纤长度2m
gas_pressure_max = 1.5; % 最大氦气气压1.5atm
n2_He = 0.4e-23;        % 氦气非线性折射率(m²/W)

% 输入脉冲参数
input_energy = 4e-3;    % 能量4mJ
input_FWHM = 25e-15;    % 脉宽25fs
wavelength = 780e-9;     % 中心波长780nm
c = 3e8;                % 光速
w0 = 2*pi*c / wavelength; % 中心角频率

% 数值计算参数
num_points = 2^14;      % 时域点数
T_window = 10e-12;      % 时域窗口10ps
dt = T_window / num_points;
t = (-num_points/2 : num_points/2-1) * dt; % 时间向量

%% 生成初始高斯脉冲
tau = input_FWHM / (2*sqrt(log(2))); % 高斯脉宽参数
A = sqrt(input_energy) * exp(-t.^2 / (2*tau^2)).'; % 高斯包络
A0 = A; % 保存初始脉冲

%% 光纤传输模拟（分步傅里叶法）
dz = 0.01;              % 步长0.01m
num_steps = fiber_length / dz;
beta2 = 1.038e-6;       % GVD参数(fs²/μm)
beta3 = 0;               % 三阶色散

% 预计算频率向量
omega = 2*pi * (-num_points/2 : num_points/2-1)/(dt*num_points); % 角频率(rad/s)
omega = fftshift(omega); % 调整频率顺序

% 主循环
for step = 1:num_steps
    % 梯度气压计算
    current_pressure = gas_pressure_max * (step*dz)/fiber_length;
    
    % 计算有效模场面积和非线性系数
    A_eff = 0.48*pi*core_radius^2; 
    gamma = (2*pi*n2_He*current_pressure) / (wavelength*A_eff);
    
    % 非线性步骤 (SPM + 自陡峭)
    A = A .* exp(1j * gamma * abs(A).^2 * dz);
    
    % 色散步骤 (频域处理)
    A_fft = fft(A);
    D = -1j*0.5*beta2*omega.^2*dz + (1j/6)*beta3*omega.^3*dz;
    A_fft = A_fft .* exp(D);
    A = ifft(A_fft);
end

%% 脉冲压缩（色散补偿）
% 寻找最优GDD补偿值
GDD_values = linspace(-2000e-30, 0, 100); % 扫描GDD范围
compressed_widths = zeros(size(GDD_values));

for i = 1:length(GDD_values)
    % 应用GDD补偿相位
    phase_comp = -0.5 * GDD_values(i) * omega.^2;
    A_comp = ifft(fft(A) .* exp(1j*phase_comp));
    
    % 计算脉宽
    intensity = abs(A_comp).^2;
    compressed_widths(i) = fwhm(t, intensity);
end

% 找到最佳压缩
[best_width, idx] = min(compressed_widths);
best_GDD = GDD_values(idx);

%% 可视化结果
% 初始和压缩后脉冲对比
subplot(2,2,1);
plot(t*1e15, abs(A0).^2/max(abs(A0).^2), 'b', t*1e15, abs(A).^2/max(abs(A).^2), 'r');
xlabel('Time (fs)'); ylabel('Normalized Intensity');
legend('Initial', 'After Fiber');
title('时域脉冲');

% 频谱展宽对比
subplot(2,2,2);
spectrum_initial = abs(fftshift(fft(A0))).^2;
spectrum_fiber = abs(fftshift(fft(A))).^2;
lambda = 2*pi*c./(omega + w0); % 转换为波长
plot(lambda*1e9, spectrum_initial/max(spectrum_initial), 'b',...
     lambda*1e9, spectrum_fiber/max(spectrum_fiber), 'r');
xlim([600 1000]); 
xlabel('Wavelength (nm)'); ylabel('Normalized Spectrum');
title('频谱展宽');

% GDD扫描结果
subplot(2,2,3);
plot(GDD_values*1e30, compressed_widths*1e15, 'k-o');
xlabel('GDD Compensation (fs²)'); 
ylabel('Pulse Width (fs)');
title('最优GDD选择');
text(best_GDD*1e30, best_width*1e15,...
    sprintf('Best: %.1f fs @ %.0f fs²', best_width*1e15, best_GDD*1e30));

% 压缩后脉冲
subplot(2,2,4);
A_final = ifft(fft(A) .* exp(1j*(-0.5*best_GDD*omega.^2)));
plot(t*1e15, abs(A_final).^2/max(abs(A_final).^2));
xlabel('Time (fs)'); ylabel('Normalized Intensity');
title(sprintf('压缩后脉宽: %.1f fs', best_width*1e15));

%% FWHM计算函数
function width = fwhm(x, y)
    half_max = max(y)/2;
    idx = find(y >= half_max, 1, 'first');
    x1 = x(idx);
    idx = find(y >= half_max, 1, 'last');
    x2 = x(idx);
    width = x2 - x1;
end