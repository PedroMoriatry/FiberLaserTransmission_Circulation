%%=====================模拟光纤中信号的传播演化 循环版====================
%主要对象>>一小段时域电场/一个时域上的脉冲
%主要循环>>在空间上轴向的循环传播
%主要输出>>腔内的电场变化 脉宽 带宽 功率等 
%备注：信号变量索引设定（时域/频域坐标,偏振分量1/分量2,轴向位置）
%===================================主程序================================
%===================================主程序================================
function FiberLaserTransmission_Circulation()
%% ============================初始设置===========================
clc,clear,close all
global n_window                %采点数量 建议为奇数 
global n_window_half           %采点数量的一半                  仅用于画图
global t_window_max            %时域积分总跨度
global dt_window               %时域积分步长 
global t_window                %时域坐标的列向量
global fre_window_fft          %频域坐标的列向量 频率fft版
global w_window_fft            %频域坐标的列向量 角频率fft版
global lambda_window_tot       %波长坐标的列向量                仅用于画图
global lambda_window_min       %波长坐标下限                    仅用于画图
global lambda_window_max       %波长坐标上限                    仅用于画图
global n_lambda_window_min     %波长坐标下限对应索引数          仅用于画图
global n_lambda_window_max     %波长坐标上限对应索引数          仅用于画图
global n_lambda                %波长坐标采点数                  仅用于画图
global lambda_window           %波长坐标的列向量 设定范围内     仅用于画图

global boolean_imagenow        %是否实时生成图表                仅用于画图
global boolean_imagefast       %是否绘图简约化                  仅用于画图
global boolean_export          %是否输出输出                    仅用于输出

global n_L                     %腔内演化采样 n                  仅用于画图
global It_L                    %腔内演化数据 It                 仅用于画图
global Ilambda_L               %腔内演化数据 Ilambda            仅用于画图
global E_L                     %腔内演化数据 E                  仅用于画图
global z_L                     %腔内演化数据 z                  仅用于画图

global boolean_HOD0            %是否不考虑三阶色散

global At2center_str

global At

global lambda0                 %中心波长
global fi                      %记录非线性相移
c=3e8;lambda0=2.78e-6;
n_window=80000+1;n_window_half=(n_window+1)/2;
t_window_max=200e-12;dfre_window=1/t_window_max;fre_window_max=n_window*dfre_window;
dt_window=t_window_max/(n_window-1);
t_window=(linspace(-t_window_max/2,+t_window_max/2,n_window))';
fre_window_fft=([linspace(0,+fre_window_max/2,(n_window+1)/2),linspace(-fre_window_max/2+dfre_window,0,(n_window-1)/2)])';
w_window_fft=fre_window_fft*2*pi;
lambda_window_tot=c./fre_window_fft;
lambda_window_min=2.625e-6;lambda_window_max=2.925e-6;%自定义信号部分再次赋值
for i=2:n_window_half
    if c/fre_window_fft(i-1)>lambda_window_min && c/fre_window_fft(i)<lambda_window_min
        n_lambda_window_min=i;n_lambda=n_lambda_window_min+1-n_lambda_window_max;end
    if c/fre_window_fft(i-1)>lambda_window_max && c/fre_window_fft(i)<lambda_window_max
        n_lambda_window_max=i;end
end;lambda_window=lambda_window_tot(n_lambda_window_max:n_lambda_window_min);
if n_lambda < 40 %判断采样点是否足够
    fprintf('>>精度过低%d\n',n_lambda);return;
else 
    fprintf('>>频域有效采样点%d个\n',n_lambda);
end

boolean_imagenow=true;
boolean_imagefast=false;
boolean_HOD0=false;
boolean_export=false;

if boolean_imagenow % 生成画布
    At2center_str='';
    f=figure;set(gcf,'Position',[550,50,700,700]);
    if boolean_imagefast==false %腔内演化数据注册
        set(f,'position',[500,80,1000,700]);
        n_L=0;n_L_pre=1e3;%预注册矩阵规模，最后再整理
        It_L=zeros(n_window,n_L_pre);
        Ilambda_L=zeros(n_lambda,n_L_pre);
        E_L=zeros(n_L_pre,1);
        z_L=zeros(n_L_pre,1);
    end
end
boolean_plotfast2end=false;
boolean_imagenow_origin=boolean_imagenow;
%% ============================初始信号===========================
if 1==0 %% 外导入信号
%待定
end
if 1==1 %% 自定义信号
fre0=c/lambda0;w0=2*pi*fre0;
lambda_window_min=2.60e-6;lambda_window_max=3.00e-6;
tao=300*1e-15;EnergyPulse=4e-9;Pmax=EnergyPulse/tao;A0=sqrt(Pmax);
Et=@(t)...  %exp(-(t/(0.5/log(2)*tao)).^2);%% (1/2/log(2)) gauss脉冲表达式E?
         sech(t/(1/(2*log(sqrt(2)+1))*tao));%%1/(2*log(sqrt(2)+-1)) sech脉冲表达式E
         %1/(1+t.^2/((sqrt(2)+1)/4*tao^2));%%(sqrt(2)+1)/4 lorentz脉冲表达式E
At=A0*Et(t_window).*exp(1i*w0*t_window);
end
if 1==1 %% 是否考虑偏振
    At0=At;
    At(:,1)=At0*sqrt(2)/2;At(:,2)=At0*sqrt(2)/2;
end
%% ============================模拟传播===========================
N=95;
for i=1:N    
OC(0.4);%默认T80%
HWP(38);%0.8 0.17Ge L3.8m Esat1.2e-9 -1 N=95>>38-180 输出端泵浦色散管理孤子
ISO();
QWP(180);
%Ge(0.17);%[输出信号]=Ge[锗棒长度]
fiber(3.0,1.2e-9,+1,1.688e-4*1); %[输出信号]=fiber[光纤长度 饱和能量 输入/输出端泵浦+1/-1 非线性参量]

At2center();fprintf('>>进度 %d/%d\n',i,N);
if boolean_imagenow %实时绘图
    if boolean_imagefast || boolean_plotfast2end==true
        plotfast();
    elseif boolean_plotfast2end==false && boolean_imagenow_origin
        plotdetail();
        if i~=N
            n_L=0;n_L_pre=1e3;%清理演化数据
            It_L=zeros(n_window,n_L_pre);
            Ilambda_L=zeros(n_lambda,n_L_pre);
            E_L=zeros(n_L_pre,1);
            z_L=zeros(n_L_pre,1);
        end
    end
end
if i==N-1
    boolean_imagenow=true;
    if boolean_imagefast;boolean_plotfast2end=true;end
    boolean_imagefast=false;
    n_L=0;n_L_pre=1e3;%预注册矩阵规模
    It_L=zeros(n_window,n_L_pre);
    Ilambda_L=zeros(n_lambda,n_L_pre);
    E_L=zeros(n_L_pre,1);
    z_L=zeros(n_L_pre,1);
end %最终绘图准备
end % END 腔循环

fprintf('>>图表生成中...\n');pause(1);plotfinal();clc;
fprintf('>>模拟结束\n');%fi;

if boolean_export
    fprintf('>>正在导出数据\n');
    
    file_path = 'C:\Users\Administrator\Desktop\校外毕设准备\空心光纤的非线性\脉冲传输数值模拟\脉冲传输数值模拟\振荡源输出数据.xlsx';
    
    writematrix(z_L', file_path, 'Sheet', '传播距离');
    writematrix(E_L', file_path, 'Sheet', '能量演化');
    writematrix(It_L', file_path, 'Sheet', '时域强度');
    writematrix(Ilambda_L', file_path, 'Sheet', '光谱强度');
    fprintf('>>数据已导出至: %s\n', file_path);
    fprintf('>>结束\n');
end

end

%=================================光纤中传播==============================
function []=fiber(L_input,Esat_input,Bpump,gamma_input)

%% 默认参数
global dt_window
global fre_window_fft
global w_window_fft 
global boolean_imagenow
global boolean_imagefast
global n_L
global It_L
global Ilambda_L
global E_L
global z_L

global boolean_HOD0

global At_fiber1
global At_fiber2
global At
global lambda0
global fi
c=3e8;
L=3.7;dL=0.005;                                              %L长度、dL步长
fre0=c/lambda0;w0=2*pi*fre0;                                 %lambda0中心波长
beta2=-0.083e-24;                                            %二阶色散 -0.083e-24(-0.083 ps^2/m)   *详见公式文档
beta3=-0.476e-39;                                             %三阶色散 0.476e-39(0.476e-3 ps^3/m)  *详见公式文档 t>-t 取负
beta4=0;%3.14e-54;                                              %四阶色散 2.97e-54（2.97e-6 ps^4/m)   *详见公式文档
gamma=1.688*1e-4;                                            %非线性参量 1.688e-04[Wm]^-1          *详见公式文档
%TR=5e-15;                                                    %拉曼响应 3-5fs                       *1.55um文献
if boolean_HOD0;beta3=0;beta4=0;end

g=8;                                                         %增益常数                             *常用值 暂无数据
a976=-1.373;                                                 %泵浦光吸收系数 976nm -1.373m-1       *测量值 
Esat=5e-9;                                                  %饱和能量                             *常用值 暂无数据
g_bandwidth=250e-9;fre_bandwidth=g_bandwidth*c/(lambda0^2);  %增益带宽 100nm量级                   *常用值 暂无数据
g_shape_fre=4*((fre_window_fft-fre0)/fre_bandwidth).^2;      %增益谱轮廓                           *近似形式

Bm=1e-6;dbetax=2*pi*Bm/lambda0;dbetay=-2*pi*Bm/lambda0;      %双折射系数 1e-6                      *常用值 暂无数据
%phase_shift_y=0.53*pi;                                       %相移                                 *来源未知
fi=0;
exphD=exp(1i*dL*(beta2/2*(w_window_fft-w0).^2+beta3/6*(w_window_fft-w0).^3+beta4/24*(w_window_fft-w0).^4 ));
exphN=@(At)  exp(1i*dL*gamma*(abs(At).^2));                                       %非线性算子
exphG=@(E,z) exp(dL*0.5*g./(1+E/(Esat*exp(a976*z))+g_shape_fre));                 %增益算子
k_N=@(At1,At2,dbeta,z) ...                                                        %非线性算子 考虑偏振
   1i*gamma*(abs(At1).^2+2/3*abs(At2).^2).*At1+1i*gamma*1/3*(conj(At1).*(At2).^2)*exp(-2*1i*dbeta*z);
%% 自定义参数
if nargin > 0 
if nargin>=1;L=L_input;end
if nargin>=2;Esat=Esat_input;end
if nargin>=3;gamma=gamma_input;end
exphN=@(At)  exp(1i*dL*gamma*(abs(At).^2));                                       %非线性算子 简易
exphG=@(E,z) exp(dL/2*0.5*g./(1+E/(Esat*exp(a976*z))+g_shape_fre));                 %增益算子
k_N=@(At1,At2,dbeta,z) ...                                                        %非线性算子 考虑偏振
   1i*gamma*(abs(At1).^2+2/3*abs(At2).^2).*At1+1i*gamma*1/3*(conj(At1).*(At2).^2)*exp(-2*1i*dbeta*z);
end

%% 模拟传播
At_fiber1=At;                                       % 保存输入光场
At=At*sqrt(0.5);                                    %耦合损耗
[~,p]=size(At);                                     % 判断偏振维度（p=1单偏振，p=2双偏振）
if p==2; Atx=At(:,1);Aty=At(:,2);end                %区分偏振分量
It=At2(At,'It');                                    % 计算时域光强（假设At2为自定义函数）
E=sum(It)*dt_window;                                %预计算能量

z=0;i=1;
while z<=L 
if p==1                                             %不考虑偏振
Aw=fft(At);Aw=exphG(E,L-z).*exphD.*Aw;              % 应用增益和色散
At=ifft(Aw);At=exphN(At).*At;                       % 应用非线性相位
It=abs(At).^2;E=sum(It)*dt_window;
end
if p==2                                             %考虑偏振
    % 频域处理增益和色散
Awx=fft(Atx);Awy=fft(Aty);
if Bpump==1;z_g=z;elseif Bpump==-1;z_g=L-z;end      % 泵浦方向;%Bpump=+1:输入端泵浦;Bpump=-1:输出端泵浦;
Awx=exphG(E,z_g).*exphD.*Awx;                       % x偏振增益+色散
Awy=exphG(E,z_g).*exphD.*Awy;                       % y偏振增益+色散
fi=max(It)*gamma*dL+fi;
Atx=ifft(Awx);Aty=ifft(Awy);                        % 回到时域

% 非线性效应（四阶龙格-库塔法或欧拉法）
if 1==0 %expN RK4
k1x=k_N(Atx,Aty,dbetax,z);  k1y=k_N(Aty,Atx,dbetay,z);  Atx1=k1x*dL/2+Atx;Aty1=k1y*dL/2+Aty;
k2x=k_N(Atx1,Aty1,dbetax,z);k2y=k_N(Aty1,Atx1,dbetay,z);Atx2=k2x*dL/2+Atx;Aty2=k2y*dL/2+Aty;
k3x=k_N(Atx2,Aty2,dbetax,z);k3y=k_N(Aty2,Atx2,dbetay,z);Atx3=k3x*dL+Atx;  Aty3=k3y*dL+Aty;
k4x=k_N(Atx3,Aty3,dbetax,z);k4y=k_N(Aty3,Atx3,dbetay,z);
Atx=dL/6*(k1x+2*k2x+2*k3x+k4x)+Atx;Aty=dL/6*(k1y+2*k2y+2*k3y+k4y)+Aty;
% 计算k1~k4并更新Atx/Aty
else 
    Atx_=k_N(Atx,Aty,dbetax,z)*dL+Atx;Aty_=k_N(Aty,Atx,dbetay,z)*dL+Aty;Atx=Atx_;Aty=Aty_;
end
%Awx=fft(Atx);Awy=fft(Aty);
%Awx=exphG(E,L-z).*exphD.*Awx;Awy=exphG(E,L-z).*exphD.*Awy;
%Atx=ifft(Awx);Aty=ifft(Awy);
It=abs(Atx).^2+abs(Aty).^2;                         % 总光强
E=sum(It)*dt_window;
end

n_L=n_L+1;
if boolean_imagenow && boolean_imagefast==false     %记录腔内演化数据
if p==2; At(:,1)=Atx;At(:,2)=Aty;end
Ilambda=At2(At,'Ilambda');                          % 计算光谱
It_L(:,n_L)=It;                                     % 记录时域光强
Ilambda_L(:,n_L)=Ilambda;                           % 记录光谱
E_L(n_L)=E;                                         % 记录能量
z_L(n_L)=max(z_L)+dL;                               % 记录位置
end
if 1==1                                             % 动态变化dL
if z<=0.01;dL=0.005;                                % 初始小步长
elseif max(It)>2e4;dL=0.01;                         % 高光强时减小步长
%elseif max(It)>5e4;dL=0.005;
else 
    dL=0.1;                                        % 低光强时增大步长
end
exphD=exp( 1i*dL*(beta2/2*(w_window_fft-w0).^2+beta3/6*(w_window_fft-w0).^3+beta4/24*(w_window_fft-w0).^4 ));   %色散相位
exphN=@(At)  exp(1i*dL*gamma*(abs(At).^2));                                                                     %非线性算子
exphG=@(E,z) exp(dL*0.5*g./(1+E/(Esat*exp(a976*z))+g_shape_fre) );                                              %增益算子
end

z=z+dL;i=i+1;
end

%% 整理输出
if p==2;At(:,1)=Atx;At(:,2)=Aty;end
At_fiber2=At;    
end

%=================================锗棒中传播==============================
function []=Ge(L_input)
%% 默认参数
global dt_window
global w_window_fft  
global boolean_imagenow
global boolean_imagefast
global n_L
global It_L
global Ilambda_L
global E_L
global z_L

global boolean_HOD0

global At_Ge1
global At_Ge2
global At
global lambda0
c=3e8;
L=0.08;dL=0.005;                                         %L长度、dL步长
fre0=c/lambda0;w0=2*pi*fre0;                             %lambda0中心波长
beta2=1.69e-24;                                          %二阶色散 1.685e-24(1.685ps^2/m)  控制脉冲展宽（正常色散） 
beta3=-3.397e-39;                                        %三阶色散 3.376e-39(3.376e-3ps^3/m) t>-t 取负  导致脉冲不对称（自陡效应）
beta4=0;%4.6e-54;                                        %四阶色散 4.51e-54(4.51e-6ps^4/m)
if boolean_HOD0;beta3=0;beta4=0;end                      % 高阶色散关闭选项
%% 色散相位补偿算子（频域）
exphD=exp(1i*dL*(beta2/2*(w_window_fft-w0).^2+beta3/6*(w_window_fft-w0).^3+beta4/24*(w_window_fft-w0).^4 ));
%% 自定义参数
if nargin > 0 
if nargin>=1;L=L_input;end
end
%% 模拟传播

At_Ge1=At;
[~,p]=size(At);if p==2; Atx=At(:,1);Aty=At(:,2);end

z=0;i=1;
while z<=L

if p==1 %不考虑偏振
Aw=fft(At);
Aw=exphD.*Aw;
At=ifft(Aw);
end
if p==2 %考虑偏振
Awx=fft(Atx);Awy=fft(Aty);
Awx=exphD.*Awx;Awy=exphD.*Awy;
Atx=ifft(Awx);Aty=ifft(Awy);
end

n_L=n_L+1;
if boolean_imagenow && boolean_imagefast==false     %记录全程数据
if p==2; At(:,1)=Atx;At(:,2)=Aty;end
It=At2(At,'It');                                    % 计算时域光强
E=sum(It)*dt_window;                                % 记录能量
Ilambda=At2(At,'Ilambda');                          % 计算光谱
It_L(:,n_L)=It;                                     % 记录时域光强
Ilambda_L(:,n_L)=Ilambda;                           % 记录光谱
E_L(n_L)=E;
z_L(n_L)=max(z_L)+dL;                               % 记录传播位置
end

z=z+dL;i=i+1;
end

%% 整理输出
if p==2;At(:,1)=Atx;At(:,2)=Aty;end
At_Ge2=At;
end

%================================数据常用处理=============================
function [output]=At2(At,type) % 输出 It Ilambda 
global n_lambda_window_min
global n_lambda_window_max

[~,p]=size(At); 

if strcmp(type,'It')
if p==1 %At不包含偏振
It=abs(At).^2;
end
if p==2 %At包含偏振
Atx=At(:,1);Aty=At(:,2);
It=abs(Atx).^2+abs(Aty).^2;
end
output=It;
end

if strcmp(type,'Ilambda')
if p==1 %At不包含偏振
Iw=abs(fft(At)).^2;
end
if p==2 %At包含偏振
Atx=At(:,1);Aty=At(:,2);
Iw=abs(fft(Atx)).^2+abs(fft(Aty)).^2;
end
output=Iw(n_lambda_window_max:n_lambda_window_min);     %截取频域窗口
end

end

%============================== 最终绘图 自定义===========================
function []=plotfinal()
%% 全局变量
global n_window
global n_lambda
global t_window
global t_window_max
global lambda_window
global n_L                      % 传播步数计数
global It_L                     % 存储传播过程中各位置的时域强度矩阵
global Ilambda_L                % 存储传播过程中各位置的频域强度矩阵
global E_L                      % 存储传播过程中各位置的脉冲能量
global z_L                      % 存储传播距离数组（米）
%% 数据处理
%矩阵整理 去除末尾为0部分
It_L=It_L(:,1:n_L);                     % 裁剪时域强度矩阵到实际传播步数
Ilambda_L=Ilambda_L(:,1:n_L);           % 裁剪频域强度矩阵
E_L=E_L(1:n_L);                         % 裁剪能量数组
z_L=z_L(1:n_L);                         % 裁剪传播距离数组
It_end=It_L(:,n_L);                     % 提取最后一步的时域强度

%时频域强度归一化
for i=1:n_L 
    It_L(:,i)=It_L(:,i)/max(It_L(:,i));
    Ilambda_L(:,i)=Ilambda_L(:,i)/max(Ilambda_L(:,i));
end
Ilambda_end=Ilambda_L(:,n_L);

%矩阵填充
z_It_L=ones(n_window,1)*z_L';                       % 时域传播距离网格（n_window行，n_L列）
z_lambda_L=ones(n_lambda,1)*z_L';                   % 频域传播距离网格（n_lambda行，n_L列）
lambda_window_L=lambda_window*ones(1,n_L);          % 频域波长网格（n_lambda行，n_L列）
t_window_L=t_window*ones(1,n_L);                    % 时域时间网格（n_window行，n_L列）

%时频域FWHM
t_FWHM=FWHM(It_L,'It');                             % 计算时域全宽半高
lambda_FWHM=FWHM(Ilambda_L,'Ilambda');              % 计算频域全宽半高

% 脉宽字符串格式化
if t_FWHM(1,n_L)~=0
    pulsewidth_str=['FWHM ' num2str(t_FWHM(1,n_L)*1e15,'%.0f'),' fs'];      % 飞秒单位
else
    pulsewidth_str='脉冲畸形';                                              % 异常处理
end
%% 绘图
close;
%f=figure;%set(f,'position',[500,80,1000,700]);

% 子图1：时域演化（3D曲面投影为2D）
figure;%subplot(3,2,1);%时域演化
surf(z_It_L,-t_window_L*1e12,It_L);                                         % 绘制时域强度曲面           
shading interp;                                                             % 平滑着色
axis([-inf,inf,-3,3,-inf,inf,0,2]);                                         % 更改设置坐标范围
view(2);                                                                    % 俯视2D投影
xlabel('Position [m]');ylabel('Time Delay [ps]');                           % xy轴标签

% 子图2：频域演化（3D曲面投影为2D）
figure;%subplot(3,2,3);%频域演化
surf(z_lambda_L,lambda_window_L*1e9,Ilambda_L);
shading interp;
axis([-inf,inf,-inf,inf,-inf,inf,0,2]);
view(2);
xlabel('Position [m]');ylabel('Wavelength [nm]');

% 子图3：时频域FWHM曲线（双Y轴）
figure; %subplot(3,2,5)                                                     % 时频域FWHM
[ax,~,h2]=plotyy(z_L,t_FWHM*1e15,z_L,lambda_FWHM*1e9,'plot');
xlabel('Position [m]');
y1_label=get(ax(1),'Ylabel');
set(y1_label,'String','Pulse Duration [fs]');                               % 左Y轴标签
y2_label=get(ax(2),'Ylabel');
set(y2_label,'String','Spectrum FWHM [nm]');                                % 右Y轴标签
set(ax,'ycolor','k');                                                       % 统一Y轴颜色为黑色
set(h2,'color','k');                                                        % 曲线颜色统一\

% 子图4：时域终态波形
figure; %subplot(3,2,2)                                                     % 时域
plot(-t_window*1e12,It_end);
axis([-t_window_max*1e12*0.04,+t_window_max*1e12*0.04,-0.05,inf]);          % 限制时间范围
xlabel('Time Delay [ps]');ylabel('Pulse Power [W]');
legend(pulsewidth_str,'Location', 'northwest');                             % 显示脉宽值
legend('boxoff');

% 子图5：频域终态波形
figure; %subplot(3,2,4)                                                     % 频域
plot(lambda_window*1e9,Ilambda_end);
axis([-inf,inf,-0.01,1.05]);
xlabel('Wavelength [nm]');ylabel('Spectrum Power [counts]');

% 子图6：脉冲能量演化
figure; %subplot(3,2,6)                                                     % 脉冲能量
plot(z_L,E_L*1e9);                                                          % 能量单位为nJ
axis([-inf,inf,-inf,inf]);
xlabel('Position [m]');ylabel('Pulse Energy [nJ]');

drawnow;                                                                    % 立即刷新图形
end
%============================== 实时绘图 详细版 ==========================
function []=plotdetail()
%% 全局变量（新增隔离器相关）
global n_window
global n_window_half
global n_lambda
global t_window
global t_window_max
global lambda_window
%global lambda_window_min
%global lambda_window_max
global n_L
global It_L
global Ilambda_L
global E_L
global z_L 
%global At_fiber1 At_fiber2
%global At_Ge1 At_Ge2
%global At_OC1  At_OC2
%global At_HWP1 At_HWP2
global At_ISO1 At_ISO2                      % 隔离器前后的脉冲数据
%global At_QWP1 At_QWP2

global T_ISO                                % 隔离器实际透过率
global At2center_str                        % 时域对齐状态字符串

%% 数据处理
%矩阵整理 去除末尾为0部分
It_L=It_L(:,1:n_L);
Ilambda_L=Ilambda_L(:,1:n_L);
E_L=E_L(1:n_L);
z_L=z_L(1:n_L);
It_end=It_L(:,n_L);
%时频域强度归一化
for i=1:n_L 
    It_L(:,i)=It_L(:,i)/max(It_L(:,i));
    Ilambda_L(:,i)=Ilambda_L(:,i)/max(Ilambda_L(:,i));
end
Ilambda_end=Ilambda_L(:,n_L);
%矩阵填充
z_It_L=ones(n_window,1)*z_L';
z_lambda_L=ones(n_lambda,1)*z_L';
lambda_window_L=lambda_window*ones(1,n_L);
t_window_L=t_window*ones(1,n_L);
%计算脉宽
pulsewidth=FWHM(It_end,'It');
if pulsewidth~=0
    pulsewidth_str=['FWHM ' num2str(pulsewidth*1e15,'%.0f'),' fs'];
else
    pulsewidth_str='脉冲畸形';
end
%隔离器通光预设
alpha=0/180*pi;                                                             % 隔离器偏振旋转角度（0弧度）
theta=linspace(0,2*pi,100)';                                                % 极角采样

% 绘制理想隔离器通光椭圆（黑色）
a1=1;x1_iso=a1*cos(theta);
b1=1/14;y1_iso=b1*sin(theta);                                               % 椭圆长短轴
temp=[x1_iso,y1_iso]*[cos(-alpha),-sin(alpha);sin(alpha),cos(alpha)];       % 旋转椭圆
x1_iso=temp(:,1);y1_iso=temp(:,2);

% 绘制实际隔离器通光椭圆（蓝色）
x2_iso=real(At_ISO1(n_window_half,1)*exp(1i*theta));                        % 提取隔离器前x偏振分量
y2_iso=real(At_ISO1(n_window_half,2)*exp(1i*theta));                        % 提取y偏振分量
xmax=max(x2_iso);ymax=max(y2_iso);Imax=sqrt(xmax^2+ymax^2);q=1;             % 计算最大强度
x2_iso=x2_iso/Imax*q;y2_iso=y2_iso/Imax*q;                                  % 归一化

% 计算隔离器透过率
It_iso1=At2(At_ISO1,'It');                                                  % 隔离器前总强度
It_iso2=At2(At_ISO2,'It');                                                  % 隔离器后总强度
T_iso=sum(It_iso2)/sum(It_iso1);                                            % 实际透过率
T_iso_str=[num2str(T_iso*100,'%.1f'),'/',num2str(T_ISO*100,'%.1f'),' %'];   % 格式化为字符串
%% 绘图
subplot(3,2,1);%时域演化
surf(z_It_L,-t_window_L*1e12,It_L);
shading interp;
axis([-inf,inf,-2,+2,-inf,inf,0,2]);
view(2);
xlabel('Position [m]');ylabel('Time Delay [ps]');

subplot(3,2,3);%频域演化
surf(z_lambda_L,lambda_window_L*1e9,Ilambda_L);
shading interp;
axis([-inf,inf,2700,2860,-inf,inf,0,2]);
view(2);
xlabel('Position [m]');ylabel('Wavelength [nm]');

subplot(3,2,5);%隔离器通光
plot(x1_iso,y1_iso,'k',x2_iso,y2_iso,'b');                                  % 绘制理想（黑）和实际（蓝）椭圆
axis([-2,2,-1,1]);
xlabel('Polarization of Beam && ISO');
legend(T_iso_str,'Location', 'northeast');                                  % 显示透过率
legend('boxoff');
set(gca,'xtick',[],'ytick',[]);                                             % 隐藏坐标刻度
box off;

subplot(3,2,2);%时域终态
plot(-t_window*1e12,It_end);
axis([-t_window_max*1e12*0.5,+t_window_max*1e12*0.5,-0.05,inf]);
xlabel(['Time Delay [ps]' At2center_str]);ylabel('Pulse Power [W]');
legend(pulsewidth_str,'Location','northwest');
legend('boxoff');

subplot(3,2,4);%频域终态
plot(lambda_window*1e9,Ilambda_end);
axis([-inf,inf,-0.01,1.05]);
xlabel('Wavelength [nm]');ylabel('Spectrum Power [counts]');

subplot(3,2,6);%脉冲能量演化
plot(z_L,E_L*1e9);axis([-inf,inf,-inf,inf]);
xlabel('Position [m]');ylabel('Pulse Energy [nJ]');
drawnow;
end
%============================== 实时绘图 快速版 ==========================
function []=plotfast()
global At                                   % 当前脉冲时域数据
global t_window                             % 时域窗口
global t_window_max                         % 时域窗口最大值
global lambda_window                        % 频域窗口
global At2center_str                        % 时域对齐状态字符串

%处理数据
It_end=At2(At,'It');                                                        % 计算当前时域强度
Ilambda_end=At2(At,'Ilambda');                                              % 计算当前频域强度
Ilambda_end=Ilambda_end/max(Ilambda_end);                                   % 频域归一化

%计算脉宽
pulsewidth=FWHM(It_end,'It');
if pulsewidth~=0
    pulsewidth_str=['FWHM ' num2str(pulsewidth*1e15,'%.0f'),' fs'];
else
    pulsewidth_str='脉冲畸形';
end

% 快速绘图（两子图）
subplot(2,1,1);                                                             % 时域
plot(-t_window*1e12,It_end);
axis([-t_window_max*1e12*0.1,+t_window_max*1e12*0.1,-0.05,inf]);            % 窄时间窗口
xlabel(['Time Delay [ps]' At2center_str]);                                  % 显示对齐状态
ylabel('Pulse Power [W]');        
legend(pulsewidth_str,'Location','northwest');
legend('boxoff');

subplot(2,1,2);                                                             % 频域
plot(lambda_window*1e9,Ilambda_end);
axis([-inf,inf,-0.01,1.05]);
xlabel('Wavelength [nm]');
ylabel('Spectrum Power [counts]');
drawnow;
end

%=============================平移数据至原点两侧==========================
function []=At2center()
%本函数将时域信号的数列顺序微调，将峰值处调整至数列中间位置
global n_window                                                             % 时域窗口点数
global n_window_half                                                        % 时域窗口半长（居中位置）
global At                                                                   % 当前脉冲时域数据
global At2center_str                                                        % 时域对齐状态字符串

[~,p]=size(At);
At_center=zeros(n_window,p);                % 初始化居中后数据
It=At2(At,'It');                            % 计算时域强度

%最高值定位法
if 1==1 
    [~,center]=max(It);% 找到时域强度最大值位置
    % 根据峰值位置循环平移数据
    if center>n_window_half
        n_left=center-n_window_half;
        At_center(1:n_window-n_left,:)=At(n_left+1:n_window,:);             % 左移数据  
        At_center(n_window-n_left+1:n_window,:)=At(1:n_left,:);
    elseif center<n_window_half
        n_right=n_window_half-center;
        At_center(n_right+1:n_window,:)=At(1:n_window-n_right,:);
        At_center(1:n_right,:)=At(n_window-n_right+1:n_window,:);
    else 
        At_center=At;                                                       % 无需调整
    end
end

%半高值定位法
if 1==0 
    It_half=max(It)/2;                                                      % 计算时域强度最大值的一半
    n1_It_half=0;
    n2_It_half=0;
    for i=1:n_window-1
        if (It(i)-It_half)*(It(i+1)-It_half)<=0                             % 检测交点（符号变化）
            n1_It_half=i;break;
        end
    end
    for i=n_window:-1:n_window_half
        if (It(i)-It_half)*(It(i-1)-It_half)<=0                             % 检测交点（符号变化）
            n2_It_half=i;break;
        end
    end
    if n1_It_half==0 && n2_It_half==0
        At_center=At;                                                       % 未找到交点，不调整
    else
        % 防止中点跨越窗口边缘（如脉冲在窗口首尾交界处）
        if abs(n1_It_half-n2_It_half)<n_window_half                         % 防止绘图时曲线在画面边缘
            center=floor((n1_It_half+n2_It_half)/2);                        % 直接取中点
        else
            center=floor((n1_It_half+n2_It_half+n_window)/2);               % 循环修正中点
        end
        if center>n_window_half
            n_left=center-n_window_half;
            At_center(1:n_window-n_left,:)=At(n_left+1:n_window,:);         % 左移
            At_center(n_window-n_left+1:n_window,:)=At(1:n_left,:);
        elseif center<n_window_half
            n_right=n_window_half-center;
            At_center(n_right+1:n_window,:)=At(1:n_window-n_right,:);       % 右移
            At_center(1:n_right,:)=At(n_window-n_right+1:n_window,:);
        else 
            At_center=At;                                                   % 中心已居中
        end
    end
end

At=At_center;                                                               % 更新数据
At2center_str='(实时调整中)';                                               % 设置状态字符串
end
%================================= 计算FWHM ==============================
function [output]=FWHM(input,type)
global t_window                 % 时域窗口
global lambda_window            % 频域窗口
global lambda0                  % 中心波长

if strcmp(type,'It')
    window=t_window;end                                                     % 时域时间轴
if strcmp(type,'Ilambda')
    window=lambda_window;end                                                % 频域波长轴

[m,n]=size(input);
output=zeros(1,n);

% 时域FWHM计算（传统方法 无法处理边带）
if strcmp(type,'It')
for j=1:n
    if max(input(:,j))~=1;input(:,j)=input(:,j)/max(input(:,j));end         % 归一化
    n_half=0;x1=0;
    % 寻找半高宽交点
    for i=2:m
        if (input(i-1,j)-0.5)*(input(i,j)-0.5)<0
            width=window(i)-x1;
            x1=window(i);
            n_half=n_half+1;
        end
    end
    if n_half==2
        output(1,j)=abs(width);                 % 有效FWHM
    end
    if n_half>2
        output(1,j)=0;                          % 无法计算
    end
end
end

% 频域FWHM计算（中心识别法 允许有边带但需轮廓规则）
if strcmp(type,'Ilambda')
    for j=1:n
        Ilambda_max=0;
        % 在中心波长附近找最大值
        for i=1:m
            if window(i)<lambda0+20e-9 && window(i)>lambda0-20e-9
                if input(i,j)>Ilambda_max
                    Ilambda_max=input(i,j);
                    n_Ilambda_max=i;
                end
            end
        end
        % 向右找半高宽点
        for i=n_Ilambda_max:n_Ilambda_max+200
            if (input(i,j)-Ilambda_max/2)*(input(i+1,j)-Ilambda_max/2)<=0
                n_Ilambda_half2=i;
                break;
            end
        end
        % 向左找半高宽点
        for i=n_Ilambda_max-200:n_Ilambda_max
            if (input(i,j)-Ilambda_max/2)*(input(i+1,j)-Ilambda_max/2)<=0
                n_Ilambda_half1=i;
                break;
            end
        end
        output(j)=window(n_Ilambda_half1)-window(n_Ilambda_half2);
    end
end

end

%=============================function 经过OC镜===========================
function []=OC(T_input)
%T 透过率
%R 反射率（损耗率）
global At                           % 当前脉冲数据
global At_OC1                       % 入射到OC前的脉冲
global At_OC2                       % 透射后的脉冲

At_OC1=At;                          % 记录入射前数据
T=0.8;                              % 默认透射率80%
if nargin>0                         % 允许自定义透射率
    T=T_input;
end
At=sqrt((1-T))*At;                  % 透射部分（反射部分被丢弃）
At_OC2=At;                          % 记录透射后数据
end
%============================function 经过1/2波片=========================
function []=HWP(theta)
%本函数模拟经过1/2波片，输入光分为水平与竖直偏振的时域信号 x与y
global At                           % 当前脉冲数据（双偏振）
global At_HWP1                      % 入射到HWP前的脉冲     
global At_HWP2                      % 透射后的脉冲

At_HWP1=At;                         % 记录入射前数据
theta=theta/180*pi;                 % 角度转弧度

% 半波片Jones矩阵作用
temp=[cos(-theta),-sin(-theta);sin(-theta)*exp(-1i*pi),cos(-theta)*exp(-1i*pi)]*At';% 逆时针旋转并施加相位延迟
At=([cos(theta),-sin(theta);sin(theta),cos(theta)]*temp)';                          % 顺时针旋转恢复方向
At_HWP2=At;                                                                         % 记录输出数据
end
%============================function 经过隔离器==========================
function []=ISO(Ty_input)
%本函数模拟经过隔离器 只通过水平方向
global At                           % 当前脉冲数据（双偏振）
global At_ISO1                      % 入射到隔离器前的脉冲
global At_ISO2                      % 透射后的脉冲
global T_ISO                        % 隔离器水平透过率

At_ISO1=At;                         % 记录入射前数据
if nargin>0                         % 允许自定义垂直偏振透过率
    Ty=Ty_input;
else 
    Ty=0.85/200;                    % 默认垂直透过率极低（模拟隔离器单向性）
end

% 水平偏振透射，垂直偏振抑制
At(:,2)=At(:,2)*sqrt(Ty);           % 垂直分量衰减
Tx=0.85;                            % 水平透过率85%
At(:,1)=At(:,1)*sqrt(Tx);
T_ISO=Tx;                           % 记录水平透过率
theta=pi/4;                         %信号偏振转45°
At=At*[cos(theta),-sin(theta);sin(theta),cos(theta)];
At_ISO2=At;
end
%============================function 经过1/4波片=========================
function []=QWP(theta)
%本函数模拟经过1/4波片，输入光分为水平与竖直偏振的时域信号 x与y
global At
global At_QWP1
global At_QWP2

At_QWP1=At;
theta=theta/180*pi;
temp=[cos(-theta),-sin(-theta);sin(-theta)*exp(-1i*pi/2),cos(-theta)*exp(-1i*pi/2)]*At';
At=([cos(theta),-sin(theta);sin(theta),cos(theta)]*temp)';
At_QWP2=At;
end
