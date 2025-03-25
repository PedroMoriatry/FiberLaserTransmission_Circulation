%%=====================ģ��������źŵĴ����ݻ� ѭ����====================
%��Ҫ����>>һС��ʱ��糡/һ��ʱ���ϵ�����
%��Ҫѭ��>>�ڿռ��������ѭ������
%��Ҫ���>>ǻ�ڵĵ糡�仯 ���� ���� ���ʵ� 
%��ע���źű��������趨��ʱ��/Ƶ������,ƫ�����1/����2,����λ�ã�
%===================================������================================
%===================================������================================
function FiberLaserTransmission_Circulation()
%% ============================��ʼ����===========================
clc,clear,close all
global n_window                %�ɵ����� ����Ϊ���� 
global n_window_half           %�ɵ�������һ��                  �����ڻ�ͼ
global t_window_max            %ʱ������ܿ��
global dt_window               %ʱ����ֲ��� 
global t_window                %ʱ�������������
global fre_window_fft          %Ƶ������������� Ƶ��fft��
global w_window_fft            %Ƶ������������� ��Ƶ��fft��
global lambda_window_tot       %���������������                �����ڻ�ͼ
global lambda_window_min       %������������                    �����ڻ�ͼ
global lambda_window_max       %������������                    �����ڻ�ͼ
global n_lambda_window_min     %�����������޶�Ӧ������          �����ڻ�ͼ
global n_lambda_window_max     %�����������޶�Ӧ������          �����ڻ�ͼ
global n_lambda                %��������ɵ���                  �����ڻ�ͼ
global lambda_window           %��������������� �趨��Χ��     �����ڻ�ͼ

global boolean_imagenow        %�Ƿ�ʵʱ����ͼ��                �����ڻ�ͼ
global boolean_imagefast       %�Ƿ��ͼ��Լ��                  �����ڻ�ͼ
global boolean_export          %�Ƿ�������                    ���������

global n_L                     %ǻ���ݻ����� n                  �����ڻ�ͼ
global It_L                    %ǻ���ݻ����� It                 �����ڻ�ͼ
global Ilambda_L               %ǻ���ݻ����� Ilambda            �����ڻ�ͼ
global E_L                     %ǻ���ݻ����� E                  �����ڻ�ͼ
global z_L                     %ǻ���ݻ����� z                  �����ڻ�ͼ

global boolean_HOD0            %�Ƿ񲻿�������ɫɢ

global At2center_str

global At

global lambda0                 %���Ĳ���
global fi                      %��¼����������
c=3e8;lambda0=2.78e-6;
n_window=80000+1;n_window_half=(n_window+1)/2;
t_window_max=200e-12;dfre_window=1/t_window_max;fre_window_max=n_window*dfre_window;
dt_window=t_window_max/(n_window-1);
t_window=(linspace(-t_window_max/2,+t_window_max/2,n_window))';
fre_window_fft=([linspace(0,+fre_window_max/2,(n_window+1)/2),linspace(-fre_window_max/2+dfre_window,0,(n_window-1)/2)])';
w_window_fft=fre_window_fft*2*pi;
lambda_window_tot=c./fre_window_fft;
lambda_window_min=2.625e-6;lambda_window_max=2.925e-6;%�Զ����źŲ����ٴθ�ֵ
for i=2:n_window_half
    if c/fre_window_fft(i-1)>lambda_window_min && c/fre_window_fft(i)<lambda_window_min
        n_lambda_window_min=i;n_lambda=n_lambda_window_min+1-n_lambda_window_max;end
    if c/fre_window_fft(i-1)>lambda_window_max && c/fre_window_fft(i)<lambda_window_max
        n_lambda_window_max=i;end
end;lambda_window=lambda_window_tot(n_lambda_window_max:n_lambda_window_min);
if n_lambda < 40 %�жϲ������Ƿ��㹻
    fprintf('>>���ȹ���%d\n',n_lambda);return;
else 
    fprintf('>>Ƶ����Ч������%d��\n',n_lambda);
end

boolean_imagenow=true;
boolean_imagefast=false;
boolean_HOD0=false;
boolean_export=false;

if boolean_imagenow % ���ɻ���
    At2center_str='';
    f=figure;set(gcf,'Position',[550,50,700,700]);
    if boolean_imagefast==false %ǻ���ݻ�����ע��
        set(f,'position',[500,80,1000,700]);
        n_L=0;n_L_pre=1e3;%Ԥע������ģ�����������
        It_L=zeros(n_window,n_L_pre);
        Ilambda_L=zeros(n_lambda,n_L_pre);
        E_L=zeros(n_L_pre,1);
        z_L=zeros(n_L_pre,1);
    end
end
boolean_plotfast2end=false;
boolean_imagenow_origin=boolean_imagenow;
%% ============================��ʼ�ź�===========================
if 1==0 %% �⵼���ź�
%����
end
if 1==1 %% �Զ����ź�
fre0=c/lambda0;w0=2*pi*fre0;
lambda_window_min=2.60e-6;lambda_window_max=3.00e-6;
tao=300*1e-15;EnergyPulse=4e-9;Pmax=EnergyPulse/tao;A0=sqrt(Pmax);
Et=@(t)...  %exp(-(t/(0.5/log(2)*tao)).^2);%% (1/2/log(2)) gauss������ʽE?
         sech(t/(1/(2*log(sqrt(2)+1))*tao));%%1/(2*log(sqrt(2)+-1)) sech������ʽE
         %1/(1+t.^2/((sqrt(2)+1)/4*tao^2));%%(sqrt(2)+1)/4 lorentz������ʽE
At=A0*Et(t_window).*exp(1i*w0*t_window);
end
if 1==1 %% �Ƿ���ƫ��
    At0=At;
    At(:,1)=At0*sqrt(2)/2;At(:,2)=At0*sqrt(2)/2;
end
%% ============================ģ�⴫��===========================
N=95;
for i=1:N    
OC(0.4);%Ĭ��T80%
HWP(38);%0.8 0.17Ge L3.8m Esat1.2e-9 -1 N=95>>38-180 ����˱���ɫɢ�������
ISO();
QWP(180);
%Ge(0.17);%[����ź�]=Ge[�������]
fiber(3.0,1.2e-9,+1,1.688e-4*1); %[����ź�]=fiber[���˳��� �������� ����/����˱���+1/-1 �����Բ���]

At2center();fprintf('>>���� %d/%d\n',i,N);
if boolean_imagenow %ʵʱ��ͼ
    if boolean_imagefast || boolean_plotfast2end==true
        plotfast();
    elseif boolean_plotfast2end==false && boolean_imagenow_origin
        plotdetail();
        if i~=N
            n_L=0;n_L_pre=1e3;%�����ݻ�����
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
    n_L=0;n_L_pre=1e3;%Ԥע������ģ
    It_L=zeros(n_window,n_L_pre);
    Ilambda_L=zeros(n_lambda,n_L_pre);
    E_L=zeros(n_L_pre,1);
    z_L=zeros(n_L_pre,1);
end %���ջ�ͼ׼��
end % END ǻѭ��

fprintf('>>ͼ��������...\n');pause(1);plotfinal();clc;
fprintf('>>ģ�����\n');%fi;

if boolean_export
    fprintf('>>���ڵ�������\n');
    
    file_path = 'C:\Users\Administrator\Desktop\У�����׼��\���Ĺ��˵ķ�����\���崫����ֵģ��\���崫����ֵģ��\��Դ�������.xlsx';
    
    writematrix(z_L', file_path, 'Sheet', '��������');
    writematrix(E_L', file_path, 'Sheet', '�����ݻ�');
    writematrix(It_L', file_path, 'Sheet', 'ʱ��ǿ��');
    writematrix(Ilambda_L', file_path, 'Sheet', '����ǿ��');
    fprintf('>>�����ѵ�����: %s\n', file_path);
    fprintf('>>����\n');
end

end

%=================================�����д���==============================
function []=fiber(L_input,Esat_input,Bpump,gamma_input)

%% Ĭ�ϲ���
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
L=3.7;dL=0.005;                                              %L���ȡ�dL����
fre0=c/lambda0;w0=2*pi*fre0;                                 %lambda0���Ĳ���
beta2=-0.083e-24;                                            %����ɫɢ -0.083e-24(-0.083 ps^2/m)   *�����ʽ�ĵ�
beta3=-0.476e-39;                                             %����ɫɢ 0.476e-39(0.476e-3 ps^3/m)  *�����ʽ�ĵ� t>-t ȡ��
beta4=0;%3.14e-54;                                              %�Ľ�ɫɢ 2.97e-54��2.97e-6 ps^4/m)   *�����ʽ�ĵ�
gamma=1.688*1e-4;                                            %�����Բ��� 1.688e-04[Wm]^-1          *�����ʽ�ĵ�
%TR=5e-15;                                                    %������Ӧ 3-5fs                       *1.55um����
if boolean_HOD0;beta3=0;beta4=0;end

g=8;                                                         %���泣��                             *����ֵ ��������
a976=-1.373;                                                 %���ֹ�����ϵ�� 976nm -1.373m-1       *����ֵ 
Esat=5e-9;                                                  %��������                             *����ֵ ��������
g_bandwidth=250e-9;fre_bandwidth=g_bandwidth*c/(lambda0^2);  %������� 100nm����                   *����ֵ ��������
g_shape_fre=4*((fre_window_fft-fre0)/fre_bandwidth).^2;      %����������                           *������ʽ

Bm=1e-6;dbetax=2*pi*Bm/lambda0;dbetay=-2*pi*Bm/lambda0;      %˫����ϵ�� 1e-6                      *����ֵ ��������
%phase_shift_y=0.53*pi;                                       %����                                 *��Դδ֪
fi=0;
exphD=exp(1i*dL*(beta2/2*(w_window_fft-w0).^2+beta3/6*(w_window_fft-w0).^3+beta4/24*(w_window_fft-w0).^4 ));
exphN=@(At)  exp(1i*dL*gamma*(abs(At).^2));                                       %����������
exphG=@(E,z) exp(dL*0.5*g./(1+E/(Esat*exp(a976*z))+g_shape_fre));                 %��������
k_N=@(At1,At2,dbeta,z) ...                                                        %���������� ����ƫ��
   1i*gamma*(abs(At1).^2+2/3*abs(At2).^2).*At1+1i*gamma*1/3*(conj(At1).*(At2).^2)*exp(-2*1i*dbeta*z);
%% �Զ������
if nargin > 0 
if nargin>=1;L=L_input;end
if nargin>=2;Esat=Esat_input;end
if nargin>=3;gamma=gamma_input;end
exphN=@(At)  exp(1i*dL*gamma*(abs(At).^2));                                       %���������� ����
exphG=@(E,z) exp(dL/2*0.5*g./(1+E/(Esat*exp(a976*z))+g_shape_fre));                 %��������
k_N=@(At1,At2,dbeta,z) ...                                                        %���������� ����ƫ��
   1i*gamma*(abs(At1).^2+2/3*abs(At2).^2).*At1+1i*gamma*1/3*(conj(At1).*(At2).^2)*exp(-2*1i*dbeta*z);
end

%% ģ�⴫��
At_fiber1=At;                                       % ��������ⳡ
At=At*sqrt(0.5);                                    %������
[~,p]=size(At);                                     % �ж�ƫ��ά�ȣ�p=1��ƫ��p=2˫ƫ��
if p==2; Atx=At(:,1);Aty=At(:,2);end                %����ƫ�����
It=At2(At,'It');                                    % ����ʱ���ǿ������At2Ϊ�Զ��庯����
E=sum(It)*dt_window;                                %Ԥ��������

z=0;i=1;
while z<=L 
if p==1                                             %������ƫ��
Aw=fft(At);Aw=exphG(E,L-z).*exphD.*Aw;              % Ӧ�������ɫɢ
At=ifft(Aw);At=exphN(At).*At;                       % Ӧ�÷�������λ
It=abs(At).^2;E=sum(It)*dt_window;
end
if p==2                                             %����ƫ��
    % Ƶ���������ɫɢ
Awx=fft(Atx);Awy=fft(Aty);
if Bpump==1;z_g=z;elseif Bpump==-1;z_g=L-z;end      % ���ַ���;%Bpump=+1:����˱���;Bpump=-1:����˱���;
Awx=exphG(E,z_g).*exphD.*Awx;                       % xƫ������+ɫɢ
Awy=exphG(E,z_g).*exphD.*Awy;                       % yƫ������+ɫɢ
fi=max(It)*gamma*dL+fi;
Atx=ifft(Awx);Aty=ifft(Awy);                        % �ص�ʱ��

% ������ЧӦ���Ľ�����-��������ŷ������
if 1==0 %expN RK4
k1x=k_N(Atx,Aty,dbetax,z);  k1y=k_N(Aty,Atx,dbetay,z);  Atx1=k1x*dL/2+Atx;Aty1=k1y*dL/2+Aty;
k2x=k_N(Atx1,Aty1,dbetax,z);k2y=k_N(Aty1,Atx1,dbetay,z);Atx2=k2x*dL/2+Atx;Aty2=k2y*dL/2+Aty;
k3x=k_N(Atx2,Aty2,dbetax,z);k3y=k_N(Aty2,Atx2,dbetay,z);Atx3=k3x*dL+Atx;  Aty3=k3y*dL+Aty;
k4x=k_N(Atx3,Aty3,dbetax,z);k4y=k_N(Aty3,Atx3,dbetay,z);
Atx=dL/6*(k1x+2*k2x+2*k3x+k4x)+Atx;Aty=dL/6*(k1y+2*k2y+2*k3y+k4y)+Aty;
% ����k1~k4������Atx/Aty
else 
    Atx_=k_N(Atx,Aty,dbetax,z)*dL+Atx;Aty_=k_N(Aty,Atx,dbetay,z)*dL+Aty;Atx=Atx_;Aty=Aty_;
end
%Awx=fft(Atx);Awy=fft(Aty);
%Awx=exphG(E,L-z).*exphD.*Awx;Awy=exphG(E,L-z).*exphD.*Awy;
%Atx=ifft(Awx);Aty=ifft(Awy);
It=abs(Atx).^2+abs(Aty).^2;                         % �ܹ�ǿ
E=sum(It)*dt_window;
end

n_L=n_L+1;
if boolean_imagenow && boolean_imagefast==false     %��¼ǻ���ݻ�����
if p==2; At(:,1)=Atx;At(:,2)=Aty;end
Ilambda=At2(At,'Ilambda');                          % �������
It_L(:,n_L)=It;                                     % ��¼ʱ���ǿ
Ilambda_L(:,n_L)=Ilambda;                           % ��¼����
E_L(n_L)=E;                                         % ��¼����
z_L(n_L)=max(z_L)+dL;                               % ��¼λ��
end
if 1==1                                             % ��̬�仯dL
if z<=0.01;dL=0.005;                                % ��ʼС����
elseif max(It)>2e4;dL=0.01;                         % �߹�ǿʱ��С����
%elseif max(It)>5e4;dL=0.005;
else 
    dL=0.1;                                        % �͹�ǿʱ���󲽳�
end
exphD=exp( 1i*dL*(beta2/2*(w_window_fft-w0).^2+beta3/6*(w_window_fft-w0).^3+beta4/24*(w_window_fft-w0).^4 ));   %ɫɢ��λ
exphN=@(At)  exp(1i*dL*gamma*(abs(At).^2));                                                                     %����������
exphG=@(E,z) exp(dL*0.5*g./(1+E/(Esat*exp(a976*z))+g_shape_fre) );                                              %��������
end

z=z+dL;i=i+1;
end

%% �������
if p==2;At(:,1)=Atx;At(:,2)=Aty;end
At_fiber2=At;    
end

%=================================����д���==============================
function []=Ge(L_input)
%% Ĭ�ϲ���
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
L=0.08;dL=0.005;                                         %L���ȡ�dL����
fre0=c/lambda0;w0=2*pi*fre0;                             %lambda0���Ĳ���
beta2=1.69e-24;                                          %����ɫɢ 1.685e-24(1.685ps^2/m)  ��������չ������ɫɢ�� 
beta3=-3.397e-39;                                        %����ɫɢ 3.376e-39(3.376e-3ps^3/m) t>-t ȡ��  �������岻�Գƣ��Զ�ЧӦ��
beta4=0;%4.6e-54;                                        %�Ľ�ɫɢ 4.51e-54(4.51e-6ps^4/m)
if boolean_HOD0;beta3=0;beta4=0;end                      % �߽�ɫɢ�ر�ѡ��
%% ɫɢ��λ�������ӣ�Ƶ��
exphD=exp(1i*dL*(beta2/2*(w_window_fft-w0).^2+beta3/6*(w_window_fft-w0).^3+beta4/24*(w_window_fft-w0).^4 ));
%% �Զ������
if nargin > 0 
if nargin>=1;L=L_input;end
end
%% ģ�⴫��

At_Ge1=At;
[~,p]=size(At);if p==2; Atx=At(:,1);Aty=At(:,2);end

z=0;i=1;
while z<=L

if p==1 %������ƫ��
Aw=fft(At);
Aw=exphD.*Aw;
At=ifft(Aw);
end
if p==2 %����ƫ��
Awx=fft(Atx);Awy=fft(Aty);
Awx=exphD.*Awx;Awy=exphD.*Awy;
Atx=ifft(Awx);Aty=ifft(Awy);
end

n_L=n_L+1;
if boolean_imagenow && boolean_imagefast==false     %��¼ȫ������
if p==2; At(:,1)=Atx;At(:,2)=Aty;end
It=At2(At,'It');                                    % ����ʱ���ǿ
E=sum(It)*dt_window;                                % ��¼����
Ilambda=At2(At,'Ilambda');                          % �������
It_L(:,n_L)=It;                                     % ��¼ʱ���ǿ
Ilambda_L(:,n_L)=Ilambda;                           % ��¼����
E_L(n_L)=E;
z_L(n_L)=max(z_L)+dL;                               % ��¼����λ��
end

z=z+dL;i=i+1;
end

%% �������
if p==2;At(:,1)=Atx;At(:,2)=Aty;end
At_Ge2=At;
end

%================================���ݳ��ô���=============================
function [output]=At2(At,type) % ��� It Ilambda 
global n_lambda_window_min
global n_lambda_window_max

[~,p]=size(At); 

if strcmp(type,'It')
if p==1 %At������ƫ��
It=abs(At).^2;
end
if p==2 %At����ƫ��
Atx=At(:,1);Aty=At(:,2);
It=abs(Atx).^2+abs(Aty).^2;
end
output=It;
end

if strcmp(type,'Ilambda')
if p==1 %At������ƫ��
Iw=abs(fft(At)).^2;
end
if p==2 %At����ƫ��
Atx=At(:,1);Aty=At(:,2);
Iw=abs(fft(Atx)).^2+abs(fft(Aty)).^2;
end
output=Iw(n_lambda_window_max:n_lambda_window_min);     %��ȡƵ�򴰿�
end

end

%============================== ���ջ�ͼ �Զ���===========================
function []=plotfinal()
%% ȫ�ֱ���
global n_window
global n_lambda
global t_window
global t_window_max
global lambda_window
global n_L                      % ������������
global It_L                     % �洢���������и�λ�õ�ʱ��ǿ�Ⱦ���
global Ilambda_L                % �洢���������и�λ�õ�Ƶ��ǿ�Ⱦ���
global E_L                      % �洢���������и�λ�õ���������
global z_L                      % �洢�����������飨�ף�
%% ���ݴ���
%�������� ȥ��ĩβΪ0����
It_L=It_L(:,1:n_L);                     % �ü�ʱ��ǿ�Ⱦ���ʵ�ʴ�������
Ilambda_L=Ilambda_L(:,1:n_L);           % �ü�Ƶ��ǿ�Ⱦ���
E_L=E_L(1:n_L);                         % �ü���������
z_L=z_L(1:n_L);                         % �ü�������������
It_end=It_L(:,n_L);                     % ��ȡ���һ����ʱ��ǿ��

%ʱƵ��ǿ�ȹ�һ��
for i=1:n_L 
    It_L(:,i)=It_L(:,i)/max(It_L(:,i));
    Ilambda_L(:,i)=Ilambda_L(:,i)/max(Ilambda_L(:,i));
end
Ilambda_end=Ilambda_L(:,n_L);

%�������
z_It_L=ones(n_window,1)*z_L';                       % ʱ�򴫲���������n_window�У�n_L�У�
z_lambda_L=ones(n_lambda,1)*z_L';                   % Ƶ�򴫲���������n_lambda�У�n_L�У�
lambda_window_L=lambda_window*ones(1,n_L);          % Ƶ�򲨳�����n_lambda�У�n_L�У�
t_window_L=t_window*ones(1,n_L);                    % ʱ��ʱ������n_window�У�n_L�У�

%ʱƵ��FWHM
t_FWHM=FWHM(It_L,'It');                             % ����ʱ��ȫ����
lambda_FWHM=FWHM(Ilambda_L,'Ilambda');              % ����Ƶ��ȫ����

% �����ַ�����ʽ��
if t_FWHM(1,n_L)~=0
    pulsewidth_str=['FWHM ' num2str(t_FWHM(1,n_L)*1e15,'%.0f'),' fs'];      % ���뵥λ
else
    pulsewidth_str='�������';                                              % �쳣����
end
%% ��ͼ
close;
%f=figure;%set(f,'position',[500,80,1000,700]);

% ��ͼ1��ʱ���ݻ���3D����ͶӰΪ2D��
figure;%subplot(3,2,1);%ʱ���ݻ�
surf(z_It_L,-t_window_L*1e12,It_L);                                         % ����ʱ��ǿ������           
shading interp;                                                             % ƽ����ɫ
axis([-inf,inf,-3,3,-inf,inf,0,2]);                                         % �����������귶Χ
view(2);                                                                    % ����2DͶӰ
xlabel('Position [m]');ylabel('Time Delay [ps]');                           % xy���ǩ

% ��ͼ2��Ƶ���ݻ���3D����ͶӰΪ2D��
figure;%subplot(3,2,3);%Ƶ���ݻ�
surf(z_lambda_L,lambda_window_L*1e9,Ilambda_L);
shading interp;
axis([-inf,inf,-inf,inf,-inf,inf,0,2]);
view(2);
xlabel('Position [m]');ylabel('Wavelength [nm]');

% ��ͼ3��ʱƵ��FWHM���ߣ�˫Y�ᣩ
figure; %subplot(3,2,5)                                                     % ʱƵ��FWHM
[ax,~,h2]=plotyy(z_L,t_FWHM*1e15,z_L,lambda_FWHM*1e9,'plot');
xlabel('Position [m]');
y1_label=get(ax(1),'Ylabel');
set(y1_label,'String','Pulse Duration [fs]');                               % ��Y���ǩ
y2_label=get(ax(2),'Ylabel');
set(y2_label,'String','Spectrum FWHM [nm]');                                % ��Y���ǩ
set(ax,'ycolor','k');                                                       % ͳһY����ɫΪ��ɫ
set(h2,'color','k');                                                        % ������ɫͳһ\

% ��ͼ4��ʱ����̬����
figure; %subplot(3,2,2)                                                     % ʱ��
plot(-t_window*1e12,It_end);
axis([-t_window_max*1e12*0.04,+t_window_max*1e12*0.04,-0.05,inf]);          % ����ʱ�䷶Χ
xlabel('Time Delay [ps]');ylabel('Pulse Power [W]');
legend(pulsewidth_str,'Location', 'northwest');                             % ��ʾ����ֵ
legend('boxoff');

% ��ͼ5��Ƶ����̬����
figure; %subplot(3,2,4)                                                     % Ƶ��
plot(lambda_window*1e9,Ilambda_end);
axis([-inf,inf,-0.01,1.05]);
xlabel('Wavelength [nm]');ylabel('Spectrum Power [counts]');

% ��ͼ6�����������ݻ�
figure; %subplot(3,2,6)                                                     % ��������
plot(z_L,E_L*1e9);                                                          % ������λΪnJ
axis([-inf,inf,-inf,inf]);
xlabel('Position [m]');ylabel('Pulse Energy [nJ]');

drawnow;                                                                    % ����ˢ��ͼ��
end
%============================== ʵʱ��ͼ ��ϸ�� ==========================
function []=plotdetail()
%% ȫ�ֱ�����������������أ�
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
global At_ISO1 At_ISO2                      % ������ǰ�����������
%global At_QWP1 At_QWP2

global T_ISO                                % ������ʵ��͸����
global At2center_str                        % ʱ�����״̬�ַ���

%% ���ݴ���
%�������� ȥ��ĩβΪ0����
It_L=It_L(:,1:n_L);
Ilambda_L=Ilambda_L(:,1:n_L);
E_L=E_L(1:n_L);
z_L=z_L(1:n_L);
It_end=It_L(:,n_L);
%ʱƵ��ǿ�ȹ�һ��
for i=1:n_L 
    It_L(:,i)=It_L(:,i)/max(It_L(:,i));
    Ilambda_L(:,i)=Ilambda_L(:,i)/max(Ilambda_L(:,i));
end
Ilambda_end=Ilambda_L(:,n_L);
%�������
z_It_L=ones(n_window,1)*z_L';
z_lambda_L=ones(n_lambda,1)*z_L';
lambda_window_L=lambda_window*ones(1,n_L);
t_window_L=t_window*ones(1,n_L);
%��������
pulsewidth=FWHM(It_end,'It');
if pulsewidth~=0
    pulsewidth_str=['FWHM ' num2str(pulsewidth*1e15,'%.0f'),' fs'];
else
    pulsewidth_str='�������';
end
%������ͨ��Ԥ��
alpha=0/180*pi;                                                             % ������ƫ����ת�Ƕȣ�0���ȣ�
theta=linspace(0,2*pi,100)';                                                % ���ǲ���

% �������������ͨ����Բ����ɫ��
a1=1;x1_iso=a1*cos(theta);
b1=1/14;y1_iso=b1*sin(theta);                                               % ��Բ������
temp=[x1_iso,y1_iso]*[cos(-alpha),-sin(alpha);sin(alpha),cos(alpha)];       % ��ת��Բ
x1_iso=temp(:,1);y1_iso=temp(:,2);

% ����ʵ�ʸ�����ͨ����Բ����ɫ��
x2_iso=real(At_ISO1(n_window_half,1)*exp(1i*theta));                        % ��ȡ������ǰxƫ�����
y2_iso=real(At_ISO1(n_window_half,2)*exp(1i*theta));                        % ��ȡyƫ�����
xmax=max(x2_iso);ymax=max(y2_iso);Imax=sqrt(xmax^2+ymax^2);q=1;             % �������ǿ��
x2_iso=x2_iso/Imax*q;y2_iso=y2_iso/Imax*q;                                  % ��һ��

% ���������͸����
It_iso1=At2(At_ISO1,'It');                                                  % ������ǰ��ǿ��
It_iso2=At2(At_ISO2,'It');                                                  % ����������ǿ��
T_iso=sum(It_iso2)/sum(It_iso1);                                            % ʵ��͸����
T_iso_str=[num2str(T_iso*100,'%.1f'),'/',num2str(T_ISO*100,'%.1f'),' %'];   % ��ʽ��Ϊ�ַ���
%% ��ͼ
subplot(3,2,1);%ʱ���ݻ�
surf(z_It_L,-t_window_L*1e12,It_L);
shading interp;
axis([-inf,inf,-2,+2,-inf,inf,0,2]);
view(2);
xlabel('Position [m]');ylabel('Time Delay [ps]');

subplot(3,2,3);%Ƶ���ݻ�
surf(z_lambda_L,lambda_window_L*1e9,Ilambda_L);
shading interp;
axis([-inf,inf,2700,2860,-inf,inf,0,2]);
view(2);
xlabel('Position [m]');ylabel('Wavelength [nm]');

subplot(3,2,5);%������ͨ��
plot(x1_iso,y1_iso,'k',x2_iso,y2_iso,'b');                                  % �������루�ڣ���ʵ�ʣ�������Բ
axis([-2,2,-1,1]);
xlabel('Polarization of Beam && ISO');
legend(T_iso_str,'Location', 'northeast');                                  % ��ʾ͸����
legend('boxoff');
set(gca,'xtick',[],'ytick',[]);                                             % ��������̶�
box off;

subplot(3,2,2);%ʱ����̬
plot(-t_window*1e12,It_end);
axis([-t_window_max*1e12*0.5,+t_window_max*1e12*0.5,-0.05,inf]);
xlabel(['Time Delay [ps]' At2center_str]);ylabel('Pulse Power [W]');
legend(pulsewidth_str,'Location','northwest');
legend('boxoff');

subplot(3,2,4);%Ƶ����̬
plot(lambda_window*1e9,Ilambda_end);
axis([-inf,inf,-0.01,1.05]);
xlabel('Wavelength [nm]');ylabel('Spectrum Power [counts]');

subplot(3,2,6);%���������ݻ�
plot(z_L,E_L*1e9);axis([-inf,inf,-inf,inf]);
xlabel('Position [m]');ylabel('Pulse Energy [nJ]');
drawnow;
end
%============================== ʵʱ��ͼ ���ٰ� ==========================
function []=plotfast()
global At                                   % ��ǰ����ʱ������
global t_window                             % ʱ�򴰿�
global t_window_max                         % ʱ�򴰿����ֵ
global lambda_window                        % Ƶ�򴰿�
global At2center_str                        % ʱ�����״̬�ַ���

%��������
It_end=At2(At,'It');                                                        % ���㵱ǰʱ��ǿ��
Ilambda_end=At2(At,'Ilambda');                                              % ���㵱ǰƵ��ǿ��
Ilambda_end=Ilambda_end/max(Ilambda_end);                                   % Ƶ���һ��

%��������
pulsewidth=FWHM(It_end,'It');
if pulsewidth~=0
    pulsewidth_str=['FWHM ' num2str(pulsewidth*1e15,'%.0f'),' fs'];
else
    pulsewidth_str='�������';
end

% ���ٻ�ͼ������ͼ��
subplot(2,1,1);                                                             % ʱ��
plot(-t_window*1e12,It_end);
axis([-t_window_max*1e12*0.1,+t_window_max*1e12*0.1,-0.05,inf]);            % խʱ�䴰��
xlabel(['Time Delay [ps]' At2center_str]);                                  % ��ʾ����״̬
ylabel('Pulse Power [W]');        
legend(pulsewidth_str,'Location','northwest');
legend('boxoff');

subplot(2,1,2);                                                             % Ƶ��
plot(lambda_window*1e9,Ilambda_end);
axis([-inf,inf,-0.01,1.05]);
xlabel('Wavelength [nm]');
ylabel('Spectrum Power [counts]');
drawnow;
end

%=============================ƽ��������ԭ������==========================
function []=At2center()
%��������ʱ���źŵ�����˳��΢��������ֵ�������������м�λ��
global n_window                                                             % ʱ�򴰿ڵ���
global n_window_half                                                        % ʱ�򴰿ڰ볤������λ�ã�
global At                                                                   % ��ǰ����ʱ������
global At2center_str                                                        % ʱ�����״̬�ַ���

[~,p]=size(At);
At_center=zeros(n_window,p);                % ��ʼ�����к�����
It=At2(At,'It');                            % ����ʱ��ǿ��

%���ֵ��λ��
if 1==1 
    [~,center]=max(It);% �ҵ�ʱ��ǿ�����ֵλ��
    % ���ݷ�ֵλ��ѭ��ƽ������
    if center>n_window_half
        n_left=center-n_window_half;
        At_center(1:n_window-n_left,:)=At(n_left+1:n_window,:);             % ��������  
        At_center(n_window-n_left+1:n_window,:)=At(1:n_left,:);
    elseif center<n_window_half
        n_right=n_window_half-center;
        At_center(n_right+1:n_window,:)=At(1:n_window-n_right,:);
        At_center(1:n_right,:)=At(n_window-n_right+1:n_window,:);
    else 
        At_center=At;                                                       % �������
    end
end

%���ֵ��λ��
if 1==0 
    It_half=max(It)/2;                                                      % ����ʱ��ǿ�����ֵ��һ��
    n1_It_half=0;
    n2_It_half=0;
    for i=1:n_window-1
        if (It(i)-It_half)*(It(i+1)-It_half)<=0                             % ��⽻�㣨���ű仯��
            n1_It_half=i;break;
        end
    end
    for i=n_window:-1:n_window_half
        if (It(i)-It_half)*(It(i-1)-It_half)<=0                             % ��⽻�㣨���ű仯��
            n2_It_half=i;break;
        end
    end
    if n1_It_half==0 && n2_It_half==0
        At_center=At;                                                       % δ�ҵ����㣬������
    else
        % ��ֹ�е��Խ���ڱ�Ե���������ڴ�����β���紦��
        if abs(n1_It_half-n2_It_half)<n_window_half                         % ��ֹ��ͼʱ�����ڻ����Ե
            center=floor((n1_It_half+n2_It_half)/2);                        % ֱ��ȡ�е�
        else
            center=floor((n1_It_half+n2_It_half+n_window)/2);               % ѭ�������е�
        end
        if center>n_window_half
            n_left=center-n_window_half;
            At_center(1:n_window-n_left,:)=At(n_left+1:n_window,:);         % ����
            At_center(n_window-n_left+1:n_window,:)=At(1:n_left,:);
        elseif center<n_window_half
            n_right=n_window_half-center;
            At_center(n_right+1:n_window,:)=At(1:n_window-n_right,:);       % ����
            At_center(1:n_right,:)=At(n_window-n_right+1:n_window,:);
        else 
            At_center=At;                                                   % �����Ѿ���
        end
    end
end

At=At_center;                                                               % ��������
At2center_str='(ʵʱ������)';                                               % ����״̬�ַ���
end
%================================= ����FWHM ==============================
function [output]=FWHM(input,type)
global t_window                 % ʱ�򴰿�
global lambda_window            % Ƶ�򴰿�
global lambda0                  % ���Ĳ���

if strcmp(type,'It')
    window=t_window;end                                                     % ʱ��ʱ����
if strcmp(type,'Ilambda')
    window=lambda_window;end                                                % Ƶ�򲨳���

[m,n]=size(input);
output=zeros(1,n);

% ʱ��FWHM���㣨��ͳ���� �޷�����ߴ���
if strcmp(type,'It')
for j=1:n
    if max(input(:,j))~=1;input(:,j)=input(:,j)/max(input(:,j));end         % ��һ��
    n_half=0;x1=0;
    % Ѱ�Ұ�߿���
    for i=2:m
        if (input(i-1,j)-0.5)*(input(i,j)-0.5)<0
            width=window(i)-x1;
            x1=window(i);
            n_half=n_half+1;
        end
    end
    if n_half==2
        output(1,j)=abs(width);                 % ��ЧFWHM
    end
    if n_half>2
        output(1,j)=0;                          % �޷�����
    end
end
end

% Ƶ��FWHM���㣨����ʶ�� �����бߴ�������������
if strcmp(type,'Ilambda')
    for j=1:n
        Ilambda_max=0;
        % �����Ĳ������������ֵ
        for i=1:m
            if window(i)<lambda0+20e-9 && window(i)>lambda0-20e-9
                if input(i,j)>Ilambda_max
                    Ilambda_max=input(i,j);
                    n_Ilambda_max=i;
                end
            end
        end
        % �����Ұ�߿��
        for i=n_Ilambda_max:n_Ilambda_max+200
            if (input(i,j)-Ilambda_max/2)*(input(i+1,j)-Ilambda_max/2)<=0
                n_Ilambda_half2=i;
                break;
            end
        end
        % �����Ұ�߿��
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

%=============================function ����OC��===========================
function []=OC(T_input)
%T ͸����
%R �����ʣ�����ʣ�
global At                           % ��ǰ��������
global At_OC1                       % ���䵽OCǰ������
global At_OC2                       % ͸��������

At_OC1=At;                          % ��¼����ǰ����
T=0.8;                              % Ĭ��͸����80%
if nargin>0                         % �����Զ���͸����
    T=T_input;
end
At=sqrt((1-T))*At;                  % ͸�䲿�֣����䲿�ֱ�������
At_OC2=At;                          % ��¼͸�������
end
%============================function ����1/2��Ƭ=========================
function []=HWP(theta)
%������ģ�⾭��1/2��Ƭ��������Ϊˮƽ����ֱƫ���ʱ���ź� x��y
global At                           % ��ǰ�������ݣ�˫ƫ��
global At_HWP1                      % ���䵽HWPǰ������     
global At_HWP2                      % ͸��������

At_HWP1=At;                         % ��¼����ǰ����
theta=theta/180*pi;                 % �Ƕ�ת����

% �벨ƬJones��������
temp=[cos(-theta),-sin(-theta);sin(-theta)*exp(-1i*pi),cos(-theta)*exp(-1i*pi)]*At';% ��ʱ����ת��ʩ����λ�ӳ�
At=([cos(theta),-sin(theta);sin(theta),cos(theta)]*temp)';                          % ˳ʱ����ת�ָ�����
At_HWP2=At;                                                                         % ��¼�������
end
%============================function ����������==========================
function []=ISO(Ty_input)
%������ģ�⾭�������� ֻͨ��ˮƽ����
global At                           % ��ǰ�������ݣ�˫ƫ��
global At_ISO1                      % ���䵽������ǰ������
global At_ISO2                      % ͸��������
global T_ISO                        % ������ˮƽ͸����

At_ISO1=At;                         % ��¼����ǰ����
if nargin>0                         % �����Զ��崹ֱƫ��͸����
    Ty=Ty_input;
else 
    Ty=0.85/200;                    % Ĭ�ϴ�ֱ͸���ʼ��ͣ�ģ������������ԣ�
end

% ˮƽƫ��͸�䣬��ֱƫ������
At(:,2)=At(:,2)*sqrt(Ty);           % ��ֱ����˥��
Tx=0.85;                            % ˮƽ͸����85%
At(:,1)=At(:,1)*sqrt(Tx);
T_ISO=Tx;                           % ��¼ˮƽ͸����
theta=pi/4;                         %�ź�ƫ��ת45��
At=At*[cos(theta),-sin(theta);sin(theta),cos(theta)];
At_ISO2=At;
end
%============================function ����1/4��Ƭ=========================
function []=QWP(theta)
%������ģ�⾭��1/4��Ƭ��������Ϊˮƽ����ֱƫ���ʱ���ź� x��y
global At
global At_QWP1
global At_QWP2

At_QWP1=At;
theta=theta/180*pi;
temp=[cos(-theta),-sin(-theta);sin(-theta)*exp(-1i*pi/2),cos(-theta)*exp(-1i*pi/2)]*At';
At=([cos(theta),-sin(theta);sin(theta),cos(theta)]*temp)';
At_QWP2=At;
end
