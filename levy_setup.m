%采用Levy法辨识传递函数
%   1.                  3 s^2 + 4 s + 2
%         G0(s)=    -----------------------
%                   4 s^3 + 3 s^2 + 2 s + 1
%   2. 获得实部虚部G(S)=Re(s)+Im(s)
%   3.U, V, S, T
%   4.θ=Φ^(-1)*Y
%   5.验证
%   分子分母阶次过大要在高次项补零
num = [3 4 2]; % 系统分子
den = [4 3 2 1]; % 系统分母
sys = tf(num, den); % 创建传递函数

% 给定角频率向量
omega = logspace(-1, 3, 1000); 
L = length(omega);


%%获得实部虚部
% [re im]=getReIm(sys,omega);
[H, f] = freqresp(sys, omega);
% 频率响应的实部
re = real(H);
% 频率响应的虚部
im = imag(H);


%%辨识模型

[~,n]=size(num);%分子阶次
[~,m]=size(den);%分母阶次

[fenzi,fenmu]=levy_fit(re,im,omega,n,m);
disp('辨识模型:');
ss=tf(fenzi,fenmu)


%%模型验证

t = 0:0.01:10;
len=length(t);
u = ones(len,1);
[y_real, ~] = lsim(sys, u, t);
[y_hat,~]= lsim(ss, u, t);

R2(y_real,y_hat);
