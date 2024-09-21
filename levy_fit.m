function [fenzi, fenmu]=levy_fit(re,im,omega,n,m)
% re频率响应的实部
% im频率响应的虚部
% Omega角速度
% n分子阶次
% m分母阶次

L=length(omega);
N=2*max(m,n)+1;
% phi=zeros(m+n+1);
V=zeros(1,N);
S=zeros(1,N);
T=zeros(1,N);
U=zeros(1,N);

for k=1:N
    for i=1:L
        V(k)=omega(i)^k+V(k);
        S(k)=re(i)*omega(i)^k+S(k);
        T(k)=im(i)*omega(i)^k+T(k);
        U(k)=(re(i)^2+im(i)^2)*omega(i)^k+U(k);
    end
end

 % 左上分块矩阵
 LU=zeros(n+1);
for k=1:m+n+1
    fuhao=1;%重新将符号置+
    for l=1:m+n+1 
        if (l<=n+1)&&(k<=n+1)   % 左上分块矩阵
            flag=k+l-2;%选取的数据
            if(flag==0)
LU(1,1)=L+1;    %V0单独填写
fuhao=-1;
            end
            if(mod(flag,2)==0)&&(flag~=0)
                LU(k,l)=fuhao*V(flag);
                fuhao=fuhao*(-1);%每填一个数将符号变号
            end
        end
    end
end




 ST=zeros(N);
 for k=1:n+1
     for l=1:n+1
         if (mod(l+k,2)==0)
             ST(k,l)=T(k+l-1);
         else
             ST(k,l)=S(k+l-1);
         end
     end

 end




% 创建一个足够大的矩阵来指定列的符号变换
signMatrix = ones(50);
% 设置每列的符号变换
for col=0:50   % 一个周期4列
signMatrix(:, col*4+1) = (-1).^(0:49)'; % 第一列: + - + - ...
signMatrix(:, col*4+3) = (-1).^(1:(50));   % 第三列: - + - + ...
signMatrix(:, col*4+4) = -ones(50, 1);      % 第四列: 全负
end

% 应用符号变换
% 右上分块矩阵
RU=ST(1:n+1,1:m).* signMatrix(1:n+1,1:m);
% 左下分块矩阵
LD=ST(1:m,1:n+1).* signMatrix(1:m,2:n+2);

%右下分块矩阵
 RD=zeros(m);
 for k=1:m
    fuhao=1;%每行重新将符号置+
    for l=1:m  
            flag=k+l;%选取的数据
            if(mod(flag,2)==0)
                RD(k,l)=fuhao*U(flag);
                fuhao=fuhao*(-1);%每填一个数将符号变号
            end
        
    end
 end

 %%indentification
 phi = [LU RU;LD RD];%将矩阵拼起来

 %构造Y
 Y=zeros(n+m+1,1);
 S0=sum(re)+1;
 Y(1)=S0;
 Y(2:n+1,1)=ST(1:n,1);
 
 Y(n+3:n+m+1,1)=RD(1:m-1);


 %find theta
 
 theta=phi \ Y;%[a0 a1 a2...b1 b2...]

 %answer
 fenzi= fliplr(theta(1:n+1)');  %[an an-1...a1 a0]
 fenmu= [fliplr(theta(n+2:n+m+1)'),1];

end