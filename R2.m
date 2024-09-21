function Rsquare=R2(y_real,y_hat)
L=length(y_real);
rss=zeros(L,1);
tss=zeros(L,1);
for i=1:L
    rss(i) = y_real(i) - y_hat(i);
    tss(i)=y_real(i)-mean(y_real);
end
 % Residual Sum of Squares
   RSS= rss' * rss; 
 % Total sum of Squares
   TSS=tss'*tss;
   Rsquare=1-RSS/TSS;

disp(['决定系数R²：', num2str(Rsquare)]);
end