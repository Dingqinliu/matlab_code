clc;clear;
%%游程测试，检验比特序列中1游程和0游程的分布是否均匀

h=0.002;t=800;n=1000000;a=10;b=8/3;c=28;r=-1;x0=1.1;y0=2.2;z0=3.3;w0=4.4;
s=zeros(1,n);

for i=1:n+t
   K11=a*(y0-x0)+w0;K12=a*(y0-(x0+K11*h/2))+w0;
    K13=a*(y0-(x0+K12*h/2))+w0;K14=a*(y0-(x0+h*K13))+w0;
    x1=x0+(K11+2*K12+2*K13+K14)*h/6;
    
    K21=c*x1-y0-x1*z0;K22=c*x1-(y0+K21*h/2)-x1*z0;
    K23=c*x1-(y0+K22*h/2)-x1*z0;K24=c*x1-(y0+h*K23)-x1*z0;
    y1=y0+(K21+2*K22+2*K23+K24)*h/6;
    
    K31=x1*y1-b*z0;K32=x1*y1-b*(z0+K31*h/2);
    K33=x1*y1-b*(z0+K32*h/2);K34=x1*y1-b*(z0+h*K33);
    z1=z0+(K31+2*K32+2*K33+K34)*h/6;
    
    K41=-y1*z1+r*w0;K42=-y1*z1+r*(w0+K41*h/2);
    K43=-y1*z1+r*(w0+K42*h/2);K44=-y1*z1+r*(w0+h*K43);
    w1=w0+(K41+2*K42+2*K43+K44)*h/6;
    
    x0=x1;y0=y1;z0=z1;w0=w1;
    if i>t
          s(i-t)=mod(floor((x1+100)*pow2(16)),2);
          if mod((i-t),3000)==0
            x0=x0+h*sin(x0);
    end
    end
end
clear K*;clear a b c h r t;clear x0 x1 y0 y1 z0 z1 w0 w1;
 
s1=2*s-1;sm=sum(s1);clear s1;sobs=abs(sm)/sqrt(n);
p_v01=erfc(sobs/sqrt(2));clear sm sobs

m=100;N=floor(n/m);f=zeros(1,N);
for i=1:N
    f(i)=sum(s((i-1)*m+1:i*m))/m;
end
chai2=4*m*sum((f-1/2).^2);
p_v02=gammainc(chai2/2,N/2,'upper');clear m N f chai2;

ps=sum(s)/n;r=bitxor(s(1:end-1),s(2:end));vn_obs=1+sum(r);%计算序列中比特1的出现频率。
p_v03=erfc(abs(vn_obs-2*n*ps*(1-ps))/(2*sqrt(2*n)*ps*(1-ps)));
clear ps r vn_obs; %清除不再使用的临时变量  