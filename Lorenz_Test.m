clc;clear;
%%FIPS 140-1随机特性测试

%产生长度为20000b的位序列 
h=0.002;n=20000;t=800;
a=10;b=8/3;c=28;r=-1;
x0=1.1;y0=2.2;z0=3.3;w0=4.4;
s=zeros(1,n);
for i=1:n+t
   K11=a*(y0-x0)+w0;K12=a*(y0-(x0+K11*h/2))+w0;
    K13=a*(y0-(x0+K12*h/2))+w0;K14=a*(y0-(x0+h*K13))+w0;
    x1=x0+(K11+K12+K13+K14)*h/6;
    
    K21=c*x1-y0-x1*z0;K22=c*x1-(y0+K21*h/2)-x1*z0;
    K23=c*x1-(y0+K22*h/2)-x1*z0;K24=c*x1-(y0+h*K23)-x1*z0;
    y1=y0+(K21+K22+K23+K24)*h/6;
    
    K31=x1*y1-b*z0;K32=x1*y1-b*(z0+K31*h/2);
    K33=x1*y1-b*(z0+K32*h/2);K34=x1*y1-b*(z0+h*K33);
    z1=z0+(K31+K32+K33+K34)*h/6;
    
    K41=-y1*z1+r*w0;K42=-y1*z1+r*(w0+K41*h/2);
    K43=-y1*z1+r*(w0+K42*h/2);K44=-y1*z1+r*(w0+h*K43);
    w1=w0+(K41+K42+K43+K44)*h/6;
    
    x0=x1;y0=y1;z0=z1;w0=w1;
    if i>t
         s(i-t)=mod(floor((x1+100)*pow2(16)),2);%floor()取整
        if mod((i-t),3000)==0
            x0=x0+h*sin(x0);
        end
    end
end

%单比特测试
k1=sum(s);k0=length(s)-k1;%计算位序列S中比特1和比特0的个数

%扑克测试
f=zeros(1,16);
for i=1:5000
    v=s(4*i-3)*8+s(4*i-2)*4+s(4*i-1)*2+s(4*i);
    f(v+1)=f(v+1)+1;
end
p=sum(f.^2)*16/5000-5000;

%求得长度为1~100的1游程和0游程的个数
run0=zeros(1,100);run1=zeros(1,100);
t0=0;t1=0;
for i=1:20000
    if i==1
        if s(i)==1
            t1=1;
        else
            t0=1;
        end
    elseif i==20000
      if s(i)==1 && t1>0
       t1=t1+1;
       run1(t1)=run1(t1)+1;
      end
      if s(i)==0 && t0>0
      t0=t0+1;
      run0(t0)=run0(t0)+1;
      end
      if s(i)==1 && t1==0
      run1(1)=run1(1)+1;
      run0(t0)=run0(t0)+1;
      end
      if s(i)==0 && t0==0
      run0(1)=run0(1)+1;
      run1(t1)=run1(t1)+1;
      end
    else
      if s(i)==1 && t1>0
      t1=t1+1;
      end
      if s(i)==0 && t0>0
      t0=t0+1;
      end  
      if s(i)==1 && t1==0
      t1=1;
      run0(t0)=run0(t0)+1;t0=0;
      end
      if s(i)==0 && t0==0
      t0=1;
      run1(t1)=run1(t1)+1;t1=0;
      end
    end
end

run0_1=run0(1);run0_2=run0(2);run0_3=run0(3);run0_4=run0(4);
run0_5=run0(5);run0_6=sum(run0(6:100));run0_26=sum(run0(26:100));

run1_1=run1(1);run1_2=run1(2);run1_3=run1(3);run1_4=run1(4);
run1_5=run1(5);run1_6=sum(run1(6:100));run1_26=sum(run1(26:100));
