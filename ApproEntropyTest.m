clc;clear;

%%近似熵测试，检测长度为m和m+1的比特序列在待测序列在待测序列中出现的频率情况
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

ps=sum(s)/n;r=bitxor(s(1:end-1),s(2:end));vn_obs=1+sum(r);
p_v03=erfc(abs(vn_obs-2*n*ps*(1-ps))/(2*sqrt(2*n)*ps*(1-ps)));
clear ps r vn_obs;

m=10^4;K=6;N=floor(n/m);
run1=zeros(1,N);
for i=1:N
    s1=s((i-1)*m+1:i*m);
    j=1;k=1;t=1;
    while j<=m
        if sum(s1(t:j))==j-t+1
            k=j-t+1;j=j+1;
        else
            t=t+1;j=t+k;
        end
    end
    run1(i)=k;
end
v=zeros(1,7);
for i=1:7
    if i==1
        v(i)=sum(run1<11);
    elseif i==7
   v(i)=sum(run1>15); 
    else
        v(i)=sum(run1==9+i);
    end
end
ps=[0.0882 0.2092 0.2483 0.1933 0.1208 0.0675 0.0727];
chai2_obs=sum((v-N*ps).^2./(N*ps));clear m N v run1 s1 i k j t ps
p_v04=gammainc(chai2_obs/2,K/2,'upper');clear chai2-obs K

M=32;Q=32;N=floor(100000/(M*Q));k0=0;k1=0;k2=0;
for i=1:N
    R=zeros(M,Q);s1=s((i-1)*M*Q+1:i*M*Q);
    for j=1:M
        R(j,:)=s1((j-1)*Q+1:j*Q);
    end
    
    for a=1:M-1
        j=a+1;
        while j<=M && R(a,a)==0
            if R(j,a)==1
                t=R(j,:);R(j,:)=R(a,:);R(a,:)=t;
            end
            j=j+1;
        end
        j=a+1;
        while j<=Q && R(a,a)==0
            if R(a,j)==1
                t=R(:,j);R(:,j)=R(:,a);R(:,a)=t;
            end
            j=j+1;
        end
        if R(a,a)==1
            for u=a+1:M
                R(u,:)=bitxor(R(u,:),R(a,:));
            end
            for u=a+1:Q
                R(:,u)=bitxor(R(:,u),R(:,a));
            end
        end
    end
    r=sum(diag(R));
    
    if r==M
        k0=k0+1;
        elseif r==M-1
            k1=k1+1;
            else
                k2=k2+1;
    end 
end
chai2_obs=(k0-0.2888*N)^2/(0.2888*N)+(k1-0.5776*N)^2/(0.5776*N)+...
    (k2-0.1336*N)^2/(0.1336*N);
p_v05=exp(-chai2_obs/2);
clear i j k a chai2_obs k* M Q r R s1 t u N

n1=10^5;T=sqrt(log(1/0.05)*n1);N0=0.95*n1/2;
X=2*s(1:n1)-1;F=fft(X);F1=abs(F(1:floor(n1/2)));
N1=sum(F1<T);d=(N1-N0)/sqrt(0.95*0.05*n1/4);
p_v06=erfc(abs(d)/sqrt(2));
clear n1 T N0 N1 X F1 F d
        
N=20;M=floor(n/N);
m=10;B=[1 0 0 1 1 0 1 0 1 1];
miyou=(M-m-1)/pow2(m);sigma2=M*(1/pow2(m)-(2*m-1)/pow2(2*m));
W=zeros(1,N);
for i=1:N
    s1=s((i-1)*M+1:i*M);
    j=1;
    while j<=M-m+1
        if sum(bitxor(B,s1(j:j+m-1)))==0
            W(i)=W(i)+1;j=j+m;
        else
            j=j+1;
        end 
    end
end
chai2_obs=sum((W-miyou).^2/sigma2);
p_v07=gammainc(chai2_obs/2,N/2,'upper');
clear Bi j m M N chai2_obs miyou sigma2 s1 W

N=968;M=1032;m=9;B=[1 1 1 1 1 1 1 1 1];
lambta=(M-m+1)/pow2(m);yita=lambta/2;
W=zeros(1,N);v=zeros(1,6);
for i=1:N
    s1=s((i-1)*M+1:i*M);
    j=1;
    while j<M-m+1
        if sum(bitxor(B,s1(j:j+m-1)))==0
            W(i)=W(i)+1;j=j+1;
        else
            j=j+1;
        end
    end
    if W(i)<5
        v(W(i)+1)=v(W(i)+1)+1;
    else
        v(6)=v(6)+1;
    end
end
ps=zeros(1,6);ps(1)=exp(-yita);ps(2)=ps(1)*yita/2;
ps(3)=ps(1)*yita/8*(yita+2);ps(4)=ps(1)*yita/8*(yita^2/6+yita+1);
ps(5)=ps(1)*yita/16*(yita^3/24+yita^2/2+3*yita/2+1);
ps(6)=1-sum(ps(1:5));chai2_obs=sum((v-N*ps).^2./(N*ps));
p_v08=gammainc(chai2_obs/2,5/2,'upper');
clear B i j m M N chai2_obs lambta yita s1 W ps v

L=7;Q=1280;K=floor(n/L)-Q;T=zeros(1,pow2(L));
for j=1:Q
    t=s((j-1)*L+1:j*L);i=sum(t.*pow2(L-1:-1:0));T(i+1)=j;
end
sm=0;
for j=Q+1:Q+K
    t=s((j-1)*L+1:j*L);i=sum(t.*pow2(L-1:-1:0));
    sm=sm+log2(j-T(i+1));T(i+1)=j;
end
c=0.7-0.8/L+(4+32/L)*power(K,-3*L)/15;
sigma=c*sqrt(3.125/K);
fn=sm/K;p_v09=erfc(abs((fn-6.1962507)/(sqrt(2)*sigma)));
clear fn i j K L Q sm t T c
M=1000;K=6;N=floor(n/M);L=zeros(1,N);
for ii=1:N
    S=s((ii-1)*M+1:ii*M);
    fx=zeros(1,N+1);lx=zeros(1,N+1);fx(1)=1;lx(1)=0;cond=0;sm=0;
    for i=1:N
        d=mod(sum(fx(1:lx(i)+1).*S(i:-1:i-lx(i))),2);
        if d==0
            lx(i+1)=lx(i);
        else
            if sum(lx(1:i))==0
                lx(i+1)=i;fx(lx(i+1)+1)=1;
            else
                for j=i-1:-1:1
                    if lx(j)<lx(j+1)
                        if j+1==i
                            cond=1;break;
                        else
                            for k=i:-1:j+1+1
                                if lx(k)==lx(k-1)
                                    sm=sm+1;
                                end
                            end
                            if sm==i-j-1
                                sm=0;cond=1;break;
                            end
                        end
                    end
                end
                if cond==1
                    cond=0;
                    fx(i-j+1:i-j+1+lx(j))=mod(fx(i-j+1:i-j+1+lx(j))+fx(1:1+lx(j)),2);
                    lx(i+1)=max(lx(i),i-lx(i));
                end
            end
        end
    end
    L(ii)=max(lx);
end
miyou=M/2+(9-(-1)^(M+1))/36-(M/3+2/9)/pow2(M);
T=(-1)^M*(L-miyou)+2/9;v=zeros(1,K+1);
for i=1:N
    if T(i)<=-2.5
        v(1)=v(1)+1;
    elseif T(i)<=-1.5
        v(2)=v(2)+1;
    elseif T(i)<=-0.5
        v(3)=v(3)+1;
    elseif T(i)<=0.5
        v(4)=v(4)+1;
    elseif T(i)<=1.5
        v(5)=v(5)+1;
    elseif T(i)<=2.5
        v(6)=v(6)+1;
    else
        v(7)=v(7)+1;
    end
end
ps=[0.010417 0.03125 0.125 0.5 0.25 0.0625 0.020833];
chai2_obs=sum((v-N*ps).^2./(N*ps));
p_v10=gammainc(chai2_obs/2,K/2,'upper');
clear cond chai2_obs d fx i ii j k K L lx M miyou N S sigma sm T v ps
m=3;
vm0=zeros(1,pow2(m));vm1=zeros(1,pow2(m-1));vm2=zeros(1,pow2(m-2));
s1=[s s(1:m-1)];
for i=1:n
    for j=1:pow2(m)
        bm0=floor(mod((j-1),2.^(m:-1:1))./(2.^(m-1:-1:0)));
        if sum(bitxor(bm0,s1(i:i+m-1)))==0
            vm0(j)=vm0(j)+1;
        end
    end
    for j=1:pow2(m-1)
        bm1=floor(mod((j-1),2.^(m-1:-1:1))./(2.^(m-2:-1:0)));
        if sum(bitxor(bm1,s1(i:i+m-2)))==0
            vm1(j)=vm1(j)+1;
        end
    end
    for j=1:pow2(m-2)
        bm2=floor(mod((j-1),2.^(m-2:-1:1))./(2.^(m-3:-1:0)));
        if sum(bitxor(bm2,s1(i:i+m-3)))==0
            vm2(j)=vm2(j)+1;
        end
    end
end
psai2m0=sum(vm0.*vm0)*pow2(m)/n-n;
psai2m1=sum(vm1.*vm1)*pow2(m-1)/n-n;
psai2m2=sum(vm2.*vm2)*pow2(m-2)/n-n;
d_ps=psai2m0-psai2m1;d_ps2=psai2m0-2*psai2m1+psai2m2;
p_v11_1=gammainc(d_ps,pow2(m-2),'upper');
p_v11_2=gammainc(d_ps2,pow2(m-3),'upper');
clear b* d* i j m sl v* ps*

m=3;
vm0=zeros(1,pow2(m));vm1=zeros(1,pow2(m+1));
s1=s(1:10000);n1=length(s1);s1=[s1 s1(1:m)];
for i=1:n1
    for j=1:pow2(m)
        bm0=floor(mod((j-1),2.^(m:-1:1))./(2.^(m-1:-1:0)));
        if sum(bitxor(bm0,s1(i:i+m-1)))==0
            vm0(j)=vm0(j)+1;
        end
    end
    for j=1:pow2(m+1)
         bm1=floor(mod((j-1),2.^(m+1:-1:1))./(2.^(m:-1:0)));
        if sum(bitxor(bm1,s1(i:i+m)))==0
            vm1(j)=vm1(j)+1;
        end
    end
end
vm0=vm0/n1;vm1=vm1/n1;fai0=sum(vm0.*log(vm0));fai1=sum(vm1.*log(vm1));
chai2=2*n1*(log(2)+fai1-fai0);
p_v12=gammainc(chai2/2,pow2(m-1),'upper');
clear b* c* f* i j m n1 s1 v*