clc;clear;
P=imread('lenna.bmp');P=rgb2gray(P);
iptsetpref('imshowborder','tight');
figure(1);subplot(1,3,1);imshow(P);title('Lenaͼ??');
[M,N]=size(P);P=double(P);
n=2*M*N;
h=0.002;t=800;
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
        s(i-t)=x1;
        if mod((i-t),3000)==0
            x0=x0+h*sin(y0);
        end
    end
end
S=mod(floor(s*pow2(16)),256);
S1=S(1:M*N);S2=S(M*N+1:2*M*N);B=zeros(M,N);C=zeros(M,N);
tic;
A=P(:);B0=0;B(1)=bitxor(bitxor(B0,S1(1)),A(1));
for i=2:M*N
    B(i)=bitxor(bitxor(B(i-1),S1(i)),A(i));
end
C0=0;C(M*N)=bitxor(bitxor(C0,S2(M*N)),B(M*N));
for i=M*N-1:-1:1
    C(i)=bitxor(bitxor(C(i+1),S2(i)),B(i));
end
C=reshape(C,M,N);toc;
subplot(1,3,2);imshow(uint8(C));title('??ɢ????ͼ??');
A=C(:);D=zeros(M,N);E=zeros(M,N);
D0=0;D(M*N)=bitxor(bitxor(D0,S2(M*N)),C(M*N));
for i=M*N-1:-1:1
    D(i)=bitxor(bitxor(C(i+1),S2(i)),C(i));
end
E0=0;E(1)=bitxor(bitxor(E0,S1(1)),D(1));
for i=2:M*N
    E(i)=bitxor(bitxor(D(i-1),S1(i)),D(i));
end
E=reshape(E,M,N);
subplot(1,3,3);imshow(uint8(E));title('??ԭ????ͼ??');