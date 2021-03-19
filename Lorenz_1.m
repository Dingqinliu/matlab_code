h=0.002;n=20000;t=800; %h�ǲ�����n���ܵ�״̬������t�ǹ���̬�ĵ���
a=10;b=8/3;c=28;r=-1;
x0=1.1;y0=2.2;z0=3.3;w0=4.4;%����ϵͳ��ʼֵ��ȡֵ��Χ�ֱ�-40��40������-40��40������1��81������-250��250��
xn=zeros(1,n);yn=zeros(1,n);zn=zeros(1,n);wn=zeros(1,n);%����1*n��ȫ����
%%��������ϵͳ����ɢ���������õĸĽ����Ľ׾�������-������������ɢ��
for i=1:n+t  %������n+t ����i���Դ�1��n+t-t?? ����ĳ�1:n ���治��t �Ƚ�һ��
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
    
    x0=x1;y0=y1;z0=z1;w0=w1;%������ʵ����ѭ����
    if i>t
        xn(i-t)=x1;yn(i-t)=y1;zn(i-t)=z1;wn(i-t)=w1;%��ȫ�����滻������ 
    end
end

subplot(3,2,1);plot(xn,yn,'-k');%��xn-yn����ͼ k��ʾ��ɫ��-k��ʾ�ú�ɫ��ֱ������
xlabel('x_n','FontName','Times New Roman');%xlabel('f','fontsize',12,'fontname','lucida calligraphy')
ylabel('y_n','FontName','Times New Roman');%frontsize��������Ǹ������С�ģ�fontname�Ǹ�������ʽ��
subplot(3,2,2);plot(xn,zn,'-r');%��xn-zn����ͼ
xlabel('x_n','FontName','Times New Roman');
ylabel('z_n','FontName','Times New Roman');
subplot(3,2,3);plot(xn,wn,'-g');%��xn-wn����ͼ
xlabel('x_n','FontName','Times New Roman');
ylabel('w_n','FontName','Times New Roman');
subplot(3,2,4);plot(yn,zn,'-b');%��yn-zn����ͼ
xlabel('y_n','FontName','Times New Roman');
ylabel('z_n','FontName','Times New Roman');
subplot(3,2,5);plot(yn,wn,'-m');%��yn-wn����ͼ
xlabel('y_n','FontName','Times New Roman');
ylabel('w_n','FontName','Times New Roman');
subplot(3,2,6);plot(zn,wn,'-c');%��zn-wn����ͼ
xlabel('z_n','FontName','Times New Roman');
ylabel('w_n','FontName','Times New Roman');%��ʾ��ɫ����ĸ��8�֣�rgb,c��ʾ���̣�m���Ϻ죬y�ƣ�k�ڣ�w��
