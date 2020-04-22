clear
N=31;%点数(选奇数点)
wc=0.5;%截止频率
n=0:1:N-1;
%wH=zeros(1,N);% 编写矩形窗
wB=0.42-0.5*cos(2*pi/N*n)+0.08*cos(2*pi/N*2*n);
hd=wc*sinc(wc*(n-(N-1)/2));%读入hd（n）函数(位移)
h1=hd.*wB;%计算h（n）
subplot(4,1,1);
stem(h1);
xlabel('n');ylabel('h(n)');
title('N=31');

N1=1024; 
H1=fft(h1,N1);%调用子程序计算H（k）
n=0:N1-1;
w=2*pi/N1*n;%频率采样点间隔
magH1=abs(H1);
subplot(4,1,2);
plot(w,abs(H1));
xlabel('w');ylabel('abs(H1)');
title('幅频');
subplot(4,1,3);
plot(w,20*log10(abs(H1)));%画幅度曲线
xlabel('w');ylabel('20*log10(abs(H1))');
title('幅频');
subplot(4,1,4);
plot(w,angle(H1));
xlabel('w');ylabel('angle(H1)');
title('相频');

sgtitle('N=31 布拉克曼窗--理想低通');
set(gcf,'position',[0,0,1000,1000]);
print(gcf,'-dbitmap','blackman31.bmp');