clear
N=31;%����(ѡ������)
wc=0.5;%��ֹƵ��
n=0:1:N-1;
%wH=zeros(1,N);% ��д���δ�
wB=0.42-0.5*cos(2*pi/N*n)+0.08*cos(2*pi/N*2*n);
hd=wc*sinc(wc*(n-(N-1)/2));%����hd��n������(λ��)
h1=hd.*wB;%����h��n��
subplot(4,1,1);
stem(h1);
xlabel('n');ylabel('h(n)');
title('N=31');

N1=1024; 
H1=fft(h1,N1);%�����ӳ������H��k��
n=0:N1-1;
w=2*pi/N1*n;%Ƶ�ʲ�������
magH1=abs(H1);
subplot(4,1,2);
plot(w,abs(H1));
xlabel('w');ylabel('abs(H1)');
title('��Ƶ');
subplot(4,1,3);
plot(w,20*log10(abs(H1)));%����������
xlabel('w');ylabel('20*log10(abs(H1))');
title('��Ƶ');
subplot(4,1,4);
plot(w,angle(H1));
xlabel('w');ylabel('angle(H1)');
title('��Ƶ');

sgtitle('N=31 ����������--�����ͨ');
set(gcf,'position',[0,0,1000,1000]);
print(gcf,'-dbitmap','blackman31.bmp');