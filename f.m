function result=f(T1,T2)
ws=0.25*pi;
wp=0.35*pi;
delta=(ws-wp); %���ɴ� 
N=63; %�������� 
% T1=0.1; %��1���������ֵ 
% T2=0.5; %��2��������ֵ 
T3=0.70;
alpha=(N-1)/2; 
l=0:N-1; %��Ƶ���ϵĲ����� 
wl=(2*pi/N)*l;
flag=0;flag1=0;flag2=0;
Hrs=zeros(1,N);         %���������Ӧ���� 
for i=1:(N+1)/2
    if wl(i)>ws&&flag==0
        Hrs(i)=T1;
        Hrs(N+1-i)=T1;
        flag=1;
        continue;
    end
    if flag==1&&flag1==0
        Hrs(i)=T2;
        Hrs(N+1-i)=T2;
        flag1=1;
        continue;
    end
    if flag1==1
        Hrs(i)=1;
        Hrs(N+1-i)=1;
    end
end
k1=0:floor((N-1)/2); %������ȡ��
k2=floor((N-1)/2)+1:N-1; 
angH=[-alpha*(2*pi)/N*k1,alpha*(2*pi)/N*(N-k2)]; %������λԼ������ --(5.235,5.236)
Hdk=Hrs.*exp(j*angH); %����Ƶ���������Hd(k) --��5-238��
h=real(ifft(Hdk,N)); %ʵ�ʵ�λ�弤��Ӧ��ȡʵ��
[c,w]=freqz(h,1);
%plot(w/pi,20*log10(abs(c)));
a=findpeaks(20*log10(abs(c)));
minmax=a(1);
for i=2:length(a)
    if(a(i)>minmax&&a(i)<-50) 
        minmax=a(i);
    end
end
result=-minmax;
end
