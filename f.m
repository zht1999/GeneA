function result=f(T1,T2)
ws=0.25*pi;
wp=0.35*pi;
delta=(ws-wp); %过渡带 
N=63; %抽样点数 
% T1=0.1; %第1个非零点点幅值 
% T2=0.5; %第2个非零点幅值 
T3=0.70;
alpha=(N-1)/2; 
l=0:N-1; %阻频带上的采样点 
wl=(2*pi/N)*l;
flag=0;flag1=0;flag2=0;
Hrs=zeros(1,N);         %理想振幅响应采样 
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
k1=0:floor((N-1)/2); %向负无穷取整
k2=floor((N-1)/2)+1:N-1; 
angH=[-alpha*(2*pi)/N*k1,alpha*(2*pi)/N*(N-k2)]; %线性相位约束条件 --(5.235,5.236)
Hdk=Hrs.*exp(j*angH); %构成频域采样向量Hd(k) --（5-238）
h=real(ifft(Hdk,N)); %实际单位冲激响应，取实部
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
