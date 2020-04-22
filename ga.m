clc;clear;close all;

%%
% f= @(a,b)(1/2*a .* sin(a) .* cos(1/2 * a) ).*(1/2*b .* sin(b) .* cos(1/2 * b) ); % �������ʽ
%f= @(a,b)(10*sin(a)+7*abs(a-5)+10).*(10*sin(b)+7*abs(b-5)+10); % �������ʽ
%%������ά����
 figure;
 [x0_1,x0_2]=meshgrid(0:.005:1);
 y0=zeros(length(x0_1),length(x0_2));
for i=1:length(x0_1)
    for j=1:length(x0_2)
        y0(i,j)=f(x0_2(j,1),x0_1(1,i));
    end
end
% y0=f(x0_1,x0_2);
 mesh(x0_1, x0_2, y0);
 hold on;

%% 

N = 50;                         % ��ʼ��Ⱥ����
chromosome = 2;                  % �ռ�ά��
ger = 200;                       % ����������    
chromlength =[16,16];
limit = [0, 1;0, 1];        %x1��x2������
pc = 0.5;                        % �������

pm = 0.06;                      % �������

all_best=0;             %���Ž�
all_x=0;


%%
% for i=1:chromosome
%     eval(['POP.chrom',num2str(i),'=round(rand(N,chromlength(i)));']);
%     for j=1:N
%         eval(['POP.x',num2str(i),'(j)=(binary2decimal(POP.chrom1(j,:)))/(2^chromlength(i)-1)*(limit(i,2)-limit(i,1)) ;'])
% %       POP.x'i'(j)=(binary2decimal(POP.chrom1(j,:)))/(2^chromlength(i)-1)*(limit(i,2)-limit(i,1)) ;
%     end
% end

POP.chrom1=round(rand(N,chromlength(1)));        %����50*16�ľ���,ÿһ��Ԫ��Ϊ0~1�������,������������
POP.chrom2=round(rand(N,chromlength(2)));

%%
% ��ʼ����
for n=1:ger
%%
    % ���滥��
for i = 1:2:N-1
    if(rand<pc)
        cpoint = round(rand*chromlength);
        POP.NEWchrom1(i,:) = [POP.chrom1(i,1:cpoint(1)),POP.chrom1(i+1,cpoint(1)+1:chromlength(1))];
        POP.NEWchrom1(i+1,:) = [POP.chrom1(i+1,1:cpoint(1)),POP.chrom1(i,cpoint(1)+1:chromlength(1))];
        POP.NEWchrom2(i,:) = [POP.chrom2(i,1:cpoint(2)),POP.chrom2(i+1,cpoint(2)+1:chromlength(2))];
        POP.NEWchrom2(i+1,:) = [POP.chrom2(i+1,1:cpoint(2)),POP.chrom2(i,cpoint(2)+1:chromlength(2))];        
    else
        POP.NEWchrom1(i,:) = POP.chrom1(i,:);
        POP.NEWchrom1(i+1,:) = POP.chrom1(i+1,:);
        POP.NEWchrom2(i,:) = POP.chrom2(i,:);
        POP.NEWchrom2(i+1,:) = POP.chrom2(i+1,:);
    end
end
    % ���ݽ��滥���Ľ�� ������Ⱥ�Ļ���
POP.chrom1=POP.NEWchrom1;
POP.chrom2=POP.NEWchrom2;

%%
    % ������첢������Ⱥ
for i=1:N
    for j1=1:chromlength(1)
        if(rand<pm)
            POP.chrom1(i,j1)=~POP.chrom1(i,j1);
        end
    end
    
    for j2=1:chromlength(2)
        if(rand<pm)
            POP.chrom2(i,j2)=~POP.chrom2(i,j2);
        end    
    end
    
end
        

%%
    % ��״�����ѡ��
    
    % ������(�����Ʊ���)ת��Ϊ�Ա�����ȡֵ��10���Ƶ�����
for i=1:N
    POP.x1(i)=(binary2decimal(POP.chrom1(i,:)))/(2^chromlength(1)-1)*(limit(1,2)-limit(1,1)) ;
end
for i=1:N
    POP.x2(i)=(binary2decimal(POP.chrom2(i,:)))/(2^chromlength(2)-1)*(limit(2,2)-limit(2,1)) ;
end
    % �����Ա�����ȡֵ�õ����������
for i=1:N
    POP.y(i)=f(POP.x1(i),POP.x2(i));
end

    %������Ľ����λ����ת��Ϊ0-1֮�����ֵ���ȣ��൱�����̶ĵĸ�������������
a=min(POP.y);
b=sum(POP.y)+N*(-a);

for i=1:N
    POP.adapt(i)=(POP.y(i)-a)/b;
end

    %��ֵ����ת��Ϊ0-1֮�������Ľڵ㣨�൱�ڰ����ת��Ϊ�����̶��ϸ�������ı߽��ߣ�
POP.NWEadapt(1)=POP.adapt(1);

for i=2:N
    POP.NWEadapt(i)=POP.adapt(i)+POP.NWEadapt(i-1);
end

    %�������̶ģ���ȡһ�������cs�����������������̶��е�λ������
    %����ĳ�����򣬾ʹ�����һ����Ⱥ�ڵ�i�������ӵ�и���������ʾ�Ļ��򣬴Ӷ��õ�����Ⱥ
for i=1:N
    cs=rand;
    [a,b]=POP_erfen(1,N,POP.NWEadapt,cs);
    POP.NEWchrom1(i,:)=POP.chrom1(b,:);
    POP.NEWchrom2(i,:)=POP.chrom2(b,:);
end

    %������Ⱥ
POP.chrom1=POP.NEWchrom1;
POP.chrom2=POP.NEWchrom2;   

%%
% �������λ��
[this_x ,this_best]=best(POP);
if all_best<this_best
    all_best=this_best;
    all_x=this_x;
end
        cla;
        mesh(x0_1, x0_2, y0);
        plot3(POP.x1,POP.x2,POP.y,'ro');
        title('״̬λ�ñ仯');
        pause(0.1);
end

figure();
mesh(x0_1, x0_2, y0);
hold on;
plot3(all_x(1),all_x(2),all_best, 'ro');title('����λ��ͼ');

disp(['���ֵ��',num2str(all_best)]);
disp(['����ȡֵ��',num2str(all_x)]);
