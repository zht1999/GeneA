function [pop2] = binary2decimal(pop)
[px,py]=size(pop);
for i = 1:py
    pop1(:,i) = 2.^(py-i).*pop(:,i);
end
%sum(.,2)������ͣ��õ�������
pop2 = sum(pop1,2); 
end
%������ת����ʮ���ƺ���
%���������
%��������Ⱥ
%�������
%ʮ������ֵ