%利用传输矩阵计算单层介质的透射率
n1=2.33;
n2=1.50;%两层介质的折射率
h1=200e-9;
h2=100e-9;%介质2厚度为0
A=1;
B=0;
wl=100:1:900;%波长wl的范围
k1=2*pi*n1./(wl*1e-9);
e1=2*pi*n1*h1./(wl*1e-9);%k1h1
k2=2*pi*n2./(wl*1e-9);
e2=2*pi*n2*h2./(wl*1e-9);%k2h2
c1=0.5*1i*(k1/k2+k2/k1);
c2=0.5*1i*(k1/k2-k2/k1);%公式部分分解
num=length(wl);
for j=1:num
    M1=[exp(1i*e1(j))*(cos(e2(j))+c1*sin(e2(j))),exp(-1i*(e1(j)))*(-c2*sin(e2(j)));exp(1i*e1(j))*c2*sin(e2(j)),exp(-1i*e1(j))*(cos(e2(j))-c1*sin(e2(j)))];%特征矩阵
    N=10;%周期数
    M=M1^N;%多个周期光子晶体的特征矩阵
    T(j)=abs((M(1,1)*M(2,2)-M(2,1)*M(1,2))./((A*M(2,2)-B*M(1,2))))^2;%透射率
end
figure(1)
plot(wl,T,'k')
box on
