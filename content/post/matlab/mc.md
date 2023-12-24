---
title: 'Monte Carlo'
subtitle: 
summary: '计算物理相关总结'
authors:
  - admin
tags: [计算物理]
categories: []
projects: []
date: 'today'
lastMod: 'today'
image:
  caption: ''
  focal_point: ''
---


### 随机数生成

**线性同余法**

若要生成 $[0,1]$上的随机数，可通过如下递推式生成伪随机数：

$$
  r_{i+1}= (ar_i+c)\mathrm{mod}M
$$

所得$\{r_i/ M\}$即想要的随机数。

```matlab
a=5;c=1;r1=3;m=2^8
N=1000;%随机数个数
r=zeros(1,N);
for i=1:N
    r(i)=mod(a*r1+c,m);
    r1=r(i);
end
figure; plot(1:N,r/m,'-*')
```


#### MC 积分


**直接抽样**：考虑 $x_i \sim \mathrm{rand}(0,1)$服从随机分布，积分等价于求 $f(x_i)$的期望
$$
 I= \int_{0}^1f(x)\mathrm{d}x =\lim_{N \to \infty}
 \frac{1}{N} \sum_{i=1}^{N}f(x_i)
$$

缺点：往往需要较大的样本量才能达到较小的误差。

**重要抽样**

考虑变量 $x_i\sim w(x)$,其中：
$$
\int_{0}^1w(x)\mathrm{d}x=1
$$

利用 ：

$$
y(x)=\int_{0}^{x}w(x')\mathrm{d} x'
$$

则原积分可以表示为：


$$
 I= \int_{0}^1 \frac{f(x(y))}{w(x(y))}\mathrm{d}y =\lim_{N \to \infty}
 \frac{1}{N} \sum_{i=1}^{N}\frac{f(x(y_i))}{w(x(y_i))}
$$

其中$y_i \sim \mathrm{rand}(0,1)$，如果 $\frac{f}{w}\sim \mathrm{const}$则抽样良好。

难点：构造合适的分布 $w(x)$与求反函数 $x(y)$。

#### 特定随机分布数构造


**Gauss 分布**

$$
w(x)=\frac{1}{\sqrt{2\pi}}e^{\frac{1}{2}x^2}
$$

由概率论：

$$
\frac{1}{12}\overline{\mathrm{rand}(0,1)}-6\sim N(0,1)
$$


```matlab
%Gauss分布的随机数；

N=100000;
Gauss=zeros(1,N);
for i=1:12
    Gauss=Gauss+rand(1,N);
end

Gauss=Gauss-6;
figure;plot(Gauss,1:N,'*')

%------判断分布是否为Gauss分布-----------------
Num=501;                                                %区间大小
Max=max(abs(Gauss))+10^-9;                       %边界
y=zeros(1,Num-1);
x=linspace(-Max,Max,Num);   %-Max,Max间分为501个点，500段
h=x(2)-x(1);
for i=1:length(Gauss)
    index=fix((Gauss(i)+Max)/h)+1;
    y(index)=y(index)+1;
end
figure;bar(x(1:end-1),y);                       %生成数绘图
```


**Von Neumanm舍选法**

考虑 $[a,b]$区间上的有界分布 $w(x)$,$w(x)\le c\equiv   w'(x)$,则可如下生成按 $w(x)$分布的随机数：

- 生成 $[a,b]$内的均匀分布随机数$x$
$$
x=(b-a)\mathrm{rand}(0,1)+a
$$
- 生成 $[0,c]$内的均匀分布随机数$y$
$$
y=c\mathrm{rand}(0,1)
$$
- 若 $y\le w(x)$，接受随机数 $x$w为所需 $w(x)$分布的随机数，否则重新抽取 $(x,y)$

```matlab
%利用Von Neumann舍选法产生随机数
a=0;b=1.0;
h=0.001;
w=@(x) 6/5*(1-0.5*x.^2);%随意选取的分布
N=5*10^6;
xx_s=Von_Neumann(w,a,b,h,N);
plot(xx_s(1:N),1:N,'*')

%---------统计随机数的分布------------
M=1000;
yy1=zeros(1,M);
h1=(b-a)/M;
xx1=a:h1:b;
for i=1:length(xx_s)
    index=fix(xx_s(i)/h1)+1;
    yy1(1,index)=yy1(1,index)+1;
end
figure;bar(xx1(1:end-1),yy1);                       %生成数绘图
S=interg(a,b,w);
yy1=yy1/N/h1*S; %因为yy1/N/h1在[a,b]上的积分为1，为了与w比较，乘以其在[a,b]上的积分S
figure;plot(xx1(1:M),yy1,'r*-');title('Von Neumann舍选法产生的随机数分布');
hold on
plot(xx1(1:M),w(xx1(1:M)),'b')
legend('Von Neumann舍选法产生的随机数分布','精确分布')

%生成步长为$h$的[a,b]区间N个按w(x)分布的随机数函数
function xx_s=Von_Neumann(w,a,b,h,N)
if nargin < 5        %nargin是用来判断输入变量个数的函数 
    h = 10^-4; %步长标准
end
x=a:h:b;
y_m=max(w(x))*1.001;%选择的w'
i=0;     
n=0;
xx_s=zeros(N,1);
while i<N
    xx=(b-a)*rand(1,1)+a;%产生[a,b]区间上的均匀分布的随机数
    yy=y_m*rand(1,1);
    n=n+1;   %统计产生了多少次随机数
    if yy<=w(xx) %舍选法的判断标准
        i=i+1;
        xx_s(i)=xx; 
    end
end
end

function S=interg(a,b,w) %w函数在[a,b]上的积分
    n=10000;
    h=(b-a)/n;
    xx_h=a:h:b;
    ww=w(xx_h);
    S=(4*sum(ww(2:2:n))+2*sum(ww(3:2:n-1))+ww(1)+ww(n+1))*h/3; %simpson积分公式
end
```

**Metropolis算法**

Metropolis算法与辐射原理一致，设变量序列 $X\sim w(X)$，应决定 $X_n$的下一点 $X_{n+1}$点生成。

- 将$w(X)$类比为能量，若$X_t$相对于$X_n$处于高能级, 即 $w(X_t)>w(X_{n})$，则 $X_t\to X_n$的自发辐射概率为1,$X_t$总是可接受的；
-  若$X_t$相对于$X_n$处于低能级, 即 $w(X_t)<w(X_{n})$，则$X_t\to X_n$的受激辐射概率为 $r=\frac{w(X_t)}{w(X_n)}$,仅有概率$r$使得$X_t$是可接受的；





