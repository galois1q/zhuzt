---
title: 'ODE求解算法'
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


### 常微分方程初值问题

**Euler法**: 基本只有教学意义的算法；

$$
\begin{array}{l}
{y_{n + 1}} = {y_n} + hf({x_n},{y_n}) + O({h^2}) \\
y({x_0}) = {y_0}
\end{array}
$$

**4阶RK法**: 实际最常用的算法，使用需要将微分方程降阶化简；

$$
\begin{array}{l}
{y_{n + 1}} = {y_n} + \frac{h}{6}({K_1} + 2{K_2} + 2{K_3} + {K_4}) + O({h^5})\\
{K_1} = f({x_n},{y_n})\\
{K_2} = f({x_n} + \frac{h}{2},{y_n} + \frac{h}{2}{K_1})\\
{K_3} = f({x_n} + \frac{h}{2},{y_n} + \frac{h}{2}{K_2})\\
{K_4} = f({x_n} + h,{y_n} + h{K_3})
\end{array}
$$

```matlab
function y = odeRK4(f,tspan,y0,n) %微分方程df/dt,f为列向量
dim = length(y0);  % 维数
h = (tspan(2) - tspan(1)) / (n-1); % 步长
y = zeros(dim, n);
y(:, 1) = y0(:); % 初值
t = linspace(tspan(1), tspan(2), n);

for i = 1:n-1
    K1 = f(t(i), y(:,i));
    K2 = f(t(i)+h/2, y(:,i)+K1*h/2);
    K3 = f(t(i)+h/2, y(:,i)+K2*h/2);
    K4 = f(t(i)+h, y(:,i)+K3*h/2);
    y(:,i+1) = y(:,i) + h/6*(K1+2*K2+2*K3+K4);
end
end

```

**Numerov算法**

对于常微分方程 

$$\frac{dy}{dx}+k^2(x)y=S(x)$$

Numerov算法给出：

$$
\begin{array}{l}
{(1+\frac{h^2}{12}k_{n+1}^2)y_{n+1}-2(1-\frac{5h^2}{12}k_n^2)y_n+(1+\frac{h^2}{12}k_{n-1}^2)y_{n-1}} \\
{=\frac{h^2}{12}(S_{{n+1}} +{10}S_{n}+S_{n-1})+O(h^6)}
\end{array}
$$

该方法不能自启动，但对于该形式的微分方程精确度和收敛性都较好。

```matlab
function y = numerov(kfun,Sfun,xspan,y0,K,n)
% numerov.m
% y'' + k(x)*y = S(x)
% 初值条件 y(0) = y0, y'(0) = K
% 时间范围 xspan=[a,b]
dim = length(y0);  % 维数
h = (xspan(2) - xspan(1)) / (n-1); % 步长
y = zeros(dim, n);
x = linspace(xspan(1), xspan(2), n);
y(:, 1) = y0(:); % 第一个初值
y(:, 2) = y0(:)+K*h; % 第二个初值(泰勒一阶近似)
k = arrayfun(kfun,x);
S = arrayfun(Sfun,x);

for i = 2:n-1
    y(:,i+1) = (2*y(:,i)*(1-5*h*h/12*k(i))-y(:,i-1)*(1+h*h/12*k(i-1))+h*h/12*(S(:,i+1)+10*S(:,i)+S(:,i-1)))/(1+h*h/12*k(i+1));
end
end
```



## 常微分方程边值问题



**迭加法**

只考虑如下二阶线性边值问题：

$$y^{\prime\prime}+p\left(x\right)y^{\prime}+q\left(x\right)y=f\left(x\right)\\ y\left(a\right)=\alpha，\quad y\left(b\right)=\beta$$

将边值问题分解为两个初值问题解的叠加，分别可四阶RK方法求解：
$$
\begin{array}{l} 
y_1^{\prime\prime}+p(x)y_1^{\prime}+q(x)y_1=f(x)\\ 
y_1(a)=\alpha\text{,}\quad y_1^{\prime}(a)=0;
\end{array}
$$

$$
\begin{array}{l} 
y_2^{\prime\prime}+p(x)y_2^{\prime}+q(x)y_2=f(x)\\
y_2(a)=0\text{,}\quad y_2^{\prime}(a)=1;
\end{array}
$$



最后原方程的解为：
$$y(x)=y_1(x)+\frac{\beta-y_1(b)}{y_2(b)}y_2(x)$$

**打靶法**
对于如下第一类边界边值问题：
$$\begin{array}{l}y^{\prime\prime}=f\left(x,y,y^{\prime}\right),\\ y\left(a\right)=\alpha,y\left(b\right)=\beta,\end{array}$$

将其转化为如下初值问题, 其中 $y'(a)=s_k$:

$$\begin{array}{l}y^{\prime\prime}=f\left(x,y,y^{\prime}\right),\\ y\left(a\right)=\alpha,\quad y^{\prime}\left(a\right)=s_k,\end{array}$$

需要通过迭代二分求根 $y(b,s_k)=0$确定$s_k$不断改进求解函数, 因而打靶法耗时较长。其他边界条件类似。


若用Numerov算法结合弦割法：

```matlab
%弦割法结合打靶法
xspan=[-1,1];
n=10001;
xx=linspace(-1,1,10001);
ya=1;
yb=1/3;
K0=-500;
K1=50;
kfun=@(x)-2./(x+2).^2;

delta=10^-12;

z0=numerov(kfun,@(x)0,xspan,ya,K0,n);
err0=z0(length(xx))-yb;
z1=numerov(kfun,@(x)0,xspan,ya,K1,n);
err1=z1(length(xx))-yb;
Dt=min(abs(err0),abs(err1));

while Dt > delta          %弦割法
    if err0*err1>0.0
       fprintf('初始斜率不在两个输入点之间，程序终止');
        break;
    end
    K2=K1-err1*(K1-K0)/(err1-err0);
    K0=K1;
    K1=K2;
    err0=err1;
    z1=numerov(kfun,@(x)0,xspan,ya,K1,n);
    err1=z1(length(xx))-yb;
    Dt=min(abs(err0),abs(err1));
    if abs(err1)>abs(err0)
        yfinal=z0;
    else
        yfinal=z1;
    end
end
```

若用4阶RK法结合简单搜索二分法：

```matlab

xx=linspace(-1,1,10001);
n=10001;
ya=1;
yb=1/3;
delta=10^-8;
K0=-12;
K1=15;

f=@(x,y)[y(2);y(1)*2./(x+2).^2];
y0=[ya;K0];
yrk0= odeRK4(f,[-1,1],y0,n);
yt0=yrk0(1,:);
err0=yt0(length(xx))-yb;

y1=[ya;K1];
yrk1= odeRK4(f,[-1,1],y1,n);
yt1=yrk1(1,:);
err1=yt1(length(xx))-yb;

Dt=min(abs(err0),abs(err1));

while Dt > delta   
    if err1*err0>0.0
       fprintf('初始斜率不在两个输入点之间，程序终止');
        break;
    end
    K2=(K0+K1)/2;
    y2=[ya;K2];
    yrk2= odeRK4(f,[-1,1],y2,n);
    yt2=yrk2(1,:);
    err2=yt2(length(xx))-yb;
    if err0*err2>0
        K0=K2;
        err0=err2; 
        yt0=yt2;
    elseif err1*err2>0
        K1=K2;
        err1=err2;
        yt1=yt2;
    end
    Dt=min(abs(err0),abs(err1));
end


```




## 常微分方程本征值问题

实际还是打靶法，从最小的猜测本征值开始依次求根。