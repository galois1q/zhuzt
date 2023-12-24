---
title: '二阶pde数值求解算法'
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




#### 二阶偏微分方程分类

 
 
$$
A{U_{xx}} + 2B{U_{xy}} + C{U_{yy}} + D{U_x} + E{U_y} + FU = 0
$$
 
 

- 若 $B^2-AC<0$,椭圆型；
- 若 $B^2-AC=0$,抛物线型；
- 若 $B^2-AC>0$,双曲型；


#### 边界条件

- 第一类边界条件，给定边界 $\Gamma$上的函数值 $\phi|_\Gamma$;
- 第二类边界条件，给定边界 $\Gamma$上的函数导数值 $\frac{\partial \phi}{\partial n}|_\Gamma$;
- 第三类边界条件，给定边界 $\Gamma$上的函数值与导数值线性组合 $a\phi+b\frac{\partial \phi}{\partial n}|_\Gamma$;

### 椭圆形偏微分方程

主要考虑Poisson方程和Laplace方程:

 
 
$$
-[\frac{{{\partial ^2}}}{{\partial {x^2}}} + \frac{{{\partial ^2}}}{{\partial {y^2}}}]\phi  = S(x,y)
$$
 
 

区域离散化：
 
 
$$
-[{\phi_{i + 1,j}} + {\phi_{i - 1,j}} + {\phi_{i,j + 1}} + {\phi_{i,j - 1}} - 4{\phi_{ij}}] = {h^2}{S_{ij}}
$$
 
 
或表示为：
 
$$ 
-[({\delta_i^2}\phi)_ {ij} + {(\delta_{j}^{2}\phi)_ {ij}}] = {h^2}{S_{ij}}
$$ 
 
 
$$
\begin{array}{l}
{(\delta_i^2\phi )_ {ij}} \equiv {\phi_{i + 1,j}} + {\phi_{i - 1,j}} - 2{\phi_{ij}}\\
{(\delta_j^2\phi )_ {ij}} \equiv {\phi_{i,j + 1}} + {\phi_{i,j - 1}} - 2{\phi_{ij}}
\end{array} 
$$ 
 

正方形边界代码示例:

```matlab
% 边值条件与初值选取
xa=0;xb=pi;
ya=0;yb=pi;
N=100;
M=100;
h1=(xb-xa)/N;
h2=(yb-ya)/M;
xinit=xa:h1:xb;
yinit=ya:h2:yb;
u=zeros(N+1,M+1);
u(:,M+1)=sin(xinit);
uold=u;
u(2:N,2:M)=1/(2*pi);% 边界平均值为1/(2*pi)以初始化内点

```

#### Jacobi迭代法

 
$$
\varphi_{ij}^{n + 1}=\frac{1}{4}[\varphi_{i - 1,j}^{n} + \varphi_{i + 1,j}^n + \varphi_{i,j - 1}^{n} + \varphi_{i,j + 1}^n + {h^2}S_{ij}^n]
$$ 
 

```matlab
% 正方形边界问题下的jacobi迭代法实现
function usolve=jacobi_iteration(u,s,h,uold,N,M,err)
if nargin <5         
    err = 10^-6;
end  
iter_num=0;
while max(max(abs(u-uold)))>=err
      uold=u;
      u(2:N,2:M)=0.25*(u(3:N+1,2:M)+u(1:N-1,2:M)+u(2:N,3:M+1)+u(2:N,1:M-1)+h^2*s(2:N,2:M));
      iter_num=iter_num+1;
end
usolve=u;
% fprintf("The error of Jacobi method is %f\n",max(max(abs(u-uold))))
fprintf("The iteration number of Jacobi method is：%d\n",iter_num)
end
```
#### Gauss-Seidal法
 
$$
\varphi_{ij}^{n + 1} = (1 - \omega )\varphi_{ij}^n + \frac{\omega }{4}[\varphi_{i - 1,j}^{n + 1} + \varphi_{i + 1,j}^n + \varphi_{i,j - 1}^{n + 1} + \varphi_{i,j + 1}^n + {h^2}S_{ij}^n]
$$ 
 

其中 $0<\omega<2$，最佳值需要依据情况确定。


```matlab
% 正方形边界问题下的Gauss_Seidel迭代法实现
function usolve=Gauss_Seidel(u,S,h,w,N,M,err)
if nargin <5         
    err = 10^-8;
end    
iter_num=0;
merr=err+1;
while merr>err
    iter_num= iter_num+1;
    merr=0;
        for i=2:N
            for j=2:M
                %Gauss-Seidel迭代法
                relax=w/4*(u(i,j+1)+u(i,j-1)+u(i+1,j)+u(i-1,j)-4*u(i,j)+h^2*S(i,j));
                u(i,j)=u(i,j)+relax;
                if abs(relax)>merr
                    merr=abs(relax);
                end
            end
        end
end
 usolve=u;       
%  fprintf("The error of Gauss-Seidel method is%f\n",abs(merr))
 fprintf("The iteration number of Gauss-Seidel method is：%d\n",iter_num)
end  

```


### 抛物线形偏微分方程

主要以扩散方程和含时薛定谔方程为主：

一维扩散方程：
 
$$\frac{{\partial \Phi }}{{\partial t}} = \frac{{{\partial ^2}\Phi }}{{\partial {x^2}}} + S(x,t)
$$ 
 

一维扩散问题边界示例

```matlab
%格点化条件
L=1;
t0=1;
a=1;%u_t=a^2u_xx
dx=0.05;x=0:dx:L;N=length(x);
dt=0.001;t=0:dt:t0;M=length(t);
[X,T]=meshgrid(x,t);
y=zeros(N,M);
y(1:(N+1)/2,1) = 2*x(1:(N+1)/2); % 边界条件
y((N+1)/2:N,1) = 2*(1-x((N+1)/2:N));
```

#### 显式差分
 
$$\Phi_i^{n + 1} = r\Phi_{i + 1}^n + (1 - 2r)\Phi_i^n + r\Phi_{i - 1}^n + S_i^n\Delta t,\quad r=\frac{\Delta t}{h^{2}}
$$ 
 

只有当 $r<\frac{1}{2}$时算法才稳定。

若记

{{< math >}}  
$$
{(\delta^{2}\Phi^n)_{i}} = \Phi_{i + 1}^{n} - 2\Phi_{i}^{n} + \Phi_{i - 1}^{n}
$$ 
{{< /math >}}  

则显式差分也可写成：
 
$$
\Phi^{n + 1} = (1 - H\Delta t)\Phi^n + S^n\Delta t
$$ 
 

其中 ：

$$
H\Phi_i \equiv-\frac{1}{{{h^2}}}{(\delta^2\Phi)_i}
$$ 


```matlab
%----显式差分-------
function  uforward=forward(y,dx,dt,M,N,a)
rr=dt/dx^2*a;
for j=1:M-1    %时间上的步数
    for i=2:N-1   %空间上的点
        y(i,j+1)=rr*y(i+1,j)+(1-2*rr)*y(i,j)+rr*y(i-1,j);
    end
end
uforward=y';
end
```

#### 直接隐式差分

{{< math >}} 
$$
\tiny
\left[ {\begin{array}{l}
{1 + 2r}&{ - r}&{}&{}&{}\\
{ - r}&{1 + 2r}&{ - r}&{}&{}\\
{}& \ddots & \ddots & \ddots &{}\\
{}&{}&{ - r}&{1 + 2r}&{ - r}\\
{}&{}&{}&{ - r}&{1 + 2r}
\end{array}} \right]\left[ {\begin{array}{l}
{\Phi_1^{n + 1}}\\
{\Phi_2^{n + 1}}\\
{\cdots }\\
{\Phi_{N - 2}^{n + 1}}\\
{\Phi_{N - 1}^{n + 1}}
\end{array}} \right]
 = \left[ {\begin{array}{l}
{\Phi_1^n + S_1^n\Delta t + r\Phi_0^{n + 1}}\\
{\Phi_2^n + S_2^n\Delta t}\\
{...}\\
{\Phi_{N - 2}^n + S_{N - 2}^n\Delta t}\\
{\Phi_{N - 1}^n + S_{N - 1}^n\Delta t + r\Phi_N^{n + 1}}
\end{array}} \right]
$$ 
{{< /math >}} 

用算符表示为：
 
$$
\Phi^{n + 1} = \frac{1}{{1 + H\Delta t}}[\Phi^n + S^n\Delta t]
$$ 
 

```matlab
%----直接隐式差分----
function  uback=back(y,dx,dt,M,N,a)
rr=dt/dx^2*a;
A=diag((1+2*rr)*ones(N-2,1),0) + diag(-rr*ones(N-3,1),-1) + diag(-rr*ones(N-3,1),1);
for j=1:M-1    %时间上的步数
    for i=2:N-1   %空间上的点
        b = y(2:N-1,j);
        b(1) = b(1) + rr*y(1,j+1);
        b(N-2) = b(N-2) + rr*y(N,j+1);
        y(2:N-1,j+1) = soldiag(A,b);
    end
end
uback=y';
end
%--soldiag---
function x = soldiag(A,d)
%求解Ax=d，A为系数矩阵
    n = length(d);
    x = zeros(n,1);
    for i = 2:n
        w = A(i,i-1)/A(i-1,i-1);
        A(i,i) = A(i,i)-w*A(i-1,i);
        d(i) = d(i)-w*d(i-1);
    end
    x(n) = d(n)/A(n,n);
    for i = (n-1):(-1):1
        x(i) = (d(i)-A(i,i+1)*x(i+1))/A(i,i);
    end
end
```

#### 平均隐式差分

{{< math >}}  
$$
\tiny
\left[ {\begin{array}{l}
{1 + r}&{ - r/2}&{}&{}&{}\\
{ - r/2}&{1 + r}&{ - r/2}&{}&{}\\
{}& \ddots & \ddots & \ddots &{}\\
{}&{}&{ - r/2}&{1 + r}&{ - r/2}\\
{}&{}&{}&{ - r/2}&{1 + r}
\end{array}} \right]\left[ {\begin{array}{l}
{\Phi_1^{n + 1}}\\
{\Phi_2^{n + 1}}\\
{\cdots }\\
{\Phi_{N - 2}^{n + 1}}\\
{\Phi_{N - 1}^{n + 1}}
\end{array}} \right] = \left[ {\begin{array}{l}
{\frac{r}{2}\Phi_2^n + (1 - r)\Phi_1^n + \frac{r}{2}\Phi_0^n + S_0^n\Delta t + \frac{r}{2}\Phi_0^{n + 1}}\\
{\frac{r}{2}\Phi_3^n + (1 - r)\Phi_2^n + \frac{r}{2}\Phi_1^n + S_1^n\Delta t}\\
{...}\\
{\frac{r}{2}\Phi_{N - 1}^n + (1 - r)\Phi_{N - 2}^n + \frac{r}{2}\Phi_{N - 3}^n + S_{N - 2}^n\Delta t}\\
{\frac{r}{2}\Phi_N^n + (1 - r)\Phi_{N - 1}^n + \frac{r}{2}\Phi_{N - 2}^n + S_{N - 1}^n\Delta t + \frac{r}{2}\Phi_N^{n + 1}}
\end{array}} \right]
$$ 
{{< /math >}}  

用算符表示为：

 
$$
\Phi^{n + 1} = \frac{1}{{1 + \frac{1}{2}H\Delta t}}[(1 - \frac{1}{2}H\Delta t)\Phi^n + S^n\Delta t]
$$ 
 

```matlab
%---------平均隐式差分-----
function uck=Crank_Nicolson(y,dx,dt,M,N,a)
rr=dt/dx^2*a;
A = diag((1+rr)*ones(N-2,1),0) + diag(-rr/2*ones(N-3,1),-1) + diag(-rr/2*ones(N-3,1),1);
for j = 1:M-1
    b = (1-rr)*y(2:N-1,j) + rr/2*y(1:N-2,j) + rr/2*y(3:N,j);
    b(1) = b(1) + rr/2*y(1,j+1);
    b(N-2) = b(N-2) + rr/2*y(N,j+1);
    y(2:N-1,j+1) = soldiag(A,b);
end
uck = y';
end
```

#### 一维薛定谔方程的平均隐式差分

 
$$
i\hbar \frac{{\partial \phi }}{{\partial t}} =  - \frac{{{\hbar ^2}}}{{2m}}{\nabla ^2}\phi  + V\phi 
$$ 
 

平均隐式差分为：

 
$$
\phi^{n + 1} = \left( {\frac{{1 - i\frac{1}{2}H\Delta t}}{{1 + i\frac{1}{2}H\Delta t}}} \right)\phi^n,\quad H=-\frac{\partial ^2}{\partial x^2}+V
$$ 
 

实际计算时：

 
$$
\phi_{}^{n + 1} = \left( {\frac{2}{{1 + i\frac{1}{2}H\Delta t}} - 1} \right)\phi_{}^n \equiv \chi  - \phi_{}^n
$$ 
 

其中 $\chi$满足：

 
$$
{\chi_{j - 1}} + [ - 2 + \frac{{2i{h^2}}}{{\Delta t}} - {h^2}{V_j}]{\chi_j} + {\chi_{j + 1}} = \frac{{4i{h^2}}}{{\Delta t}}\phi_j^n,j = 1, \cdot  \cdot  \cdot ,N - 1
$$ 
 
或矩阵形式：

令 ：

 
$$
r=- 2 + \frac{{2i{h^2}}}{{\Delta t}} - {h^2}{V_j} 
$$ 
 

{{< math >}}  
$$
\left[ {\begin{array}{l}
r&1&{}&{}&{}\\
1&r&1&{}&{}\\
{}& \ddots & \ddots & \ddots &{}\\
{}&{}&1&r&1\\
{}&{}&{}&1&r
\end{array}} \right]\left[ {\begin{array}{l}
{{\chi_1}}\\
{{\chi_2}}\\
{\cdots }\\
{{\chi_{N - 2}}}\\
{{\chi_{N - 1}}}
\end{array}} \right] = \left[ {\begin{array}{l}
{\frac{{4i{h^2}}}{{\Delta t}}\phi_1^n{\rm{ - }}{\chi_0}}\\
{\frac{{4i{h^2}}}{{\Delta t}}\phi_2^n}\\
{...}\\
{\frac{{4i{h^2}}}{{\Delta t}}\phi_{N - 2}^n}\\
{\frac{{4i{h^2}}}{{\Delta t}}\phi_{N - 1}^n + {\chi_{N}}}
\end{array}} \right]
$$ 
{{< /math >}}  

一维高斯波包初值边界值示例

```matlab
%设定条件
xa=-30;xb=30;dx=0.1;
x=xa:dx:xb;
ta=0;tb=5;dt=0.01;t=ta:dt:tb;
N = length(x);
M = length(t);
V = transpose(arrayfun(@v,x));
sigma=3;k0=2.5;%Gauss波包参数
x0=-10;              %Gauss波包参数

psi0 = exp(k0*1i.*x).*exp(-(x-x0).^2*log10(2)/(2*sigma^2)); 
psi0(1) = 0; psi0(N) = 0; % 波函数边界条件
psi0 = transpose(psi0);
```

平均隐式差分

```matlab
%平均隐式差分
psi=zeros(N,M);
psi(:,1)=psi0;
chi = zeros(N,M-1);
r = 1i*dt/(dx^2);
A = diag((1+r)*ones(N-2,1)+1i/2*dt*V(2:N-1)) + diag(-r/2*ones(N-3,1),-1) + diag(-r/2*ones(N-3,1),1);
for j = 1:M-1
    b = 2*psi(2:N-1,j);
    b(1) = b(1) + r/2*chi(1,j);
    b(N-2) = b(N-2) + r/2*chi(N,j);
    chi(2:N-1,j) = inv(A)b;
    psi(2:N-1,j+1) = chi(2:N-1,j)-psi(2:N-1,j);
end
```

