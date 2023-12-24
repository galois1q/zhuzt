---
title: '数值积分与数值求根'
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


### 插值

插值的目的是尽可能用多项式曲线拟合原函数的两点间曲线，然而值得注意的是并非越高阶的多项式拟合越好，过高的多项式次数会带来振荡，即Runge现象。

常见的插值：
- 线性插值：
$$\color{red}{\varphi_1(x)=y_0l_0(x)+y_1l_1(x)}$$
$$l_0(x)=\frac{x-x_1}{x_0-x_1},\quad l_1(x)=\frac{x-x_0}{x_1-x_0}$$
- 二次插值：
$$\color{red}{\varphi_2(x)=y_0l_0(x)+y_1l_1(x)+y_2l_2(x)}$$

$$l_0(x)=\frac{(x-x_1)(x-x_2)}{(x_0-x_1)(x_0-x_2)}\quad l_1(x)=\frac{(x-x_0)(x-x_2)}{(x_1-x_0)(x_1-x_2)}\quad l_2(x)=\frac{(x-x_0)(x-x_1)}{(x_2-x_0)(x_2-x_1)}$$

- n次插值：


$$\phi (x) = \sum\limits_{i = 0}^n {{y_i}{l_i}(x)} $$

{{< math >}} 
$$
\tiny
\begin{aligned}
{l_i}(x) &= \frac{{(x - {x_0}) \cdots (x - {x_{i - 1}})(x - {x_{i + 1}}) \cdots (x - {x_n})}}{{({x_i} - {x_0}) \cdots ({x_i} - {x_{i - 1}})({x_i} - {x_{i + 1}}) \cdots ({x_i} - {x_n})}}\\
& = \frac{{\omega (x)}}{{(x - {x_i})\omega '({x_i})}}    
\end{aligned}
$$
{{< /math >}} 

### 数值积分公式


**复化梯形公式**：对 $f(x)$线性插值，其中 $f(x_0)=a,f(x_n)=b,x_i=a+ih$，可得：

$$\begin{array}{c}
\int_a^b {f(x)dx}  = \frac{h}{2}[f(a) + 2\sum\limits_{i = 1}^{n - 1} {f(a + ih)}  + f(b)] + O({h^2})\\
 = \frac{h}{2}[2\sum\limits_{i = 0}^n {f(a + ih)}  - f(a) - f(b)] + O({h^2})
\end{array}$$

```matlab
function F=echelon(y,xa,xb,n)
    h=(xb-xa)/n;
    xi=xa:h:xb;
    y=y(xi);
    F=h/2*(-y(1)-y(n+1)+2*sum(y));
end    
```

**复化Simpson公式**: 对 $f(x)$在 $[a,b]$分 $2n$区间每个区间 $[x_i,x_{i+1}]$进行二次插值，可得到如下公式：

{{< math >}} 
$$
\tiny
\begin{aligned}
    \int_a^b {f(x)dx} & = \frac{h}{3}[f(a) + 4\sum\limits_{k = 0}^{n - 1} {f(a + (2k + 1)h)}  + 2\sum\limits_{k = 1}^{n - 1} {f(a + 2kh)}  + f(b)],\\
    &\quad {{h = (b - a)/2n,\quad }}{x_{k}= a + kh (k = 0,1,...2n)}
\end{aligned}
$$
{{< /math >}} 

```matlab
function F_s=my_simpson(f,x_min,x_max,n)
h=(x_max-x_min)/n;
x=x_min:h:x_max;
yy=f(x);
F_s=(4*sum(yy(2:2:n))+2*sum(yy(3:2:n-1))+yy(1)+yy(n+1))*h/3;
end        
```


**3/8 Simpson公式**:采取三阶插值的积分算法,将 $[a,b]$区间 $3n$等分，$h = (b - a)/3n,\quad {x_{k}= a + kh (k = 0,1,...3n)}$,则可得：

$$\int_a^b {f(x)dx}  = \frac{{3h}}{8}[f(a) + 3\sum\limits_{k = 0}^{n - 1} {[f(a + (3k + 1)h)}  + f(a + (3k + 2)h)] + 2\sum\limits_{k = 1}^{n - 1} {f(a + 3kh)}  + f(b)]$$

```matlab
function F_s=my_simpson(f,x_min,x_max,n)
h=(x_max-x_min)/n;
x=x_min:h:x_max;
y1=f(x);
F_s = (3*sum(y1(2:3:n-1))+3*sum(y1(3:3:n))+2*sum(y1(4:3:n-2))+y1(1)+y1(n+1))*h*3/8;
end        
```


对于无穷区间或边界带奇点的积分，应先将积分变形再进行积分。


### 数值求根

**二分法**: 最常见最简单的求根方式，利用满足 $f(a)f(b)<0$必在 $[a,b]$内有根进行二分求根。优点是可以找到所有根，但收敛速度慢。

```matlab
function [x,n]=bisection_root(f,a,b,delta)
 n=0; %搜索次数
 if f(a)*f(b)> 0
    fprintf('方程在[%d,%d]区间上无根',a,b)
    return
 else    
    while abs(b-a)>delta    %判断标准
        x=(a+b)/2;    %二分法
       	n=n+1;
        if f(x)*f(a)>0;%
            a=x;
        else
            b=x;
        end
    end
    x=(a+b)/2;
 end

```

**简单搜索法结合二分法**:二分法的改进

```matlab
%简单搜索法结合二分法
function root=Solve(f, x0, h, delta) 
%求根函数
if nargin < 4        %nargin是用来判断输入变量个数的函数 
    delta = 10^-7; %误差标准
end
 while h > delta          %简单搜索结合二分法
    if f(x0+h)*f(x0)>0
         x0=x0+h;
    else
         a=x0;
         b=x0+h;
         while abs(a-b)>delta
            x0=(a+b)/2;     
            if f(x0)*f(a)<0
                b=x0;
            else
                a=x0;
            end
         end
         x0=(a+b)/2;    
         h=b-a;
    end 
 end
root=x0;
end
```

**Newton法**:利用如下公式逼近根，在 $\lvert f(x_k) \rvert \le \delta, \lvert x_k-x_{k+1} \rvert \le\epsilon $停止收敛。
$${x_{k + 1}} = {x_k} - \frac{{f({x_k})}}{{f'({x_k})}},\quad k = 0,1,2, \cdots $$
Newton法收敛速度快，但依赖于给定函数的导数表达式且未必能找到所有根。

```matlab
%用Newton法求根
function [x1,k]=newton(y,dy,x0,delta,delta1)
k=0; %迭代次数
x1=x0-y(x0)/dy(x0); %Newton迭代公式
while abs(y(x0))>delta ||abs(x1-x0)>delta1%判断标准
    x0=x1;
    x1=x0-y(x0)/dy(x0); %Newton迭代公式
    k=k+1;
     if k>10^8 %如果迭代次数超过10^8
         fprintf('Newton法迭代次数已达10^8,很可能无法找到根\n'); %显示提示信息 
         break; %跳出循环 
     end
end   
end
```


**弦割法**: 利用数值微分估计实际导数的改进Newton法：
$${x_{k + 1}} = {x_k} - \frac{{f({x_k})({x_k} - {x_{k - 1}})}}{{f({x_k}) - f({x_{k - 1}})}}$$


但选取不好起始点仍然有可能找不到根。

```matlab
%弦割法求根
function [z,m]=secant(y,x0,delta)
a=-0.5+x0;b=0.5+x0;%两个启动点
m=1;
while abs(a-b)>delta | abs(y(b))>delta %判断标准
    z=b-y(b)*(b-a)/(y(b)-y(a)); %弦割法的迭代公式
    m=m+1;
    a=b;
    b=z;
end
end
```
