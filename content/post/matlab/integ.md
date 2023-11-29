---
title: '数值积分求解算法'
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


#### 插值

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

$${l_i}(x) = \frac{{(x - {x_0}) \cdots (x - {x_{i - 1}})(x - {x_{i + 1}}) \cdots (x - {x_n})}}{{({x_i} - {x_0}) \cdots ({x_i} - {x_{i - 1}})({x_i} - {x_{i + 1}}) \cdots ({x_i} - {x_n})}} = \frac{{\omega (x)}}{{(x - {x_i})\omega '({x_i})}}$$

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

$$\int_a^b {f(x)dx}  = \frac{h}{3}[f(a) + 4\sum\limits_{k = 0}^{n - 1} {f(a + (2k + 1)h)}  + 2\sum\limits_{k = 1}^{n - 1} {f(a + 2kh)}  + f(b)],\quad {{h = (b - a)/2n,\quad }}{x_{k}= a + kh (k = 0,1,...2n)}$$

**3/8 Simpson公式**:采取三阶插值的积分算法,将 $[a,b]$区间 $3n$等分，$h = (b - a)/3n,\quad {x_{k}= a + kh (k = 0,1,...3n)}$,则可得：

$$\int_a^b {f(x)dx}  = \frac{{3h}}{8}[f(a) + 3\sum\limits_{k = 0}^{n - 1} {[f(a + (3k + 1)h)}  + f(a + (3k + 2)h)] + 2\sum\limits_{k = 1}^{n - 1} {f(a + 3kh)}  + f(b)]$$

对于无穷区间或边界带奇点的积分，应先将积分变形再进行积分。
