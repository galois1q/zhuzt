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

#### 线性同余法

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
