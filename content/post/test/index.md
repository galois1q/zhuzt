---
title: 网站测试
subtitle: 网站测试

# Summary for listings and search engines
summary: 对相关功能测试

# Link this post with a project
projects: []

# Date published
date: 'today'

# Date updated
lastmod: 'today'

# Is this an unpublished draft?
draft: false

# Show this page in the Featured widget?
featured: false

# Featured image
# Place an image named `featured.jpg/png` in this page's folder and customize its options here.
image:
  caption: 'Image credit: [**Unsplash**](https://unsplash.com/photos/CpkOjOcXdUY)'
  focal_point: ''
  placement: 2
  preview_only: false

authors:
  - alice
  - bob

tags:
  - 测试
  - 学习

categories:
  - Demo
  - 教程
---
## 列表测试

> 你说得对，但是《原神》是由米哈游自主研发的一款全新开放世界冒险游戏。
- 游戏发生在一个被称作「提瓦特」的幻想世界，在这里，被神选中的人将被授予「神之眼」，导引元素之力。
- 你将扮演一位名为「旅行者」的神秘角色。
- 在自由的旅行中邂逅性格各异、能力独特的同伴们，和他们一起击败强敌，找回失散的亲人
- 同时，逐步发掘「原神」的真相。

1. 你是谁
2. 你从哪来

## 公式测试

$a^{2}+b^{2}=c^{2}$

$$a^{2}+b^{2}=c^{2}$$

\begin{equation}
  a^{2}+b^{2}=c^{2}
\end{equation}

$\bra{x}\ket{\psi}$


## 表格测试

| Tables   |      Are      |  Cool |
|----------|:-------------:|------:|
| col 1 is |  left-aligned | $1600 |
| col 2 is |    centered   |   $12 |
| col 3 is | right-aligned |    $1 |


## pdf

```
{{ < embed-pdf url="./path/to/your/file.pdf"  > }}
```
#### This shortcode loads the first page and displays pagination

#### Hide pagination

If you wish to hide the pagination use the following attribute
```
{{ < embed-pdf url="./path/to/your/file.pdf" hidePaginator="true"  > }}
```

#### Hide loading spinner

If you wish to hide the loading spinner use the following attribute
```
{{ < embed-pdf url="./path/to/your/file.pdf" hideLoader="true"  > }}
```

{{< embed-pdf url="./temple.pdf" hideLoader="false" hidePaginator="false" >}}