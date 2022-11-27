# 矩阵分析大作业

本项目为关于矩阵分解的LU、QR(Gram-Schmidt)、Orthogonal Reduction (Householder reduction 和Givens reduction)和URV程序实现。并包含一个综合程序，根据选择参数的不同，实现不同的矩阵分解；并实现Ax=b方程组的求解，以及计算A的行列式。


## 环境
* python 3.8
* numpy 1.23.4

## 文件描述:
```
main.py                       # 主程序
MatrixExample.py              # 存储了一些矩阵样本
LUFactorization.py            # LU 分解
QRFactorization.py            # QR(Gram-Schmidt) 分解
HousholderReduction.py        # Housholder 约简
GivensReduction.py            # Givens 约简
URVFactorization.py           # URV 分解
```
## 运行结果
配置好环境以后，直接运行py文件。

### 主程序
主程序运行后会进入循环，首先根据提示输入1-5选择分解方法，输入6退出程序；其次输入1-7选择矩阵，矩阵在MatrixExample.py文件中定义。然后程序会计算矩阵的行列式以及使用选定分解方法的结果。

在预设的矩阵中，1-2包含A和b矩阵，会返回A的分解结果和Ax=b的解；3-7只包含A矩阵，不包含b矩阵，故只返回A的分解结果

注意：本程序中的Ax=b仅支持含有唯一解的情况，即rank(A)=rank(A|b)
```
> python .\main.py 
```

```
----------------------------------------------
Please choose one method:
 1. LU Factorization
 2. QR Factorization
 3. Householder Reduction
 4. Givens Reduction
 5. URV Reduction
 6. Exit the program!

2
Please type in the index (1-7) of the matrix defined in the MatrixExample.py!

2
The determinant of the input matrix is:  25.0
The input matrix A is: 
[[ 2. -3.  1.]
 [ 1. -1. -2.]
 [ 3.  1. -1.]]
The input matrix b is: 
[ 7. -2.  0.]
The input matrix is: 
[[ 2. -3.  1.]
 [ 1. -1. -2.]
 [ 3.  1. -1.]]
The result of the QR factorization is shown as follows: 
-----The Q matrix is: 
[[ 0.535 -0.774  0.341]
 [ 0.267 -0.228 -0.936]
 [ 0.802  0.592  0.085]]
-----The R matrix is: 
[[ 3.742 -1.069 -0.802]
 [ 0.     3.14  -0.91 ]
 [-0.    -0.     2.128]]
The solution of Qy = b is:
[ 3.207 -4.96   4.256]
The solution of Rx = y is:
[[ 1.]
 [-1.]
 [ 2.]]
```

### LU分解
LU分解要求被分解矩阵为方阵，且为非奇异矩阵。测试样例中，包含两个3x3非奇异矩阵、一个3x4矩阵和一个3x3奇异矩阵

```
> python .\LUFactorization.py 
```
结果如下：

```
The input matrix is: 
[[ 1.  2.  5.]
 [ 4. 20.  6.]
 [ 7.  1.  9.]]
The result of the LU factorization is shown as follows: 
-----The L matrix is: 
[[1.    0.    0.   ]
 [0.571 1.    0.   ]
 [0.143 0.096 1.   ]]
-----The U matrix is: 
[[ 7.     1.     9.   ]
 [ 0.    19.429  0.857]
 [ 0.     0.     3.632]]
-----The P matrix is: 
[[0. 0. 1.]
 [0. 1. 0.]
 [1. 0. 0.]]
----------------------------------------------
The input matrix is: 
[[10.  2.  4.]
 [ 6.  5.  2.]
 [ 7. 15.  9.]]
The result of the LU factorization is shown as follows: 
-----The L matrix is: 
[[1.    0.    0.   ]
 [0.7   1.    0.   ]
 [0.6   0.279 1.   ]]
-----The U matrix is: 
[[10.     2.     4.   ]
 [ 0.    13.6    6.2  ]
 [ 0.    -0.    -2.132]]
-----The P matrix is: 
[[1. 0. 0.]
 [0. 0. 1.]
 [0. 1. 0.]]
----------------------------------------------
The input matrix is: 
[[10.  2.  4.  8.]
 [ 6.  5.  2.  1.]
 [ 7. 15.  9.  6.]]
ERROR: The input matrix is not a square matrix.
----------------------------------------------
The input matrix is: 
[[10.  2.  4.]
 [10.  2.  4.]
 [ 7. 15.  9.]]
ERROR: The matrix is not invertible.
```

### QR分解
QR分解要求被分解矩阵为列向量无关的mxn矩阵，测试样例中，包含两个3x3非奇异矩阵、一个3x4矩阵和一个3x3奇异矩阵

```
> python .\QRFactorization.py 
```
结果如下：
```
The input matrix is: 
[[ 1.  2.  5.]
 [ 4. 20.  6.]
 [ 7.  1.  9.]]
The result of the QR factorization is shown as follows: 
-----The Q matrix is: 
[[ 0.123  0.039  0.992]
 [ 0.492  0.865 -0.095]
 [ 0.862 -0.5   -0.087]]
-----The R matrix is: 
[[ 8.124 10.955 11.324]
 [ 0.    16.881  0.885]
 [-0.    -0.     3.602]]
----------------------------------------------
The input matrix is: 
[[10.  2.  4.]
 [ 6.  5.  2.]
 [ 7. 15.  9.]]
The result of the QR factorization is shown as follows: 
-----The Q matrix is: 
[[ 0.735 -0.572  0.363]
 [ 0.441 -0.002 -0.897]
 [ 0.515  0.82   0.251]]
-----The R matrix is: 
[[13.601 11.396  8.455]
 [ 0.    11.142  5.084]
 [ 0.     0.     1.914]]
----------------------------------------------
The input matrix is: 
[[10.  2.  4.  8.]
 [ 6.  5.  2.  1.]
 [ 7. 15.  9.  6.]]
ERROR: The column vertor is not linearly independent.
----------------------------------------------
The input matrix is: 
[[10.  2.  4.]
 [10.  2.  4.]
 [ 7. 15.  9.]]
ERROR: The column vertor is not linearly independent.
```

### Housholder 约简
Housholder 约简的简单实现，测试样例中，包含两个3*3非奇异矩阵、一个3x3奇异矩阵、一个3x4矩阵和一个4x3矩阵

```
> python .\HousholderReduction.py 
```
结果如下：
```
The input matrix is: 
[[1. 2. 3.]
 [4. 5. 6.]
 [7. 8. 9.]]
The result of the Housholder Reduction is shown as follows: 
-----The P matrix is: 
[[ 0.123  0.492  0.862]
 [ 0.905  0.302 -0.302]
 [-0.408  0.816 -0.408]]
-----The T matrix is: 
[[ 8.124  9.601 11.078]
 [-0.     0.905  1.809]
 [ 0.     0.     0.   ]]
----------------------------------------------
The input matrix is: 
[[10.  2.  4.]
 [ 6.  5.  2.]
 [ 7. 15.  9.]]
The result of the Housholder Reduction is shown as follows: 
-----The P matrix is: 
[[ 0.735  0.441  0.515]
 [-0.572 -0.002  0.82 ]
 [ 0.363 -0.897  0.251]]
-----The T matrix is: 
[[13.601 11.396  8.455]
 [-0.    11.142  5.084]
 [ 0.    -0.     1.914]]
----------------------------------------------
The input matrix is: 
[[10.  2.  4.]
 [10.  2.  4.]
 [ 7. 15.  9.]]
The result of the Housholder Reduction is shown as follows: 
-----The P matrix is: 
[[ 0.634  0.634  0.444]
 [-0.314 -0.314  0.896]
 [ 0.707 -0.707  0.   ]]
-----The T matrix is: 
[[15.78   9.189  9.062]
 [ 0.    12.189  5.557]
 [-0.     0.    -0.   ]]
----------------------------------------------
The input matrix is: 
[[10.  2.  4.  8.]
 [ 6.  5.  2.  1.]
 [ 7. 15.  9.  6.]]
The result of the Housholder Reduction is shown as follows: 
-----The P matrix is: 
[[ 0.735  0.441  0.515]
 [-0.572 -0.002  0.82 ]
 [ 0.363 -0.897  0.251]]
-----The T matrix is: 
[[13.601 11.396  8.455  9.411]
 [-0.    11.142  5.084  0.337]
 [ 0.    -0.     1.914  3.511]]
----------------------------------------------
The input matrix is: 
[[10.  6.  7.]
 [ 2.  5. 15.]
 [ 4.  2.  9.]
 [ 8.  1.  6.]
 [10.  2.  4.]]
The result of the Housholder Reduction is shown as follows: 
-----The P matrix is: 
[[ 0.593  0.119  0.237  0.475  0.593]
 [ 0.411  0.771  0.092 -0.36  -0.314]
 [-0.608  0.452  0.502  0.418 -0.017]
 [ 0.127  0.116 -0.531  0.677 -0.48 ]
 [ 0.306 -0.417  0.634  0.111 -0.564]]
-----The T matrix is: 
[[16.852  6.29  13.292]
 [ 0.     5.517 11.853]
 [-0.    -0.     9.477]
 [ 0.    -0.    -0.   ]
 [ 0.     0.     0.   ]]
```

### Givens 约简
Givens 约简的简单实现，测试样例中，包含两个3*3非奇异矩阵、一个3x3奇异矩阵、一个3x4矩阵和一个4x3矩阵

```
> python .\GivensReduction.py 
```
结果如下：
```
The input matrix is: 
[[1. 2. 3.]
 [4. 5. 6.]
 [7. 8. 9.]]
3
The result of the Givens Reduction is shown as follows: 
-----The P matrix is: 
[[ 0.123  0.492  0.862]
 [ 0.905  0.302 -0.302]
 [-0.408  0.816 -0.408]]
-----The T matrix is: 
[[ 8.124  9.601 11.078]
 [-0.     0.905  1.809]
 [ 0.    -0.     0.   ]]
----------------------------------------------
The input matrix is: 
[[10.  2.  4.]
 [ 6.  5.  2.]
 [ 7. 15.  9.]]
3
The result of the Givens Reduction is shown as follows: 
-----The P matrix is: 
[[ 0.735  0.441  0.515]
 [-0.572 -0.002  0.82 ]
 [ 0.363 -0.897  0.251]]
-----The T matrix is: 
[[13.601 11.396  8.455]
 [ 0.    11.142  5.084]
 [-0.     0.     1.914]]
----------------------------------------------
The input matrix is: 
[[10.  2.  4.]
 [10.  2.  4.]
 [ 7. 15.  9.]]
3
The result of the Givens Reduction is shown as follows: 
-----The P matrix is: 
[[ 0.634  0.634  0.444]
 [-0.314 -0.314  0.896]
 [ 0.707 -0.707 -0.   ]]
-----The T matrix is: 
[[15.78   9.189  9.062]
 [-0.    12.189  5.557]
 [-0.    -0.    -0.   ]]
----------------------------------------------
The input matrix is: 
[[10.  2.  4.  8.]
 [ 6.  5.  2.  1.]
 [ 7. 15.  9.  6.]]
4
The result of the Givens Reduction is shown as follows: 
-----The P matrix is: 
[[ 0.735  0.441  0.515]
 [-0.572 -0.002  0.82 ]
 [ 0.363 -0.897  0.251]]
-----The T matrix is: 
[[13.601 11.396  8.455  9.411]
 [ 0.    11.142  5.084  0.337]
 [-0.     0.     1.914  3.511]]
----------------------------------------------
The input matrix is: 
[[10.  6.  7.]
 [ 2.  5. 15.]
 [ 4.  2.  9.]
 [ 8.  1.  6.]]
3
The result of the Givens Reduction is shown as follows: 
-----The P matrix is: 
[[ 0.737  0.147  0.295  0.59 ]
 [ 0.261  0.8    0.026 -0.539]
 [-0.615  0.458  0.498  0.405]
 [-0.101  0.358 -0.815  0.444]]
-----The T matrix is: 
[[13.565  6.34  13.565]
 [-0.     5.08  10.827]
 [-0.    -0.     9.475]
```

### URV 分解
URV 分解的简单实现，测试样例中，包含两个3*3非奇异矩阵

```
> python .\URVFactorization.py 
```
结果如下：
```
The input matrix is: 
[[1. 2. 3.]
 [4. 5. 6.]
 [7. 8. 9.]]
The result of the URV Factorization is shown as follows: 
-----The U matrix is: 
[[ 0.123  0.905 -0.408]
 [ 0.492  0.302  0.816]
 [ 0.862 -0.302 -0.408]]
-----The R matrix is: 
[[16.76  -0.     0.   ]
 [ 1.714  1.074 -0.   ]
 [ 0.     0.     0.   ]]
-----The V matrix is: 
[[ 0.485  0.573  0.661]
 [-0.774 -0.072  0.63 ]
 [ 0.408 -0.816  0.408]]
----------------------------------------------
The input matrix is: 
[[10.  2.  4.]
 [10.  2.  4.]
 [ 7. 15.  9.]]
The result of the URV Factorization is shown as follows: 
-----The U matrix is: 
[[ 0.634 -0.314  0.707]
 [ 0.634 -0.314 -0.707]
 [ 0.444  0.896  0.   ]]
-----The R matrix is: 
[[20.385  0.     0.   ]
 [ 7.964 10.771  0.   ]
 [-0.     0.    -0.   ]]
-----The V matrix is: 
[[ 0.774  0.451  0.445]
 [-0.572  0.798  0.187]
 [-0.271 -0.399  0.876]]
```











