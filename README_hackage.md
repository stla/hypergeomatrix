# hypergeomatrix

## Evaluation of the hypergeometric function of a matrix argument (Koev & Edelman's algorithm)

Let $(a_1, \ldots, a_p)$ and $(b_1, \ldots, b_q)$ be two vectors of real or 
complex numbers, possibly empty, $\alpha > 0$ and $X$ a real symmetric or a 
complex Hermitian matrix. 
The corresponding *hypergeometric function of a matrix argument* is defined by 

$${}_pF_q^{(\alpha)} \left(\begin{matrix} a_1, \ldots, a_p \\\\ b_1, \ldots, b_q\end{matrix}; X\right) = \sum_{k=0}^{\infty}\sum_{\kappa \vdash k} \frac{{(a_1)}_{\kappa}^{(\alpha)} \cdots {(a_p)}_{\kappa}^{(\alpha)}} {{(b_1)}_{\kappa}^{(\alpha)} \cdots {(b_q)}_{\kappa}^{(\alpha)}} \frac{C_{\kappa}^{(\alpha)}(X)}{k!}.$$

The inner sum is over the integer partitions $\kappa$ of $k$ (which we also 
denote by $|\kappa| = k$). The symbol ${(\cdot)}_{\kappa}^{(\alpha)}$ is the 
*generalized Pochhammer symbol*, defined by

$${(c)}^{(\alpha)}_{\kappa} = \prod_{i=1}^{\ell}\prod_{j=1}^{\kappa_i} \left(c - \frac{i-1}{\alpha} + j-1\right)$$

when $\kappa = (\kappa_1, \ldots, \kappa_\ell)$. 
Finally, $C_{\kappa}^{(\alpha)}$ is a *Jack function*. 
Given an integer partition $\kappa$ and $\alpha > 0$, and a 
real symmetric or complex Hermitian matrix $X$ of order $n$, 
the Jack function 

$$C_{\kappa}^{(\alpha)}(X) = C_{\kappa}^{(\alpha)}(x_1, \ldots, x_n)$$

is a symmetric homogeneous polynomial of degree $|\kappa|$ in the 
eigen values $x_1$, $\ldots$, $x_n$ of $X$. 

The series defining the hypergeometric function does not always converge. 
See the references for a discussion about the convergence. 

The inner sum in the definition of the hypergeometric function is over 
all partitions $\kappa \vdash k$ but actually 
$C_{\kappa}^{(\alpha)}(X) = 0$ when $\ell(\kappa)$, the number of non-zero 
entries of $\kappa$, is strictly greater than $n$.

For $\alpha=1$, $C_{\kappa}^{(\alpha)}$ is a *Schur polynomial* and it is 
a *zonal polynomial* for $\alpha = 2$. 
In random matrix theory, the hypergeometric function appears for $\alpha=2$ 
and $\alpha$ is omitted from the notation, implicitely assumed to be $2$. 

Koev and Edelman (2006) provided an efficient algorithm for the evaluation 
of the truncated series 

$$\sideset{_p^m}{_q^{(\alpha)}}F \left(\begin{matrix} a_1, \ldots, a_p \\\\ b_1, \ldots, b_q\end{matrix}; X\right) = \sum_{k=0}^{m}\sum_{\kappa \vdash k} \frac{{(a_1)}_{\kappa}^{(\alpha)} \cdots {(a_p)}_{\kappa}^{(\alpha)}} {{(b_1)}_{\kappa}^{(\alpha)} \cdots {(b_q)}_{\kappa}^{(\alpha)}} 
\frac{C_{\kappa}^{(\alpha)}(X)}{k!}.$$

Hereafter, $m$ is called the *truncation weight of the summation* 
(because $|\kappa|$ is called the weight of $\kappa$), the vector 
$(a_1, \ldots, a_p)$ is called the vector of *upper parameters* while 
the vector $(b_1, \ldots, b_q)$ is called the vector of *lower parameters*. 
The user has to supply the vector $(x_1, \ldots, x_n)$ of the eigenvalues 
of $X$. 

For example, to compute

$$\sideset{_2^{15}}{_3^{(2)}}F \left(\begin{matrix} 3, 4 \\\\ 5, 6, 7\end{matrix}; 0.1, 0.4\right)$$

you have to enter 

```haskell
hypergeomat 15 2 [3.0, 4.0], [5.0, 6.0, 7.0] [0.1, 0.4]
```

We said that the hypergeometric function is defined for a real symmetric 
matrix or a complex Hermitian matrix $X$. Thus the eigenvalues of $X$ 
are real. However we do not impose this restriction in `hypergeomatrix`. 
The user can enter any list of real or complex numbers for the eigenvalues. 

### Gaussian rational numbers

The library allows to use **Gaussian rational numbers**, i.e. complex numbers 
with a rational real part and a rational imaginary part. The Gaussian rational 
number $a + ib$ is obtained with `a +: b`, e.g. `(2%3) +: (5%2)`. The imaginary 
unit usually denoted by $i$ is represented by `e(4)`:

```haskell
ghci> import Math.HypergeoMatrix
ghci> import Data.Ratio
ghci> alpha = 2%1
ghci> a = (2%7) +: (1%2)
ghci> b = (1%2) +: (0%1)
ghci> c = (2%1) +: (3%1)
ghci> x1 = (1%3) +: (1%4)
ghci> x2 = (1%5) +: (1%6)
ghci> hypergeomat 3 alpha [a, b] [c] [x1, x2]
26266543409/25159680000 + 155806638989/3698472960000*e(4)
```

### Univariate case

For $n = 1$, the hypergeometric function of a matrix argument is known as the 
[generalized hypergeometric function](https://mathworld.wolfram.com/HypergeometricFunction.html). 
It does not depend on $\alpha$. The case of $\sideset{_{2\thinspace}^{}}{_1^{}}F$ is the most known, 
this is the Gauss hypergeometric function. Let's check a value. It is known that

$$\sideset{_{2\thinspace}^{}}{_1^{}}F \left(\begin{matrix} 1/4, 1/2 \\\\ 3/4\end{matrix}; 80/81\right) = 1.8.$$

Since $80/81$ is close to $1$, the convergence is slow. We compute the truncated series below 
for $m = 300$.

```haskell
ghci> h <- hypergeomat 300 2 [1/4, 1/2] [3/4] [80/81]
ghci> h
1.7990026528192298
```


## References

- Plamen Koev and Alan Edelman. 
*The efficient evaluation of the hypergeometric function of a matrix argument*.
Mathematics of computation, vol. 75, n. 254, 833-846, 2006.

- Robb Muirhead. 
*Aspects of multivariate statistical theory*. 
Wiley series in probability and mathematical statistics. 
Probability and mathematical statistics. 
John Wiley & Sons, New York, 1982.

- A. K. Gupta and D. K. Nagar. 
*Matrix variate distributions*. 
Chapman and Hall, 1999.
