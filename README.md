# Marceau 
# Overview
This module provide a fast and efficient way to compute the [Panjer's Algorithm][panjer] in a Python shell.

## Usage

In the following paragraphs, I am going to describe how you can get and use Marceau for your own projects.

###  Getting it

To download Marceau, either from this [Github][git-repo-url] repository or simply use Pypi via pip.
```sh
pip install Marceau
```

### Module

Marceau uses two modules to work properly, you need to make sure to have the following on your computer:

- [Scipy.stats] - Used to generate probability mass function from discrete distributions.
- [Numpy] - For usefull calculations.

You are then ready to use it:
```sh
import Marceau
```

### Using it

The class Cossette built in the Marceau module calculate the Probability Density Function (PDF) and the Cumulative Distribution Function (CDF) of a Compound Distribution.

```sh
from Marceau import Cossette
```

The command
```sh
Cossette.help()
```
and
```sh
Cossette.example()
```

provide respectivly an brief help and two example of the following algorithm.

## Panjer's Algorithm

We are interested in the compound random variable: $$X=\sum_{i=1}^{N}B_{i}$$

where:
* $M$ is a frequence random variable from [Panjer-Katz] probability distribution family, otherwise known as (a,b, $0$)[class of distributions]. For $M=0$ we have $X=0$.
* $\underline{B}={B_{k},k\in\mathbb{N}^{+}}$ are positive i.i.d random variable defined on $\mathbb{N}$.
* $\underline{B}$ and $M$ are independant.

Therefore, the random variable $X$ has value in $\mathbb{N}$. And the Panjer's recursive method works as follow:
* If $B_{i}$ are distributed on a lattice $h\mathbb{N}$ with latticewidth $h>0$. $B\in$\{ $0$, $1h$, $2h$,....\}
* We have $X\in$ $A_{h}$=\{ $0$, $1h$, $2h$,....\}
* With $W_{M}$ beeing the probability generating function of M, we compute $f_{X}(0)=W_{M}(f_{B}(0))$
* The Panjer's recursive relation states for $k>0$: $$f_{X}(kh)=\frac{1}{1-af_{B}(0)}\sum_{i=1}^{k}(a+b\frac{jh}{kh})f_{B}(jh)\times f_{X}((k-j)h)$$


### Implementation

In order to compute the Panjer's Algorithm, we need to enter the following feature to our class Cossette.

| Arguments | Data Type| Description| 
| ------ | ------ | ------ |
| k | a positive integer| the epoch of recursion to find X distribution |
| h | a strictly positive integer|  the latticewidth of $B_{i}$ distribution |
| parameters | a list of length $1$ (poisson or geometric) or $2$ (binomial or negative binomial)| the parameters for the $X$ compound distribution |
| method | a string with value 'Binomial', 'NegBinomial', 'Geometric' or 'Poisson'| the law of $X$ compound distribution |
| fb | a list of length $k+1$ | this correspond to the $f_{B}$ values when those are given, default value is an empty list |
| generator\_param| a list of length $1$ (poisson or geometric) or $2$ (binomial or negative binomial) | the parameters of the $B$ distribution, only needed if $f_{B}$ is empty, default value is empty|
| generator\_method|   a string with value 'Binomial', 'NegBinomial', 'Geometric' or 'Poisson' | the law of $B$ distribution, only needed if $f_{B}$ is empty, default value is empty |


## Example
#### Example 1
Let $X\sim PComp(\lambda=2,F_{B}),$ with $B\sim Bin(10,0.4)$.

We implement the following

```sh
model= Marceau.Cossette(k=10,parameters=[2],method='Poisson',generator_method='Binomial',generator_param=[10,0.4]) 
```
And we get our output with the call model.panjer():
```sh
model.panjer()
 >>> f(10*1)=0.05434563071580669 
 F(10*1)=0.6980136730471336 
 ```
#### Example 2 

Let $X\sim PComp(\lambda=2,F_{B}),$ with $B \in$ \{ $1000$, $2000$, ... , $6000$ \} and the following values for $f_{B}(hk)$ with $h=1000$:

| $k$ | $0$ | $1$   | $2$   | $3$   | $4$    | $5$   | $6$   |
|---|---|-----|-----|-----|------|-----|------|
| $f_{B}(hk)$ |$0$ | $0.2$ | $0.3$ | $0.2$ | $0.15$ | $0.1$ | $0.05$ |

We implement the following:
```sh
fb=np.zeros(30*1000+1)
fb[0]=0
fb[1000]=0.2
fb[2000]=0.3
fb[3000]=0.2
fb[4000]=0.15
fb[5000]=0.1
fb[6000]=0.05 
model= Marceau.Cossette(k=10,h=1000,parameters=[1.25],method='Poisson',fb=fb) 
```

And we get our output with the call model.panjer():
```sh
model.panjer()
 >>> f(10*1000)=0.02089842353538644 
F(10*1000)=0.9536818666811318  
```


## Aknowledgement
This module was built with the help of [Marceau] lecture of Risk Theory.


## License

MIT
Copyright (c) 2022 Rayane Vigneron




[//]: # (These are reference links used in the body of this note and get stripped out when the markdown processor does its job)

   [git-repo-url]: <https://github.com/despervita/Marceau>
   [panjer]: <https://www.casact.org/sites/default/files/database/astin_vol12no1_22.pdf>
   [scipy.stats]: <https://docs.scipy.org/doc/scipy/reference/stats.html>
   [numpy]: <https://numpy.org/doc/stable/index.html>
   [Panjer-Katz]: <https://doi.org/10.1016/j.insmatheco.2010.03.010>
   [class of distributions]: <https://www.actuaries.org/ASTIN/Colloquia/Helsinki/Papers/S7_13_Fackler.pdf>
   [jQuery]: <http://jquery.com>
   [Marceau]: <https://www.act.ulaval.ca/departement-et-professeurs/professeurs-et-personnel/professeurs/fiche-de-professeur/etienne-marceau-138>
   
