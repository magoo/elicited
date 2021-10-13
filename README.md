# Elicited

Helper tools to construct probability distributions built from expert elicited data for use in monte carlo simulations. 

Credit to Brett Hoover, packaging by @magoo

## Usage

```
pip install elicited
```

```
import elicited as e
```

`elicited` is just a helper tool when using numpy and scipy, so you'll need these in your code.


```
import numpy as np
import scipy
```

### Lognormal

See [Occurance and Applications](https://en.wikipedia.org/wiki/Log-normal_distribution#Occurrence_and_applications) for examples of lognormal distributions in nature. 

> **Expert**: Most customers hold around \$20K (`val_mod`) but I could imagine a customer with $2.5M (`val_max`)

``` python
logN_mean, logN_stdv = e.elicitLogNormal(val_mod, val_max)
```


### Pareto

The 80/20 rule. See [Occurance and Applications](https://en.wikipedia.org/wiki/Pareto_distribution#Occurrence_and_applications)

> Expert: The legal costs of an incident could be devastating. Typically costs are almost zero (`val_min`) but a black swan could be $100M (`val_max`). 

``` python
b = e.elicitPareto(val_min, val_max)
p = pareto(b, loc=val_min-1., scale=1.))
```

### PERT

See [PERT Distribution](https://en.wikipedia.org/wiki/PERT_distribution)

> Expert: Our customers have anywhere from \$500-\$6000 (`val_min` / `val_max`), but it's most typically around $4500 (`val_mod`)



``` python
PERT_a, PERT_b = e.elicitPERT(val_min, val_mod, val_max)
pert = beta(PERT_a, PERT_b, loc=val_min, scale=val_max-val_min)
```


### Zipf's

See [Applications](https://en.wikipedia.org/wiki/Zipf%27s_law#Applications)

> Expert: If we get sued, there will only be a few litigants (`nMin`). Very rarely it could be 30 or more litigants (`nMax`), maybe once every thousand cases (`pMax`) it would be more.


``` python
Zs = e.elicitZipf(nMin, nMax, pMax, report=True)

pd = zipf(Zs, nMin-1)
```


