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


``` python
import numpy as np
from scipy.stats import poisson, zipf, beta, pareto, lognorm
```

### Lognormal

See [Occurance and Applications](https://en.wikipedia.org/wiki/Log-normal_distribution#Occurrence_and_applications) for examples of lognormal distributions in nature. 

> **Expert**: Most customers hold around \$20K (`mode`) but I could imagine a customer with $2.5M (`max`)

``` python

mode = 20000
max = 2500000

mean, stdv = e.elicitLogNormal(mode, max)
asset_values = lognorm(s=stdv, scale=np.exp(mean))
asset_values.rvs(100)

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
nMin = 1
nMax = 30
pMax = 1/1000

Zs = e.elicitZipf(nMin, nMax, pMax, report=True)

litigants = zipf(Zs, nMin-1)

litigants.rvs(100)
```

## Reference: Other Useful Elicitations

Listed as a courtesy, these distributions are simple enough to elicit data into directly without a helper function.

### Uniform

A "zero knowledge" distribution where all values within the range have equal probability of appearing. Similar to `random.randint(a, b)`

> Expert: The crowd will be between 50 (`min`) and 500 (`max`) due to fire code restrictions and the existing residents in the building.

``` python

from scipy.stats import uniform

min = 50
max = 500

range = max - min

crowd_size = uniform(min, range)
crowd_size.rvs(100)
```

### Poisson

> Expert: About 3000 Customers (`average`) add a credit card to their account every quarter.

``` python
from scipy.stats import poisson
average = 3000
upsells = poisson(average)
upsells.rvs(100)

```
