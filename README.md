# Elicited

Helper tools to construct probability distributions built from expert elicited data for use in monte carlo simulations. 

Credit to Brett Hoover, packaging by @magoo. This is early code, happy to look at issues and feedback.

## Usage

```bash
pip install elicited
```

```python
import elicited as e

```
`elicited` is just a helper tool when using numpy and scipy, so you'll need them in your code to pass parameters to distributions you're looking for.


```python
import numpy as np
from scipy.stats import poisson, zipf, beta, pareto, lognorm
```

Build distributions for monte carlo trials.

```python
# Elicted values from an expert
mode = 20000
max = 2500000

# Parameters for lognormal distribution
mean, stdv = e.elicitLogNormal(mode, max)

# Freeze distribution w/ parameters
asset_values = lognorm(s=stdv, scale=np.exp(mean))

# Draw values from lognormal distribution
asset_values.rvs(100)
```





### Lognormal

`elicitLogNormal` takes a `mode`, and a value at quantile `max`. The default `max` quantile is a 95% percentile estimate, and can be changed with `quantP`.

See [Occurance and Applications](https://en.wikipedia.org/wiki/Log-normal_distribution#Occurrence_and_applications) for examples of lognormal distributions in nature, and [elicited technical docs](docs/lognormal.md).

> **Expert**: Most customers hold around \$20K (`mode`) but I could imagine a customer with $2.5M (`max`)

``` python

mode = 20000
max = 2500000

mean, stdv = e.elicitLogNormal(mode, max, quantP=0.95)
asset_values = lognorm(s=stdv, scale=np.exp(mean))
asset_values.rvs(100)

```


### PERT

See [PERT Distribution](https://en.wikipedia.org/wiki/PERT_distribution), and [elicited technical docs](docs/pert.md).

> Expert: Our customers have anywhere from \$500-\$6000 (`val_min` / `val_max`), but it's most typically around $4500 (`val_mod`)


``` python
PERT_a, PERT_b = e.elicitPERT(val_min, val_mod, val_max)
pert = beta(PERT_a, PERT_b, loc=val_min, scale=val_max-val_min)
```

### Pareto

`elicitPareto` takes a `min`, and a value at quantile `max`. The default `max` quantile is a 95% percentile estimate, and can be changed with `quantP`.

The 80/20 rule. See [Occurance and Applications](https://en.wikipedia.org/wiki/Pareto_distribution#Occurrence_and_applications), and [elicited technical docs](docs/pareto.md).

> Expert: The legal costs of an incident could be devastating. Typically costs are almost zero (`min`) but a black swan could be $100M (`max`). 

``` python
b = e.elicitPareto(min, max, quantP=0.95)
p = pareto(b, loc=val_min-1., scale=1.))
```

### Zipf's

`elicitZipfs` takes a `min`, and a value at quantile `max`. The default `max` quantile is a 95% percentile estimate, and can be changed with `quantP`. This function also has a `report` argument to report error. 

See [Applications](https://en.wikipedia.org/wiki/Zipf%27s_law#Applications), and [elicited technical docs](docs/zipf.md).

> Expert: If we get sued, there will only be a few litigants (`min`). Very rarely it could be 30 or more litigants (`max`), maybe once every thousand cases (`p`) it would be more.


``` python
min = 1
max = 30
p = 1/1000

Zs = e.elicitZipf(min, max, quantP=p)

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
