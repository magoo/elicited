# Elicit Distributions

Assists with the construction of probability distributions built from expert elicited data for use in monte carlo simulations.

## Usage

Until this is packaged for pip, copy `elicit_distibutions.py` in your code. Then:

```
import elicited
```

`elicited` is just a helper tool when using numpy and scipy, so you'll need these too. 

```
import numpy as np
import scipy
```



### Lognormal

See [Occurance and Applications](https://en.wikipedia.org/wiki/Log-normal_distribution#Occurrence_and_applications) for examples of lognormal distributions in nature. 

> Expert: I have assets at risk that would generate a wide range of losses.
> 
> Elicitor: What is the most common value of these assets?
>
> Expert: About $ 20K (`val_mod`)
> 
> Elicitor: What's the largest asset value you can imagine?
> 
> Expert: I suppose it could go as high as $2.5M (`val_max`)

Lognormal requires mean and standard deviation.

```
logN_mean, logN_stdv = elicitLogNormal(val_mod, val_max)
```


### Pareto

The 80/20 rule. See [Occurance and Applications](https://en.wikipedia.org/wiki/Pareto_distribution#Occurrence_and_applications)

> Expert: The legal costs of an incident could be devastating.
> 
> Elicitor: How devastating are we talking?
> 
> Expert: Well, typically costs are zero (`val_min`), but a black swan could be $100M (`val_max`). 
> 
> Elicitor: So we can assume yoru minimum legal costs for an incident are zero, and your maximum costs are $100M?
> 
> Expert: Sure.

```
b = elicitPareto(val_min, val_max)
p = pareto(b, loc=val_min-1., scale=1.))
```

### PERT

See [PERT Distribution](https://en.wikipedia.org/wiki/PERT_distribution)

> Expert: We have accounts that could be lost and result in losses.
>
> Elicitor: What is the dollar value of these accounts?
>
> Expert: About $500-$6000 (`val_min` / `val_max`).
>
> Elicitor: What's the most common account? (`val_mod`)
>
> Expert: Probably around $4500.

```
PERT_a, PERT_b = elicitPERT(val_min, val_mod, val_max)
pert = beta(PERT_a, PERT_b, loc=val_min, scale=val_max-val_min)
```


### Poisson

See [Occurance and Application](https://en.wikipedia.org/wiki/Poisson_distribution#Occurrence_and_applications)

This is done in numpy without assistance from elicitor. As a courtesy for those looking to use it, here's an example. 

`Example Code`


### Zipf's

See [Applications](https://en.wikipedia.org/wiki/Zipf%27s_law#Applications)

> Expert: We are concerned about lawsuits relatd to a breach.
>
> Elicitor: Assuming a breach happens, how many litigants will there be?
>
> Expert: One or a few. We could also see an Equifax-like situation. (`nMin`)
>
> Elicitor: So most likely only a handful of litigants. What's a nightmare situation?
>
> Expert: I'd guess maybe 30 or more litigants? (`nMax`)
>
> Elicitor: How likely would it be to have more than 30 litigants?
>
> Expert: Very unlikely, most cases would only have a few, as I said.
>
> Elicitor: Let's give it a number. Is it one-in-a-thousand, or million cases?
>
> Expert: I'd say one in a million cases. (`pMax`)

```
Zs = elicitZipf(nMin, nMax, pMax, report=True)

pd = zipf(Zs, nMin-1)
```


