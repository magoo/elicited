import elicited as e
import numpy as np
from scipy.stats import zipf, beta, pareto, lognorm
import time

min = 10000
mode = 20000
max = 2500000

# Lognormal
start = time.time()
mean, stdv = e.elicitLogNormal(mode, max, quantP=0.95)
asset_values = lognorm(s=stdv, scale=np.exp(mean))
#print(asset_values.rvs(10))
end = time.time()
print("Lognormal Elapsed: ", end - start)

# Pareto
start = time.time()
b = e.elicitPareto(min,max,quantP=0.95)
pareto_values = pareto(b, loc=min-1., scale=1.)
#print(pareto_values.rvs(10))
end = time.time()
print("Pareto Elapsed: ", end - start)

# PERT
start = time.time()
PERT_a, PERT_b = e.elicitPERT(min, mode, max)
pert_values = beta(PERT_a, PERT_b, loc=min, scale=max-min)
pert_values.rvs(100)
#print(pert_values.rvs(10))
end = time.time()
print("PERT Elapsed: ", end - start)

# Zipf's
start = time.time()
Zs = e.elicitZipf(min, max, quantP=0.95, report=True)
litigants_values = zipf(Zs, min-1)
#print(litigants_values.rvs(10))
end = time.time()
print("Zipfs Elapsed: ", end - start)