"""Independent SciPy reference values for the HTML's teaching calculations."""

import json
from pathlib import Path

import scipy
from scipy.stats import beta, nbinom, poisson


result = {
    "counts": [],
    "alleles": [],
    "scipy_version": scipy.__version__,
    "reference": "NB n=r, p=r/(r+mu); beta(HB+1, n-HB+1)",
}
for model in ["poisson", "nb", "zinb"]:
    for mu, r, pi in [(0.1, 0.2, 0), (0.5, 1, 0.2), (1, 2, 0.2), (12, 0.2, 0.8), (12, 40, 0.8)]:
        dist = poisson(mu) if model == "poisson" else nbinom(r, r / (r + mu))
        probs = [float(dist.pmf(k)) for k in range(16)] + [float(dist.sf(15))]
        if model == "zinb":
            probs = [p * (1 - pi) for p in probs]
            probs[0] += pi
        result["counts"].append(dict(model=model, mu=mu, r=r, pi=pi, probabilities=probs))
for n, hb in [(0, 0), (4, 3), (100, 75), (200, 0), (200, 200), (200, 100)]:
    dist = beta(hb + 1, n - hb + 1)
    result["alleles"].append(dict(n=n, hb=hb, mean=float(dist.mean()), low=float(dist.ppf(.025)), high=float(dist.ppf(.975))))
Path(__file__).with_name("math_reference.json").write_text(json.dumps(result, indent=2) + "\n")
print("Saved 15 count distributions and 6 Beta posterior references.")
