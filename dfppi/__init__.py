"""dfppi: prediction-powered fusion of a randomised trial (RCT) and observational data (OBS) for the treatment
effect of the trial population.

    fuse_ate(trial, obs, trial_propensity)    Algorithm 1 (AIPWF): ATE, standard error, Wald interval
    fuse_cate(trial, obs, trial_propensity)   Algorithm 2 (DRF or RF): a fitted CATE function

Read SPEC.md first (inputs, outputs, assumptions, when not to use).  Modules: aipwf (Algorithm 1), drf and rf
(Algorithm 2), method (nuisances, density ratio, scores), crossfit (rotations and the cross-fitted variance), api
(the two functions above), cli (python -m dfppi).
"""
__version__ = "0.1.0"

from .api import CATEFit, fuse_ate, fuse_cate  # noqa: E402

__all__ = ["fuse_ate", "fuse_cate", "CATEFit", "__version__"]
