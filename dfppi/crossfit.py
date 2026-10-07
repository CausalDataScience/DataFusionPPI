"""Cross-fitting of Algorithm 1: the three trial blocks rotate through the roles (nuisance, tuning, evaluation),
the OBS roles stay fixed, and the ATE is the mean of the three estimates.  crossfit_ate gives that mean and the
first-order standard error of the mean (the variance of Thm. 2 with the cross-rotation terms).
"""
import numpy as np

from .method import var

ROTATIONS = ((0, 1, 2), (1, 2, 0), (2, 0, 1))    # trial blocks in the roles (nuisance, tuning, evaluation)
EVALUATED_IN = {roles[2]: k for k, roles in enumerate(ROTATIONS)}    # block -> the rotation that evaluates it
NUISANCE_IN = {roles[0]: k for k, roles in enumerate(ROTATIONS)}     # block -> the rotation that fits on it


def crossfit_ate(parts, g_blocks, g_On, calibrated_ratio):
    """The ATE averaged over the rotations and the first-order standard error of the average.

    parts[k]          the aipwf output of rotation k, and r_k at the OBS nuisance rows
    g_blocks[b]       g at the rows of trial block b;  g_On: g at the OBS nuisance rows
    calibrated_ratio  the ratio is estimated and calibrated to balance g on the nuisance samples.  Then every
                      trial unit enters twice, through the score of the rotation that evaluates its block and
                      through the calibration of the rotation that fits on its block:
                          phi_i = psi_R(V_i) + omega g(X_i),
                      and the OBS nuisance rows enter through chi = mean_k omega_k r_k g.
    """
    phi = []
    for b in range(3):
        value = parts[EVALUATED_IN[b]][0]["psi_R"]
        if calibrated_ratio:
            value = value + parts[NUISANCE_IN[b]][0]["omega"] * g_blocks[b]
        phi.append(value - value.mean())                 # centred within the block
    phi = np.concatenate(phi)
    V = phi @ phi / (len(phi) - 3) / len(phi)
    psi_O = np.mean([out["psi_O"] for out, _ in parts], axis=0)
    V += var(psi_O) / len(psi_O)
    if calibrated_ratio:
        chi = np.mean([out["omega"] * r_On * g_On for out, r_On in parts], axis=0)
        V += var(chi) / len(chi)
    return float(np.mean([out["estimate"] for out, _ in parts])), float(np.sqrt(V))
