"""Fit gene-wise standardization on discovery samples and apply it unchanged."""
import numpy as np


def discovery_standardize(expression, discovery_samples):
    samples = expression.columns.intersection(discovery_samples)
    if len(samples) < 2:
        raise ValueError("At least two discovery samples are required for gene-wise standardization")
    train = expression.loc[:, samples]
    center = train.mean(axis=1)
    scale = train.std(axis=1).replace(0, np.nan)
    return expression.sub(center, axis=0).div(scale, axis=0), center, scale
