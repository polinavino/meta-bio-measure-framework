"""Estimate signature direction on discovery only; freeze it for validation."""
from scipy.stats import spearmanr


def discovery_signs(frame, names):
    train = frame[(frame["cohort"] == "discovery") & frame["mortality"].notna()]
    return {name: (1 if spearmanr(train[name], train["mortality"]).statistic >= 0 else -1)
            for name in names}
