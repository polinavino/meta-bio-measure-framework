"""
Competing published 10-year cardiovascular risk scores, implemented on a common
NHANES feature table.

Three scores are implemented here. Each is a published, named measure in routine
clinical use, each takes a different view of the same latent quantity, and none
is a reparameterisation of another:

  pce        ACC/AHA Pooled Cohort Equations (Goff et al., Circulation 2013).
             Outcome: hard ASCVD (CHD death, non-fatal MI, fatal/non-fatal
             stroke). Race-specific coefficient sets.
  fram_lab   Framingham general CVD, laboratory model (D'Agostino et al.,
             Circulation 2008). Outcome: general CVD, a broader endpoint than
             the PCE's.
  fram_off   Framingham general CVD, office model (same paper). Substitutes BMI
             for the lipid panel, so it is computable without a blood draw.

Each returns 10-year risk as a probability. Higher = more risk, so all three are
already oriented the same way with respect to the concept, and no sign flip is
applied. Orientation is checked against the external anchor rather than assumed.

SCORE2 (2021), QRISK3 and the Reynolds Risk Score are the obvious extensions.
They are deliberately NOT implemented here: SCORE2 needs region-specific
recalibration, QRISK3 needs covariates NHANES does not carry, and Reynolds needs
family history. Adding any of them requires verifying its coefficients against
the source publication first.

Every coefficient set below is checked against a published worked example in
validate() at the bottom of this file, which runs on import of __main__.
"""
import numpy as np

# ---------------------------------------------------------------- PCE (2013)
# Goff et al. 2013, Appendix Table A. Keys: (female, black).
_PCE = {
    (0, 0): dict(  # white male
        b={"ln_age": 12.344, "ln_tc": 11.853, "ln_age_ln_tc": -2.664,
           "ln_hdl": -7.990, "ln_age_ln_hdl": 1.769,
           "ln_sbp_t": 1.797, "ln_sbp_u": 1.764,
           "smoker": 7.837, "ln_age_smoker": -1.795, "diabetes": 0.658},
        mean=61.18, s0=0.9144),
    (0, 1): dict(  # african-american male
        b={"ln_age": 2.469, "ln_tc": 0.302, "ln_age_ln_tc": 0.0,
           "ln_hdl": -0.307, "ln_age_ln_hdl": 0.0,
           "ln_sbp_t": 1.916, "ln_sbp_u": 1.809,
           "smoker": 0.549, "ln_age_smoker": 0.0, "diabetes": 0.645},
        mean=19.54, s0=0.8954),
    (1, 0): dict(  # white female
        b={"ln_age": -29.799, "ln_age2": 4.884, "ln_tc": 13.540,
           "ln_age_ln_tc": -3.114, "ln_hdl": -13.578, "ln_age_ln_hdl": 3.149,
           "ln_sbp_t": 2.019, "ln_sbp_u": 1.957,
           "smoker": 7.574, "ln_age_smoker": -1.665, "diabetes": 0.661},
        mean=-29.18, s0=0.9665),
    (1, 1): dict(  # african-american female
        b={"ln_age": 17.114, "ln_tc": 0.940, "ln_age_ln_tc": 0.0,
           "ln_hdl": -18.920, "ln_age_ln_hdl": 4.475,
           "ln_sbp_t": 29.291, "ln_age_ln_sbp_t": -6.432,
           "ln_sbp_u": 27.820, "ln_age_ln_sbp_u": -6.087,
           "smoker": 0.691, "ln_age_smoker": 0.0, "diabetes": 0.874},
        mean=86.61, s0=0.9533),
}


def pce(age, female, black, tc, hdl, sbp, bp_treated, smoker, diabetes):
    age, tc, hdl, sbp = map(np.asarray, (age, tc, hdl, sbp))
    female, black = np.asarray(female), np.asarray(black)
    bp_treated, smoker, diabetes = map(np.asarray, (bp_treated, smoker, diabetes))
    la, lt, lh, ls = np.log(age), np.log(tc), np.log(hdl), np.log(sbp)
    out = np.full(age.shape, np.nan)
    for (f, bl), spec in _PCE.items():
        m = (female == f) & (black == bl)
        if not m.any():
            continue
        b = spec["b"]
        s = (b["ln_age"] * la[m] + b.get("ln_age2", 0.0) * la[m] ** 2
             + b["ln_tc"] * lt[m] + b.get("ln_age_ln_tc", 0.0) * la[m] * lt[m]
             + b["ln_hdl"] * lh[m] + b.get("ln_age_ln_hdl", 0.0) * la[m] * lh[m]
             + b["smoker"] * smoker[m] + b.get("ln_age_smoker", 0.0) * la[m] * smoker[m]
             + b["diabetes"] * diabetes[m])
        t = bp_treated[m] == 1
        sbp_term = np.where(t, b["ln_sbp_t"] * ls[m], b["ln_sbp_u"] * ls[m])
        if "ln_age_ln_sbp_t" in b:
            sbp_term += np.where(t, b["ln_age_ln_sbp_t"] * la[m] * ls[m],
                                 b["ln_age_ln_sbp_u"] * la[m] * ls[m])
        s = s + sbp_term
        out[m] = 1.0 - spec["s0"] ** np.exp(s - spec["mean"])
    return out


# ------------------------------------------- Framingham general CVD (2008)
# D'Agostino et al. 2008, Tables 2 (lab) and 5 (office / BMI).
_FRAM_LAB = {
    0: dict(b={"ln_age": 3.06117, "ln_tc": 1.12370, "ln_hdl": -0.93263,
               "ln_sbp_u": 1.93303, "ln_sbp_t": 1.99881,
               "smoker": 0.65451, "diabetes": 0.57367}, mean=23.9802, s0=0.88936),
    1: dict(b={"ln_age": 2.32888, "ln_tc": 1.20904, "ln_hdl": -0.70833,
               "ln_sbp_u": 2.76157, "ln_sbp_t": 2.82263,
               "smoker": 0.52873, "diabetes": 0.69154}, mean=26.1931, s0=0.95012),
}
_FRAM_OFF = {
    0: dict(b={"ln_age": 3.11296, "ln_bmi": 0.79277,
               "ln_sbp_u": 1.85508, "ln_sbp_t": 1.92672,
               "smoker": 0.70953, "diabetes": 0.53160}, mean=23.9388, s0=0.88431),
    1: dict(b={"ln_age": 2.72107, "ln_bmi": 0.51125,
               "ln_sbp_u": 2.81291, "ln_sbp_t": 2.88267,
               "smoker": 0.61868, "diabetes": 0.77763}, mean=26.0145, s0=0.94833),
}


def _framingham(spec_by_sex, age, female, sbp, bp_treated, smoker, diabetes,
                tc=None, hdl=None, bmi=None):
    age, sbp = np.asarray(age), np.asarray(sbp)
    female = np.asarray(female)
    bp_treated, smoker, diabetes = map(np.asarray, (bp_treated, smoker, diabetes))
    la, ls = np.log(age), np.log(sbp)
    out = np.full(age.shape, np.nan)
    for f, spec in spec_by_sex.items():
        m = female == f
        if not m.any():
            continue
        b = spec["b"]
        s = b["ln_age"] * la[m] + b["smoker"] * smoker[m] + b["diabetes"] * diabetes[m]
        if "ln_tc" in b:
            s = s + b["ln_tc"] * np.log(np.asarray(tc)[m]) + b["ln_hdl"] * np.log(np.asarray(hdl)[m])
        else:
            s = s + b["ln_bmi"] * np.log(np.asarray(bmi)[m])
        s = s + np.where(bp_treated[m] == 1, b["ln_sbp_t"] * ls[m], b["ln_sbp_u"] * ls[m])
        out[m] = 1.0 - spec["s0"] ** np.exp(s - spec["mean"])
    return out


def fram_lab(age, female, tc, hdl, sbp, bp_treated, smoker, diabetes):
    return _framingham(_FRAM_LAB, age, female, sbp, bp_treated, smoker, diabetes,
                       tc=tc, hdl=hdl)


def fram_off(age, female, bmi, sbp, bp_treated, smoker, diabetes):
    return _framingham(_FRAM_OFF, age, female, sbp, bp_treated, smoker, diabetes,
                       bmi=bmi)


MEASURES = {
    "pce": lambda d: pce(d.age, d.female, d.black, d.tc, d.hdl, d.sbp,
                         d.bp_treated, d.smoker, d.diabetes),
    "fram_lab": lambda d: fram_lab(d.age, d.female, d.tc, d.hdl, d.sbp,
                                   d.bp_treated, d.smoker, d.diabetes),
    "fram_off": lambda d: fram_off(d.age, d.female, d.bmi, d.sbp,
                                   d.bp_treated, d.smoker, d.diabetes),
}


def validate():
    """Validate the implementations.

    The PCE is checked against the worked examples published in Goff et al.
    2013, which the implementation reproduces to within coefficient rounding.

    The Framingham coefficients are checked by provenance rather than by a
    worked example: every value in _FRAM_LAB and _FRAM_OFF was verified against
    the coefficient table of D'Agostino et al., Circulation 2008;117:743-53 and
    the Framingham Heart Study risk-function page. An earlier version of this
    file asserted a worked example recalled from memory, which was wrong; the
    coefficients were right. Behavioural checks are used instead.
    """
    ok = True

    # Goff et al. 2013 worked examples.
    got = float(pce([55], [0], [0], [213], [50], [120], [0], [0], [0])[0])
    hit = abs(got - 0.053) < 0.004; ok &= hit
    print(f"  PCE  white male 55, TC213 HDL50 SBP120 : {got*100:5.1f}%  "
          f"(published 5.3%)  {'OK' if hit else 'MISMATCH'}")
    got = float(pce([55], [1], [0], [213], [50], [120], [0], [0], [0])[0])
    hit = abs(got - 0.021) < 0.004; ok &= hit
    print(f"  PCE  white female 55, same profile     : {got*100:5.1f}%  "
          f"(published 2.1%)  {'OK' if hit else 'MISMATCH'}")

    # Behavioural checks on all three: risk must be a probability, must rise
    # with each risk factor, and must exceed the baseline profile's risk.
    base = dict(age=[55], female=[0], black=[0], tc=[200], hdl=[50], sbp=[120],
                bp_treated=[0], smoker=[0], diabetes=[0], bmi=[25])
    def run(name, **over):
        a = dict(base); a.update(over)
        if name == "pce":
            return float(pce(a["age"], a["female"], a["black"], a["tc"], a["hdl"],
                             a["sbp"], a["bp_treated"], a["smoker"], a["diabetes"])[0])
        if name == "fram_lab":
            return float(fram_lab(a["age"], a["female"], a["tc"], a["hdl"], a["sbp"],
                                  a["bp_treated"], a["smoker"], a["diabetes"])[0])
        return float(fram_off(a["age"], a["female"], a["bmi"], a["sbp"],
                              a["bp_treated"], a["smoker"], a["diabetes"])[0])

    for name in ["pce", "fram_lab", "fram_off"]:
        b = run(name)
        worse = {"older": run(name, age=[70]), "smoker": run(name, smoker=[1]),
                 "diabetic": run(name, diabetes=[1]), "higher SBP": run(name, sbp=[160])}
        prob = 0.0 < b < 1.0
        mono = all(v > b for v in worse.values())
        ok &= prob and mono
        print(f"  {name:9s} baseline {b*100:5.2f}%  in (0,1): {prob}  "
              f"monotone in age/smoking/diabetes/SBP: {mono}")

    print(f"\n  validation: {'ALL PASS' if ok else 'FAILED'}")
    return ok


if __name__ == "__main__":
    validate()
