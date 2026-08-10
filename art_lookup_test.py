"""
Compare on-ART HIV mortality rates between the two competing implementations
under review right now:

  - OLD: the ARTMortalityTable-style lookup table (age x ART-duration x CD4
    x sex x adherence), from branch `art_adherence`, `HIV.get_art_mortality_hazard`.
  - NEW: PR #561 (`feat/simplify-art-mortality`)'s closed-form function --
    a single baseline rate scaled by flat sex/adherence multipliers and
    piecewise age/CD4 multipliers, with the ART-duration axis dropped
    (optional, off by default).

Both implementations' parameters are embedded directly below (pulled via
`git show <branch>:stisim/diseases/hiv.py`), rather than importing `stisim`,
since only one branch's version of hiv.py can be checked out/importable at a
time -- this way the comparison doesn't depend on which branch happens to be
on disk.

Produces a 2 (method) x 4 (sex x adherence) grid of CD4-by-age heatmaps of
the annual on-ART mortality RATE (not converted to a per-timestep
probability -- we're comparing the underlying rate functions, not running a
sim), all on a shared log color scale for direct visual comparison, plus a
printed new/old ratio summary per panel as a numeric "ballpark" check.

Caveat: the OLD table has an ART-duration-since-initiation axis that the NEW
function doesn't use by default. Per her commit message and the PR's own docs
diff, this isn't an average over duration -- she argues STIsim's CD4-
reconstitution curve already reproduces the *direction* of the duration
effect (mortality falls as CD4 recovers over time on ART) but only "~15-30%"
of its *magnitude*, and explicitly documents that the NEW function understates
excess mortality in the first 1-2 years on ART relative to EMOD. So OLD is
evaluated here at a single, explicitly-chosen duration bin (not averaged) --
set DURATION_BIN below (see OLD_DUR_BIN_LABELS for the options) -- and
comparing against the early bins (`<6mo`, `6-12mo`) is precisely where her own
caveat says to expect the biggest gap.

Also worth knowing: the docs table in this PR lists `art_death_rate = 0.0186`,
but the actual code default is `0.00554` -- a ~3.4x discrepancy between the
two. NEW_ART_DEATH_RATE below uses the code value (what actually runs), so if
that's the stale one, every new/old ratio in this script is ~3.4x too low.
Worth confirming with her which value is correct.
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

# ----------------------------------------------------------------------
# OLD: ARTMortalityTable-style lookup (branch `art_adherence`, hiv.py)
# ----------------------------------------------------------------------
OLD_AGE_BINS = np.array([25, 35, 45, 125])  # years; upper edges of <25/25-35/35-45/45+
OLD_DUR_BINS = np.array([182, 365, 730, 1095, 45625])  # days since ART initiation
OLD_CD4_BINS = np.array([0, 25, 74.5, 149.5, 274.5, 424.5, 624.5])  # CD4 count

OLD_EFFECTIVE_MALE = np.array([
    [[0.2015, 0.2015, 0.1128, 0.0625, 0.0312, 0.0206, 0.0162],
     [0.2176, 0.2176, 0.1219, 0.0675, 0.0337, 0.0223, 0.0175],
     [0.2350, 0.2350, 0.1316, 0.0729, 0.0364, 0.0240, 0.0189],
     [0.2538, 0.2538, 0.1421, 0.0787, 0.0393, 0.0260, 0.0205]],
    [[0.0875, 0.0875, 0.0490, 0.0271, 0.0136, 0.0062, 0.0041],
     [0.0945, 0.0945, 0.0529, 0.0293, 0.0146, 0.0067, 0.0044],
     [0.1021, 0.1021, 0.0572, 0.0316, 0.0158, 0.0073, 0.0047],
     [0.1102, 0.1102, 0.0617, 0.0342, 0.0171, 0.0079, 0.0051]],
    [[0.0255, 0.0255, 0.0181, 0.0128, 0.0085, 0.0058, 0.0038],
     [0.0288, 0.0288, 0.0204, 0.0145, 0.0096, 0.0065, 0.0043],
     [0.0326, 0.0326, 0.0231, 0.0164, 0.0108, 0.0074, 0.0049],
     [0.0368, 0.0368, 0.0261, 0.0185, 0.0123, 0.0084, 0.0055]],
    [[0.0164, 0.0164, 0.0116, 0.0083, 0.0055, 0.0037, 0.0025],
     [0.0186, 0.0186, 0.0131, 0.0093, 0.0062, 0.0042, 0.0042],
     [0.0210, 0.0210, 0.0148, 0.0106, 0.0070, 0.0048, 0.0048],
     [0.0237, 0.0237, 0.0168, 0.0119, 0.0079, 0.0054, 0.0054]],
    [[0.0119, 0.0119, 0.0081, 0.0066, 0.0033, 0.0033, 0.0033],
     [0.0135, 0.0135, 0.0092, 0.0074, 0.0037, 0.0037, 0.0037],
     [0.0152, 0.0152, 0.0103, 0.0084, 0.0042, 0.0042, 0.0042],
     [0.0172, 0.0172, 0.0117, 0.0095, 0.0047, 0.0047, 0.0047]],
])

OLD_NONSUPP_MALE = np.array([
    [[0.2015, 0.2015, 0.1128, 0.0625, 0.0312, 0.0206, 0.0162],
     [0.2176, 0.2176, 0.1219, 0.0675, 0.0337, 0.0223, 0.0175],
     [0.2350, 0.2350, 0.1316, 0.0729, 0.0364, 0.0240, 0.0189],
     [0.2538, 0.2538, 0.1421, 0.0787, 0.0393, 0.0260, 0.0205]],
    [[0.1715, 0.1715, 0.5600, 0.3100, 0.1550, 0.0713, 0.0465],
     [0.1852, 0.1852, 0.6048, 0.3348, 0.1674, 0.0770, 0.0502],
     [0.2000, 0.2000, 0.6532, 0.3616, 0.1808, 0.0832, 0.0542],
     [0.2160, 0.2160, 0.7054, 0.3905, 0.1953, 0.0898, 0.0586]],
    [[0.0532, 0.0532, 0.0362, 0.0293, 0.0171, 0.0116, 0.0095],
     [0.0601, 0.0601, 0.0409, 0.0331, 0.0193, 0.0131, 0.0107],
     [0.0679, 0.0679, 0.0462, 0.0374, 0.0218, 0.0148, 0.0121],
     [0.0768, 0.0768, 0.0522, 0.0422, 0.0246, 0.0168, 0.0137]],
    [[0.0335, 0.0335, 0.0228, 0.0184, 0.0108, 0.0073, 0.0060],
     [0.0379, 0.0379, 0.0258, 0.0208, 0.0122, 0.0083, 0.0068],
     [0.0428, 0.0428, 0.0291, 0.0235, 0.0137, 0.0094, 0.0076],
     [0.0484, 0.0484, 0.0329, 0.0266, 0.0155, 0.0106, 0.0086]],
    [[0.0234, 0.0234, 0.0159, 0.0129, 0.0091, 0.0069, 0.0064],
     [0.0265, 0.0265, 0.0180, 0.0145, 0.0103, 0.0077, 0.0073],
     [0.0299, 0.0299, 0.0203, 0.0164, 0.0116, 0.0088, 0.0082],
     [0.0338, 0.0338, 0.0230, 0.0186, 0.0131, 0.0099, 0.0093]],
])

OLD_EFFECTIVE_FEMALE = np.array([
    [[0.2015, 0.2015, 0.0993, 0.0518, 0.0259, 0.0171, 0.0135],
     [0.2156, 0.2156, 0.1062, 0.0554, 0.0277, 0.0183, 0.0144],
     [0.2307, 0.2307, 0.1137, 0.0593, 0.0296, 0.0196, 0.0154],
     [0.2468, 0.2468, 0.1216, 0.0634, 0.0317, 0.0209, 0.0165]],
    [[0.0875, 0.0875, 0.0431, 0.0225, 0.0112, 0.0052, 0.0034],
     [0.0936, 0.0936, 0.0461, 0.0241, 0.0120, 0.0055, 0.0036],
     [0.1002, 0.1002, 0.0494, 0.0257, 0.0129, 0.0059, 0.0039],
     [0.1072, 0.1072, 0.0528, 0.0276, 0.0138, 0.0063, 0.0041]],
    [[0.0241, 0.0241, 0.0166, 0.0135, 0.0067, 0.0044, 0.0044],
     [0.0262, 0.0262, 0.0181, 0.0147, 0.0073, 0.0048, 0.0048],
     [0.0286, 0.0286, 0.0197, 0.0160, 0.0080, 0.0052, 0.0052],
     [0.0312, 0.0312, 0.0215, 0.0175, 0.0087, 0.0057, 0.0057]],
    [[0.0149, 0.0149, 0.0103, 0.0084, 0.0042, 0.0042, 0.0042],
     [0.0163, 0.0163, 0.0112, 0.0091, 0.0046, 0.0046, 0.0046],
     [0.0177, 0.0177, 0.0122, 0.0099, 0.0050, 0.0050, 0.0050],
     [0.0193, 0.0193, 0.0133, 0.0108, 0.0054, 0.0054, 0.0054]],
    [[0.0084, 0.0084, 0.0057, 0.0046, 0.0023, 0.0023, 0.0023],
     [0.0092, 0.0092, 0.0062, 0.0051, 0.0025, 0.0025, 0.0025],
     [0.0100, 0.0100, 0.0068, 0.0055, 0.0028, 0.0028, 0.0028],
     [0.0109, 0.0109, 0.0074, 0.0060, 0.0030, 0.0030, 0.0030]],
])

OLD_NONSUPP_FEMALE = np.array([
    [[0.2015, 0.2015, 0.1128, 0.0625, 0.0312, 0.0206, 0.0162],
     [0.2176, 0.2176, 0.1219, 0.0675, 0.0337, 0.0223, 0.0175],
     [0.2350, 0.2350, 0.1316, 0.0729, 0.0364, 0.0240, 0.0189],
     [0.2538, 0.2538, 0.1421, 0.0787, 0.0393, 0.0260, 0.0205]],
    [[0.1837, 0.1837, 0.0845, 0.0441, 0.0220, 0.0101, 0.0066],
     [0.1965, 0.1965, 0.0904, 0.0472, 0.0236, 0.0108, 0.0071],
     [0.2103, 0.2103, 0.0967, 0.0505, 0.0252, 0.0116, 0.0076],
     [0.2250, 0.2250, 0.1035, 0.0540, 0.0270, 0.0124, 0.0081]],
    [[0.0461, 0.0461, 0.0318, 0.0258, 0.0151, 0.0103, 0.0084],
     [0.0502, 0.0502, 0.0346, 0.0281, 0.0164, 0.0113, 0.0091],
     [0.0547, 0.0547, 0.0378, 0.0306, 0.0179, 0.0123, 0.0100],
     [0.0596, 0.0596, 0.0412, 0.0334, 0.0195, 0.0134, 0.0109]],
    [[0.0286, 0.0286, 0.0197, 0.0160, 0.0113, 0.0086, 0.0080],
     [0.0311, 0.0311, 0.0215, 0.0174, 0.0124, 0.0094, 0.0087],
     [0.0339, 0.0339, 0.0234, 0.0190, 0.0135, 0.0102, 0.0095],
     [0.0370, 0.0370, 0.0255, 0.0207, 0.0147, 0.0111, 0.0104]],
    [[0.0161, 0.0161, 0.0110, 0.0089, 0.0063, 0.0047, 0.0044],
     [0.0176, 0.0176, 0.0119, 0.0097, 0.0068, 0.0052, 0.0048],
     [0.0192, 0.0192, 0.0130, 0.0105, 0.0075, 0.0056, 0.0053],
     [0.0209, 0.0209, 0.0142, 0.0115, 0.0081, 0.0061, 0.0057]],
])

OLD_TABLES = {
    ('m', True): OLD_EFFECTIVE_MALE, ('m', False): OLD_NONSUPP_MALE,
    ('f', True): OLD_EFFECTIVE_FEMALE, ('f', False): OLD_NONSUPP_FEMALE,
}

# Labels for OLD_DUR_BINS' 5 duration-since-ART-initiation bins, in order.
OLD_DUR_BIN_LABELS = ['<6mo', '6-12mo', '1-2yr', '2-3yr', '>3yr']


def _digitize_clip(vals, bins):
    return np.clip(np.digitize(vals, bins), 0, len(bins) - 1)


def old_hazard(age, cd4, sex, effective, duration_bin='<6mo'):
    """
    Annual on-ART mortality rate per the OLD lookup table, at a single,
    explicitly-chosen ART-duration bin (no averaging across bins).

    Args:
        age, cd4: broadcastable arrays (e.g. same-shape 2D grids from meshgrid)
        sex: 'm' or 'f'
        effective: True for effective (suppressive) ART, False for non-suppressive
        duration_bin: which ART-duration-since-initiation bin to pull from --
            one of OLD_DUR_BIN_LABELS: '<6mo', '6-12mo', '1-2yr', '2-3yr', '>3yr'.
            Defaults to '<6mo', the first-6-months-on-ART bin.
    """
    table = OLD_TABLES[(sex, effective)]  # shape (5, 4, 7): [duration_bin, age_bin, cd4_bin]
    age_idx = _digitize_clip(age, OLD_AGE_BINS)
    cd4_idx = _digitize_clip(cd4, OLD_CD4_BINS)
    dur_idx = OLD_DUR_BIN_LABELS.index(duration_bin)
    return table[dur_idx, age_idx, cd4_idx]


# ----------------------------------------------------------------------
# NEW: closed-form function (branch `feat/simplify-art-mortality`, PR #561, hiv.py)
# ----------------------------------------------------------------------
NEW_ART_DEATH_RATE = 0.00554  # Annual rate, male, suppressed, age<25, healthiest CD4
NEW_REL_DEATH_UNSUPP = 2.03   # Multiplier for non-suppressive ART
NEW_REL_DEATH_F = 0.74        # Multiplier for females, on ART
NEW_ART_DEATH_AGE = [  # (age_lo, age_hi, mult)
    (0, 25, 1.0),
    (25, 35, 1.10),
    (35, 45, 1.21),
    (45, 125, 1.32),
]
NEW_ART_DEATH_CD4 = [  # (cd4_lo, cd4_hi, mult)
    (0, 25, 7.12),
    (25, 74.5, 4.81),
    (74.5, 149.5, 3.28),
    (149.5, 274.5, 1.82),
    (274.5, 424.5, 1.22),
    (424.5, np.inf, 1.0),
]
# NB: art_death_dur is None (off) by default on the PR branch, so it's
# intentionally not included here -- this comparison uses the PR's own default.


def new_hazard(age, cd4, sex, effective):
    """ Annual on-ART mortality rate per the NEW closed-form function. """
    age = np.asarray(age, dtype=float)
    cd4 = np.asarray(cd4, dtype=float)
    rate = np.full(age.shape, NEW_ART_DEATH_RATE, dtype=float)
    if not effective:
        rate = rate * NEW_REL_DEATH_UNSUPP
    if sex == 'f':
        rate = rate * NEW_REL_DEATH_F

    age_mult = np.ones_like(rate)
    for lo, hi, mult in NEW_ART_DEATH_AGE:
        age_mult[(age >= lo) & (age < hi)] = mult

    cd4_mult = np.ones_like(rate)
    for lo, hi, mult in NEW_ART_DEATH_CD4:
        cd4_mult[(cd4 >= lo) & (cd4 < hi)] = mult

    return rate * age_mult * cd4_mult


def compute_grids(age_vals, cd4_vals, duration_bin='<6mo'):
    """
    Returns dict[(sex, effective)] -> dict('old'/'new' -> 2D array), shape
    (len(age_vals), len(cd4_vals)).
    """
    age_grid, cd4_grid = np.meshgrid(age_vals, cd4_vals, indexing='ij')
    grids = {}
    for sex in ('m', 'f'):
        for effective in (True, False):
            grids[(sex, effective)] = dict(
                old=old_hazard(age_grid, cd4_grid, sex, effective, duration_bin=duration_bin),
                new=new_hazard(age_grid, cd4_grid, sex, effective),
            )
    return grids


def compute_global_range(age_vals, cd4_vals, duration_bins=None):
    """
    Compute the (vmin, vmax) of hazard values across ALL of the given OLD
    duration bins (default: all 5, OLD_DUR_BIN_LABELS) and the NEW function,
    for the same age/CD4 grid used elsewhere. Pass the result to
    plot_heatmaps(..., vmin=, vmax=) to freeze the color scale across
    multiple compute_grids()/plot_heatmaps() calls for different
    duration_bin values, so the same rate always maps to the same color.
    """
    if duration_bins is None:
        duration_bins = OLD_DUR_BIN_LABELS
    all_vals = []
    for duration_bin in duration_bins:
        grids = compute_grids(age_vals, cd4_vals, duration_bin=duration_bin)
        for g in grids.values():
            all_vals.append(g['old'].ravel())
            all_vals.append(g['new'].ravel())
    all_vals = np.concatenate(all_vals)
    return all_vals.min(), all_vals.max()


def print_ratio_summary(grids):
    """ Print new/old ratio summary stats (median, 5th/95th pctile) per panel. """
    print('New-vs-old ratio summary (new/old; 1.0 = identical):')
    for (sex, effective), g in grids.items():
        ratio = g['new'] / g['old']
        label = f"sex={sex}, {'effective' if effective else 'nonsuppressive'}"
        print(f'  {label:28s}: median={np.median(ratio):.2f}, '
              f'p5={np.percentile(ratio, 5):.2f}, p95={np.percentile(ratio, 95):.2f}, '
              f'min={ratio.min():.2f}, max={ratio.max():.2f}')


def plot_heatmaps(age_vals, cd4_vals, grids, duration_bin='<6mo', vmin=None, vmax=None):
    """
    2 (old/new) x 4 (sex x adherence) grid of CD4-by-age heatmaps, shared
    log color scale.

    By default (vmin/vmax=None), the color scale is auto-ranged from this
    call's own grids -- fine for a single plot, but it means the same rate
    value can get a different color in separate calls with different
    duration_bin selections. Pass explicit vmin/vmax (e.g. from
    compute_global_range()) to freeze the scale across multiple calls.
    """
    panels = [('m', True), ('m', False), ('f', True), ('f', False)]
    if vmin is None or vmax is None:
        all_vals = np.concatenate([g[m].ravel() for g in grids.values() for m in ('old', 'new')])
        vmin = all_vals.min() if vmin is None else vmin
        vmax = all_vals.max() if vmax is None else vmax
    norm = LogNorm(vmin=vmin, vmax=vmax)

    fig, axes = plt.subplots(2, 4, figsize=(16, 7), sharex=True, sharey=True)
    im = None
    for col, key in enumerate(panels):
        sex, effective = key
        label = f"{'M' if sex == 'm' else 'F'}, {'effective' if effective else 'nonsuppressive'}"
        for row, method in enumerate(('old', 'new')):
            ax = axes[row, col]
            im = ax.pcolormesh(cd4_vals, age_vals, grids[key][method], norm=norm, cmap='viridis', shading='auto')
            if row == 0:
                ax.set_title(label)
            if col == 0:
                ax.set_ylabel(f"{'OLD (table)' if method == 'old' else 'NEW (function)'}\nAge (years)")
            if row == 1:
                ax.set_xlabel('CD4 count')
    fig.colorbar(im, ax=axes, label='Annual on-ART mortality rate', shrink=0.8)
    fig.suptitle(f'OLD (lookup table, duration bin={duration_bin}) vs. NEW (closed-form function) '
                 f'on-ART mortality: age x CD4, by sex and adherence')
    return fig


def compute_adherence_ratio_grids(grids):
    """
    Ratio of non-suppressive to effective ART mortality (nonsuppressive /
    effective) at every age/CD4 grid point, for each (sex, method)
    combination -- 2 sexes x 2 methods (old/new) = 4 grids in total, each
    shape (len(age_vals), len(cd4_vals)). Takes the output of compute_grids().
    """
    ratios = {}
    for sex in ('m', 'f'):
        for method in ('old', 'new'):
            ratios[(sex, method)] = grids[(sex, False)][method] / grids[(sex, True)][method]
    return ratios


def print_adherence_ratio_tables(age_vals, cd4_vals, ratios):
    """ Print each of the 4 (sex, method) ratio grids as a readable age x CD4 table. """
    for (sex, method), ratio in ratios.items():
        df = pd.DataFrame(ratio, index=age_vals, columns=cd4_vals)
        df.index.name = 'age'
        df.columns.name = 'cd4'
        print(f'--- Non-suppressive/effective mortality ratio: sex={sex}, method={method.upper()} ---')
        print(df.round(2).to_string())
        print()


def plot_adherence_ratio_heatmaps(age_vals, cd4_vals, ratios, duration_bin='<6mo', vmin=None, vmax=None):
    """ 2 (method: old/new) x 2 (sex) grid of CD4-by-age heatmaps of the non-suppressive/effective ratio. """
    if vmin is None or vmax is None:
        all_vals = np.concatenate([r.ravel() for r in ratios.values()])
        vmin = all_vals.min() if vmin is None else vmin
        vmax = all_vals.max() if vmax is None else vmax
    norm = LogNorm(vmin=vmin, vmax=vmax)

    fig, axes = plt.subplots(2, 2, figsize=(9, 7), sharex=True, sharey=True)
    im = None
    for row, method in enumerate(('old', 'new')):
        for col, sex in enumerate(('m', 'f')):
            ax = axes[row, col]
            im = ax.pcolormesh(cd4_vals, age_vals, ratios[(sex, method)], norm=norm, cmap='magma', shading='auto')
            if row == 0:
                ax.set_title(f'sex={sex}')
            if col == 0:
                ax.set_ylabel(f"{'OLD' if method == 'old' else 'NEW'}\nAge (years)")
            if row == 1:
                ax.set_xlabel('CD4 count')
    fig.colorbar(im, ax=axes, label='Non-suppressive / effective mortality ratio', shrink=0.8)
    fig.suptitle(f'Non-suppressive vs. effective ART mortality ratio (OLD duration bin={duration_bin})')
    return fig


if __name__ == '__main__':

    # Which OLD ART-duration bin to compare against -- one of OLD_DUR_BIN_LABELS:
    # '<6mo', '6-12mo', '1-2yr', '2-3yr', '>3yr'.
    DURATION_BIN = '<6mo'

    age_vals = np.arange(15, 66, 2)
    cd4_vals = np.arange(0, 810, 20)

    # Frozen so the same rate value always gets the same color, whichever
    # DURATION_BIN you plot -- computed once across all 5 duration bins.
    VMIN, VMAX = compute_global_range(age_vals, cd4_vals)

    grids = compute_grids(age_vals, cd4_vals, duration_bin=DURATION_BIN)
    print_ratio_summary(grids)
    plot_heatmaps(age_vals, cd4_vals, grids, duration_bin=DURATION_BIN, vmin=VMIN, vmax=VMAX)

    # Non-suppressive vs. effective ART mortality ratio: 4 tables (2 sexes x
    # 2 methods), each age x CD4. Printed on a coarser grid than the fine
    # heatmap grid above, since a 26x41 table isn't readable as text.
    print()
    table_age_vals = np.array([20, 30, 40, 50, 60])
    table_cd4_vals = np.array([0, 100, 200, 300, 400, 500, 600, 700, 800])
    table_grids = compute_grids(table_age_vals, table_cd4_vals, duration_bin=DURATION_BIN)
    adherence_ratios = compute_adherence_ratio_grids(table_grids)
    print_adherence_ratio_tables(table_age_vals, table_cd4_vals, adherence_ratios)

    # Same ratio, as a heatmap on the finer grid for a visual read
    fine_ratios = compute_adherence_ratio_grids(grids)
    plot_adherence_ratio_heatmaps(age_vals, cd4_vals, fine_ratios, duration_bin=DURATION_BIN)

    plt.show()
