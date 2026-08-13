"""
ART mortality test/validation script -- tracks and plots HIV-related
mortality both via STIsim simulation and via closed-form (no-simulation)
numerical evaluation of the exact same formulas.

Consolidates what were previously three separate files (art_adherence_test.py,
art_adherence_test_2.py, art_adherence_test_3.py) into one, updated for the
CURRENT `HIV.get_art_mortality_hazard` (anchored to the off-ART CD4 hazard;
see stisim/diseases/hiv.py and art_implementation_notes.md). The earlier
files' "old" (`use_art_mortality_table=False`) vs. "new" (PR #561's
`art_death_rate`) comparison is gone -- both were superseded when
get_art_mortality_hazard was rewritten to fix a real bug (on-ART mortality
could exceed off-ART mortality at high CD4; see below), so there's only one
current implementation to validate now, not two to compare.

`get_art_mortality_hazard` computes:

    rate = off_art_rate(cd4) * rel_art_mortality[effective? / sex] * age_mult(age) * rel_death_f?

The non-suppressive relative-mortality factor is sex-specific
(`rel_art_mortality_unsupp_m`/`_f`), so the non-suppressive/effective
mortality RATIO can differ by sex -- currently 2x higher for men
(`0.7/0.25=2.8`) than for women (`0.35/0.25=1.4`); `rel_death_f` cancels out
of that ratio since it applies equally to both adherence categories.
Anchoring to `off_art_rate(cd4)` (the same CD4 table `make_p_hiv_death` uses
off-ART) guarantees, BY CONSTRUCTION, that on-ART mortality can never exceed
off-ART mortality at the same CD4 count -- as long as
`rel_art_mortality_unsupp_m * (largest age/duration multiplier) <= 1`
(true for the shipped defaults: `0.7 * 1.32 = 0.924`; `_f` is smaller, so
it's covered too).

Provides:
  1. `AgeCD4ARTMortalityTracker` -- a single tracker that records deaths and
     person-years per timestep, stratified by age bracket x CD4 bracket x
     sex x ART status (off_art / on_effective_art / on_nonsuppressive_art).
     One simulation run feeds BOTH simulation-based views below.
  2a. Yearly mortality RATE time series, pooled across age/CD4/sex, by
      category ('any' / 'effective' / 'nonsuppressive' / 'off_art').
  2b. Age x CD4 faceted mortality rate (for a chosen sex), aggregated over
      the full sim, plus a pass/fail check (`check_off_art_higher`) that
      off-ART mortality is never lower than on-ART mortality at the same
      age/CD4/sex -- the actual regression test for the invariant above.
  3. Closed-form (no-simulation) views: the exact same off-ART/on-ART rate
     formulas evaluated directly against live `sti.HIVPars()` values --
     numeric tables (2, one per sex) and heatmaps -- so simulation noise and
     cell sample-size limitations can be ruled out when interpreting 2a/2b.

For the ART retention-*duration* analysis (how long people stay on ART, as
opposed to how likely they are to die while on it), see art_duration_test.py.
For a frozen historical comparison against the older, now-superseded
ARTMortalityTable lookup-table implementation, see art_lookup_test.py (not
live code -- a snapshot from before this fix).
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import starsim as ss
import stisim as sti


def _bin_label(lo, hi, open_ended_at):
    return f'{lo}-{hi}' if hi < open_ended_at else f'{lo}+'


def _bin_lookup(value, bins, open_ended_at):
    for lo, hi in bins:
        if lo <= value < hi:
            return _bin_label(lo, hi, open_ended_at)
    return None


class AgeCD4ARTMortalityTracker(ss.Analyzer):
    """
    Tracks deaths and person-years, stratified by (user-specified) age
    bracket, CD4 bracket, sex, and ART status: 'off_art' (infected but not
    on ART, i.e. ART-naive or discontinued), 'on_effective_art', or
    'on_nonsuppressive_art'.

    Deaths are attributed using each agent's (age bracket, CD4 bracket, sex,
    status) as of the END of the PREVIOUS timestep, cached across steps --
    same pattern as ARTStatusTracker in art_transmission_test.py, and for
    the same reason: HIV.step_die() clears on_effective_art/
    on_nonsuppressive_art (and zeroes cd4) for agents who die this step, so
    the live state can't be read after the fact.
    """
    CATEGORIES = ['off_art', 'on_effective_art', 'on_nonsuppressive_art']
    AGE_OPEN_ENDED_AT = 200
    CD4_OPEN_ENDED_AT = 100_000

    def __init__(self, age_bins, cd4_bins):
        super().__init__()
        self.age_bins = age_bins
        self.cd4_bins = cd4_bins
        self.age_labels = [_bin_label(lo, hi, self.AGE_OPEN_ENDED_AT) for lo, hi in age_bins]
        self.cd4_labels = [_bin_label(lo, hi, self.CD4_OPEN_ENDED_AT) for lo, hi in cd4_bins]
        self._cache = {}  # uid -> (age_label, cd4_label, sex, category), as of end of previous timestep

    def init_results(self):
        super().init_results()
        results = []
        for age_label in self.age_labels:
            for cd4_label in self.cd4_labels:
                for sex in ('m', 'f'):
                    for cat in self.CATEGORIES:
                        key = f'{age_label}_{cd4_label}_{sex}_{cat}'
                        results += [
                            ss.Result(f'deaths_{key}', dtype=int),
                            ss.Result(f'person_years_{key}', dtype=float),
                        ]
        self.define_results(*results)

    def step(self):
        hiv = self.sim.diseases.hiv
        ti = self.ti
        ppl = self.sim.people
        dt_years = self.dt.years

        # Deaths this step, attributed via the pre-death cache
        died_uids = (hiv.ti_dead == ti).uids
        for uid in died_uids:
            cached = self._cache.pop(int(uid), None)
            if cached is None:
                continue
            age_label, cd4_label, sex, cat = cached
            self.results[f'deaths_{age_label}_{cd4_label}_{sex}_{cat}'][ti] += 1

        # Person-time this step, by CURRENT age/CD4/sex/status
        cache = {}
        cat_uids = [
            ('off_art', (hiv.infected & ~hiv.on_art).uids),
            ('on_effective_art', hiv.on_effective_art.uids),
            ('on_nonsuppressive_art', hiv.on_nonsuppressive_art.uids),
        ]
        for cat, uids in cat_uids:
            for uid in uids:
                u = ss.uids([uid])
                age_label = _bin_lookup(ppl.age[u][0], self.age_bins, self.AGE_OPEN_ENDED_AT)
                cd4_label = _bin_lookup(hiv.cd4[u][0], self.cd4_bins, self.CD4_OPEN_ENDED_AT)
                if age_label is None or cd4_label is None:
                    continue
                sex = 'm' if bool(ppl.male[u][0]) else 'f'
                key = f'{age_label}_{cd4_label}_{sex}_{cat}'
                self.results[f'person_years_{key}'][ti] += dt_years
                cache[int(uid)] = (age_label, cd4_label, sex, cat)
        self._cache = cache


def run_sim(age_bins, cd4_bins, n_agents=20_000, start=2000, stop=2030, coverage_target=0.6):
    """
    Run a single sim with realistic ART testing/coverage and the tracker
    attached. Returns (tracker, yearvec) -- one run feeds both the yearly
    time-series view and the age x CD4 faceted view below.
    """
    hiv_test = sti.HIVTest(test_prob_data=0.3, start=start)
    ramp_end = min(start + 4, stop)
    art = sti.ART(coverage={'year': [start, ramp_end, stop], 'value': [0, coverage_target, coverage_target]}, p_effective_art=0.5)
    sim = sti.Sim(
        diseases=sti.HIV(init_prev=ss.bernoulli(p=0.2)),
        n_agents=n_agents, start=start, stop=stop,
        interventions=[hiv_test, art],
        analyzers=[AgeCD4ARTMortalityTracker(age_bins, cd4_bins)],
    )
    sim.run(verbose=0)
    return sim.analyzers.agecd4artmortalitytracker, sim.t.yearvec


def yearly_mortality_by_category(tracker, yearvec, age_bins, cd4_bins):
    """
    Total mortality per calendar year, pooled across all age/CD4 brackets
    and both sexes, broken out into four categories: 'any' (on ART of any
    kind, i.e. effective + non-suppressive combined), 'effective',
    'nonsuppressive', and 'off_art'. One row per (year, category).
    """
    age_labels = [_bin_label(lo, hi, AgeCD4ARTMortalityTracker.AGE_OPEN_ENDED_AT) for lo, hi in age_bins]
    cd4_labels = [_bin_label(lo, hi, AgeCD4ARTMortalityTracker.CD4_OPEN_ENDED_AT) for lo, hi in cd4_bins]
    years = np.floor(yearvec).astype(int)
    unique_years = sorted(set(years))

    deaths_by_cat = {cat: np.zeros(len(yearvec)) for cat in AgeCD4ARTMortalityTracker.CATEGORIES}
    py_by_cat = {cat: np.zeros(len(yearvec)) for cat in AgeCD4ARTMortalityTracker.CATEGORIES}
    for cat in AgeCD4ARTMortalityTracker.CATEGORIES:
        for age_label in age_labels:
            for cd4_label in cd4_labels:
                for sex in ('m', 'f'):
                    key = f'{age_label}_{cd4_label}_{sex}_{cat}'
                    deaths_by_cat[cat] += tracker.results[f'deaths_{key}'].values
                    py_by_cat[cat] += tracker.results[f'person_years_{key}'].values

    deaths_any = deaths_by_cat['on_effective_art'] + deaths_by_cat['on_nonsuppressive_art']
    py_any = py_by_cat['on_effective_art'] + py_by_cat['on_nonsuppressive_art']

    by_category = [
        ('any', deaths_any, py_any),
        ('effective', deaths_by_cat['on_effective_art'], py_by_cat['on_effective_art']),
        ('nonsuppressive', deaths_by_cat['on_nonsuppressive_art'], py_by_cat['on_nonsuppressive_art']),
        ('off_art', deaths_by_cat['off_art'], py_by_cat['off_art']),
    ]

    rows = []
    for year in unique_years:
        mask = years == year
        for category, deaths_arr, py_arr in by_category:
            deaths = deaths_arr[mask].sum()
            person_years = py_arr[mask].sum()
            rate = deaths / person_years if person_years > 0 else np.nan
            rows.append(dict(year=year, category=category, deaths=deaths,
                              person_years=person_years, annual_death_rate=rate))
    return pd.DataFrame(rows)


def plot_yearly_mortality(yearly):
    """ Annual death rate (deaths / person-years) per year, one panel per category. """
    categories = ['any', 'effective', 'nonsuppressive', 'off_art']
    fig, axes = plt.subplots(1, len(categories), figsize=(5 * len(categories), 4), sharey=True)
    for ax, category in zip(axes, categories):
        sub = yearly.query('category == @category').sort_values('year')
        ax.plot(sub['year'], sub['annual_death_rate'], marker='o', color='tab:blue')
        ax.set_title(category)
        ax.set_xlabel('Year')
    axes[0].set_ylabel('Annual death rate (deaths / person-years)')
    fig.suptitle('Annual HIV death rate per year, by ART status')
    fig.tight_layout()
    plt.close(fig)  # avoid double-display in notebooks: inline backend auto-shows open figures, AND returning one triggers rich-display
    return fig


def summarize_mortality_by_status(tracker, age_bins, cd4_bins, sex):
    """
    Tidy DataFrame, one row per (age_band, cd4_band, category), for a given
    sex ('m' or 'f'), aggregated over the FULL simulation (not per-year).
    category is one of 'total' (all three ART-status categories combined),
    'off_art', 'on_effective_art', 'on_nonsuppressive_art'.
    """
    age_labels = [_bin_label(lo, hi, AgeCD4ARTMortalityTracker.AGE_OPEN_ENDED_AT) for lo, hi in age_bins]
    cd4_labels = [_bin_label(lo, hi, AgeCD4ARTMortalityTracker.CD4_OPEN_ENDED_AT) for lo, hi in cd4_bins]

    rows = []
    for age_label in age_labels:
        for cd4_label in cd4_labels:
            cat_deaths, cat_py = {}, {}
            for cat in AgeCD4ARTMortalityTracker.CATEGORIES:
                key = f'{age_label}_{cd4_label}_{sex}_{cat}'
                cat_deaths[cat] = tracker.results[f'deaths_{key}'].values.sum()
                cat_py[cat] = tracker.results[f'person_years_{key}'].values.sum()

            entries = [('total', sum(cat_deaths.values()), sum(cat_py.values()))]
            entries += [(cat, cat_deaths[cat], cat_py[cat]) for cat in AgeCD4ARTMortalityTracker.CATEGORIES]

            for category, deaths, person_years in entries:
                rate = deaths / person_years if person_years > 0 else np.nan
                rows.append(dict(age_band=age_label, cd4_band=cd4_label, category=category,
                                  deaths=deaths, person_years=person_years, annual_death_rate=rate))
    return pd.DataFrame(rows)


def check_off_art_higher(summary):
    """
    Print a pass/fail check: is off_art's death rate higher than both
    on-ART categories, in every age x CD4 cell with enough data to compare?

    Since get_art_mortality_hazard anchors on-ART mortality to the off-ART
    CD4 hazard (see module docstring), this is expected to ALWAYS pass -- a
    violation here (beyond stochastic noise in thin cells) means the
    invariant has broken, e.g. via an age_mult/art_death_dur override large
    enough to push rel_art_mortality_unsupp_m/f * (largest multiplier) above 1.
    """
    categories = ['off_art', 'on_effective_art', 'on_nonsuppressive_art']
    wide = summary[summary['category'].isin(categories)].pivot_table(
        index=['age_band', 'cd4_band'], columns='category', values='annual_death_rate')
    wide = wide.dropna()  # drop cells with no person-time in one of the categories

    violations = wide[(wide['off_art'] <= wide['on_effective_art']) | (wide['off_art'] <= wide['on_nonsuppressive_art'])]
    print(f'Checked {len(wide)} age x CD4 cells with data in all three categories.')
    if len(violations):
        print(f'{len(violations)} cell(s) where off_art is NOT higher than both on-ART categories:')
        print(violations.round(4).to_string())
    else:
        print('off_art death rate is higher than both on-ART categories in every comparable cell.')
    return wide


def plot_mortality_by_status(summary, age_labels, cd4_labels, sex):
    """
    Grid of bar-chart panels faceted by age (rows) x CD4 bracket (columns).
    Each panel compares annual death rate for total / off_art /
    on_effective_art / on_nonsuppressive_art, side by side. Y axes are
    independent per panel (not shared) since CD4 drives ~100x differences
    in scale -- a shared axis would flatten the high-CD4 panels unreadably.
    """
    categories = ['total', 'off_art', 'on_effective_art', 'on_nonsuppressive_art']
    colors = {'total': 'tab:gray', 'off_art': 'tab:red', 'on_effective_art': 'tab:green', 'on_nonsuppressive_art': 'tab:orange'}
    short_labels = {'total': 'total', 'off_art': 'off', 'on_effective_art': 'eff', 'on_nonsuppressive_art': 'nonsupp'}

    fig, axes = plt.subplots(len(age_labels), len(cd4_labels),
                              figsize=(3.2 * len(cd4_labels), 2.6 * len(age_labels)),
                              sharey=False, squeeze=False)
    for i, age_label in enumerate(age_labels):
        for j, cd4_label in enumerate(cd4_labels):
            ax = axes[i, j]
            sub = summary.query('age_band == @age_label and cd4_band == @cd4_label')
            rates = [sub.loc[sub['category'] == cat, 'annual_death_rate'].iloc[0] for cat in categories]
            ax.bar(range(len(categories)), rates, color=[colors[c] for c in categories])
            ax.set_xticks(range(len(categories)))
            ax.set_xticklabels([short_labels[c] for c in categories], rotation=45, ha='right', fontsize=7)
            if i == 0:
                ax.set_title(f'CD4 {cd4_label}')
            if j == 0:
                ax.set_ylabel(f'Age {age_label}\nAnnual death rate')
    fig.suptitle(f'HIV mortality rate by ART status, age x CD4 (sex={sex})')
    fig.tight_layout()
    plt.close(fig)  # avoid double-display in notebooks: inline backend auto-shows open figures, AND returning one triggers rich-display
    return fig


def pull_hiv_pars():
    """
    Fresh sti.HIVPars() instance -- the single source of truth for every
    constant used below. Deliberately NOT hard-coded/copy-pasted (unlike
    art_lookup_test.py's OLD-branch tables, which had to be embedded because
    two different branches' code can't be imported in the same process) --
    this script only ever needs the currently-checked-out branch's own
    values, so pulling them live avoids the exact stale-duplicate problem
    that produced an earlier `art_death_rate` docs/code mismatch.
    """
    return sti.HIVPars()


def off_art_rate_analytic(cd4, pars=None):
    """ Closed-form annual off-ART mortality rate. Matches HIV.make_p_hiv_death exactly. """
    pars = pars or pull_hiv_pars()
    cd4 = np.asarray(cd4, dtype=float)
    idx = np.digitize(cd4, pars.cd4_death_bins)
    return pars.cd4_death_rates[idx] * pars.rel_death


def on_art_rate_analytic(age, cd4, sex, effective, pars=None):
    """
    Closed-form annual on-ART mortality rate. Matches HIV.get_art_mortality_hazard
    exactly (the art_death_dur term is omitted since it's None/off by default
    on this branch -- pass a HIVPars with art_death_dur set if you want it).

    Anchored to off_art_rate_analytic(cd4) -- same CD4 table used off-ART --
    scaled by a relative-mortality factor per adherence category (sex-specific
    for non-suppressive ART, so the non-suppressive/effective ratio can differ
    by sex), then by age/sex multipliers. This structurally guarantees on-ART
    <= off-ART at the same CD4 (the fix for the upper-bound violation found earlier).
    """
    pars = pars or pull_hiv_pars()
    age = np.asarray(age, dtype=float)
    cd4 = np.asarray(cd4, dtype=float)

    rate = off_art_rate_analytic(cd4, pars=pars) * np.ones_like(age)
    if effective:
        rel_art_mortality = pars.rel_art_mortality_effective
    else:
        rel_art_mortality = pars.rel_art_mortality_unsupp_m if sex == 'm' else pars.rel_art_mortality_unsupp_f
    rate = rate * rel_art_mortality
    if sex == 'f':
        rate = rate * pars.rel_death_f

    age_mult = np.ones_like(rate)
    for lo, hi, mult in pars.art_death_age:
        age_mult[(age >= lo) & (age < hi)] = mult

    return rate * age_mult


def compute_analytic_grids(age_vals, cd4_vals, pars=None):
    """
    Closed-form (no simulation) annual mortality rate grids, shape
    (len(age_vals), len(cd4_vals)), for both sexes and all three categories.
    Returns dict[sex] -> dict['off_art'/'on_effective_art'/'on_nonsuppressive_art'] -> 2D array.
    """
    pars = pars or pull_hiv_pars()
    age_grid, cd4_grid = np.meshgrid(age_vals, cd4_vals, indexing='ij')
    grids = {}
    for sex in ('m', 'f'):
        grids[sex] = dict(
            off_art=off_art_rate_analytic(cd4_grid, pars=pars) * np.ones_like(age_grid),
            on_effective_art=on_art_rate_analytic(age_grid, cd4_grid, sex, True, pars=pars),
            on_nonsuppressive_art=on_art_rate_analytic(age_grid, cd4_grid, sex, False, pars=pars),
        )
    return grids


def print_analytic_tables(age_vals, cd4_vals, grids):
    """
    Print 2 tables (one per sex): for each sex, one age x CD4 table per
    category (off_art / on_effective_art / on_nonsuppressive_art), grouped
    together under that sex's heading.
    """
    categories = ['off_art', 'on_effective_art', 'on_nonsuppressive_art']
    for sex in ('m', 'f'):
        print(f'===== Closed-form annual mortality rate, sex={sex} =====')
        for cat in categories:
            df = pd.DataFrame(grids[sex][cat], index=age_vals, columns=cd4_vals)
            df.index.name, df.columns.name = 'age', 'cd4'
            print(f'--- {cat} ---')
            print(df.round(4).to_string())
        print()


def print_multiplier_stages(pars=None):
    """
    Print the raw, hard-coded values that feed into on_art_rate_analytic/
    off_art_rate_analytic in isolation -- the shared off-ART CD4 table, the
    relative-mortality factors, and the age multipliers -- rather than just
    the final combined rate, plus the worst-case-multiplier arithmetic that
    confirms *why* on-ART mortality can never exceed off-ART mortality at
    the same CD4 (see module docstring).
    """
    pars = pars or pull_hiv_pars()
    print('--- Off-ART CD4 bins/rates (HIVPars.cd4_death_bins/cd4_death_rates), shared with on-ART ---')
    print(pd.DataFrame({'cd4_bin_upper': pars.cd4_death_bins, 'rate': pars.cd4_death_rates}).to_string(index=False))
    print(f'  (rel_death = {pars.rel_death} multiplies these, for both off-ART and on-ART)')
    print('--- Relative-mortality factors (fraction of off-ART rate retained on ART) ---')
    print(f'  rel_art_mortality_effective: {pars.rel_art_mortality_effective}')
    print(f'  rel_art_mortality_unsupp_m (non-suppressive ART, males): {pars.rel_art_mortality_unsupp_m}')
    print(f'  rel_art_mortality_unsupp_f (non-suppressive ART, females): {pars.rel_art_mortality_unsupp_f}')
    print(f'  rel_death_f (female, on ART, both adherence categories): {pars.rel_death_f}')
    print(f'  --> non-suppressive/effective mortality ratio: male={pars.rel_art_mortality_unsupp_m / pars.rel_art_mortality_effective:.2f}, '
          f'female={pars.rel_art_mortality_unsupp_f / pars.rel_art_mortality_effective:.2f} '
          f'(rel_death_f cancels out of this ratio, since it applies equally to both categories)')
    print('--- On-ART age_mult (HIVPars.art_death_age) ---')
    print(pd.DataFrame(pars.art_death_age, columns=['age_lo', 'age_hi', 'mult']).to_string(index=False))
    print()
    # Invariant check: on-ART rate = off_art_rate(cd4) * rel_art_mortality * age_mult
    # [* rel_death_f]. Since off_art_rate(cd4) is the SAME shared factor on both
    # sides, on-ART <= off-ART at every CD4 count iff
    # rel_art_mortality * (largest age multiplier) <= 1 -- checked here for the
    # worst case (non-suppressive males, since rel_art_mortality_unsupp_m is the
    # larger of the two sex-specific values, and rel_death_f < 1 only helps females).
    max_age_mult = max(mult for _, _, mult in pars.art_death_age)
    worst_case_factor = pars.rel_art_mortality_unsupp_m * max_age_mult
    print(f'Worst-case combined multiplier (non-suppressive, oldest age band, male): '
          f'{pars.rel_art_mortality_unsupp_m} * {max_age_mult} = {worst_case_factor:.4f}')
    if worst_case_factor <= 1:
        print('  --> <= 1: on-ART mortality can never exceed off-ART mortality at the same CD4, for any age/sex/adherence.')
    else:
        print('  --> > 1: the on-ART-never-exceeds-off-ART invariant is VIOLATED for this combination -- '
              'lower rel_art_mortality_unsupp_m or the age multipliers.')


def plot_analytic_heatmaps(age_vals, cd4_vals, grids, vmin=None, vmax=None):
    """ 2 (sex) x 3 (category) grid of CD4-by-age heatmaps of the closed-form annual mortality rate. """
    categories = ['off_art', 'on_effective_art', 'on_nonsuppressive_art']
    if vmin is None or vmax is None:
        all_vals = np.concatenate([grids[sex][cat].ravel() for sex in ('m', 'f') for cat in categories])
        vmin = all_vals.min() if vmin is None else vmin
        vmax = all_vals.max() if vmax is None else vmax
    norm = LogNorm(vmin=vmin, vmax=vmax)

    fig, axes = plt.subplots(2, 3, figsize=(13, 7), sharex=True, sharey=True)
    im = None
    for row, sex in enumerate(('m', 'f')):
        for col, cat in enumerate(categories):
            ax = axes[row, col]
            im = ax.pcolormesh(cd4_vals, age_vals, grids[sex][cat], norm=norm, cmap='viridis', shading='auto')
            if row == 0:
                ax.set_title(cat)
            if col == 0:
                ax.set_ylabel(f'sex={sex}\nAge (years)')
            if row == 1:
                ax.set_xlabel('CD4 count')
    fig.suptitle('Closed-form (no simulation) annual mortality rate: age x CD4, by sex and ART status')
    # NB: fig.colorbar(im, ax=axes) followed by tight_layout() causes the colorbar
    # to overlap the rightmost column -- tight_layout() re-flows the grid without
    # knowing space was already reserved for the colorbar. Instead, explicitly
    # reserve the right margin via subplots_adjust, then place the colorbar in its
    # own axis within that margin.
    fig.subplots_adjust(left=0.08, right=0.88, top=0.90, bottom=0.08, wspace=0.15, hspace=0.15)
    cbar_ax = fig.add_axes([0.90, 0.15, 0.02, 0.7])  # [left, bottom, width, height], figure-fraction
    fig.colorbar(im, cax=cbar_ax, label='Annual mortality rate')
    plt.close(fig)  # avoid double-display in notebooks: inline backend auto-shows open figures, AND returning one triggers rich-display
    return fig


if __name__ == '__main__':

    # Specify your own age/CD4 brackets, and which sex to plot (for the age x
    # CD4 faceted view), here
    age_bins = [(0, 25), (25, 35), (35, 45), (45, 200)]
    cd4_bins = [(0, 200), (200, 350), (350, 500), (500, 100_000)]
    SEX = 'f'  # 'm' or 'f'

    age_labels = [_bin_label(lo, hi, AgeCD4ARTMortalityTracker.AGE_OPEN_ENDED_AT) for lo, hi in age_bins]
    cd4_labels = [_bin_label(lo, hi, AgeCD4ARTMortalityTracker.CD4_OPEN_ENDED_AT) for lo, hi in cd4_bins]

    # ==================================================================
    # PART 1: simulation-based. One run, two derived views.
    # ==================================================================
    tracker, yearvec = run_sim(age_bins, cd4_bins, n_agents=20_000, start=2000, stop=2030)

    # 1a. Yearly mortality rate, pooled across age/CD4/sex, by ART status
    yearly = yearly_mortality_by_category(tracker, yearvec, age_bins, cd4_bins)
    print('Annual HIV death rate per year, by ART status:')
    print(yearly.pivot_table(index='year', columns='category', values='annual_death_rate').to_string(float_format='%.4f'))
    plot_yearly_mortality(yearly)

    # 1b. Age x CD4 faceted mortality rate, for SEX, aggregated over the full sim
    print()
    summary = summarize_mortality_by_status(tracker, age_bins, cd4_bins, sex=SEX)
    print(f'Annual HIV death rate by ART status, age x CD4, sex={SEX}:')
    print(summary.pivot_table(index=['age_band', 'cd4_band'], columns='category',
                               values='annual_death_rate').to_string(float_format='%.4f'))
    print()
    check_off_art_higher(summary)
    plot_mortality_by_status(summary, age_labels, cd4_labels, sex=SEX)

    # ==================================================================
    # PART 2: closed-form (no simulation) -- exact same math as
    # make_p_hiv_death/get_art_mortality_hazard, evaluated directly on an
    # age x CD4 grid, pulling live parameter values from sti.HIVPars().
    # ==================================================================
    print()
    print_multiplier_stages()

    analytic_age_vals = np.array([20, 30, 40, 50, 60])
    analytic_cd4_vals = np.array([0, 100, 200, 300, 400, 500, 600, 700, 800])
    analytic_grids = compute_analytic_grids(analytic_age_vals, analytic_cd4_vals)
    print_analytic_tables(analytic_age_vals, analytic_cd4_vals, analytic_grids)

    fine_age_vals = np.arange(15, 66, 2)
    fine_cd4_vals = np.arange(0, 810, 20)
    fine_grids = compute_analytic_grids(fine_age_vals, fine_cd4_vals)
    plot_analytic_heatmaps(fine_age_vals, fine_cd4_vals, fine_grids)

    plt.show()
