"""
ART mortality-table vs. simple-mortality comparison (art.md step 4.5).

Compares STIsim's two ART-mortality modes:
  - `use_art_mortality_table=False` (default/"old"): zero HIV mortality
    while on ART, regardless of age, sex, or adherence.
  - `use_art_mortality_table=True` ("new"): ARTMortalityTable-style hazard,
    stratified by age, sex, ART adherence category (effective vs.
    non-suppressive), duration since ART initiation, and current CD4 (see
    HIV.get_art_mortality_hazard in stisim/diseases/hiv.py).

Produces, for a user-specified set of age brackets and calendar-year date
range:
  - A printed table of annual on-ART death rate (deaths / person-years),
    broken out by (version, age bracket, sex, adherence category).
  - A grid of bar-chart plots (rows = sex, columns = adherence category)
    comparing old vs. new death rate per age bracket.

For the ART retention-*duration* analysis (how long people stay on ART, as
opposed to how likely they are to die while on it), see art_duration_test.py.
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import starsim as ss
import stisim as sti


def _bin_label(lo, hi, open_ended_at=200):
    return f'{lo}-{hi}' if hi < open_ended_at else f'{lo}+'


def _age_band(age, age_bins):
    for lo, hi in age_bins:
        if lo <= age < hi:
            return _bin_label(lo, hi)
    return None


class AgeSexAdherenceMortalityTracker(ss.Analyzer):
    """
    Tracks deaths and person-years, stratified by (user-specified) age
    bracket, sex, and treatment status: ART adherence category (effective
    vs. non-suppressive) for people on ART, plus a separate 'off_art'
    category (infected but not currently on ART, i.e. ART-naive or
    discontinued) for comparison.

    Deaths are attributed using each agent's (age bracket, sex, status) as
    of the END of the PREVIOUS timestep, cached across steps -- the same
    pattern ARTStatusTracker uses in art_transmission_test.py, and for the
    same reason: HIV.step_die() clears on_effective_art/on_nonsuppressive_art
    for agents who die this step, so the live state can't be read after
    the fact.
    """
    ADHERENCE = ['on_effective_art', 'on_nonsuppressive_art']  # on-ART categories only (used for old-vs-new comparisons)
    OFF_ART = 'off_art'  # infected & not on ART (art_naive or art_discontinued); mortality here is identical old vs. new

    def __init__(self, age_bins):
        super().__init__()
        self.age_bins = age_bins
        self.age_labels = [_bin_label(lo, hi) for lo, hi in age_bins]
        self._cache = {}  # uid -> (age_label, sex, status), as of end of previous timestep

    def init_results(self):
        super().init_results()
        results = []
        for age_label in self.age_labels:
            for sex in ('m', 'f'):
                for cat in self.ADHERENCE + [self.OFF_ART]:
                    key = f'{age_label}_{sex}_{cat}'
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
            label, sex, cat = cached
            self.results[f'deaths_{label}_{sex}_{cat}'][ti] += 1

        # Person-time this step, by CURRENT age band/sex/status
        cache = {}
        for cat, uids in [*((cat, hiv[cat].uids) for cat in self.ADHERENCE),
                           (self.OFF_ART, (hiv.infected & ~hiv.on_art).uids)]:
            for uid in uids:
                u = ss.uids([uid])
                label = _age_band(ppl.age[u][0], self.age_bins)
                if label is None:
                    continue
                sex = 'm' if bool(ppl.male[u][0]) else 'f'
                self.results[f'person_years_{label}_{sex}_{cat}'][ti] += dt_years
                cache[int(uid)] = (label, sex, cat)
        self._cache = cache


def run_mortality_comparison(age_bins, n_agents=10_000, start=2000, stop=2020, coverage_target=0.6):
    """
    Run identical scenarios with the old (zero on-ART mortality) and new
    (ARTMortalityTable-style) behavior, and return each run's tracker plus
    the shared calendar-year time vector.
    """
    trackers = {}
    sim = None
    for version, use_table in [('old', False), ('new', True)]:
        hiv_test = sti.HIVTest(test_prob_data=0.3, start=start)
        ramp_end = min(start + 4, stop)
        art = sti.ART(coverage={'year': [start, ramp_end, stop], 'value': [0, coverage_target, coverage_target]}, p_effective_art=0.5)
        sim = sti.Sim(
            diseases=sti.HIV(init_prev=ss.bernoulli(p=0.2), use_art_mortality_table=use_table),
            n_agents=n_agents, start=start, stop=stop,
            interventions=[hiv_test, art],
            analyzers=[AgeSexAdherenceMortalityTracker(age_bins)],
        )
        sim.run(verbose=0)
        trackers[version] = sim.analyzers.agesexadherencemortalitytracker
    return trackers, sim.t.yearvec  # yearvec is identical for both runs (same start/stop/dt)


def summarize_mortality_comparison(trackers, yearvec, age_bins, date_range=None):
    """
    Build a tidy DataFrame of annual death rate (deaths / person-years), one
    row per (version, age band, sex, adherence category), aggregated only
    over timesteps whose calendar year falls in date_range.

    Args:
        date_range (tuple, optional): (start_year, end_year), inclusive of
            start_year and exclusive of end_year. Defaults to the full run.
    """
    if date_range is None:
        mask = np.ones(len(yearvec), dtype=bool)
    else:
        lo, hi = date_range
        mask = (yearvec >= lo) & (yearvec < hi)

    age_labels = [_bin_label(lo, hi) for lo, hi in age_bins]
    rows = []
    for version, tracker in trackers.items():
        for age_label in age_labels:
            for sex in ('m', 'f'):
                for cat in AgeSexAdherenceMortalityTracker.ADHERENCE:
                    key = f'{age_label}_{sex}_{cat}'
                    deaths = tracker.results[f'deaths_{key}'].values[mask].sum()
                    person_years = tracker.results[f'person_years_{key}'].values[mask].sum()
                    rate = deaths / person_years if person_years > 0 else np.nan
                    rows.append(dict(version=version, age_band=age_label, sex=sex, adherence=cat,
                                      deaths=deaths, person_years=person_years, annual_death_rate=rate))
    return pd.DataFrame(rows)


def yearly_mortality_by_adherence(trackers, yearvec, age_bins):
    """
    Total mortality per calendar year, pooled across all age brackets and
    both sexes, broken out into four categories: 'any' (on ART of any kind,
    i.e. effective + non-suppressive combined), 'effective', 'nonsuppressive',
    and 'off_art' (infected but not currently on ART). One row per
    (version, year, category).

    Note 'off_art' mortality is computed by the same code path regardless of
    `use_art_mortality_table`, so old vs. new should only differ there by
    simulation noise -- it's included for a full picture, not because the
    ARTMortalityTable toggle is expected to change it.
    """
    all_cats = AgeSexAdherenceMortalityTracker.ADHERENCE + [AgeSexAdherenceMortalityTracker.OFF_ART]
    age_labels = [_bin_label(lo, hi) for lo, hi in age_bins]
    years = np.floor(yearvec).astype(int)
    unique_years = sorted(set(years))

    rows = []
    for version, tracker in trackers.items():
        # Sum deaths/person-years across all age brackets and sexes, per timestep, per category
        deaths_by_cat = {cat: np.zeros(len(yearvec)) for cat in all_cats}
        py_by_cat = {cat: np.zeros(len(yearvec)) for cat in all_cats}
        for cat in all_cats:
            for age_label in age_labels:
                for sex in ('m', 'f'):
                    key = f'{age_label}_{sex}_{cat}'
                    deaths_by_cat[cat] += tracker.results[f'deaths_{key}'].values
                    py_by_cat[cat] += tracker.results[f'person_years_{key}'].values

        deaths_any = deaths_by_cat['on_effective_art'] + deaths_by_cat['on_nonsuppressive_art']
        py_any = py_by_cat['on_effective_art'] + py_by_cat['on_nonsuppressive_art']

        by_category = [
            ('any', deaths_any, py_any),
            ('effective', deaths_by_cat['on_effective_art'], py_by_cat['on_effective_art']),
            ('nonsuppressive', deaths_by_cat['on_nonsuppressive_art'], py_by_cat['on_nonsuppressive_art']),
            ('off_art', deaths_by_cat[AgeSexAdherenceMortalityTracker.OFF_ART], py_by_cat[AgeSexAdherenceMortalityTracker.OFF_ART]),
        ]
        for year in unique_years:
            mask = years == year
            for category, deaths_arr, py_arr in by_category:
                deaths = deaths_arr[mask].sum()
                person_years = py_arr[mask].sum()
                rate = deaths / person_years if person_years > 0 else np.nan
                rows.append(dict(version=version, year=year, category=category,
                                  deaths=deaths, person_years=person_years, annual_death_rate=rate))
    return pd.DataFrame(rows)


def plot_yearly_mortality(yearly):
    """ Total deaths per year, one panel per category, old vs. new overlaid. """
    categories = ['any', 'effective', 'nonsuppressive', 'off_art']
    fig, axes = plt.subplots(1, len(categories), figsize=(5 * len(categories), 4), sharey=True)
    for ax, category in zip(axes, categories):
        for version, color in [('old', 'tab:gray'), ('new', 'tab:blue')]:
            sub = yearly.query('version == @version and category == @category').sort_values('year')
            ax.plot(sub['year'], sub['deaths'], marker='o', label=version, color=color)
        ax.set_title(category)
        ax.set_xlabel('Year')
        ax.legend()
    axes[0].set_ylabel('Total deaths')
    fig.suptitle('Old vs. new: total deaths per year, by ART status')
    fig.tight_layout()
    return fig


def plot_mortality_comparison(comparison, age_labels):
    """
    Grid of bar-chart panels: rows = sex, columns = ART adherence category.
    Each panel shows old vs. new annual on-ART death rate per age bracket.
    """
    adherence_cats = AgeSexAdherenceMortalityTracker.ADHERENCE
    fig, axes = plt.subplots(2, len(adherence_cats), figsize=(5.5 * len(adherence_cats), 7), sharey=True)
    x = np.arange(len(age_labels))
    width = 0.35
    for row, sex in enumerate(('m', 'f')):
        for col, cat in enumerate(adherence_cats):
            ax = axes[row, col]
            for offset, version, color in [(-1, 'old', 'tab:gray'), (1, 'new', 'tab:blue')]:
                rates = [
                    comparison.query('version == @version and sex == @sex and adherence == @cat and age_band == @label')['annual_death_rate'].iloc[0]
                    for label in age_labels
                ]
                ax.bar(x + offset * width / 2, rates, width, label=version, color=color)
            ax.set_xticks(x)
            ax.set_xticklabels(age_labels)
            ax.set_title(f'{cat}, sex={sex}')
            if col == 0:
                ax.set_ylabel('Annual death rate while on ART')
            ax.legend()
    fig.suptitle('Old (simple) vs. new (ARTMortalityTable) on-ART mortality, by age and sex')
    fig.tight_layout()
    return fig


if __name__ == '__main__':

    # Specify your own age brackets and calendar-year window here
    age_bins = [(0, 25), (25, 35), (35, 45), (45, 200)]
    date_range = (2010, 2020)

    trackers, yearvec = run_mortality_comparison(age_bins, n_agents=10_000, start=2000, stop=2020)
    comparison = summarize_mortality_comparison(trackers, yearvec, age_bins, date_range=date_range)

    print(f'Old vs. new ARTMortality comparison, {date_range[0]}-{date_range[1]} (annual death rate while on ART):')
    print(comparison.pivot_table(index=['adherence', 'sex', 'age_band'], columns='version',
                                  values='annual_death_rate').to_string(float_format='%.4f'))

    age_labels = [_bin_label(lo, hi) for lo, hi in age_bins]
    plot_mortality_comparison(comparison, age_labels)

    # Total (age/sex-pooled) on-ART mortality per calendar year, for people
    # on ART of any kind, on effective ART only, and on non-suppressive ART only.
    yearly = yearly_mortality_by_adherence(trackers, yearvec, age_bins)
    print()
    print('Total on-ART deaths per year, by adherence category:')
    print(yearly.pivot_table(index='year', columns=['category', 'version'], values='deaths').to_string(float_format='%.0f'))
    plot_yearly_mortality(yearly)

    plt.show()
