"""
ART retention-duration test script (art.md, step 5.2), stratified by
user-specified age brackets.

Estimates the EMPIRICAL distribution of time spent on ART per treatment
"spell" (ART initiation -> discontinuation, death, or end of sim), broken
out by the agent's age at the START of that spell (age at ART initiation).

This is distinct from just reading off `HIVPars.dur_on_art` (the lognormal
draw at initiation, mean 3yr / std 1.5yr, see hiv.py) because the REALIZED
time on ART also depends on:
  1. Death while on ART, which truncates a spell early.
  2. `ART`'s coverage-correction logic, which can pull people off ART ahead
     of their drawn duration if the on-ART population exceeds a coverage
     target (see art_implementation_notes.md section 3b) -- or hold people
     on longer than drawn if the sim never reaches capacity.
  3. Right-censoring: agents still on ART when the sim ends have an
     unknown, only partially-observed spell length.

Produces, both overall and per age bracket:
  - Summary stats (n completed vs. censored spells, mean/median completed
    duration, breakdown by how the spell ended, Kaplan-Meier restricted
    mean).
  - A histogram of completed-spell durations, one panel per age bracket.
  - An overlaid Kaplan-Meier "fraction still on ART" curve vs. years since
    initiation, one line per age bracket plus an overall line, which
    correctly accounts for censored spells (naive mean/median of completed
    spells alone is biased short, since it ignores that censored people have
    already survived at least that long).
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


class ARTDurationTracker(ss.Analyzer):
    """
    Record start/stop calendar-years (and age at start) for every ART
    "spell" -- a contiguous period an agent spends on_art (effective or
    non-suppressive), ending in voluntary/program discontinuation, death, or
    (if ongoing at sim end) right-censoring.
    """
    def __init__(self, age_bins):
        super().__init__()
        self.age_bins = age_bins
        self._active = {}  # uid -> (spell start year, age at spell start)
        self.spells = []    # list of dicts: uid, start, stop, duration, end_reason, age_at_start, age_band

    def step(self):
        hiv = self.sim.diseases.hiv
        now = self.now.years  # numeric calendar year, e.g. 2013.5
        on_art_now = set(int(u) for u in hiv.on_art.uids)

        # New spells: agents who (re-)started ART this step
        starting = (hiv.ti_art == self.ti) & hiv.on_art
        ages = self.sim.people.age
        for uid in starting.uids:
            if int(uid) not in self._active:
                age_at_start = float(ages[ss.uids([uid])][0])
                self._active[int(uid)] = (now, age_at_start)

        # Ended spells: previously-tracked agents no longer on ART. Covers
        # both stop_art() discontinuation and death (HIV.step_die clears
        # on_art for agents who die, whatever the cause of death).
        ended_uids = [uid for uid in self._active if uid not in on_art_now]
        alive = self.sim.people.alive
        for uid in ended_uids:
            start, age_at_start = self._active.pop(uid)
            is_alive = bool(alive[ss.uids([uid])][0])
            reason = 'stopped' if is_alive else 'death'
            self.spells.append(dict(
                uid=uid, start=start, stop=now, duration=now - start, censored=False,
                end_reason=reason, age_at_start=age_at_start, age_band=_age_band(age_at_start, self.age_bins),
            ))


def kaplan_meier(durations, censored):
    """
    Manual Kaplan-Meier estimator of the survival function (fraction of
    spells still ongoing) as a function of time since ART initiation.

    Args:
        durations (array-like): observed spell length for each spell
            (time-to-event for completed spells, time-to-censoring for
            ongoing ones)
        censored (array-like of bool): True if the spell was still ongoing
            at the end of the observation window (right-censored)

    Returns:
        (times, survival): step-function coordinates suitable for plotting
    """
    df = pd.DataFrame({'t': durations, 'censored': censored}).sort_values('t')
    times, surv = [0.0], [1.0]
    survival = 1.0
    remaining = df
    for t in sorted(df['t'].unique()):
        n_at_risk = len(remaining)
        if n_at_risk == 0:
            break
        n_events = ((remaining['t'] == t) & (~remaining['censored'])).sum()
        survival *= (1 - n_events / n_at_risk)
        times.append(t)
        surv.append(survival)
        remaining = remaining[remaining['t'] > t]
    return np.array(times), np.array(surv)


def run_duration_test(age_bins, n_agents=10_000, start=2000, stop=2040, coverage_target=0.6):
    """
    Run a sim with a ramping-then-held ART coverage target and return a
    DataFrame of every observed ART spell (completed or censored), tagged
    with the agent's age band at spell start.
    """
    hiv_test = sti.HIVTest(test_prob_data=0.3, start=start)
    # Ramping coverage target with a long hold, similar to real usage
    # (hivsim_examples/zimbabwe/sim.py) -- this activates ART's
    # coverage-correction logic, so realized spell lengths reflect BOTH the
    # random dur_on_art draw AND program-capacity-driven early stops.
    ramp_end = min(start + 10, stop)
    art = sti.ART(coverage={'year': [start, ramp_end, stop], 'value': [0, coverage_target, coverage_target]})

    sim = sti.Sim(
        diseases=sti.HIV(init_prev=ss.bernoulli(p=0.2)),
        n_agents=n_agents, start=start, stop=stop,
        interventions=[hiv_test, art],
        analyzers=[ARTDurationTracker(age_bins)],
    )
    sim.run(verbose=0)

    # NB: Sim deep-copies analyzers on init, so the tracker actually
    # populated during the run is `sim.analyzers.artdurationtracker`, not
    # the object constructed above.
    tracker = sim.analyzers.artdurationtracker

    # Anyone still on ART at sim end is a right-censored spell
    final_year = sim.t.yearvec[-1]
    censored_spells = [
        dict(uid=uid, start=s, stop=final_year, duration=final_year - s, censored=True,
             end_reason='ongoing', age_at_start=a, age_band=_age_band(a, age_bins))
        for uid, (s, a) in tracker._active.items()
    ]
    return pd.DataFrame(tracker.spells + censored_spells)


class AgeDependentARTDropout(ss.Intervention):
    """
    Standalone custom intervention: overrides HIV.ti_stop_art for
    newly-initiated ART agents with an age-dependent EXPONENTIAL retention
    duration -- i.e. a constant annual dropout HAZARD per age band -- instead
    of HIV's built-in single lognormal `dur_on_art`, which has no age
    dependence at all (see art_implementation_notes.md section 8, "Option
    A/B", for why exponential is the distribution that exactly reproduces a
    flat annual dropout rate: it's the constant-hazard distribution).

    Must be placed AFTER `ART` in `interventions=[...]`: it runs later in
    the same timestep, after ART.step()/HIV.start_art() has already set
    ti_art/ti_stop_art for anyone starting ART this step, and simply
    overwrites ti_stop_art for those same agents with an age-conditioned
    draw. This doesn't disturb ART's handling of agents *stopping* ART this
    step, since that check (against durations set on a prior timestep)
    already happened earlier in the same ART.step() call.

    Agents outside every band in `age_dropout` (e.g. under-15 agents
    diagnosed at sim initialization, bypassing the usual age-restricted
    HIVTest eligibility) are left with whatever ti_stop_art HIV's own
    start_art() already set -- i.e. they silently fall back to the
    unmodified lognormal dur_on_art.

    Args:
        age_dropout: list of (age_lo, age_hi, mean_years) tuples. mean_years
            is the mean retention duration under a constant-hazard
            (exponential) dropout model, i.e. mean_years = 1 / annual_dropout_rate.
            Example: [(15, 25, 20/3), (25, 200, 20)] reproduces a 15%/yr
            dropout hazard under 25 and a 5%/yr dropout hazard 25+.
    """
    def __init__(self, age_dropout, **kwargs):
        super().__init__()
        self.age_dropout = age_dropout
        self.define_pars(**{
            f'dur_{i}': ss.expon(scale=ss.years(mean_years))
            for i, (_, _, mean_years) in enumerate(age_dropout)
        })
        self.update_pars(**kwargs)

    def step(self):
        hiv = self.sim.diseases.hiv
        ti = self.ti
        just_started = ((hiv.ti_art == ti) & hiv.on_art).uids
        if len(just_started) == 0:
            return
        ages = self.sim.people.age[just_started]
        for i, (lo, hi, _) in enumerate(self.age_dropout):
            in_band = just_started[(ages >= lo) & (ages < hi)]
            if len(in_band) == 0:
                continue
            dur = self.pars[f'dur_{i}'].rvs(in_band)
            hiv.ti_stop_art[in_band] = ti + dur.astype(int)


def run_duration_test_age_dropout(age_dropout, min_age=15, n_agents=10_000, start=2000, stop=2040, coverage_target=0.6):
    """
    Run a sim where ART is restricted to agents aged >= min_age, and
    retention duration is drawn from an age-dependent exponential dropout
    hazard (AgeDependentARTDropout) instead of HIV's built-in single
    lognormal dur_on_art. Returns a spells DataFrame in the same shape as
    run_duration_test(), tagged by age band at initiation using the SAME
    age bands as `age_dropout`, so the resulting plots line up with the
    hazard groups actually used.
    """
    age_bins = [(lo, hi) for lo, hi, _ in age_dropout]

    hiv_test = sti.HIVTest(
        test_prob_data=0.3, start=start,
        eligibility=lambda sim: ~sim.diseases.hiv.diagnosed & (sim.people.age >= min_age),
    )
    ramp_end = min(start + 10, stop)
    art = sti.ART(coverage={'year': [start, ramp_end, stop], 'value': [0, coverage_target, coverage_target]})
    dropout = AgeDependentARTDropout(age_dropout)

    sim = sti.Sim(
        diseases=sti.HIV(init_prev=ss.bernoulli(p=0.2)),
        n_agents=n_agents, start=start, stop=stop,
        interventions=[hiv_test, art, dropout],  # dropout must come AFTER art
        analyzers=[ARTDurationTracker(age_bins)],
    )
    sim.run(verbose=0)
    tracker = sim.analyzers.artdurationtracker

    final_year = sim.t.yearvec[-1]
    censored_spells = [
        dict(uid=uid, start=s, stop=final_year, duration=final_year - s, censored=True,
             end_reason='ongoing', age_at_start=a, age_band=_age_band(a, age_bins))
        for uid, (s, a) in tracker._active.items()
    ]
    return pd.DataFrame(tracker.spells + censored_spells)


def print_duration_summary(df, age_labels):
    """ Print n spells, completed-spell mean/median, and KM restricted mean, overall and per age band. """
    for label in list(age_labels) + ['overall']:
        sub = df if label == 'overall' else df[df['age_band'] == label]
        if len(sub) == 0:
            print(f'{label}: no spells observed')
            continue
        completed = sub[~sub['censored']]
        times, surv = kaplan_meier(sub['duration'].values, sub['censored'].values)
        rmst = np.sum(np.diff(times) * (surv[:-1] + surv[1:]) / 2)
        mean_c = completed['duration'].mean() if len(completed) else np.nan
        median_c = completed['duration'].median() if len(completed) else np.nan
        print(f'{label}: n={len(sub)} ({len(completed)} completed, {len(sub) - len(completed)} censored), '
              f'completed mean/median={mean_c:.2f}/{median_c:.2f} yr, '
              f'KM restricted mean (lower bound)={rmst:.2f} yr')


def plot_duration_histograms(df, age_labels):
    """ Grid of histograms of completed-spell durations, one panel per age bracket. """
    completed = df[~df['censored']]
    fig, axes = plt.subplots(1, len(age_labels), figsize=(4 * len(age_labels), 4), sharey=True)
    if len(age_labels) == 1:
        axes = [axes]
    for ax, label in zip(axes, age_labels):
        sub = completed[completed['age_band'] == label]
        ax.hist(sub['duration'], bins=20)
        ax.set_title(label)
        ax.set_xlabel('Time on ART (years)')
    axes[0].set_ylabel('Number of completed spells')
    fig.suptitle('Distribution of completed ART-spell durations, by age at initiation')
    fig.tight_layout()
    return fig


def plot_km_by_age(df, age_labels):
    """ Overlaid Kaplan-Meier retention curves, one line per age bracket plus an overall line. """
    fig, ax = plt.subplots()
    for label in age_labels:
        sub = df[df['age_band'] == label]
        if len(sub) == 0:
            continue
        times, surv = kaplan_meier(sub['duration'].values, sub['censored'].values)
        ax.step(times, surv, where='post', label=label)
    times, surv = kaplan_meier(df['duration'].values, df['censored'].values)
    ax.step(times, surv, where='post', label='overall', color='k', linestyle='--')
    ax.set_xlabel('Years since ART initiation')
    ax.set_ylabel('Fraction still on ART')
    ax.set_title('Kaplan-Meier ART retention curve, by age at initiation')
    ax.set_ylim(0, 1.05)
    ax.legend()
    return fig


if __name__ == '__main__':

    # ------------------------------------------------------------------
    # Check 1: baseline duration distribution using HIV's built-in
    # dur_on_art (lognormal, no age dependence), stratified by age at
    # initiation purely for reporting -- age has no causal effect here.
    # ------------------------------------------------------------------
    age_bins = [(0, 25), (25, 35), (35, 45), (45, 200)]
    age_labels = [_bin_label(lo, hi) for lo, hi in age_bins]

    df = run_duration_test(age_bins, n_agents=10_000, start=2000, stop=2040)

    print('=== Check 1: baseline (HIV.dur_on_art, no age dependence) ===')
    print_duration_summary(df, age_labels)
    plot_duration_histograms(df, age_labels)
    plot_km_by_age(df, age_labels)

    # ------------------------------------------------------------------
    # Check 2: age-dependent EXPONENTIAL dropout hazard via the standalone
    # AgeDependentARTDropout intervention -- 15%/yr (mean 20/3 yr) for ages
    # 15-25, 5%/yr (mean 20 yr) for ages 25+. ART is restricted to agents
    # 15+ (via HIVTest's eligibility). This should show visibly shorter
    # retention and faster KM decay for the 15-25 band than Check 1.
    # ------------------------------------------------------------------
    age_dropout = [(15, 25, 20 / 3), (25, 200, 20)]  # (age_lo, age_hi, mean retention years = 1/annual_dropout_rate)
    age_dropout_labels = [_bin_label(lo, hi) for lo, hi, _ in age_dropout]

    df_custom = run_duration_test_age_dropout(age_dropout, min_age=15, n_agents=10_000, start=2000, stop=2040)

    print()
    print('=== Check 2: age-dependent exponential dropout hazard (AgeDependentARTDropout) ===')
    print_duration_summary(df_custom, age_dropout_labels)
    plot_duration_histograms(df_custom, age_dropout_labels)
    plot_km_by_age(df_custom, age_dropout_labels)

    plt.show()
