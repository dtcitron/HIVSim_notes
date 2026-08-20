"""
Test-drive the redesigned PrEP implementation (branch feat/prep_parameters).

Two parts:
  1. Three-scenario comparison (no PrEP / oral PrEP / long-acting PrEP),
     plotting prevalence and new infections over time.
  2. A small controlled sim verifying PrEP courses expire exactly on
     schedule (pure duration mechanics, no epidemic dynamics).

How the new PrEP implementation works
--------------------------------------
State lives on HIV, not on the Prep intervention, mirroring ART's on_art/diagnosed pattern -- this lets multiple PrEP products share one "currently protected" status instead of each tracking its own.
  - hiv.on_prep / hiv.prep_naive / hiv.prep_discontinued: cascade states, same convention as ART's never_art/on_art/post_art.
  - hiv.prep_eff: the REALIZED efficacy for a person's current course (base efficacy x adherence, computed once at enrollment), so update_transmission() doesn't need to know which product someone is on.
  - hiv.prep_source: which Prep instance enrolled this person, used for that instance's own coverage-target bookkeeping.
  - hiv.ti_prep_start / hiv.ti_prep_stop: when the current course started and when it expires.

There are no named "varieties" -- any product is just a sti.Prep instance with its own prep_eff/prep_dur/prep_adh. Adherence is a continuous multiplier, not a bernoulli split: realized efficacy = prep_eff * prep_adh, baked in once at enrollment. Long-acting products just leave prep_adh=1.0.

HIV.start_prep()/stop_prep() are the only things that touch these states. start_prep() is idempotent -- enrolling someone already on ANY product is a no-op, so products can't override each other. Coverage is a prevalence (stock) target like ART/VMMC, not a per-step hazard: scalar, time-varying dict, or age/sex-stratified DataFrame (n_prep/p_prep columns), and each Prep instance tops itself up to target every step.

Expiry isn't automatic -- each Prep instance calls hiv.stop_prep() on its own enrollees (via prep_source) once ti_prep_stop has passed (see Prep.step()). Anything enrolling people via hiv.start_prep() directly, like the one-off intervention in part 2 below, must replicate that check itself or courses never end.

How to use it
-------------
Construct one sti.Prep(...) instance per product, each with its own
prep_eff/prep_dur/prep_adh/coverage/eligibility/name, and add them all to
the sim's interventions list alongside the usual testing/ART interventions.
Multiple instances naturally compose: a person eligible for two products is
enrolled in whichever tops up its target first, and is automatically
excluded from the other until their current course ends. See the
scenario_oral/scenario_la examples below.
"""

import hivsim
import stisim as sti
import starsim as ss
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

n_agents = 10_000
dur = 15


def make_sim(interventions, label):
    sim = hivsim.demo(
        'simple', run=False, plot=False, n_agents=n_agents, dur=dur,
        networks=[sti.StructuredSexual(), ss.MaternalNet(), ss.BreastfeedingNet()],
    )
    hiv_test = sti.HIVTest(name='hiv_test', test_prob_data=0.3)
    art = sti.ART(coverage=0.6)
    sim.pars.interventions = [hiv_test, art] + list(interventions)
    sim.label = label
    return sim


# 1. No PrEP (baseline)
scenario_none = []

# 2. Oral PrEP: 10% of women 15-35, 3-month courses, 85% efficacy, 75% adherence.
#    Flat coverage + an eligibility filter (Prep defaults to FSW-only otherwise).
prep_oral = sti.Prep(
    prep_eff=0.85, prep_dur=ss.months(3), prep_adh=0.75,
    coverage=0.1,
    eligibility=lambda sim: sim.people.female & (sim.people.age >= 15) & (sim.people.age < 35),
    name='prep_oral',
)
scenario_oral = [prep_oral]

# 3. Long-acting PrEP: 500 men + 1000 women aged 15-35, 12-month courses,
#    100% efficacy/adherence. Numbers differ by sex, so this needs an
#    age/sex-stratified absolute-count (n_prep) coverage target rather than
#    a flat fraction.
la_coverage = pd.DataFrame({
    'Year':   [2000, 2000],
    'AgeBin': ['[15,35)', '[15,35)'],
    'Gender': ['m', 'f'],
    'n_prep': [500, 1000],
})
prep_la = sti.Prep(
    prep_eff=1.0, prep_dur=ss.months(12), prep_adh=1.0,
    coverage=la_coverage,
    eligibility=lambda sim: sim.people.alive,  # age/sex restriction comes from la_coverage's strata
    name='prep_la',
)
scenario_la = [prep_la]

scenarios = {
    'No PrEP': scenario_none,
    'Oral PrEP (10% women 15-35)': scenario_oral,
    'Long-acting PrEP (500M/1000F, 15-35)': scenario_la,
}

sims = [make_sim(intvs, label) for label, intvs in scenarios.items()]
for sim in sims:
    sim.run(verbose=0)

def bin_timeseries(timevec, values, bin_size=6, agg='mean'):
    """
    Downsample a time series into non-overlapping chunks of bin_size steps.
    Use agg='mean' for rates (e.g. prevalence) and agg='sum' for per-step
    counts (e.g. new_infections), so a count bin means "total over the
    window" rather than being diluted by averaging. Drops any leftover
    steps that don't fill a full bin.
    """
    values = np.asarray(values)
    n = len(values) - (len(values) % bin_size)
    values = values[:n].reshape(-1, bin_size)
    timevec = np.asarray(timevec)[:n].reshape(-1, bin_size)
    binned_time = timevec.mean(axis=1)
    binned_values = values.mean(axis=1) if agg == 'mean' else values.sum(axis=1)
    return binned_time, binned_values


def total_new_prep(sim):
    """Sum new_prep across every Prep instance in the sim (all-zero if none,
    e.g. the no-PrEP scenario) -- a scenario can have more than one product."""
    total = np.zeros(len(sim.t.yearvec))
    for intv in sim.interventions():
        if isinstance(intv, sti.Prep):
            total += np.asarray(intv.results.new_prep)
    return total


fig, axes = plt.subplots(1, 3, figsize=(18, 5))
for sim in sims:
    hiv = sim.results.hiv
    t, prev = bin_timeseries(sim.t.yearvec, hiv.prevalence, bin_size=6, agg='mean')
    _, new_inf = bin_timeseries(sim.t.yearvec, hiv.new_infections, bin_size=6, agg='sum')
    t_annual, new_prep = bin_timeseries(sim.t.yearvec, total_new_prep(sim), bin_size=12, agg='sum')
    axes[0].plot(t, prev, label=sim.label)
    axes[1].plot(t, new_inf, label=sim.label)
    axes[2].plot(t_annual, new_prep, label=sim.label)

axes[0].set_title('HIV prevalence')
axes[0].set_xlabel('Year')
axes[0].set_ylabel('Prevalence')
axes[1].set_title('New infections')
axes[1].set_xlabel('Year')
axes[1].set_ylabel('New infections per 6 months')
axes[2].set_title('PrEP distributions')
axes[2].set_xlabel('Year')
axes[2].set_ylabel('New PrEP starts per year')
for ax in axes:
    ax.legend(fontsize=8)
fig.suptitle('PrEP scenario comparison')
fig.tight_layout()
plt.savefig('prep_test_scenarios.png', dpi=150)
print('Saved figure to prep_test_scenarios.png')

for sim in sims:
    hiv = sim.results.hiv
    print(f'{sim.label:40s} final prevalence={hiv.prevalence[-1]:.4f}  '
          f'cum_infections={int(hiv.cum_infections[-1])}')


# %% Duration check: 100 people on PrEP for 12 months, verify expiry timing
print('\n--- PrEP duration check ---')

class OneTimePrepEnrollment(ss.Intervention):
    """Enrolls everyone once at t=0, then only ever expires courses (never
    tops back up) -- unlike sti.Prep's continuous stock-target coverage,
    which would immediately re-enroll people the same step they expire,
    masking the drop in n_on_prep. Expiry isn't automatic on HIV's side --
    each Prep-like instance is responsible for calling stop_prep() on its
    own enrollees every step (mirrors Prep.step())."""
    def __init__(self, eff, dur, adh=1.0):
        super().__init__()
        self.eff, self.dur, self.adh = eff, dur, adh
        self._source_id = -1
        self._enrolled = False

    def step(self):
        hiv = self.sim.diseases.hiv
        if not self._enrolled:
            hiv.start_prep(self.sim.people.alive.uids, eff=self.eff, dur=self.dur,
                            source_id=self._source_id, adh=self.adh)
            self._enrolled = True
        else:
            expiring = hiv.on_prep & (hiv.prep_source == self._source_id) & (hiv.ti_prep_stop <= self.ti)
            if expiring.any():
                hiv.stop_prep(expiring.uids)
        return


sim_dur = hivsim.demo(
    'simple', run=False, plot=False, n_agents=100, dur=3,
    networks=[sti.StructuredSexual(), ss.MaternalNet(), ss.BreastfeedingNet()],
)
sim_dur.pars.interventions = [OneTimePrepEnrollment(eff=0.9, dur=ss.months(12))]
sim_dur.run(verbose=0)

timevec = sim_dur.t.yearvec  # plain float years -- hiv.timevec is a DateArray, awkward for arithmetic
n_on_prep = np.asarray(sim_dur.results.hiv.n_on_prep)
enrolled_idx = 1  # first step after the one-time enrollment fires

fig2, ax = plt.subplots(figsize=(8, 4))
ax.plot(timevec, n_on_prep)
ax.axvline(timevec[enrolled_idx] + 1.0, color='red', linestyle='--', label='expected expiry (12 months)')
ax.set_xlabel('Year')
ax.set_ylabel('Number on PrEP')
ax.set_title('PrEP duration check: 100 agents, 12-month course, single enrollment')
ax.legend()
fig2.tight_layout()
plt.savefig('prep_test_duration_check.png', dpi=150)
print('Saved figure to prep_test_duration_check.png')

# Look for the bulk expiry (everyone shares the same ti_prep_stop) after
# enrollment, not the pre-enrollment zero baseline or incidental single-agent
# dips from background mortality along the way.
after = n_on_prep[enrolled_idx:]
expiry_idx = enrolled_idx + np.argmax(after < n_on_prep[enrolled_idx] / 2)
expiry_year = timevec[expiry_idx] - timevec[enrolled_idx]
print(f'Enrolled: {int(n_on_prep[enrolled_idx])} agents at year {timevec[enrolled_idx]:.3f}')
print(f'Dropped to {int(n_on_prep[expiry_idx])} agents at year {timevec[expiry_idx]:.3f} '
      f'({expiry_year:.3f} years after enrollment; expected 1.000)')
assert abs(expiry_year - 1.0) < sim_dur.t.dt_year * 1.5, 'PrEP did not expire on schedule'
print('PASS: PrEP expired on schedule.')
