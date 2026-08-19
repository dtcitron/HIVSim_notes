"""
Tutorial: specifying VLS (viral load suppression) coverage on sti.ART

By default, ART assumes every newly-initiated agent achieves viral
suppression (effective ART). This demo shows how to use the new
`vls_coverage` parameter to specify what fraction of ART initiators
achieve suppression instead — as a flat rate, a time-varying rate, or
a rate stratified by age/sex within the eligible population.

`vls_coverage` accepts the same input formats as `coverage` (see
sti.ART's docstring), except values must always be proportions (0-1),
since VLS coverage is a probability, not a count of people.
"""

import hivsim
import stisim as sti
import starsim as ss
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

n_agents = 3000
dur = 20


class StratifyCascade(ss.Analyzer):
    """
    Tracks ART cascade counts (diagnosed, on ART, on effective ART) by sex.

    hiv natively stratifies prevalence/n_infected by sex, but not the
    downstream cascade steps — this analyzer fills that gap.
    """
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.name = 'stratify_cascade'
        return

    def init_results(self):
        super().init_results()
        self.define_results(
            ss.Result('n_diagnosed_m', dtype=int, auto_plot=False),
            ss.Result('n_diagnosed_f', dtype=int, auto_plot=False),
            ss.Result('n_on_art_m', dtype=int, auto_plot=False),
            ss.Result('n_on_art_f', dtype=int, auto_plot=False),
            ss.Result('n_on_effective_art_m', dtype=int, auto_plot=False),
            ss.Result('n_on_effective_art_f', dtype=int, auto_plot=False),
        )
        return

    def step(self):
        ti = self.ti
        hiv = self.sim.diseases.hiv
        ppl = self.sim.people
        self.results.n_diagnosed_m[ti] = (hiv.diagnosed & ppl.male).count()
        self.results.n_diagnosed_f[ti] = (hiv.diagnosed & ppl.female).count()
        self.results.n_on_art_m[ti] = (hiv.on_art & ppl.male).count()
        self.results.n_on_art_f[ti] = (hiv.on_art & ppl.female).count()
        self.results.n_on_effective_art_m[ti] = (hiv.on_effective_art & ppl.male).count()
        self.results.n_on_effective_art_f[ti] = (hiv.on_effective_art & ppl.female).count()
        return


def make_sim(art, label):
    sim = hivsim.demo(
        'simple', run=False, plot=False, n_agents=n_agents, dur=dur,
        networks=[sti.StructuredSexual(), ss.MaternalNet(), ss.BreastfeedingNet()],
        analyzers=[StratifyCascade()],
    )
    hiv_test = sti.HIVTest(name='hiv_test', test_prob_data=0.3)
    sim.pars.interventions = [hiv_test, art]
    sim.label = label
    return sim


# 1. Default behavior: 100% of ART initiators achieve viral suppression
art_default = sti.ART(coverage=0.8)

# 2. Flat VLS coverage: 60% of ART initiators are virally suppressed,
#    the remaining 40% are on ART but non-suppressive
art_flat = sti.ART(coverage=0.8, vls_coverage=0.6)

# 3. Time-varying VLS coverage: suppression rates improve over time,
#    e.g. as better regimens/adherence support roll out
art_time_varying = sti.ART(
    coverage=0.8,
    vls_coverage={'year': [2000, 2010, 2025], 'value': [0.4, 0.6, 0.85]},
)

# 4. Stratified VLS coverage: specify suppression fractions per age/sex
#    stratum within the eligible population. Any stratum not listed here
#    defaults to 100% suppression. Column names follow the same flexible
#    aliasing as `coverage` (Year/AgeBin/Gender, case-insensitive).
vls_strat = pd.DataFrame({
    'Year':   [2000, 2000, 2000, 2000],
    'AgeBin': ['[15,25)', '[15,25)', '[25,100)', '[25,100)'],
    'Gender': ['m',       'f',       'm',        'f'],
    'p_vls':  [0.35,      0.55,      0.65,       0.80],
})
art_stratified = sti.ART(coverage=0.8, vls_coverage=vls_strat)

sims = [
    make_sim(art_default, 'default (100% VLS)'),
    make_sim(art_flat, 'flat 60% VLS'),
    make_sim(art_time_varying, 'time-varying VLS'),
    make_sim(art_stratified, 'age/sex-stratified VLS'),
]

# Run each scenario in series (no ss.parallel/multiprocessing)
for sim in sims:
    sim.run(verbose=0)

for sim in sims:
    n_eff = int(sim.results.hiv.n_on_effective_art[-1])
    n_nonsupp = int(sim.results.hiv.n_on_nonsuppressive_art[-1])
    n_art = n_eff + n_nonsupp
    frac_eff = n_eff / n_art if n_art else float('nan')
    print(f'{sim.label:28s}  on_art={n_art:4d}  effective={n_eff:4d}  '
          f'nonsuppressive={n_nonsupp:4d}  frac_effective={frac_eff:.2f}')


# %% 4x4 grid: rows = scenario, columns = prevalence + 90-90-90 cascade,
# each panel showing men (blue) and women (red) separately.
#
# 90-90-90 (UNAIDS convention, cascade computed relative to the previous
# stage, not relative to total PLHIV each time):
#   1st-90: % of PLHIV diagnosed        = n_diagnosed / n_infected
#   2nd-90: % of diagnosed on ART       = n_on_art / n_diagnosed
#   3rd-90: % of those on ART who are
#           virally suppressed          = n_on_effective_art / n_on_art
#
# prevalence/n_infected are native hiv results (already split by sex);
# n_diagnosed/n_on_art/n_on_effective_art come from the StratifyCascade analyzer.

def safe_ratio(num, denom):
    num, denom = np.asarray(num), np.asarray(denom)
    return np.divide(num, denom, out=np.full(num.shape, np.nan), where=denom > 0)


def get_result(sim, name):
    """ Look up a (possibly sex-suffixed) result by name, from hiv or the StratifyCascade analyzer """
    hiv = sim.results.hiv
    if name in hiv:
        return np.asarray(hiv[name])
    return np.asarray(sim.results.stratify_cascade[name])


# (label, numerator template, denominator template or None for a direct value)
# '{s}' is filled in with 'm' or 'f'
metrics = [
    ('Prevalence', 'prevalence_{s}', None),
    ('1st-90: % diagnosed', 'n_diagnosed_{s}', 'n_infected_{s}'),
    ('2nd-90: % of diagnosed on ART', 'n_on_art_{s}', 'n_diagnosed_{s}'),
    ('3rd-90: % of on-ART suppressed', 'n_on_effective_art_{s}', 'n_on_art_{s}'),
]

sexes = [('Male', 'm', 'tab:blue'), ('Female', 'f', 'tab:red')]

fig, axes = plt.subplots(len(sims), len(metrics), figsize=(4 * len(metrics), 3 * len(sims)),
                          sharex='col', sharey='col')

for row, sim in enumerate(sims):
    timevec = sim.results.hiv.timevec
    for col, (metric_label, num_tmpl, denom_tmpl) in enumerate(metrics):
        ax = axes[row, col]
        for sex_label, s, color in sexes:
            num = get_result(sim, num_tmpl.format(s=s))
            y = num if denom_tmpl is None else safe_ratio(num, get_result(sim, denom_tmpl.format(s=s)))
            ax.plot(timevec, y, color=color, label=sex_label)
        if col > 0:
            ax.set_ylim(0, 1)
        if row == 0:
            ax.set_title(metric_label)
        if col == 0:
            ax.set_ylabel(sim.label, fontsize=9)
        if row == len(sims) - 1:
            ax.set_xlabel('Year')

axes[0, 0].legend(loc='lower right', fontsize=8)
fig.suptitle('HIV prevalence and 90-90-90 cascade by VLS scenario and sex')
fig.tight_layout()
plt.savefig('art_vls_intervention_demo.png', dpi=150)
print('Saved figure to art_vls_intervention_demo.png')
plt.show()
