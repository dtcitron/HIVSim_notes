"""
Demo: program VMMC vs. non-program traditional male circumcision.

--------------------------------------------------------------------------
WHAT CHANGED IN STISIM, AND WHY
--------------------------------------------------------------------------
Before this change, everything about circumcision (who's circumcised, when,
and how protective it is) lived entirely on the `VMMC` intervention: it
defined its own `circumcised`/`ti_circumcised` states and its own `eff_circ`
par, and applied the `rel_sus` reduction itself inside `VMMC.step()`.

That's fine as long as there is exactly one way to become circumcised. It
breaks down as soon as you want a second pathway — e.g. traditional
(non-program, non-medical) circumcision of boys before sexual debut — because
that second pathway has no way to see "has this person already been
circumcised by the other mechanism?" without reaching into a specific
intervention instance that might not even be present in a given sim.

The fix follows a pattern already established elsewhere in stisim for ART:
`on_art`, `diagnosed`, etc. live on the `HIV` disease module itself, not on
the `ART` intervention — `ART.step()` just calls `hiv.start_art(uids)`. We
did the same thing for circumcision:

  * `HIV` (stisim/diseases/hiv.py) now owns:
      - States: `circumcised` (bool), `circ_traditional` (bool — True if
        traditional/non-program, False if program VMMC; only meaningful when
        `circumcised` is True), `ti_circumcised`. A plain bool rather than an
        int code, since every read already gates on `circumcised` first, so
        there's no separate "none" value to encode — it's a genuine binary
        split (program vs. non-program), not an open-ended category set.
      - Pars: `eff_circ` (program VMMC efficacy, default 0.6 — unchanged
        value) and `eff_circ_traditional` (new, default 0.6, independently
        configurable).
      - `HIV.circumcise(uids, traditional=False)`: the single write path.
        It's idempotent — uids that are already circumcised are silently
        skipped, so a person's circ_traditional/ti_circumcised reflect
        whichever mechanism got to them FIRST, and a second mechanism can
        never overwrite it.
      - The `rel_sus` reduction is applied inside `HIV.update_transmission`
        (same place `rel_sus_age` and ART's `rel_trans` effect are applied),
        keyed off `circ_traditional`, so it works automatically for whichever
        mechanism did the circumcising — no `rel_sus` math lives in `VMMC`
        anymore.

  * `VMMC` (stisim/interventions/hiv_interventions.py) now only decides WHO
    and WHEN, via two independent pathways that both terminate in
    `hiv.circumcise()`:
      1. Program VMMC (unchanged behavior): coverage-target logic
         (`coverage=`), prioritized by a per-agent `willingness` score,
         restricted by the `eligibility` callable. Calls `hiv.circumcise(uids)`.
      2. Traditional circumcision (new): a one-time probabilistic
         assessment — `traditional_prob` (probability, default 0 = off) —
         applied to males once they reach `traditional_age` (default 15)
         IF they have not yet had sexual debut (checked via the sexual
         network's `over_debut`). This pathway is independent of
         `eligibility` (which only scopes program VMMC) and independent of
         `coverage` (it can run entirely on its own, with `coverage=None`).
         Calls `hiv.circumcise(uids, traditional=True)`.

Because both pathways write through the same `hiv.circumcise()` sink, a boy
who got traditional circumcision at 15 is automatically excluded from
program VMMC's "who still needs circumcising" pool as an adult, and vice
versa — with no special-case code required in either pathway.

--------------------------------------------------------------------------
THE FOUR SCENARIOS BELOW
--------------------------------------------------------------------------
1. baseline      - no circumcision of any kind.
2. program VMMC  - maintains 25% of uncircumcised men aged 25+ circumcised
                   (eff_circ, HIV default = 0.6).
3. traditional   - 10% of pre-sexual-debut boys are traditionally
                   circumcised at age 15 (eff_circ_traditional = 0.5, less
                   protective than the program's medical procedure).
4. combined      - both (2) and (3) running together.

Coverage note: VMMC's `coverage` is a *stock* (prevalence) target, not a
flow — see the comment in `VMMC.step()`. A plain proportion (0.25 here) is
interpreted directly as that target: VMMC tops up new circumcisions each
step to hold the eligible pool at 25% circumcised, rather than circumcising
a fixed count.
"""

import numpy as np
import matplotlib.pyplot as plt
import starsim as ss
import stisim as sti

n_agents = 10000
start = 2000
dur = 25  # years
eligible_age = 25       # program VMMC target population: men 25+
traditional_age = 15    # traditional circumcision assessed at this age
traditional_prob = 0.1  # 10% of pre-debut boys get traditional circumcision
eff_circ_traditional = 0.5  # less protective than program VMMC (default eff_circ=0.6)


class CircumcisionTracker(ss.Analyzer):
    """
    Records alive-male population and circumcision-by-type counts EVERY STEP.

    Both quantities change over time (births/deaths grow and shrink the
    population; `HIV.step_die()` resets `circumcised`/`circ_traditional` to
    False when someone dies, matching how `on_art`/`diagnosed` are
    cleared on death). That means neither can be recovered correctly from a
    single post-run snapshot:
      - `people.male.raw` at the end of the run includes preallocated
        capacity for not-yet-born agents, so a single final count wildly
        overstates the population at any earlier point in the sim.
      - Reconstructing "who was circumcised when" from the final
        `ti_circumcised` snapshot silently drops anyone circumcised early who
        has since died -- and early cohorts have had more time to die than
        recent ones, which biases the reconstructed curve to look like
        uptake ramped up gradually even when (as here) it was actually
        immediate. Recording the counts live, as the sim progresses, avoids
        both problems.
    """
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.name = 'circ_track'
        return

    def init_results(self):
        super().init_results()
        self.define_results(
            ss.Result('n_alive_male', dtype=int, auto_plot=False),
            ss.Result('n_circ_program', dtype=int, auto_plot=False),
            ss.Result('n_circ_traditional', dtype=int, auto_plot=False),
        )
        return

    def step(self):
        ti = self.ti
        ppl = self.sim.people
        hiv = self.sim.diseases.hiv
        self.results.n_alive_male[ti] = (ppl.male & ppl.alive).count()
        self.results.n_circ_program[ti] = (hiv.circumcised & ~hiv.circ_traditional).count()
        self.results.n_circ_traditional[ti] = (hiv.circumcised & hiv.circ_traditional).count()
        return


def make_sim(vmmc, label):
    """A fresh HIV+VMMC sim; only `vmmc` (and, via it, HIV's eff_circ_traditional) varies."""
    hiv = sti.HIV(beta_m2f=0.05, beta_m2c=0.1, init_prev=0.05,
                  eff_circ_traditional=eff_circ_traditional)
    sim = sti.Sim(
        start=start, dur=dur, n_agents=n_agents,
        diseases=[hiv],
        networks=[sti.StructuredSexual(recall_prior=True), sti.PriorPartners()],
        demographics=[ss.Pregnancy(fertility_rate=10), ss.Deaths(death_rate=10)],
        interventions=[vmmc] if vmmc is not None else [],
        analyzers=[CircumcisionTracker()],
    )
    sim.label = label
    return sim


# A plain proportion is a *stock* (prevalence) target: VMMC tops up to keep
# 25% of the eligible pool circumcised, not a count quota (see the coverage
# note in VMMC.step()'s docstring).
program_coverage = 0.25

program_eligibility = lambda sim: sim.people.age >= eligible_age  # noqa: E731

scenarios = [
    make_sim(None, 'baseline'),
    make_sim(sti.VMMC(coverage=program_coverage, eligibility=program_eligibility),
              'program VMMC'),
    make_sim(sti.VMMC(traditional_prob=traditional_prob, traditional_age=traditional_age),
              'traditional MC'),
    make_sim(sti.VMMC(coverage=program_coverage, eligibility=program_eligibility,
                       traditional_prob=traditional_prob, traditional_age=traditional_age),
              'program + traditional'),
]

for sim in scenarios:
    sim.run(verbose=0)

# --------------------------------------------------------------------------
# Summary
# --------------------------------------------------------------------------
print(f'{"scenario":24s}  {"cum. male infections":>22s}  {"n circ (program)":>17s}  {"n circ (traditional)":>21s}')
for sim in scenarios:
    hiv = sim.diseases.hiv
    male_infections = int(sum(hiv.results.new_infections_m))
    # Final counts among people CURRENTLY alive (matches step_die's reset-on-death
    # semantics); this is a live count, not a reconstruction, so no survivorship bias.
    n_program = (hiv.circumcised & ~hiv.circ_traditional).count()
    n_traditional = (hiv.circumcised & hiv.circ_traditional).count()
    print(f'{sim.label:24s}  {male_infections:22d}  {n_program:17d}  {n_traditional:21d}')

# --------------------------------------------------------------------------
# Plot: circumcision uptake by type, and cumulative male infections, over time
# --------------------------------------------------------------------------
fig, axes = plt.subplots(1, 2, figsize=(12, 5))
colors = dict(zip([s.label for s in scenarios], plt.cm.tab10.colors))

ax = axes[0]
for sim in scenarios:
    track = sim.results.circ_track
    timevec = track.timevec
    # Both numerator and denominator are live per-step counts (see
    # CircumcisionTracker) -- no post-hoc reconstruction, so no survivorship bias.
    n_alive_male = np.maximum(np.asarray(track.n_alive_male), 1)
    frac_program = np.asarray(track.n_circ_program) / n_alive_male
    frac_traditional = np.asarray(track.n_circ_traditional) / n_alive_male
    ax.plot(timevec, frac_program, color=colors[sim.label], linestyle='-', label=f'{sim.label} (program)')
    ax.plot(timevec, frac_traditional, color=colors[sim.label], linestyle='--', label=f'{sim.label} (traditional)')
ax.set_title('Circumcised fraction of all males, by type')
ax.set_xlabel('Year')
ax.set_ylabel('Fraction circumcised')
ax.legend(fontsize=7)

ax = axes[1]
for sim in scenarios:
    hiv = sim.diseases.hiv
    timevec = hiv.results.timevec
    cum_inf = np.cumsum(hiv.results.new_infections_m)
    ax.plot(timevec, cum_inf, color=colors[sim.label], label=sim.label)
ax.set_title('Cumulative new HIV infections in men')
ax.set_xlabel('Year')
ax.set_ylabel('Cumulative infections')
ax.legend(fontsize=8)

fig.suptitle('Program VMMC vs. traditional circumcision')
fig.tight_layout()
plt.savefig('vmmc_traditional_circumcision_demo.png', dpi=150)
print('Saved figure to vmmc_traditional_circumcision_demo.png')
plt.show()
