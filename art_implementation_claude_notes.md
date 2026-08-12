# ART implementation notes (for the `art` branch rewrite)

Companion to `art.md`. Covers: how ART currently works in STIsim/HIVsim, what EMOD's
`ARTMortalityTable` does (as a reference design), and what that implies for adding
an adherence/suppression-category state. All citations are `file:line` in
`/Users/daniel/Documents/IDM/stisim` and `/Users/daniel/Documents/IDM/EMOD` unless noted.

---

## 1. How ART distribution works in HIVsim today

Two-stage pipeline: `HIVTest` schedules *when* someone starts ART; `ART` decides
*whether* they actually do, subject to an optional coverage target.

1. **Testing** (`stisim/interventions/hiv_interventions.py:125-153`, base logic in
   `base_interventions.py:181-220`): eligible agents (default = undiagnosed) are tested
   each timestep at a probability derived from `test_prob_data`.
2. **Diagnosis**: positive → `hiv.diagnosed = True`, `hiv.ti_diagnosed = ti` (`hiv_interventions.py:130-131`).
3. **Scheduling**: `HIVTest` (not `ART`) draws a `dur_dx2tx` delay (default 0) and sets
   `hiv.ti_art = ti + delay` (`hiv_interventions.py:132-136`). `ART.init_pre` warns if no
   `HIVTest` is present, or if it's ordered after `ART` in the intervention list
   (`hiv_interventions.py:270-288`) — diagnosis must happen first in the same step.
4. **`ART.step()`** (`hiv_interventions.py:314-382`):
   - Stop anyone whose `ti_stop_art <= ti` (§3) via `hiv.stop_art(...)`.
   - `ready = (ti_art == ti) & ~on_art`, filtered through `art_initiation`
     (default `ss.bernoulli(p=0.9)` — flat 90%, **not** age/sex/CD4-dependent).
   - **No coverage target** (`coverage=None`, default): everyone who passes
     `art_initiation` starts unconditionally, no capacity constraint.
   - **Coverage target set**: only `n_to_treat - len(on_art)` slots are filled, via
     `prioritize_art` (`hiv_interventions.py:384-405`) — weights =
     `cd4_counts * (1/care_seeking)`, i.e. sicker + more care-seeking agents go first.
   - **`art_coverage_correction`** (`hiv_interventions.py:407-489`) actively pushes the
     on-ART population toward the target every timestep — can both add people and
     **actively remove** people (stopping ART for highest-CD4/lowest-care-seeking first)
     if the on-ART count exceeds target. Global or per-age/sex-stratum variants.

Real usage (`hivsim_examples/zimbabwe/sim.py:86,92` + `art_coverage.csv`): absolute-count
coverage target (`n_art` column), no age/sex stratification, no adherence modeling.

`docs/examples/art_interruptions.qmd` — the existing "ART interruption" example is purely
an **aggregate coverage-curve perturbation** (edit the `n_art`/`p_art` time series, let
`art_coverage_correction` chase the lower target). There is no individual-level
interruption/adherence state today.

## 2. Current effects of ART on individuals

**Mortality** — not a relative-risk multiplier. `hiv.start_art()` (`diseases/hiv.py:586-634`)
snapshots `cd4_preart`, draws `dur_on_art` (see §3) to set `ti_stop_art`, and **nullifies
the previously-scheduled natural-history death time** (`ti_latent`/`ti_falling`/`ti_zero`
set to NaN, acute/latent/falling states cleared). Each step, on-ART agents get CD4
reconstitution via a logistic curve (`cd4_increase`, `hiv.py:273-295`). Crucially: the
probabilistic CD4-stratified death draw (`make_p_hiv_death`) is computed **only for
`off_art` agents** (`hiv.py:400-404`) — on-ART agents are excluded from that hazard
entirely, not down-weighted. `HIVPars.art_efficacy` (default 0.96) is **not** used for
mortality — naming trap, it's transmission-only (below).

**Transmission — confirmed to exist**, contrary to the branch's stated assumption.
`HIV.update_transmission` (`hiv.py:456-498`):
```python
if self.on_art.any():
    full_eff = self.pars.art_efficacy                 # default 0.96
    time_to_full_eff = self.pars.time_to_art_efficacy  # default 6 months
    timesteps_on_art = ti - self.ti_art[art_uids]
    new_on_art = timesteps_on_art < time_to_full_eff/self.dt
    efficacy_to_date = full_eff  # ramped linearly to full_eff over new_on_art period
    self.rel_trans[art_uids] *= 1 - efficacy_to_date
```
`rel_trans` ramps linearly 0% → 96% suppression over 6 months post-`ti_art`, recomputed
fresh every step from `on_art`/`ti_art` (reset to 1 each step, `hiv.py:468-469`).

**PMTCT** — separate effect: mothers `on_art & pregnant`/`breastfeeding` reduce **infant**
`rel_sus` by `pmtct_efficacy` (default 0.96) on top of the mother's own `rel_trans`
reduction (`hiv_interventions.py:358-380`).

## 3. What affects ART retention time

**(a) Random duration draw at initiation** — the actual mechanism. `dur_on_art =
pars.dur_on_art.rvs(uids) * rel_dur_on_art`, `dur_on_art = ss.lognorm_ex(3yr, 1.5yr)`
(`hiv.py:65-66`), optionally time-trended (`dur_on_art_trend`). Drawn **once at start**,
**not conditioned on anything about the agent** (no CD4, age, or adherence) except the
flat `rel_dur_on_art` scalar and the trend. Checked every step: `ti_stop_art <= ti` →
`hiv.stop_art()`.

**(b) Coverage-driven forced discontinuation** — independent of (a); `art_coverage_correction`
pulls people off early if on-ART count exceeds target, prioritizing removal of
highest-CD4/lowest-care-seeking agents.

**On stop** (`hiv.py:636-660`): new post-ART CD4-decline duration drawn from `cd4_preart`,
current `cd4`, `cd4_potential` → new `ti_zero`. Restart cycles are fully supported at the
state-machine level.

**Nothing that models graded adherence exists today** — no adherence probability, no
suppression-failure branch, no drug resistance, no distinction between "stopped by choice"
vs "stopped by program capacity" downstream of `stop_art`. `care_seeking` (`hiv.py:331-339`,
`~N(1, 0.5)` floored at 0.01) is the closest individual trait, but only acts as a
*prioritization weight* for who starts/stops under capacity constraints — it doesn't gate
initiation probability or retention duration directly.

## 4. Implementing a 4-state ART-adherence property in STIsim

**No true categorical/enum `Arr` type exists** in `starsim/starsim/arrays.py`. Options are
`FloatArr`, `IntArr` (int-backed, `nan=int_nan` sentinel, docstring recommends avoiding it,
used exactly once in all of stisim as a lifetime *counter*, not a category), `BoolArr`,
`BoolState` (auto-registers for `n_<name>` result generation via `_auto_states`,
`modules.py:576-578`). "Individual property" is EMOD/DTK terminology — the STIsim-native
concept is just **state**, registered via `self.define_states(...)`.

**Idiomatic pattern used everywhere else in this codebase (SEIS: `exposed`/`asymptomatic`/
`symptomatic`/`pid`, HIV: `acute`/`latent`/`falling`/`on_art`/`never_art`): one `BoolState`
per named category**, with mutual exclusivity hand-maintained by the update code (not
framework-enforced — a real gap; e.g. `hiv.py:373-376`, `609-611` each duplicate the
toggle-boilerplate inline rather than centralizing it).

**Recommendation**: follow that convention —
```python
self.define_states(
    ss.BoolState('art_naive', default=True),
    ss.BoolState('on_effective_art'),
    ss.BoolState('on_nonsuppressive_art'),
    ss.BoolState('art_discontinued'),
)
```
over a single `IntArr` code, because: (1) matches the existing idiom exactly, (2) each
state gets an automatic `n_<name>` result for free, (3) composable with `&`/`|`/`~` exactly
like `on_art` is used throughout `ART.step()`/`hiv.py` today. Worth adding a
`set_art_state(uids, new_state)` helper to centralize the mutual-exclusivity toggling that's
currently duplicated inline. `on_effective_art | on_nonsuppressive_art` would become the new
definition of "on ART" for compatibility with existing `rel_trans`/mortality code, with the
two sub-states then driving *differential* transmission/CD4-growth parameters.

## 5. EMOD's `ARTMortalityTable` — reference design

Class hierarchy: `AntiretroviralTherapy` (base — Cox/IeDEA parametric mortality model, handles transmission suppression + viral-suppression ramp-up) → `ARTMortalityTable` (overrides only the mortality-duration sampler with an empirical lookup table) and `AntiretroviralTherapyFull` (sibling — adds random dropout timer).

**What it does**: draws a single random *survival time from ART enrollment to AIDS death* **once**, at ART (re-)distribution — not a per-timestep hazard. That duration replaces the individual's natural-history death timer entirely (`InfectionHIV::SetupSuppressedDiseaseTimers`). `ARTMortalityTable` can be **redistributed to someone already on ART**, which fully resets this death timer (unlike the base class, which refuses to redistribute to already-on-ART
individuals).

**Table axes**: `ART_Duration_Days_Bins` × `Age_Years_Bins` × `CD4_Count_Bins` only (`ARTMortalityTable.h:30-33`). **Sex and adherence are NOT table axes** — both are handled by campaign-config targeting (distributing differently-parameterized intervention instances to different sub-populations by `Target_Gender` / `Property_Restrictions`), not inside the C++ class. One real example (`Regression/HIV/61_HIV_ARTMortalityTable/campaign.json`) routes `Adherence:HIGH` individuals to `ARTMortalityTable` and `Adherence:LOW` individuals to the plain Cox-model `AntiretroviralTherapy` — i.e. adherence selects *between model classes*, not axes of one table. **If the STIsim port wants age/sex/duration + adherence as first-class axes of one unified table, that's a deliberate generalization beyond what EMOD's C++ actually implements** — flag as an explicit design decision, not an EMOD port.

**Sampling algorithm** (`ARTMortalityTable.cpp:111-142`): a piecewise-exponential
survival-time sampler exploiting the memoryless property — draw `Exp(rate)` using the current ART-duration bin's rate; if the draw exceeds the bin boundary, discard and redraw with the next bin's rate; repeat. Table values are **annual hazard rates**. Age and CD4 are captured **once** at the top of the function and reused unchanged across the whole loop — no re-aging or CD4-reconstitution mid-loop (a known EMOD simplification, worth deciding whether to replicate or improve on).

**Adherence in EMOD** is a `viral_suppression` boolean (campaign-configured, not
individually stochastic) feeding an `ARTStatus` enum
(`WITHOUT_ART_* / ON_PRE_ART / ON_ART_BUT_NOT_VL_SUPPRESSED / ON_VL_SUPPRESSED /
ON_BUT_FAILING / ON_BUT_ADHERENCE_POOR / OFF_BY_DROPOUT`, `HIVEnums.h:22-33`):
- **`viral_suppression=true`**: table-drawn death timer installed, transmission ramps to
  the suppression multiplier (default 0.08) — *but* a **second, independent "failing"
  timer** carves a terminal window (default 9 months, `AIDS_duration_in_months`) out of
  the *end* of that same table-drawn duration, during which the person is fully infectious
  again (`ON_BUT_FAILING` — no transmission suppression, plus a 10x AIDS-stage infectivity
  multiplier) even though mortality still follows the original table-drawn schedule.
- **`viral_suppression=false`** (`ON_BUT_ADHERENCE_POOR`, i.e. non-suppressive ART): no
  mortality benefit at all (natural-history timer untouched) and no transmission benefit.

| ARTStatus | Mortality | Transmission |
|---|---|---|
| ramping up (VL not yet suppressed) | table-drawn timer | ramping 1.0 → multiplier |
| VL suppressed | table-drawn timer | fixed at multiplier (0.08 default) |
| on but failing (terminal window) | same table-drawn timer | **no suppression** + 10x AIDS multiplier |
| non-suppressive ART | **no ART effect** | **no suppression** |
| off by dropout | natural-history timer recomputed from CD4-at-dropout | no suppression |

**Table shape example** (`componentTests/testdata/ARTMortalityTableTest/DocExample.json`):
3 duration bins × 4 age bins × 5 CD4 bins, `MortalityTable[artDurationIndex][ageIndex][cd4Index]`
= annual hazard rate (e.g. 0.6 down to 0.01). A degenerate 1×1×1 table is legal. All-zero
table is legal and means "no AIDS mortality for that cell" (`daysUntilDeath = FLT_MAX`).

## 6. Open questions to resolve before implementing (step 3 in `art.md`)

- Do we want age/sex/duration/adherence as axes of **one** mortality table (a genuine
  generalization of EMOD), or keep EMOD's split of "adherence selects model class,
  table handles age/duration/CD4 only"? Affects whether `ARTMortalityTable`-equivalent
  data structure needs an adherence dimension at all.
- Do we want EMOD's "terminal failing window with restored infectiousness" behavior
  (a real, deliberate EMOD design choice) or a simpler two-state (effective/non-suppressive)
  model without a third transitional "failing" state?
- EMOD's mortality-table draw is a one-time point-value; captures age/CD4 once and never
  updates until redistribution. STIsim's existing mortality mechanism (excluding on-ART
  agents from the stochastic death draw + reconstituting CD4 continuously) is already a
  *continuous* per-timestep model, not a point-value draw. Porting EMOD's table-and-draw
  approach wholesale vs. keeping the continuous CD4-driven approach but branching it by
  adherence category are architecturally different choices — decide before starting §3 of
  `art.md`.
- Typical ART retention/discontinuation duration by age (checklist item in `art.md`) is
  still unanswered — current `dur_on_art` (lognormal, mean 3yr/std 1.5yr) has no age
  dependence; needs a data source.

## 7. `ARTMortalityTable` vs. current HIVsim mortality — direct comparison (art.md step 4.1)

The campaign-file excerpt in `art.md` gives four `ARTMortalityTable` instances — males/females × effective/non-suppressive ART — each an annual-mortality-rate lookup table with three axes: age (4 bins), CD4 count at initiation (7 bins), and time already spent on ART (5 bins, from `ART_Duration_Days_Bins: [182, 365, 730, 1095, 45625]`, i.e. `<6mo / 6-12mo / 1-2yr / 2-3yr / >3yr`). The numbers below all come from the lowest-CD4 column, for males aged 25-35, unless stated otherwise — verified by direct computation, not by eye.

**Finding 1: non-suppressive ART is far deadlier than effective ART, at the same age/CD4/duration.** Comparing the two tables column-by-column at the 6-12-month duration bin: except for the two lowest-CD4 columns (where both tables are already near-ceiling and only differ by ~2x, 0.185 vs. 0.095), every other CD4 column shows non-suppressive ART at a remarkably consistent ~11.4x higher mortality than effective ART (e.g. 0.335 vs. 0.029; 0.167 vs. 0.015 — the ratio holds to two significant figures across columns). That near-constant ratio suggests the non-suppressive table was built as a simple scalar multiple of the effective table over most of the CD4 range, rather than independently fit.

**Finding 2: mortality drops steeply the longer someone has already been on ART.** Within a fixed age/CD4/adherence cell, effective-ART mortality for males aged 25-35 goes 0.2176 → 0.0945 → 0.0288 → 0.0186 → 0.0135 across the five duration bins (`<6mo` through `>3yr`) — roughly a 16-fold drop from the first bin to the last. This encodes "early mortality risk": people who just started ART are often still sick (opportunistic infections, not yet immune-reconstituted), and that risk fades the longer they stay on treatment.

**Finding 3: mortality rises mildly with age, at fixed CD4/duration/adherence.** For the 6-12-month duration bin, lowest-CD4 column: 0.0875 / 0.0945 / 0.1021 / 0.1102 across the four age bins `<25 / 25-35 / 35-45 / 45+`. Don't confuse this with the AYA-dropout-rate discussion in step 5 below — that's about how often young people *leave* ART; this is about mortality risk for people who *stay* on ART, which the table says rises slightly with age. They're two different effects, both real, easy to conflate when porting.

**Comparison table:**

| Dimension | EMOD `ARTMortalityTable` | Current HIVsim |
|---|---|---|
| When is mortality risk assigned? | Once, at ART (re-)initiation: draws a full survival-time-to-AIDS-death via piecewise-exponential sampling over the duration bins, replacing the natural-history death timer (`ARTMortalityTable.cpp:111-142`) | Continuously: every timestep, `make_p_hiv_death` (`hiv.py:310-319`) computes a fresh CD4-stratified hazard — but **only for agents with `on_art == False`** (`hiv.py:405-406`) |
| Mortality while on effective ART | Non-zero, table-driven, declining with duration (see above) | **Exactly zero** — `on_art` agents are excluded from the hazard entirely, regardless of duration or CD4 |
| Mortality while on non-suppressive ART | Non-zero, table-driven, markedly *higher* than effective ART at the same age/CD4/duration (this campaign routes both to `ARTMortalityTable`, just with different tables — a stricter design than the general split described in §5, where non-suppressive ART elsewhere gets *no* mortality benefit at all) | **Also exactly zero** — `on_nonsuppressive_art` is not distinguished from `on_effective_art` anywhere in `make_p_hiv_death`/`step_state`; the *only* place adherence category currently matters is `rel_trans` (transmission), not mortality |
| Age dependence | Explicit axis (4 bins), modest effect at fixed CD4/duration (see above) | None at all — CD4 stratification only |
| Sex dependence | Handled by routing to sex-specific campaign events/tables, not a table axis | None |
| Duration-since-initiation dependence | Explicit axis (5 bins), large effect (~16x from `<6mo` to `>3yr` in the example above) | None — the only "duration" dependence is indirect, via the logistic CD4 reconstitution curve (`cd4_increase`, `hiv.py:273-295`), which does not directly gate a mortality hazard while on ART |
| CD4 dependence | Captured **once** at ART initiation, frozen for the rest of the table-drawn spell | Continuously updated every step via `cd4_increase`/`post_art_decline`, but the resulting mortality-relevant CD4 value only matters for agents who are **off** ART |

**Bottom line**: this is the single largest behavioral gap relative to `art.md`'s stated goal 1.2 ("mortality should also depend on whether an individual adheres to ART"). Today, adherence category has zero mortality consequence — someone on non-suppressive ART has identical (zero) HIV mortality risk to someone on fully-suppressive ART, differing only in transmission efficacy (35% vs. 96%, `hiv.py:62-63`). EMOD's reference design says the opposite: failing/non-suppressive treatment should carry mortality risk close to (or, in the early duration bins, exceeding) untreated risk, not zero excess risk. Any mortality-model rewrite needs to at minimum re-enable *some* CD4-stratified hazard for `on_nonsuppressive_art` agents, whether or not the full age/duration table structure is ported.

## 8. Plan: demographic-dependent ART retention duration (art.md step 5.3)

Empirical grounding comes from `art_adherence_test.py` (step 5.2): it tracks every agent's on-ART "spells" (start at `ti_art`, end at `stop_art`/death/censoring) and produces both a raw histogram of completed-spell durations and a Kaplan-Meier retention curve that correctly accounts for right-censored (still-on-ART-at-sim-end) spells — a plain mean/median of completed spells alone is biased short, since it implicitly ignores the people who are still retained and simply haven't left yet. In a 10k-agent test run with a coverage target ramping to 60% and holding (2000–2040), completed spells averaged ~3.0yr (median 2.6yr), but the KM-based restricted-mean estimate was ~3.4yr and rising at the last observed time point, already above the "nominal" `dur_on_art` mean of 3yr — i.e., in this scenario, coverage correction is (if anything) partially offsetting the raw draw's tail rather than aggressively truncating it. That balance will shift with different coverage trajectories, which is exactly why an empirical check matters more than reasoning about `dur_on_art`'s parameters in isolation.

Two architecturally different ways to introduce age/sex-dependence, given the target numbers in `art.md` step 5.1 (~15%/yr dropout for AYA, ~5%/yr for adults 25+):

**Option A — age/sex multiplier on the existing point-draw** (small change). Add a
`rel_dur_on_art_age` parameter following the exact convention already used for
`rel_sus_age` (`hiv.py:42-46`): a list of `(age_lo, age_hi, sex, multiplier)` tuples, applied
multiplicatively to the `dur_on_art` draw in `start_art` (`hiv.py:637-638`) alongside the
existing scalar `rel_dur_on_art` and `dur_on_art_trend`. Pros: minimal code change, reuses an
existing pattern, keeps `ti_stop_art` precomputed (no touch to `ART.step`/coverage-correction
logic). Cons: annual dropout probabilities don't map cleanly onto "multiply the mean of a
lognormal" — because the lognormal's *hazard* isn't constant over time, a single multiplier
can't exactly reproduce "15%/yr" as a flat rate; it can only be tuned to approximately match
one reference time point (e.g., 1-year retention) via search, using `art_adherence_test.py`
itself as the fitting target: run the sim with a candidate multiplier, read the empirical
1-year discontinuation rate for that age band off the KM curve, and adjust until it matches
15% (AYA) / 5% (adult).

**Option B — per-timestep age/sex-dependent dropout hazard** (bigger change, recommended).
Replace (or optionally override) the single point-draw with a per-timestep discontinuation
probability `p_stop_art(age, sex)`, evaluated every step the same way `make_p_hiv_death`
already evaluates CD4-based mortality (`ss.peryear(rate).to_prob(dt)`), filtered against
`on_art` agents in `ART.step`/`HIV.step_state`. This maps *directly* onto the given data
(`0.15/yr` under some age cutoff, `0.05/yr` above it — a lookup table keyed by age bin, easily
extended to sex or `on_effective_art`/`on_nonsuppressive_art` if adherence categories turn out
to have different dropout rates too) with no reparametrization/fitting step required. It's
also more consistent with the rest of `hiv.py`'s style (continuous per-timestep hazards) than
the current point-draw-at-initiation design. Cost: `ti_stop_art` can no longer be precomputed
at `start_art` time; discontinuation becomes a live Bernoulli filter each step (removing the
"redraw scaled by `rel_dur_on_art`/trend" simplicity, or reimplementing those as scalars
multiplying the per-step rate instead of the drawn duration) and touches the
`art_coverage_correction` interaction more directly, since now two independent processes
(age-based dropout hazard + coverage-driven forced removal) compete to end a spell each step.

**Recommendation**: Option B for the eventual implementation — it encodes the target rates
directly, generalizes cleanly to more covariates, and matches the existing mortality-hazard
idiom — but it's a bigger refactor than Option A, so it's reasonable to prototype with Option
A first if a quick approximate answer is more valuable than getting the exact target rates.
Either way, the validation loop is the same: after implementing, rerun
`art_adherence_test.py` with the tracker's spells stratified by age at ART initiation (the
script would need a small extension — stash `hiv.age[uid]` at spell start alongside
`start`/`stop`) and confirm the age-specific 1-year KM discontinuation rate lands near 15%
(AYA) and 5% (adults 25+).

## 9. Distance from current state to an `ARTMortalityTable`-equivalent (art.md step 3, mortality half)

Step 3 of `art.md` was implemented for the *adherence-state and transmission* half only (`on_effective_art`/`on_nonsuppressive_art` states, differential `rel_trans`, `art_transmission_test.py`) — the *mortality* half of goal 1.2 ("mortality should also depend on whether an individual adheres to ART") was never wired up, which is exactly the gap §7 quantifies. The good news:
most of the state this needs already exists, for reasons unrelated to mortality. What's
actually missing is small and localized.

**Already in place (no new work needed):**
- `on_effective_art` / `on_nonsuppressive_art` — the adherence selector itself (`hiv.py:122-123`)
- `ti_art` — time of initiation, so duration-since-initiation (`ti - ti_art`) is a one-line
  computation, already used for the transmission ramp-up (`hiv.py:501`)
- `cd4` (continuously updated) *and* `cd4_preart` (frozen at initiation) — both EMOD-faithful
  (frozen) and continuous CD4-lookup variants are already tracked, no new state needed either way
- `people.age`, `people.female`/`people.male` — standard Starsim `People` attributes, already
  used elsewhere in the module (e.g. `rel_sus_age`)
- The **per-timestep hazard idiom itself** — `make_p_hiv_death` (`hiv.py:310-319`) already shows
  the exact pattern (`np.digitize` into bins → `ss.peryear(rate).to_prob(dt)`) that a table
  lookup would reuse

**What's actually missing:**
1. **A table data structure + parameter.** Something like `HIVPars.art_mortality_table`, keyed
   by adherence category, holding age/duration/CD4 bin edges and a hazard array — analogous in
   spirit to `rel_sus_age`'s `(age_lo, age_hi, sex, multiplier)` tuple-list convention, or a
   direct 3-D (or 4-D, if sex is folded in as an axis rather than routed via separate tables
   the way EMOD does) array + bin-edges dict. Default `None` so existing sims are unaffected
   unless a user opts in — same convention as `dur_on_art_trend`/`rel_sus_age`.
2. **A lookup method**, e.g. `get_art_mortality_hazard(uids)`: digitize age, `ti - ti_art`, and
   CD4 (continuous `cd4` recommended over frozen `cd4_preart` — see decision below) into their
   respective bins, select the effective- or non-suppressive-ART table via the existing
   adherence booleans, and convert the annual rate to a per-step probability the same way
   `make_p_hiv_death` already does.
3. **A ~10-15 line change to `step_state`'s death block** (`hiv.py:403-417`): today it computes
   `p_death` only for `off_art = infected & ~on_art` and filters just that group. It needs a
   second hazard computed for on-ART agents via (2), with both groups' filtered UIDs feeding
   the same `sim.people.request_death(...)` call.
4. **Sex handling.** EMOD routes sex via separate campaign-targeted intervention instances, not
   a table axis. STIsim has no equivalent per-population routing mechanism, so the simplest port
   is folding sex in as a 4th lookup axis on one table object rather than maintaining two nearly-
   identical parameter blocks.
5. **Default table values.** art.md's campaign excerpt gives one real (if narrow) example. Two
   options: ship those exact numbers as STIsim's built-in default, or leave the parameter unset
   by default and use those numbers only as a worked example/test fixture (e.g. a new
   `art_mortality_test.py`, mirroring how `art_transmission_test.py` already exercises the
   adherence split). Recommend the latter — baking in someone else's calibrated numbers as a
   silent default is the kind of thing that quietly breaks a downstream user's calibration.

**Two decisions still open (carried over from §6, now concrete enough to actually make):**

*CD4 lookup: continuous vs. frozen-at-initiation.* When the new mortality table looks up someone's CD4 bin, should it use their CD4 count *right now, this timestep* (`cd4`), or a snapshot of what their CD4 was *on the day they started ART* (`cd4_preart`, never updated again after that)? EMOD uses the frozen snapshot, but only because of how EMOD's mechanism works: it draws one death-timer value a single time, at initiation, so age/CD4 only ever get read once anyway — freezing them costs nothing extra in that design. STIsim's mortality hazard, by contrast, is already re-evaluated fresh every timestep (see the point above about `make_p_hiv_death` reusing this idiom), and STIsim already recomputes a person's live CD4 every timestep regardless, for the reconstitution curve (`cd4_increase`). So there's no engineering reason to freeze it here — recommend using the live `cd4` value at each lookup. Concretely: someone who started ART two years ago at a low CD4 count has likely reconstituted significantly since then; looking up their *current* CD4 places them in a healthier (lower-mortality) bin, which better reflects reality, whereas freezing it at the two-year-old low value would overstate their ongoing mortality risk indefinitely.

*Terminal "failing window" behavior.* This is a specific EMOD design detail worth walking through concretely, because it's easy to gloss over. When EMOD puts someone on suppressive ("effective") ART, it schedules two separate clocks at once, not one: (1) a death timer, drawn once from the mortality table, saying roughly when this person might eventually die of AIDS despite treatment; and (2) a "failing" timer, which carves out a fixed window (9 months by default) immediately before that scheduled death. During most of someone's time on effective ART, they're genuinely suppressed — low transmission risk, as expected. But once they enter that final pre-death window, EMOD flips them into a distinct `ON_BUT_FAILING` state: even though they are nominally still "on effective ART," their transmission suppression is switched off and they become fully infectious again (in fact 10x more infectious, matching late-stage AIDS), while their mortality countdown keeps following the original schedule unaffected. The intent is to model a real clinical pattern — ART eventually stops controlling the virus for some people even while they keep taking it (e.g. resistance, advanced disease), and during that terminal decline they become transmissible again. Reproducing this in STIsim would mean adding a third adherence sub-state on top of the two we already built (`on_effective_art`, `on_nonsuppressive_art`) — something like "on effective ART but currently failing" — with its own transmission rule and its own entry/exit timing logic. Since this isn't anything `art.md` actually asked for, and it's a nontrivial addition on top of everything else in this plan, recommend leaving it out unless you specifically want to model this late-stage-failure-and-renewed-infectiousness behavior.

**Why this is a small patch, not a rewrite:** the alternative — porting EMOD's actual mechanism
(a one-time piecewise-exponential survival-time draw at initiation, replacing the natural-
history death timer) — would be the big-effort path, since it conflicts with how `ti_zero`/
natural-history timers are already used elsewhere in this class and would need those to be
suppressed/replaced for on-ART agents rather than just adding a parallel hazard. Staying with
STIsim's existing continuous-hazard architecture and merely *extending its stratification* (age
× duration × CD4 × adherence, instead of CD4 only) is architecturally consistent with the rest
of `hiv.py` and confined to items 1-3 above.

**Validation:** no new test script is strictly required — `art_transmission_test.py`'s existing
`ARTStatusTracker` already plots deaths and death *rate* stratified by ART status
(`art_naive`/`on_effective_art`/`on_nonsuppressive_art`/`art_discontinued`), and today those
on-ART lines sit at exactly zero. Rerunning that unmodified script after implementing (1)-(3)
should show non-zero, adherence-differentiated, age/duration-varying mortality among the

## 10. Summary of substantive code changes: ART adherence & differential transmission (art.md step 3)

This section is a plain inventory of what actually changed in `stisim` to add the adherence split and its transmission effect — the mortality-side changes are §7/§9, this is transmission/adherence only. All citations are `file:line` in `/Users/daniel/Documents/IDM/stisim` unless noted.

**New states** (`hiv.py:243-247`): `art_naive` (default `True`) replaces the old `never_art` — confirmed no remaining references to `never_art` anywhere in `stisim`. Two new sub-states, `on_effective_art` and `on_nonsuppressive_art`, sit alongside the existing `on_art` (which now means "on ART of either kind," true whenever either sub-state is true) and `art_discontinued`. The four states together form the intended `ARTNaive`/`OnEffectiveART`/`OnNonSuppressiveART`/`ARTDiscontinued` categorization from `art.md` step 3.1, implemented as four separate `BoolState`s rather than one enum/int code — matching the existing codebase idiom (see §4 above) rather than introducing a new pattern.

**New parameters** (`hiv.py:173-176`): `effective_art_efficacy=0.96` and `nonsupp_art_efficacy=0.35` replace the old single `art_efficacy`. `p_effective_art=ss.bernoulli(p=1.0)` is new — the probability that a newly-initiated agent achieves viral suppression rather than landing on non-suppressive ART; default preserves old behavior (always effective) unless overridden.

**`start_art()`** (`hiv.py:775-809`): gained an optional `p_effective_art` argument (float, or an already-initialized `ss.Dist`). For each newly-treated cohort, agents are split into `effective_uids`/`nonsupp_uids` via `self.pars.p_effective_art.filter(uids)` (or the passed-in override), and `on_effective_art`/`on_nonsuppressive_art` are set mutually exclusively (`hiv.py:806-809`). `art_naive` is cleared only for agents who were actually naive before this call (`newly_treated = uids[self.art_naive[uids]]`, `hiv.py:792-793`), so re-initiation after a discontinuation doesn't incorrectly re-trigger naive-only logic elsewhere.

**`update_transmission()`** (`hiv.py:675`): the flat `full_eff = self.pars.art_efficacy` lookup was replaced with `full_eff = np.where(self.on_effective_art[art_uids], self.pars.effective_art_efficacy, self.pars.nonsupp_art_efficacy)` — a per-agent efficacy selected by adherence category. The surrounding mechanism (linear ramp from 0 to `full_eff` over `time_to_art_efficacy`, reset to 1 and recomputed fresh every step) is unchanged from before the adherence split; only the efficacy *value* being ramped toward became adherence-dependent.

**State cleanup**: both new sub-states are cleared alongside the rest of the ART state machine wherever the others already were — on death (`step_die`, `hiv.py:605-608`) and on stopping ART (`stop_art`, `hiv.py:858-859`) — so there's no dangling `on_effective_art=True` on an agent who's `art_discontinued` or dead.

**`ART` intervention** (`stisim/interventions/hiv_interventions.py`): gained its own `p_effective_art` par (default `ss.bernoulli(p=1.0)`, `hiv_interventions.py:266`), documented and forwarded to `hiv.start_art()`.

**Known gotcha — `p_effective_art` is only forwarded on the coverage-target path.** `ART.step()` has two call sites for `hiv.start_art()`: the coverage-target path, via `prioritize_art()`, correctly forwards it (`hiv.start_art(start_uids, p_effective_art=self.pars.p_effective_art)`, `hiv_interventions.py:412`). The no-coverage-target path does not: `hiv.start_art(dx_to_treat)` (`hiv_interventions.py:353`) omits the argument entirely, so it silently falls back to `HIVPars.p_effective_art` (module-level default, always-effective) rather than whatever the `ART` intervention instance was configured with. Concretely: `sti.ART(p_effective_art=0.5)` with no `coverage=` argument will NOT produce any non-suppressive ART at all — everyone ends up on effective ART, silently, because the 50% only ever gets applied on the code path that requires a coverage target. `art_transmission_test.py`'s own working example sidesteps this by always setting `coverage=` alongside `p_effective_art=` (`art_transmission_test.py:123-126`), so this gap wouldn't have been caught by that test. Worth a follow-up fix (`hiv_interventions.py:353` → `hiv.start_art(dx_to_treat, p_effective_art=self.pars.p_effective_art)`) — flagging rather than fixing here since it's outside what was asked for this pass.

**Test coverage**: `art_transmission_test.py`'s `ARTStatusTracker` analyzer verifies the adherence split end-to-end — it plots population counts on effective vs. non-suppressive ART over time, and (per §7/§9) deaths/death-rate/onward-transmission-count by ART status. It's also the script that first surfaced, empirically, that on-ART mortality was flat zero regardless of adherence category (the finding that motivated §7).
on-ART groups — the same script becomes its own before/after check.

## 11. Where the numbers actually came from, and reconciling the Weibull literature against continuous CD4 (art.md step 6)

Ingested everything in `/Users/daniel/Documents/IDM/ART_Mortality_Table_Docs/`: Monisha Sharma's summary email plus three paper abstracts (`artmorttable_abstracts.md`), two internal supplemental write-ups (`Supplement ART Mortality and Transmission DRAFT 2019-05-22.docx`, `Supplemental ART mortality and transmission write up 9 6 19.docx` — read via `textutil -convert txt`, since the Read tool can't open `.docx` directly), and the three source PDFs (Gupta 2011, May 2010 Lancet, Johnson 2013 PLoS Medicine; a fourth, "May mortality by CD4 count.pdf", is almost certainly May et al. 2016 CID — cited by role and exact numbers in the write-up below, not opened directly, since `pdftoppm`/poppler isn't installed in this environment and the write-up already quotes everything needed from it).

**There are two distinct EMOD-lineage mortality models, not one**, and Monisha's email ("Both the Weibull mortality function and ART mortality table are based on the IeDEA ART program dataset") describes both:

1. **The Weibull mortality function** (`Supplement ART Mortality and Transmission DRAFT 2019-05-22.docx`, dated May 2019 — the earlier of the two write-ups). A **pure time-since-ART-initiation** hazard, with **no CD4 term at all** — the model's own variable list is shape/scale parameters, viral-load-stratum weights, the VL hazard ratio, and `x` = time since ART initiation. Fit by grid search (λ∈[100,400], k∈[0.01,0.50], minimizing squared error against IeDEA survival) to `λ=302, k=0.32` for the virally-suppressed group; the non-suppressed group's hazard is that same curve scaled by a **constant** hazard ratio `HRns=1.96` "at all time points" — i.e. a flat multiplier, not itself CD4- or age-varying. `k=0.32 < 1` means a Weibull hazard that is *monotonically decreasing* in time since ART initiation and never flattens to a plateau — a much sharper, longer-tailed early-mortality signature than an exponential decay would produce.

2. **The ART mortality table** (`Supplemental ART mortality and transmission write up 9 6 19.docx`, September 2019 — the later, more elaborated write-up; this is the one EMOD's `ARTMortalityTable` JSON actually implements). Four **relative survival models** (Johnson et al. 2013, ref [1] in that write-up) fit separately by sex and by duration-since-ART (`<12mo` / `≥12mo`), stratified by **age, baseline CD4 category, and regional program**, background-mortality-subtracted, adjusted by a `0.5` WHO-stage correction factor (their high-CD4 patients were disproportionately symptomatic, since the source guidelines only allowed CD4>200 ART initiation for WHO stage IV), with CD4≥350 risk ratios **extrapolated from a separate high-income cohort** (May et al. 2016, since IeDEA itself has almost no empirical data above CD4 350 — this is a real evidentiary gap, not a modeling choice, and it's directly relevant below since it's exactly the CD4 range where the current code's invariant-violation bug was found). Table Xa (baseline rates, CD4<50) is further split into 5 duration bins (`1st 6mo / 6-11mo / 12-23mo / 24-35mo / ≥36mo` — identical bin edges to EMOD's `ART_Duration_Days_Bins`), reflecting the *same* `HRns=1.96` VL/adherence hazard ratio as the Weibull write-up (verified: `0.1715/0.0875 = 0.183671/0.093710 = 1.96` exactly, at every duration past the first 6 months, for both sexes) — the two write-ups agree on this one number even though they're different model families. Table Xb's age effect is a **continuous per-decade** risk ratio (`1.07`-`1.13` depending on sex/duration), not a step table.

**On the duration-6-months detail**: the write-up states plainly *"we assumed it takes 6 months to achieve viral suppression; therefore individuals experienced the same mortality rate during the first 6 months after ART initiation regardless of treatment adherence."* This is exactly what I found by direct computation of the EMOD JSON back in §7 — male effective/non-suppressive tables are byte-identical in the `<6mo` duration bin — confirming it's a deliberate documented modeling choice (mirroring `Days_To_Achieve_Viral_Suppression: 183.0` in the EMOD campaign JSON itself), not a coincidence or data-entry artifact.

**Now, the actual question (art.md 6.2): is the Weibull's declining-hazard-over-time shape already expressed in how STIsim's on-ART mortality changes with CD4 as CD4 reconstitutes?** Only in *direction*, not in *mechanism*, *shape*, or *magnitude*:

- **Mechanism**: the literature's Weibull decline is a hazard that's an explicit function of *time since ART initiation* (`x` in equation X.X.1) — CD4 doesn't appear in that equation at all. STIsim's decline is instead *entirely mediated through CD4*: `get_art_mortality_hazard()` has no time-since-`ti_art` term by default (`art_death_dur=None`); mortality falls over time only as a side effect of `cd4_increase()`'s logistic CD4-recovery curve pushing agents into lower-mortality `off_art_rate(cd4)` bins. These are two different causal stories that happen to both predict "mortality falls after starting ART" — one directly, one by proxy through a *different* biological process (immune reconstitution) that happens to correlate with it.
- **Shape**: a Weibull hazard with `k=0.32` is a smooth power-law decay. STIsim's implicit decline is the composition of a smooth logistic CD4-growth curve (`art_cd4_growth=0.1`) with a *discrete, piecewise-constant* CD4→mortality lookup (`cd4_death_bins`/`cd4_death_rates`, 6 bins) — the resulting mortality-vs-time trajectory is a step function that jumps down each time an agent's reconstituting CD4 crosses a bin edge, not a smooth curve of any particular family, let alone a Weibull one.
- **Magnitude**: see §12 immediately below — the gap is large and concentrated exactly where you'd expect given the mechanism mismatch (early duration, high initial CD4, where the literature's dedicated early-mortality signal has no analog in a CD4-mediated-only model).

So: no, the current continuous-CD4 mechanism is not a faithful stand-in for either literature model's duration effect. It's a plausible-looking proxy — CD4-reconstituting agents do become less likely to die over time under both — that happens to reuse machinery already needed for the reconstitution curve, but it isn't quantitatively equivalent to the fitted duration effect in either source, and (per §12) understates it substantially in exactly the regime (early duration, decent CD4) both papers identify as the highest-differential window.

## 12. Quantitative reconciliation: current branch vs. EMOD JSON table vs. new literature (art.md step 6.3)

Three numbers matter most for reconciling the branch's `rel_art_mortality_unsupp_m/f` against real data, and none of them agree:

| Source | Non-suppressive/effective mortality ratio |
|---|---|
| Literature (`HRns`, both write-ups, Lee et al. 2017 AIDS journal) | **1.96**, constant — same value for both sexes, all CD4, all durations past 6 months |
| EMOD JSON table, lowest CD4 column | **1.96 (male)**, **~1.94-2.10 (female)** — matches the literature almost exactly |
| EMOD JSON table, mid-to-high CD4 columns | **~11.4 (male)**, **~1.94-2.10 (female, unchanged)** |
| Current branch (`rel_art_mortality_unsupp_m/f`) | **2.8 (male)**, **1.4 (female)** — explicitly flagged in the commit message as "not based on anything empirical" |

**The EMOD JSON table's own internal inconsistency, found by direct computation**: the male non-suppressive/effective ratio is *exactly* `1.96` at the lowest CD4 column (matching the literature precisely) but jumps to `~11.4` at every other CD4 column — verified programmatically, not by eye (`art_lookup_test.py`'s embedded tables: `0.1852/0.0945=1.96` at column 0-1, `0.6048/0.0529=11.43` at column 2, consistently ~11.4 through column 6). The female table, by contrast, stays close to `1.94`-`2.10` across every CD4 column — i.e. **consistent with the literature's flat 1.96 everywhere**, unlike the male table. Since neither write-up describes a CD4-dependent adherence penalty anywhere (Table Xb's CD4/age multipliers are applied identically regardless of VL/adherence status, by construction — they're not indexed by VL stratum at all), the male table's CD4-dependent blow-up to 11.4x has no support in either source document I have access to. Two explanations are plausible: either EMOD's implementers made an additional, undocumented male-specific adjustment not captured in either draft I've read, or there's a data-entry/transcription issue specific to the male non-suppressive table. I can't distinguish between these without either the actual EMOD source code that generated this JSON or a version of the write-up postdating September 2019 — worth asking Monisha directly whether a male-specific high-CD4 adjustment exists in her records.

**Given that, the current branch's `2.8`/`1.4` split is defensible as "closer to the literature than a single flat EMOD-derived multiplier would be" (it's *sex-differentiated*, which the literature's `HRns` explicitly is not, but the female table empirically is) but is not itself literature-derived, and understates the male ratio dramatically relative to what's actually hard-coded in the EMOD JSON's non-lowest-CD4 cells (2.8 vs. up to 11.4). The most defensible number to adopt, if the goal is fidelity to the literature specifically (rather than the EMOD JSON's possibly-anomalous male table): a flat `HRns=1.96` for both sexes, with the caveat that this contradicts the *observed* EMOD JSON male behavior at non-trivial CD4, which itself may or may not be intentional.**

**Absolute magnitude, not just the ratio, also diverges substantially — and does so in a duration/CD4-dependent way that traces directly back to §11's mechanism gap.** Computed both sides for two matched cells (male, effective/suppressed ART):

- **High CD4 (≥500), early (0-5 months on ART)**: literature = `0.2015 × 0.0806 = 0.01624`/yr (Table Xa baseline × Table Xb CD4 multiplier). Current branch = `off_art_rate(cd4≥500)=0.003 × rel_art_mortality_effective(0.25) × age_mult(25-35: 1.10) = 0.000825`/yr. **Literature is ~20x higher** — because the current branch has no duration term at all by default, and its CD4-scaling borrows the *off-ART* table's shape (which bottoms out at a low `0.003`/yr ceiling) rather than the literature's own steeper, duration-and-CD4-interacted on-ART gradient.
- **Low CD4 (<50), late (≥36 months on ART)**: literature = `0.011909`/yr (this CD4 band *is* Table Xa's own reference baseline, no Xb multiplier needed). Current branch = `off_art_rate(cd4<50)=0.05 × 0.25 × 1.10 = 0.01375`/yr. **These are close** — within 15% of each other.

So the branch's fidelity to the literature isn't uniformly bad — it's *specifically* worst in exactly the region §11 identified mechanistically: early duration and higher CD4, where the literature has a dedicated, large, well-evidenced early-mortality signal (partly Weibull-derived, partly clinical-consensus-derived per Monisha's email: "mortality is highest within the first 3-6 months of ART initiation... more pronounced at lower CD4 counts") that the current CD4-reconstitution-only mechanism has no way to reproduce, by construction. Reasonable next step (art.md step 7): re-derive `rel_art_mortality_effective`/`rel_art_mortality_unsupp_m`/`_f` and `art_death_age` directly from Table Xa/Xb rather than from the EMOD JSON's possibly-anomalous male column, and separately decide whether to re-enable `art_death_dur` (already wired, currently unused) with duration multipliers derived from Table Xa's own 5-bin structure — which would close most of the ~20x early/high-CD4 gap directly, since that's precisely the axis the current default is missing.

## 13. Is the current on-ART mortality process history-dependent, and does it resemble the EMOD/literature Weibull?

**Short answer: yes it's history-dependent (not exponential/memoryless), but its shape doesn't resemble the Weibull's.** `get_art_mortality_hazard()` recomputes the hazard fresh every timestep from "current state" (age, sex, current CD4), which looks Markovian/memoryless at first glance. But CD4 isn't a freely-fluctuating state — it's a **deterministic function of τ = time since `ti_art`**, anchored to two values (`cd4_preart`, `cd4_potential`) fixed once at ART initiation (`cd4_increase()`, `hiv.py:309-330`):

```
cd4(τ) = cd4_preart + (cd4_potential − cd4_preart) · tanh(0.05·τ)          [τ in months; 0.05 = art_cd4_growth/2]
```

Since `off_art_rate(cd4)` is just a lookup on that, the hazard for a *given individual* is a genuine, well-defined function of τ alone: `h(τ) = rate_bin(cd4(τ)) × rel_art_mortality × age_mult × rel_death_f?`. A hazard that varies with elapsed duration is by definition not exponential — the dependence on "time since ART initiation" is real, it's just smuggled in entirely through CD4 rather than referenced explicitly (`art_death_dur=None` by default, so there's no *additional* explicit duration term on top of this implicit one).

**But the shape is structurally different from the literature's Weibull (`λ=302` days, `k=0.32`, §11), in two concrete ways, verified numerically (see the CSV described below):**

1. **STIsim's hazard hits a hard floor and stays there.** CD4 saturates near `cd4_potential` (a `tanh` curve, bounded) and once it lands in the lowest `cd4_death_rates` bin, the hazard is *exactly* flat forever after — e.g. for a starting CD4 of 500, the hazard is already bit-identical from month 12 onward for the rest of the CSV's 60-month window. The Weibull, with `k=0.32 < 1`, **never plateaus** — it keeps declining indefinitely (approaching, but never reaching, zero) and is mathematically unbounded (infinite instantaneous hazard) at the exact instant of ART initiation, `τ=0` — a real feature of `k<1` Weibull shapes, not a computational artifact (shown as `inf` in the CSV rather than silently truncated).
2. **Different functional family for the decline itself.** STIsim's decline is a smooth `tanh` growth curve composed with a **discrete 6-bin step function** (`cd4_death_bins`) — the actual hazard-vs-τ trajectory is a staircase, flat between CD4-bin crossings, not the continuous power-law `h(τ) ∝ τ^(k−1)` shape a Weibull produces.

### `weibull_hazard_comparison.csv`

Generated by `weibull_hazard_comparison.py`, which pulls STIsim's actual formula live from `sti.HIVPars()` (not hard-coded) and computes the literature Weibull hazard from the transcribed `λ=302, k=0.32, HRns=1.96` parameters (§11). 768 rows: full factorial over

- **age**: `20, 30, 40, 50` — one representative age per `art_death_age` bin (`<25`/`25-35`/`35-45`/`45+`)
- **sex**: `m, f`
- **starting_cd4** (`cd4_preart`, i.e. CD4 *at ART initiation*): `50, 200, 350, 500`
- **adherence**: `effective` / `nonsuppressive`
- **tau_months** (time since ART initiation): `0, 1, 2, 3, 6, 9, 12, 18, 24, 36, 48, 60`

Columns: `cd4_at_tau`/`cd4_potential` (STIsim's reconstituted CD4 trajectory, for inspection), `stisim_hazard_per_year`, `weibull_hazard_per_year`, and `ratio_weibull_over_stisim`.

**How to read it:**

- **The Weibull column is identical across `age`/`sex`/`starting_cd4` for a fixed `(adherence, tau_months)`.** This isn't a bug in the CSV — the literature Weibull model (§11) genuinely has no age, sex, or CD4 covariates; it's purely a function of time-since-ART and suppression status. So every block of 16 rows (4 ages × 2 sexes × 4 starting CD4s) sharing the same `adherence`/`tau_months` has the same `weibull_hazard_per_year` value by construction. This is itself the point: it shows one axis where STIsim's approach is *richer* than the base Weibull (covariate stratification) even as it's *poorer* on the duration axis the Weibull actually specializes in.
- **`ratio_weibull_over_stisim` is large and highly variable** (order 40-2000x across the sampled grid, and `inf` at `τ=0`) — read this as confirming §12's magnitude finding at far more grid points than the two spot-checks there, not as a literal "STIsim should be scaled up by this factor" — the ratio itself varies by orders of magnitude across the grid (largest at high starting CD4 + early τ, smallest at low starting CD4 + late τ), so no single rescaling constant would fix it; the fix has to be structural (a real duration term), per §12's recommendation.
- **Age is held fixed, not aged forward with τ.** A row with `age=20` represents an agent who was 20 at evaluation, for every `tau_months` value in that row's trajectory — it does *not* model an agent who started ART at 20 and is now 20+τ/12 years old. This isolates each covariate's effect independently (clean comparison design) at the cost of not being a literal single-agent trajectory; if you want the latter, add `age + tau_months/12` to the `age_mult` lookup for that row.
  
### One last comment about how EMOD's mortality seems to be consistently higher

Prompt: "From the table it appears that the mortality hazard in STISim is consistently much lower than the mortality hazards in EMOD. Do I understand this correctly? Is this because everything is bounded by the off-art mortality which is also somewhat low?"

Confirmed on both counts — and the second part goes further than "somewhat low," it's structurally impossible for the current design to close the gap in the early months, regardless of what multiplier you'd choose.

Yes, you're reading it correctly. Across all 704 finite-ratio rows in the CSV (everything except τ=0), STIsim's hazard is lower than the Weibull's, every single time — smallest gap is 17.5x, and it goes up from there (into the hundreds/thousands at high starting CD4 + early τ). Zero rows where STIsim ≥ Weibull.

And yes, the off-ART table is the reason — but it's a harder ceiling than "somewhat low" suggests. Two things compound:

The invariant we built forces on-ART ≤ off-ART, by construction. rel_art_mortality_effective=0.25 and rel_art_mortality_unsupp_m/f=0.7/0.35 are fractions less than 1, specifically chosen (§9-10 of this doc) so on-ART mortality can never exceed the off-ART rate at the same CD4. That's deliberate — it's the fix for the earlier bug where non-suppressive ART looked worse than no treatment at all.

But even the off-ART table's own ceiling — 0.30/yr, its rate for the sickest possible CD4 (near zero) — is lower than what the Weibull says early on-ART mortality should be. Concretely:

τ (months)	Weibull hazard/yr	Off-ART table's worst possible rate
1	1.84	0.30
3	0.87	0.30
6	0.55	0.30
9	0.41	0.30
12	0.34	0.30
For roughly the first 9 months, the Weibull's predicted hazard exceeds even the off-ART table's sickest-CD4 rate — the ceiling the on-ART formula is anchored under. That means no choice of rel_art_mortality_effective/rel_art_mortality_unsupp_m/f (any fraction ≤1) could ever reproduce the Weibull's early spike, no matter how you tune those numbers — the off-ART table itself doesn't reach high enough, structurally, regardless of multiplier.

That's a stronger statement than "the current numbers are miscalibrated" — it means the invariant and the Weibull's early-mortality signal are fundamentally incompatible during the first ~9 months on ART, since they were never trying to describe the same thing: the off-ART table describes ongoing untreated natural history (calibrated independently, predates this whole effort), while the Weibull describes a specific, acute, transient window (recent initiation — possibly IRIS, unresolved opportunistic infections, not-yet-effective treatment) that has no analog in "what's your CD4 right now," even at CD4≈0. Closing this gap for real would need either relaxing the ≤off-ART invariant specifically for an early window, or accepting that the invariant only makes sense past some minimum duration on treatment (e.g. reintroducing art_death_dur as an unbounded multiplier applied before the ceiling check, rather than folding duration entirely through the bounded CD4 channel as it is now).