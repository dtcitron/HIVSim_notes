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
