# ART Implementation Notes (detailed companion)

Companion to `art_implementation_notes.md` — same section order, more detail and code references. All file:line references are in `/Users/daniel/Documents/IDM/stisim` unless noted.

## How ART Distribution Works in HIVSim

1. **Testing**: `HIVTest` schedules *when* someone starts ART — sets `hiv.ti_art` on diagnosis
2. **Diagnosis**: positive test sets `hiv.diagnosed=True`, `hiv.ti_diagnosed=ti`
3. **Scheduling**: delay from diagnosis to treatment start (`dur_dx2tx`, default 0)
4. **`ART` intervention** decides *whether* they actually initiate, each step:
   - `ready = (ti_art == ti) & ~on_art`, filtered through `art_initiation` (default 90%, flat — not age/sex/CD4-dependent)
   - No coverage target set: everyone who passes `art_initiation` starts, uncapped
   - Coverage target set: `prioritize_art()` fills only the available slots, prioritizing by CD4 count and care-seeking
   - `art_coverage_correction` pushes the on-ART population toward the target every step — can also actively remove people early if the target shrinks

`HIVTest` must appear before `ART` in the interventions list, or newly-diagnosed agents won't start ART until the next step.

## Effects of ART on Individuals

- **CD4**: reconstitutes via a logistic growth curve while on ART (`cd4_increase`, `hiv.py:309-330`); declines linearly after stopping (`post_art_decline`)
- **Mortality**: see "Current ART Mortality Implementation" below
- **Transmission**: see "Effect of ART on Transmission" below
- **Retention/duration**: `dur_on_art ~ LogNormal(3yr, 1.5yr) × rel_dur_on_art`, drawn once at initiation, no age/sex dependence. Can also end early via `art_coverage_correction` if coverage targets tighten (rare in practice — care-seeking behavior makes already-on-ART agents unlikely to be the ones removed)
- **Time to full ART efficacy**: 6 months (`time_to_art_efficacy`), linear ramp from 0

## ART Adherence States

Four `BoolState`s, not a single enum — matches the existing codebase idiom (e.g. `acute`/`latent`/`falling`), and each gets a free `n_<name>` result:

- `art_naive` (default `True`)
- `on_effective_art`
- `on_nonsuppressive_art`
- `art_discontinued`

`on_art` = `on_effective_art | on_nonsuppressive_art`. Split at initiation via `p_effective_art` (`ss.bernoulli`, default: always effective). `art_naive` only clears for agents who were actually naive before that call, so re-initiation after discontinuation doesn't misfire naive-only logic elsewhere.

**Known gotcha**: `ART.step()` has two call sites for `hiv.start_art()`. Only the coverage-target path forwards `p_effective_art` correctly. The no-coverage-target path (`hiv_interventions.py:353`) omits it, so `sti.ART(p_effective_art=0.5)` with no `coverage=` set silently produces zero non-suppressive ART. Not yet fixed.

## Effect of ART on Transmission

`rel_trans` ramps linearly (over 6 months from `ti_art`) down to `1 - efficacy`:

- `effective_art_efficacy = 0.96` (96% reduction — suppressed, "U=U")
- `nonsupp_art_efficacy = 0.35` (35% reduction — not suppressed, still substantially infectious: ~16x more transmissible than an effective-ART agent once both are fully ramped)

## EMOD's `ARTMortalityTable` — what it did differently

EMOD's original table drew a single survival time *once*, at ART initiation (not re-evaluated afterward), stratified by age (4 bins) × CD4-at-initiation (7 bins, frozen) × time-since-initiation (5 bins) × sex (routed via separate campaign events) × adherence (separate tables for effective/non-suppressive).

Compared to HIVsim's mortality *before this branch*:

1. Non-suppressive ART was far deadlier than effective ART, same age/sex/CD4 — HIVsim had zero difference
2. Mortality dropped steeply the longer someone had been on ART — HIVsim had no duration axis at all
3. Mortality rose mildly with age, same sex/CD4/adherence — HIVsim had no age dependence
4. Excess mortality while on ART was non-zero — HIVsim's on-ART mortality was exactly zero, always

EMOD's own JSON table turned out to be internally inconsistent on point 1 once actually computed — see "Numerical mismatch vs. literature" below.

## Current ART Mortality Implementation

Branch: `pr-561-simplify-art-mortality` (local) / `feat/simplify-art-mortality` (remote).

```
rate = off_art_rate(cd4) × rel_art_mortality[effective?] × age_mult(age) × rel_death_f?
```

- `off_art_rate(cd4)`: the *same* CD4-stratified table used for off-ART mortality (`cd4_death_bins`/`cd4_death_rates`), looked up at the agent's **current** (not frozen-at-initiation) CD4
- `rel_art_mortality_effective = 0.25` — fraction of the off-ART rate retained on effective ART, both sexes
- `rel_art_mortality_unsupp_m = 0.7` / `rel_art_mortality_unsupp_f = 0.35` — non-suppressive fraction, sex-specific. Ratio to `rel_art_mortality_effective`: 2.8x for men, 1.4x for women (men's ratio set to exactly 2x women's, by design — not fit to data)
- `age_mult`: 1.0 (<25) / 1.10 (25-35) / 1.21 (35-45) / 1.32 (45+)
- `rel_death_f = 0.74` — flat additional multiplier for females, applied equally to both adherence categories (so it doesn't affect the 2.8/1.4 ratio above — it's a separate sex offset)

Recalculated fresh every `step_state()` from current CD4. Since CD4 rises over time on ART (`cd4_increase`), mortality implicitly declines over time too — but only via this CD4 channel; there's no explicit duration term (`art_death_dur=None` by default — wired in the function but unused).

**Guarantee, by construction**: on-ART mortality can never exceed off-ART mortality at the same CD4, since every `rel_art_mortality_*` value is ≤1 and `rel_art_mortality_unsupp_m × max(age_mult) = 0.7 × 1.32 = 0.924 ≤ 1`. This was a deliberate fix — the function used to compute on-ART mortality from an independent baseline + CD4 table, unrelated to the off-ART table, which let non-suppressive (and sometimes effective) ART mortality exceed off-ART mortality at high CD4.

**Is this history-dependent (not exponential)?** Yes. CD4 isn't a free-floating state — it's a deterministic function of τ = time since `ti_art`: `cd4(τ) = cd4_preart + (cd4_potential − cd4_preart) · tanh(0.05·τ)`. So the hazard genuinely varies with duration-on-ART, just implicitly, through CD4, rather than via an explicit term. It does *not* resemble EMOD/literature's Weibull shape, though — it hits a hard floor and stays flat once CD4 saturates, where a Weibull (`k=0.32<1`) keeps declining indefinitely and never plateaus.

## What May Need to Be Fixed

### Numerical values are much lower than EMOD/literature

`weibull_hazard_comparison.csv` (`weibull_hazard_comparison.py` to regenerate) compares current HIVsim mortality against the literature Weibull (`λ=302` days, `k=0.32`, `HR=1.96` for non-suppressed — from Monisha Sharma's source material in `/Users/daniel/Documents/IDM/ART_Mortality_Table_Docs/`). Current HIVsim is lower everywhere in the comparison grid (17.5x-2000x lower, depending on age/CD4/duration). For roughly the first 9 months on ART, the literature's hazard exceeds even the off-ART table's own worst-case rate (0.30/yr) — meaning **no choice of `rel_art_mortality_*` can close this gap**; it's structural, not a calibration problem. The off-ART table describes ongoing untreated natural history; the literature Weibull describes a specific, acute, transient early-ART risk window that has no analog in "what's your CD4 right now," even at CD4≈0.

`rel_art_mortality_unsupp_m/f` ratio (2.8/1.4, i.e. men's is 2x women's) is not empirically derived — flagged as a placeholder in the commit that introduced it. For comparison: EMOD's actual JSON table gives ~11.4x for men at most CD4 levels, but only ~1.96x at the lowest CD4 column — and that 1.96x matches Monisha's source literature (`HRns=1.96`) exactly. EMOD's female table stays close to ~2.0x across all CD4 levels, also matching the literature. The male 11.4x has no support in the literature reviewed so far — possibly an undocumented EMOD-specific adjustment, possibly a data issue in that one table. Either way, it can't be reproduced under the current off-ART-bounded design without breaking the guarantee above (11.4x at high CD4 pushes on-ART mortality above off-ART).

### Check: does CD4 increase faster on effective vs. non-suppressive ART?

Not currently modeled — `cd4_increase()` doesn't condition on adherence category at all, only on `cd4_preart`/`cd4_potential`. Worth checking whether this is realistic (clinically, poor adherence is associated with slower/incomplete CD4 recovery, not just a mortality penalty).

### Terminal "failing window" behavior (not implemented)

EMOD models a window of time before death where an agent is nominally still on effective ART but has entered `ON_BUT_FAILING` — transmission suppression turns off (10x infectious, matching late AIDS) even though they're still counted as on effective ART, in the ~9 months before their scheduled death. Not implemented here; would need a third adherence sub-state with its own transmission rule and entry/exit timing. Out of scope unless specifically requested.
