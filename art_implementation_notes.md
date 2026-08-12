# ART Implementation Notes

## How ART Distribution Works in HIVSim

1. Testing: `HIVTest` schedules when someone starts ART
2. Diagnosis
3. Scheduling
4. `ART` decides whether they actually initiate ART

## Previous effects of ART on individuals

1. CD4 count: Each step, on-ART agents get CD4 reconstitution via a logistic curve (`cd4_increase`, `hiv.py:273-295`)
2. Mortality - HIV-related mortality is zero for people on ART
3. Adherence - no differentiation between adherence/non-adherence, or VLS/Non-VLS
4. Transmission - reduced by 96% while on ART (`rel_trans = 0.04` regardless of adherence)
5. Retention/duration - duration while on ART drawn from lognormal distribution, 
   1. Also possible people are dropped from ART as the distribution of ART and ART coverage are shuffled around while trying to match distribution targets. (Unlikely, due to care-seeking behavior from individuals already on ART, likely to continue being on ART.)
6. Time to full ART efficacy = 6 monhts

## Implementing a 4-state ART-adherence property in STISim

1. ART Naive
2. On Effective ART
3. On Nonsuppressive ART
4. ART Discontinued

## Effect of ART on Transmission

`effective_art_efficacy=0.96` and `nonsupp_art_efficacy=0.35`

## How to upgrade mortality while on ART

Mortality while on ART should depend on the following things:

1. Age
2. Sex
3. CD4 count
4. ART adherence/non-adherence (ie, virally suppressed or no)

EMOD's ARTMortality table also incorporates a "time since initiating ART" to deal with the history dependence of mortality.

How is it different?

1. Non-suppressive ART is far deadlier than effective ART, for same age/sex/CD4 stratum
2. Mortality drops steeply the longer someone has been on ART
3. Mortality rises mildly with age, for same sex/CD4/adherence stratum
4. Excess mortality while on ART is now non-zero

## Current ART Mortality Implementation

ART mortality rate is bounded from above by the `off_art_rate`, stratified by CD4 count. The mortality is then modified by effective/nonsuppressive ART status, age, and gender. Mortality while on effective ART is always lower than on nonsuppressive ART; mortality while on nonsuppressive ART is always lower than mortality off art for all strata.

Mortality while on nonsuppressive ART is higher among men than among women.

`rate = off_art_rate(cd4) × rel_art_mortality[effective?] × age_mult(age) × rel_death_f?`

ART mortality is re-calculated at each time step (`step_state()`) based on CD4 count (remember, CD4 count is steadily increasing for individuals who are on ART over time).

## What may need to be fixed in the future about ART Mortality

### Check: does CD4 count increase faster for effective ART than nonsuppressive ART?

### Check: numerical values for mortality while on ART

Comparing the mortality hazard rates for EMOD-HIV's ARTMortality Table - see `weibull_hazard_comparison.csv` - it appears that mortality in EMOD was elevated compared to HIVSim. Given time, we will need to revisit this, and reconcile the numerical values, if model calibration proves too difficult with lower mortality rates.

Currently `rel_art_mortality_unsupp_m/f` ratio is 2, which is not based on anything empirical. In original ARTMortalityTable this was closer to 11.4, but that cannot be reconciled with the upper bound mortality while off ART.

### Terminal "failing window" behavior

EMOD includes a window of time, pre-death, when they are in an `ON_BUT_FAILING` ART state. Nominally they may still be on ART, but transmission blocking is off. This is meant to represent ART failure.
