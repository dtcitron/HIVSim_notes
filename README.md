# HIVSim_notes

I will use this repository to keep track of notebooks and scripts which I use as I start to use HIVSim

* stisim_tutorial.ipynb - [STISim Tutorials](https://docs.stisim.org/tutorials/tut_intro.html)
* Testing_new_ART_transmission_and_mortality.ipynb - ART transmission, mortality, duration spot-checks

## Changes made to HIVSim

* [x] Add EffectiveART/NonSuppressiveART tracking, to enable 90-90-90 metrics tracking
* [x] Allow transmission to vary by EffectiveART/NonSuppressiveART status
* [x] Allow mortality to vary by EffectiveART/NonSuppressiveART status, as well as age, sex, and duration on ART, qualitatively following ARTMortalityTable from EMOD-HIV
* [x] Adopting a calculated art mortality that doesn't quite use the ARTMortality Table hard-coded in.
  * [x] Verified that mortality while on ART is now bounded by mortality while off ART for all demographic groups

## Future change ideas

* [ ] Make duration on ART adjustable, so we can adjust it by age
* [ ] Fully adopt ARTMortalityTable, which is based on Weibull-distribution fits to survival curves while on ART