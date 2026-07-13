# Planning New ART Implementation

I have just checked out a new branch called ART. The purpose of this branch will be changing the way ART is implemented in HIVSim and STISim.

When we implement ART, it should do a few things

1. Reduce mortality of the individuals who are taking it
   1. Mortality should depend on age, sex, and time since ART initiation
   2. Mortality should also depend on whether individual adheres to ART (ie, is on effective art vs non suppressive art)
2. Reduce transmission of the individuals who are taking it (undetectable = untransmissible)
   1. Transmission should depend on whether individual adheres to ART

Currently, ART is implemented in a way that does reduce mortality but does not do any of the other things.

# Claude instructions:

1. Give me a new document in HIVSim_notes which explains the following things
   1. How does ART distribution work in HIVSim?
   2. What are the current effects of ART on individuals who take it?
   3. What processes affect ART retention time?
   4. How would I implement a new "individual property" which tracks ART adherence and whether or not someone has achieved viral load suppression using STIsim (is it even called an individual property in stisim, or something else)
2. Look through the EMOD repo, locate ARTMortalityTable intervention, and understand how it works and what it does - the next steps will involve
3. Implement changes to ART:
   1. Add a new property which tracks ART adherence - whether or not an individual is currently taking effective or nonsuppressive ART
      1. new property should have possible four states: ARTNaive, OnEffectiveART, OnNonSuppressiveART, ARTDiscontinued
      2. new property should have default to artnaive when individual is initiated
      3. In hiv.py, in the care and treatment states, art_naive replaces never_art; on_art should stay there and turn on and off according to art status. on_effective_art and on_nonsuppressive art can be there in addition to on_art
   2. When ART is initiated, need the function for distributing ART to take as an argument the probability of receiving EffectiveART or NonSuppressiveART
      1. default to always giving effectiveart `P_EffectiveART = 1`
   3. Need now to add two more attributes similar to `self.pars.art_efficacy`
      1. `self.pars.effective_art_efficacy` = 0.96
      2. `self.pars.nonsupp_art_efficacy` = 0.35
      3. And change it such that `self.pars.art_efficacy` gets replaced by effective or nonsupp art efficacy everywhere, based on a check of whether someone is on effectiveart or nonsuppressive art
   4. Write me a test which distributes 50% effectiveart and 50% nonsuppressiveart to a total of 1000 people in an HIVSim simulation with 10k people in it. 
      1. Give me scripts for catching the output and plotting the total number of people on either type of ART.
      2. give me scripts for catching the output and plotting the mortality of different people depending on art status - art_naive, art_discontinued, on_effective, on_nonsuppressive, and on_art
      3. put the test in HIVSim_notes/art_transmission_test.py


Some future changes that will need to be made eventually

   5. Change ART distribution to be able label individuals in different age categories with whether they are on effectiveart or onnonsuppressiveart
   6. Change effects of ART on mortality according to oneffectiveart or onnonsuppressiveart, following what ARTMortalityTable does


Some future checks that will need to be made eventually

* [ ] What is the typical duration of time spent on ART (that is to say, retention), and is there a way to let this vary by age?