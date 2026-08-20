# Notes for how to implement changes to PrEP

## Current implementation of PrEP

* How is PrEP currently implemented? When someone receives PrEP, what happens? Where is the information about the intervention stored (ie part of the HIV or the PrEP module or somewhere else) and what information is stored?

* What are the differences between the current implementation and the desired changes listed below?

## Desired changes to PrEP

1. I would like PrEP coverage to work the same way as ART and VMMC, where the user can specify numbers of people, fractions of people, and have those things change over time.
2. Just as with ART, I would like the following states added to the PrEP module to track people with a history of taking PrEP; initiate everyone with prep_naive = True, when someone receives ART they go prep_naive = False and on_prep = True, when prep is discontinued they go prep_discontinued = True
   1. prep_naive
   2. on_prep
   3. prep_discontinued
3. We will want to be able to specify multiple varieties of PrEP - oral PrEP, or long-acting cabotegravir or lenacapivir. For this reason, I want to be able to specify the duration and efficacy of the intervention
4. We will also want to specify adherence, especially for oral PrEP, since variable adherence reduces the efficacy
5. This means we want to be able to parameterize PrEP with 3 parameters - efficacy, adherence, and duration
6. Make sure that prep duration corresponds to the time when prep is due to expire, and that prep expires at the correct time.
7. Give me a script for test-driving these new changes to PrEP. 
   1. I would like to run a simulation with 10k people and no prep; a scenario where we put 10% of women between the ages of 15 and 35 on oral PrEP (lasts 3 months, 85% efficacy, 75% adherence) at the start of every year; a scenario where we put 500 men and 1000 women between the ages of 15 and 35 on long acting prep (lasts 12 months, 100% efficacy, 100% adherence)
   2. I would like plots which show prevalence and the number of new infections with each of the scenarios over time
   3. I also want to run a simulation where 100 people are put on prep for 12 months, then it expires, and I want the simulation to show the discrete duration of PrEP - that it expires on time