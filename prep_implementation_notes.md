# Notes for how to implement changes to PrEP

## Current implementation of PrEP

* How is PrEP currently implemented? When someone receives PrEP, what happens? Where is the information about the intervention stored (ie part of the HIV or the PrEP module or somewhere else) and what information is stored?

* What are the differences between the current implementation and the desired changes listed below?

## Desired changes to PrEP

1. We will want to be able to specify multiple varieties of PrEP - oral PrEP, or long-acting cabotegravir or lenacapivir. For this reason, I want to be able to specify the duration and efficacy of the intervention
2. We will also want to specify adherence, especially for oral PrEP, since variable adherence reduces the efficacy
3. This means we want to be able to specify PrEP with three interventions
4. I would also like to add a flag which specifies whether or not someone is currently taking PrEP (ie, nobody double-takes prep or takes multiple kinds of prep) and also a flag which says oral or long-acting prep