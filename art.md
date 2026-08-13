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
4.  Second round of changes to ART model within HIVSim: ART Mortality variation
    1.  First parse the sample ARTMortalityTable intervention excerpted from an EMOD-HIV campaign file below, and summarize for me how different this approach to mortality calculation while on ART is from what is currently in HIVSim
    2.  Implement ARTMortalityTable-style (more complex) differential art mortality
        1.  As in EMOD ARTMortalityTable, continuously compute new cd4-stratified hazard for people who are on ART based on continuous (not frozen at initiation). Use the lookup table shared in this document as a model for how to compute the hazard I need.
        2.  Mortality while on ART should depend on Age, sex, EffectiveART/NonsuppressiveART, and duration since initiation, as in ARTMortality tables below
        3.  Allow for nonzero mortality while on ART - this is especially important during the early months on ART for people who have low CD4 count
        4.  I would like the artmortalitytable numbers to be the default, but also make it possible for a user to input their own tables if they want to change this.
        5.  Ignore the terminal failing window behavior for now; will consult with team
    3.  Allow the user to choose whether or not to implement the simpler version with no mortality variance or dependency while on ART or ARTMortalityTable when constructing the model
    4.  Need you to add a new third diagram to the docs under: user_guide/diseases/hiv.html which reflect the new changes to ART mortality
    5.  Write me a test which allows me to compare the old and new versions of ARTMortality across different. Add this to the art_adherence_test script

5. ARTMortalityTable- second round: Preferred style would be to not include the hard-coded ARTmortalityTable used in EMOD, and instead to include a single function. I am going to work with some comments from colleague Monisha Sharma to implement this a little more cleanly
   1. Need to adjust the ratio of mortality for onnonsuppressiveART/oneffectiveART such that it is different for men and for women. For men it should be 11.4, for women it should be 2.3.
6. ARTMortalityTable - third round: The original ARTMortalityTable was based on Weibull curves - it remains unclear to me exactly where those tables came from, but I do have a lot of material from Monisha sharma about the sources for everything, including some citations which reference survival models fit to weibull curves;
   1. Ingest the documents in this folder: /Users/daniel/Documents/IDM/ART_Mortality_Table_Docs/
      1. The artmorttable_abstracts.md contains an email describing where artmortality table came from as well as the abstracts from associated papers used to build artmortality table
      2. This folder also contains three supporting manuscripts and some supplemental documents for building artmortality table.
   2. Give me a written answer in art_implementation_notes.md describing in detail how to reconcile the weibull distribution in the literature against the current way of representing/implementing mortality on ART. Is this something already expressed in how mortailty rates shift with cd4 counts (which are changing over time)?
   3. Also give me a written answer describing in detail whether the numerical parameters for calculating artmortality implemented in the current branch, the hard-coded ARTmortalityTable from EMOD, and this new background material can be reconciled quantitatively.
7. ARTMortalityTable - third round, implementation
   1. for the current version, re-implement the scaling factors for oneffectiveart and on nonsuppressive art using quantitative support from this literature - describe in detail for me in art_implementation_notes how you propose to incorporate the new information in recalculating the mortality values
   2. Implement the weibull. By the end of this I should have 3 different versions of ART mortality:
      1. the current version, which has exponential mortality while on art and  enforces off art mortality always above on art mortality (commit 3d1b25603643524de007c34a8f8e2212275c5a9e (HEAD -> pr-561-simplify-art-mortality)) 
      2. updated version which has exponential mortality while on art, accounting for the data ingested
      3. updated version which has weibull mortality

8.  Third round of changes to ART model within HIVSim: ART retention duration
    1.  Currently for EMOD we have about 15% annual dropout for AYA and 5% annual dropout for adults older than 25 years of age
    2.  Unclear to me exactly how long people remain on ART for on average in a given simulation - give me a python script called art_adherence_test.py which lets me estimate the distribution of time spent on ART
    3.  Give me a plan for how we might specify average time spent on ART depending on demographics (age or sex or other property)

## Sample EMOD-HIV campaign file code which we use for specifying ART mortality depending on Age, Sex, ART Adherence, and time spent on ART

```       {
            "Event_Name": "Initiate effective ART: MALES",
            "class": "CampaignEventByYear",
            "Nodeset_Config": { "class": "NodeSetAll" },
            "Start_Year": 1980,
            "Event_Coordinator_Config": {
                "class": "StandardInterventionDistributionEventCoordinator",
                "Demographic_Coverage": 1,
                "Intervention_Config": {
                    "class": "NodeLevelHealthTriggeredIV",
                    "Trigger_Condition_List": [ "EffectiveART" ],
                    "Target_Demographic": "ExplicitGender",
                    "Target_Gender": "Male",
                    "Actual_IndividualIntervention_Config": {
                        "class": "ARTMortalityTable",
                        "Intervention_Name": "SuppressedOnART",
                        "Cost_To_Consumer" : 1,
                        "ART_Multiplier_On_Transmission_Prob_Per_Act"             : 0.04,
                        "ART_Is_Active_Against_Mortality_And_Transmission"        : 1,
                        "Days_To_Achieve_Viral_Suppression"                       : 183.0,
                        "ART_Duration_Days_Bins": [182, 365, 730, 1095, 45625],
                        "Age_Years_Bins": [25, 35, 45, 125],
                        "CD4_Count_Bins": [0, 25, 74.5, 149.5, 274.5, 424.5, 624.5],
                        "MortalityTable": [
                            [
                                [0.2015, 0.2015, 0.1128, 0.0625, 0.0312, 0.0206, 0.0162],
                                [0.2176, 0.2176, 0.1219, 0.0675, 0.0337, 0.0223, 0.0175],
                                [0.2350, 0.2350, 0.1316, 0.0729, 0.0364, 0.0240, 0.0189],
                                [0.2538, 0.2538, 0.1421, 0.0787, 0.0393, 0.0260, 0.0205]
                            ],
                            [
                                [0.0875, 0.0875, 0.0490, 0.0271, 0.0136, 0.0062, 0.0041],
                                [0.0945, 0.0945, 0.0529, 0.0293, 0.0146, 0.0067, 0.0044],
                                [0.1021, 0.1021, 0.0572, 0.0316, 0.0158, 0.0073, 0.0047],
                                [0.1102, 0.1102, 0.0617, 0.0342, 0.0171, 0.0079, 0.0051]
                            ],
                            [
                                [0.0255, 0.0255, 0.0181, 0.0128, 0.0085, 0.0058, 0.0038],
                                [0.0288, 0.0288, 0.0204, 0.0145, 0.0096, 0.0065, 0.0043],
                                [0.0326, 0.0326, 0.0231, 0.0164, 0.0108, 0.0074, 0.0049],
                                [0.0368, 0.0368, 0.0261, 0.0185, 0.0123, 0.0084, 0.0055]
                            ],
                            [
                                [0.0164, 0.0164, 0.0116, 0.0083, 0.0055, 0.0037, 0.0025],
                                [0.0186, 0.0186, 0.0131, 0.0093, 0.0062, 0.0042, 0.0042],
                                [0.0210, 0.0210, 0.0148, 0.0106, 0.0070, 0.0048, 0.0048],
                                [0.0237, 0.0237, 0.0168, 0.0119, 0.0079, 0.0054, 0.0054]
                            ],
                            [
                                [0.0119, 0.0119, 0.0081, 0.0066, 0.0033, 0.0033, 0.0033],
                                [0.0135, 0.0135, 0.0092, 0.0074, 0.0037, 0.0037, 0.0037],
                                [0.0152, 0.0152, 0.0103, 0.0084, 0.0042, 0.0042, 0.0042],
                                [0.0172, 0.0172, 0.0117, 0.0095, 0.0047, 0.0047, 0.0047]
    
                            ]
                        ]
                    }
                }
            }
        },
        {
            "Event_Name": "Initiate non-suppressive ART: MALES",
            "class": "CampaignEventByYear",
            "Nodeset_Config": { "class": "NodeSetAll" },
            "Start_Year": 1980,
            "Event_Coordinator_Config": {
                "class": "StandardInterventionDistributionEventCoordinator",
                "Demographic_Coverage": 1,
                "Intervention_Config": {
                    "class": "NodeLevelHealthTriggeredIV",
                    "Trigger_Condition_List": [ 
                        "NonSuppressiveART", 
                        "NonSuppressiveARTDepressed" 
                    ],
                    "Target_Demographic": "ExplicitGender",
                    "Target_Gender": "Male",
                    "Actual_IndividualIntervention_Config": {
                        "class": "ARTMortalityTable",
                        "Intervention_Name": "NotSuppressedOnART",
                        "Cost_To_Consumer" : 1,
                        "ART_Multiplier_On_Transmission_Prob_Per_Act"             : 0.65,
                        "ART_Is_Active_Against_Mortality_And_Transmission"        : 1,
                        "Days_To_Achieve_Viral_Suppression"                       : 183.0,
                        "ART_Duration_Days_Bins": [182, 365, 730, 1095, 45625],
                        "Age_Years_Bins": [25, 35, 45, 125],
                        "CD4_Count_Bins": [0, 25, 74.5, 149.5, 274.5, 424.5, 624.5],
                        "MortalityTable": [
                            [
                                [0.2015, 0.2015, 0.1128, 0.0625, 0.0312, 0.0206, 0.0162],
                                [0.2176, 0.2176, 0.1219, 0.0675, 0.0337, 0.0223, 0.0175],
                                [0.2350, 0.2350, 0.1316, 0.0729, 0.0364, 0.0240, 0.0189],
                                [0.2538, 0.2538, 0.1421, 0.0787, 0.0393, 0.0260, 0.0205]
                            ],
                            [
                                [ 0.1715, 0.1715, 0.5600, 0.3100, 0.1550, 0.0713, 0.0465 ],
                                [ 0.1852, 0.1852, 0.6048, 0.3348, 0.1674, 0.0770, 0.0502 ],
                                [ 0.2000, 0.2000, 0.6532, 0.3616, 0.1808, 0.0832, 0.0542 ],
                                [ 0.2160, 0.2160, 0.7054, 0.3905, 0.1953, 0.0898, 0.0586 ]
                            ],
                            [
                                [ 0.0532, 0.0532, 0.0362, 0.0293, 0.0171, 0.0116, 0.0095 ],
                                [ 0.0601, 0.0601, 0.0409, 0.0331, 0.0193, 0.0131, 0.0107 ],
                                [ 0.0679, 0.0679, 0.0462, 0.0374, 0.0218, 0.0148, 0.0121 ],
                                [ 0.0768, 0.0768 ,0.0522, 0.0422, 0.0246, 0.0168, 0.0137 ]
                            ],
                            [
                                [0.0335, 0.0335, 0.0228, 0.0184, 0.0108, 0.0073, 0.0060 ],
                                [0.0379, 0.0379, 0.0258, 0.0208, 0.0122, 0.0083, 0.0068 ],
                                [0.0428, 0.0428, 0.0291, 0.0235, 0.0137, 0.0094, 0.0076 ],
                                [0.0484, 0.0484, 0.0329, 0.0266, 0.0155, 0.0106 ,0.0086 ]
                            ],
                            [
                                [0.0234, 0.0234, 0.0159, 0.0129, 0.0091, 0.0069, 0.0064 ],
                                [0.0265, 0.0265, 0.0180, 0.0145, 0.0103, 0.0077, 0.0073 ],
                                [0.0299, 0.0299, 0.0203, 0.0164, 0.0116, 0.0088, 0.0082 ],
                                [0.0338, 0.0338, 0.0230, 0.0186, 0.0131, 0.0099, 0.0093 ]
    
                            ]
                        ]
                    }
                }
            }
        },
        {
            "Event_Name": "Initiate effective ART: FEMALES",
            "class": "CampaignEventByYear",
            "Nodeset_Config": { "class": "NodeSetAll" },
            "Start_Year": 1980,
            "Event_Coordinator_Config": {
                "class": "StandardInterventionDistributionEventCoordinator",
                "Demographic_Coverage": 1,
                "Intervention_Config": {
                    "class": "NodeLevelHealthTriggeredIV",
                    "Trigger_Condition_List": [ "EffectiveART" ],
                    "Target_Demographic": "ExplicitGender",
                    "Target_Gender": "Female",
                    "Actual_IndividualIntervention_Config": {
                        "class": "ARTMortalityTable",
                        "Intervention_Name": "SuppressedOnART",
                        "Cost_To_Consumer" : 1,
                        "ART_Multiplier_On_Transmission_Prob_Per_Act"             : 0.04,
                        "ART_Is_Active_Against_Mortality_And_Transmission"        : 1,
                        "Days_To_Achieve_Viral_Suppression"                       : 183.0,
                        "ART_Duration_Days_Bins": [182, 365, 730, 1095, 45625],
                        "Age_Years_Bins": [25, 35, 45, 125],
                        "CD4_Count_Bins": [0, 25, 74.5, 149.5, 274.5, 424.5, 624.5],
                        "MortalityTable": [
                            [
                                    [ 0.2015, 0.2015, 0.0993, 0.0518, 0.0259, 0.0171, 0.0135 ],
                                    [ 0.2156, 0.2156, 0.1062, 0.0554, 0.0277, 0.0183, 0.0144 ],
                                    [ 0.2307, 0.2307, 0.1137, 0.0593, 0.0296, 0.0196, 0.0154 ],
                                    [ 0.2468, 0.2468, 0.1216, 0.0634, 0.0317, 0.0209, 0.0165 ]
                                ],
                                [
                                    [ 0.0875, 0.0875, 0.0431, 0.0225, 0.0112, 0.0052, 0.0034 ],
                                    [ 0.0936, 0.0936, 0.0461, 0.0241, 0.0120, 0.0055, 0.0036 ],
                                    [ 0.1002, 0.1002, 0.0494, 0.0257, 0.0129, 0.0059, 0.0039 ],
                                    [ 0.1072, 0.1072, 0.0528, 0.0276, 0.0138, 0.0063, 0.0041 ]
                                ],
                                [
                                    [0.0241, 0.0241, 0.0166, 0.0135, 0.0067, 0.0044, 0.0044 ],
                                    [0.0262, 0.0262, 0.0181, 0.0147, 0.0073, 0.0048, 0.0048 ],
                                    [0.0286, 0.0286, 0.0197, 0.0160, 0.0080, 0.0052, 0.0052 ],
                                    [0.0312, 0.0312, 0.0215, 0.0175, 0.0087, 0.0057, 0.0057 ]
                                ],
                                [
                                    [0.0149, 0.0149, 0.0103, 0.0084, 0.0042, 0.0042, 0.0042 ],
                                    [0.0163, 0.0163, 0.0112, 0.0091, 0.0046, 0.0046, 0.0046 ],
                                    [0.0177, 0.0177, 0.0122, 0.0099, 0.0050, 0.0050, 0.0050 ],
                                    [0.0193, 0.0193, 0.0133, 0.0108, 0.0054, 0.0054, 0.0054 ]
                                ],
                                [
                                    [0.0084, 0.0084, 0.0057, 0.0046, 0.0023, 0.0023, 0.0023 ],
                                    [0.0092, 0.0092, 0.0062, 0.0051, 0.0025, 0.0025, 0.0025 ],
                                    [0.0100, 0.0100, 0.0068, 0.0055, 0.0028, 0.0028, 0.0028 ],
                                    [0.0109, 0.0109, 0.0074, 0.0060, 0.0030, 0.0030, 0.0030 ]
                                ]
                        ]
                    }
                }
            }
        },
        {
            "Event_Name": "Initiate non-suppressive ART: FEMALES",
            "class": "CampaignEventByYear",
            "Nodeset_Config": { "class": "NodeSetAll" },
            "Start_Year": 1980,
            "Event_Coordinator_Config": {
                "class": "StandardInterventionDistributionEventCoordinator",
                "Demographic_Coverage": 1,
                "Intervention_Config": {
                    "class": "NodeLevelHealthTriggeredIV",
                    "Trigger_Condition_List": [ 
                        "NonSuppressiveART", 
                        "NonSuppressiveARTDepressed" 
                    ],
                    "Target_Demographic": "ExplicitGender",
                    "Target_Gender": "Female",
                    "Actual_IndividualIntervention_Config": {
                        "class": "ARTMortalityTable",
                        "Intervention_Name": "NotSuppressedOnART",
                        "Cost_To_Consumer" : 1,
                        "ART_Multiplier_On_Transmission_Prob_Per_Act"             : 0.65,
                        "ART_Is_Active_Against_Mortality_And_Transmission"        : 1,
                        "Days_To_Achieve_Viral_Suppression"                       : 183.0,
                        "ART_Duration_Days_Bins": [182, 365, 730, 1095, 45625],
                        "Age_Years_Bins": [25, 35, 45, 125],
                        "CD4_Count_Bins": [0, 25, 74.5, 149.5, 274.5, 424.5, 624.5],
                        "MortalityTable": [
                            [
                                    [ 0.2015, 0.2015, 0.1128, 0.0625, 0.0312, 0.0206, 0.0162 ],
                                    [ 0.2176, 0.2176, 0.1219, 0.0675, 0.0337, 0.0223, 0.0175 ],
                                    [ 0.2350, 0.2350, 0.1316, 0.0729, 0.0364, 0.0240, 0.0189 ],
                                    [ 0.2538, 0.2538, 0.1421, 0.0787, 0.0393, 0.0260, 0.0205 ]
                                ],
                                [
                                    [ 0.1837, 0.1837, 0.0845, 0.0441, 0.0220, 0.0101, 0.0066 ],
                                    [ 0.1965, 0.1965, 0.0904, 0.0472, 0.0236, 0.0108, 0.0071 ],
                                    [ 0.2103, 0.2103, 0.0967, 0.0505, 0.0252, 0.0116, 0.0076 ],
                                    [ 0.2250, 0.2250, 0.1035, 0.0540, 0.0270, 0.0124, 0.0081 ]
                                ],
                                [
                                    [ 0.0461, 0.0461, 0.0318, 0.0258, 0.0151, 0.0103, 0.0084 ],
                                    [ 0.0502, 0.0502, 0.0346, 0.0281, 0.0164, 0.0113, 0.0091 ],
                                    [ 0.0547, 0.0547, 0.0378, 0.0306, 0.0179, 0.0123, 0.0100 ],
                                    [ 0.0596, 0.0596, 0.0412, 0.0334, 0.0195, 0.0134, 0.0109 ]
                                ],
                                [
                                    [ 0.0286, 0.0286, 0.0197, 0.0160, 0.0113, 0.0086, 0.0080 ],
                                    [ 0.0311, 0.0311, 0.0215, 0.0174, 0.0124, 0.0094, 0.0087 ],
                                    [ 0.0339, 0.0339, 0.0234, 0.0190, 0.0135, 0.0102, 0.0095 ],
                                    [ 0.0370, 0.0370, 0.0255, 0.0207, 0.0147, 0.0111, 0.0104 ]
                                ],
                                [
                                    [ 0.0161, 0.0161, 0.0110, 0.0089, 0.0063, 0.0047, 0.0044 ],
                                    [ 0.0176, 0.0176, 0.0119, 0.0097, 0.0068, 0.0052, 0.0048 ],
                                    [ 0.0192, 0.0192, 0.0130, 0.0105, 0.0075, 0.0056, 0.0053 ],
                                    [ 0.0209, 0.0209, 0.0142, 0.0115, 0.0081, 0.0061, 0.0057 ]
                                ]
                        ]
                    }
                }
            }
        },
```