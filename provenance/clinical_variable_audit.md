# Clinical-variable availability audit

## GSE163256 / primary cohort

The primary cohort was recruited through ICHseq, a substudy of the MISTIE III surgical arm. The source report states that a catheter was placed 23-90 hours after ICH onset, followed by recombinant tissue plasminogen activator (rtPA) irrigation approximately every 8 hours to facilitate hematoma drainage. Blood and hematoma effluent were collected longitudinally while the catheter remained in place.

The source study directly examined rtPA as a technical/clinical exposure. Samples collected at initial surgery before rtPA administration followed the same main temporal expression trends as post-rtPA samples, and ex vivo rtPA exposure did not significantly alter the highlighted pathways. Table S8 of the source report contains the rtPA differential-expression analysis.

The public log-FPKM matrix used here encodes patient, compartment, and day. It does not encode observation-linked rtPA dose/time, catheter manipulation, rebleeding, antibiotics, corticosteroids, infection status, or other medication histories. Public source materials describe cohort-level intervention and selected clinical characteristics, but do not supply a complete observation-level medication/exposure table that can be joined unambiguously to the 136 locked profiles. These variables therefore were not entered as covariates.

## GSE166638 / exploratory case study

The one-patient report documents catheter drainage, an asymptomatic rebleeding event estimated from CT, drainage, and sequencing changes, and cessation of tPA treatment. These data are clinically informative within that patient but cannot estimate population-level treatment or rebleeding effects.

## Analysis decision

- Compartment and collection day are included in the primary model.
- Patient is modeled as the repeated-measures unit.
- rtPA/drainage are treated as cohort design features and discussed using the source study's direct sensitivity evidence.
- Rebleeding, infection, antibiotics, corticosteroids, and other medications are acknowledged as unavailable observation-level covariates, not coded as absent.

## Primary sources

1. Askenase MH et al. Sci Immunol. 2021;6:eabd6279. https://doi.org/10.1126/sciimmunol.abd6279
2. Goods BA et al. JCI Insight. 2021;6:e145857. https://doi.org/10.1172/jci.insight.145857
3. MISTIE III registration, NCT01827046. https://clinicaltrials.gov/study/NCT01827046

