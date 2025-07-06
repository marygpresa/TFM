# TFM
Bioinformatics thesis project

Option A: AML
First approach was using HiPathia to perform a mechanistic analysis on a cohort of more than 3000 RNAseq samples to study the effect of the inhibition of METTL3 in cellular pathways in AML.

Some opportunities/obstacles:

    -no paired data: AML is in blood so no healthy samples can be taken in parallel for the same patient, a cohort of AML patients and a cohort of healthy patients has to be taken separately. I used TARGET-AML cohort for AML and GTex for healthy.

    - METTL3 gene is not part of HiPathia so can't directly inhibitd it in silico and simulate the effect: my alternative is to use a new cohort were inhibition has been carried out in vivo https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE94613 from study : https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE94613.


Option B: analyse lung cancer
