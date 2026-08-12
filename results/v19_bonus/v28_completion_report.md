# v28 completion report

## Priority additions completed

- Four-subtype resampling stability: compared k = 3, 4 and 5 using PAC, silhouette width, within-cluster consensus, adjusted Rand index, sample-label retention and cluster-wise Jaccard similarity across 200 independent 80% subsamples.
- Full adjusted immune phenotype: tested 50 ESTIMATE, marker-score, ssGSEA, CIBERSORT and checkpoint features in separate ER-adjusted and PAM50-adjusted nested linear models, with stage and outcome-appropriate ESTIMATEScore adjustment.

## Main findings

- k = 4: PAC 0.0967; median ARI 0.930; median held-out label retention 97.4%; median mean cluster Jaccard 0.947. k = 3 remained a more parsimonious competitor with slightly stronger resampling indices, which is stated explicitly in the manuscript.
- Immune features passing model-wide BH-FDR < 0.05: 38/50 after ER adjustment and 19/50 after PAM50 adjustment. The immune phenotype is partly, but not wholly, explained by intrinsic subtype and tumour context.

## QA

- TCGA sample mappings verified for all 952 samples; locked k = 4 labels matched the manuscript labels exactly after optimal relabelling.
- Figure-source audit: 0 failures. PDF minimum text size: S13 6.0 pt; S14 5.3 pt.
- v28 render: 35 pages, 16 inline main figures and 5 main tables retained; full-page contact-sheet inspection passed.
