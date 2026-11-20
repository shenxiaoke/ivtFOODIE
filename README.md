# ivtFOODIE

Code accompanying the ivtFOODIE study.

| Workflow | Purpose |
| --- | --- |
| [Mapping](1-ivtFOODIE-mapping-pipeline/) | Process sequencing reads and generate conversion ratios. |
| [Kd fitting](2-ivtFOODIE-fit-Kd-pipeline/) | Fit Kd. |
| [Seq2Kd–eQTL](3-ivtFOODIE-Seq2Kd-eQTL/) | Relate predicted affinity changes to GTEx eQTL effects. |
| [Seq2Kd–ClinVar](4-ivtFOODIE-Seq2Kd-ClinVar/) | Compare coding and promoter predictions between pathogenic and benign variants. |

Each workflow documents its required inputs, execution order and output locations. The association workflows reference frozen data releases and the corresponding Seq2Kd model inputs.
