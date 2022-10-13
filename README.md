# GAN-WGCNA

## Publication
### GAN-WGCNA: calculating gene modules as a way to find key intermediate regulators in cocaine addiction
#### Workflow
```mermaid
graph TD;
    GAN_training -- Latent_space_interploation --> Time-series_gene_expression_profile;
    Time-series_gene_expression_profile -- A --> GAN-WGCNA;
    Time-series_gene_expression_profile -- A --> rescued_DEG;

```
## Reference
- Park J, Kim H, Kim J, Cheon M (2020) A practical application of generative adversarial networks for RNA-seq analysis to predict the molecular progress of Alzheimer's disease. PLOS Computational Biology 16(7): e1008099. https://doi.org/10.1371/journal.pcbi.1008099
