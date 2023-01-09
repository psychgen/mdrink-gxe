Genotype-environment interplay in associations between maternal drinking
and offspring emotional and behavioral problems: analysis code
================
Laurie Hannigan (<laurie.hannigan@bristol.ac.uk>)

- <a href="#introduction" id="toc-introduction">Introduction</a>
- <a href="#data-preparation" id="toc-data-preparation">Data
  preparation</a>
- <a href="#running-linear-models" id="toc-running-linear-models">Running
  linear models</a>
- <a href="#running-multilevel-models"
  id="toc-running-multilevel-models">Running multilevel models</a>
- <a href="#running-snp-wise-models"
  id="toc-running-snp-wise-models">Running SNP-wise models</a>
- <a href="#creating-xpgs-and-running-xpgs-models"
  id="toc-creating-xpgs-and-running-xpgs-models">Creating xPGS and running
  xPGS models</a>
- <a href="#collating-results-and-plotting"
  id="toc-collating-results-and-plotting">Collating results and
  plotting</a>
- <a href="#used-packages" id="toc-used-packages">Used packages</a>

------------------------------------------------------------------------

## Introduction

This is the analytic code repository for the project entitled
[Genotype-environment interplay in associations between maternal
drinking and offspring emotional and behavioral
problems](add_preprint_link). In it, we walk through the scripts used to
complete the analyses for this project.

## Data preparation

The data preparation script
([./scripts/00_data_prep.R](https://github.com/psychgen/mdrink-gxe/blob/master/scripts/00_data_prep.R))
contains code to prepare the variables for analyses. It relies primarily
on the [MoBatools](add%20link) packages `phenotools` and `genotools`.

## Running linear models

The linear models (i.e., unadjusted maternal drinking by PGS moderation
analyses, plus models with adjustment for covariates in different tiers)
are run by the
[./scripts/01_run_linear_models.R](https://github.com/psychgen/mdrink-gxe/blob/master/scripts/01_run_linear_models.R)
script, with reference to a script containing helper functions for this
part of the analyses
([01a…funs.R](https://github.com/psychgen/mdrink-gxe/blob/master/scripts/01_linear_model_funs.R)).

## Running multilevel models

Models with the final level of adjustment - that is, adjustment for
unmeasured familial confounding - are run using the
[./scripts/02_run_multilevel_models.R](https://github.com/psychgen/mdrink-gxe/blob/master/scripts/02_run_multilevel_models.R)
script, with reference to a script containing helper functions for this
part of the analyses
([02a…funs.R](https://github.com/psychgen/mdrink-gxe/blob/master/scripts/02a_multilevel_model_funs.R)).
Since these models are actuall run in Mplus, these scripts refer to
scripts in the
[./scripts/mplus](https://github.com/psychgen/mdrink-gxe/blob/master/scripts/mplus)
folder.

## Running SNP-wise models

In order to create the interaction effect-based “xPGS”, we first run
SNP-wise interaction models. These are set up using the
[./scripts/03_run_SNP_models.R](https://github.com/psychgen/mdrink-gxe/blob/master/scripts/03_run_SNP_models.R)
code, which creates job scripts that are manually submitted to a high
performance computing cluster (see
[./scripts/bash](https://github.com/psychgen/mdrink-gxe/blob/master/scripts/bash)
folder.)

## Creating xPGS and running xPGS models

The scripts for creating the xPGS using PRSice (again on an HPC) are
made using the
[./scripts/04_make_xpgs.R](https://github.com/psychgen/mdrink-gxe/blob/master/scripts/04_make_xpgs.R)
code, supported by a helper functions script
([./scripts/04a_xpgs_funs.R](https://github.com/psychgen/mdrink-gxe/blob/master/scripts/04a_xpgs_funs.R)).
Then models using these xPGS are specified and run by scripts that
essentially mirror earlier ones (scripts 05-06 ≈ scripts 01-02).

## Collating results and plotting

Finally, all results are collated and plotted in the
[./scripts/07_collate_plot_res.R](https://github.com/psychgen/mdrink-gxe/blob/master/scripts/07_collate_plot_res.R)
code, supported by a helper functions script
([./scripts/07a_collate_funs.R](https://github.com/psychgen/mdrink-gxe/blob/master/scripts/07a_collate_funs.R)).

## Used packages

The following packages were used in this project, and we are grateful to
their developers and maintainers:

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-ggthemes" class="csl-entry">

Arnold, Jeffrey B. 2021. *Ggthemes: Extra Themes, Scales and Geoms for
’Ggplot2’*. <https://CRAN.R-project.org/package=ggthemes>.

</div>

<div id="ref-MplusAutomation" class="csl-entry">

Hallquist, Michael N., and Joshua F. Wiley. 2018. “MplusAutomation: An R
Package for Facilitating Large-Scale Latent Variable Analyses in Mplus.”
*Structural Equation Modeling*, 1–18.
<https://doi.org/10.1080/10705511.2017.1402334>.

</div>

<div id="ref-phenotools" class="csl-entry">

Hannigan, Laurie. 2022. *Phenotools: Facilitates Reproducible Workflows
with Data from MoBa and Linked Registry Sources*.
<https://github.com/psychgen/phenotools>.

</div>

<div id="ref-tictoc" class="csl-entry">

Izrailev, Sergei. 2021. *Tictoc: Functions for Timing r Scripts, as Well
as Implementations of Stack and List Structures*.
<https://CRAN.R-project.org/package=tictoc>.

</div>

<div id="ref-patchwork" class="csl-entry">

Pedersen, Thomas Lin. 2020. *Patchwork: The Composer of Plots*.
<https://CRAN.R-project.org/package=patchwork>.

</div>

<div id="ref-grateful" class="csl-entry">

Rodríguez-Sánchez, Francisco, Connor P. Jackson, and Shaurita D.
Hutchins. 2022. *Grateful: Facilitate Citation of r Packages*.
<https://github.com/Pakillo/grateful>.

</div>

<div id="ref-lavaan" class="csl-entry">

Rosseel, Yves. 2012. “<span class="nocase">lavaan</span>: An R Package
for Structural Equation Modeling.” *Journal of Statistical Software* 48
(2): 1–36. <https://www.jstatsoft.org/v48/i02/>.

</div>

<div id="ref-tidyverse" class="csl-entry">

Wickham, Hadley, Mara Averick, Jennifer Bryan, Winston Chang, Lucy
D’Agostino McGowan, Romain François, Garrett Grolemund, et al. 2019.
“Welcome to the <span class="nocase">tidyverse</span>.” *Journal of Open
Source Software* 4 (43): 1686. <https://doi.org/10.21105/joss.01686>.

</div>

</div>
