# Integrated data state-space models for estimating Northwest raptor trends
This repo contains supplementary data, R code, JAGS code, and a website for 
C. J. W McClure, B. W. Rolek, J. Fleischer. Composite population trends reveal status of wintering raptors in the Northwestern USA. 2023. Biological Conservation.  

Zenodo repo : [![DOI](https://zenodo.org/badge/566442152.svg)](https://zenodo.org/badge/latestdoi/566442152)

Corresponding author for repo: rolek.brian at peregrinefund.org  

Materials are permanently archived online at (INSERT ZOTERO LINK after acceptance).

We assessed population trends of raptors in the Northwestern United States using state-space models that integrate data from disparate monitoring programs (Christmas Bird Counts and Winter Raptor Surveys) to create composite population trends.

Examples of code for integrated models are also available as a [workflow HTML website here.](https://the-peregrine-fund.github.io/northwest-trends/workflow.html "Workflow for integrated data state-space models")

# METADATA 
Data are included in the "data" folder as an .Rdata file. Once loaded into R, a list object named "dall" contains data for 11 raptor species. Here we provide metadata for "dall". Each species has it's own list within "dall". For example, dall[["RTHA"]] will return all Red-tailed Hawk data necessary for the integrated model for JAGS.

Each species list has the following fields:
* count- WRS counts of that species.
* dist- total distance driven during a WRS survey.
* time- The winter when the WRS survey was conducted. E.g. 2004-2005=1, ..., 2020-2021=X
* rt- route number of WRS survey.
* frst- first time a WRS route was surveyed.
* strat1- strata index of WRS routes
* sp- average speed during WRS surveys. Note this variable was only available for WRS data.
* ncounts- total number of WRS counts.
* nroutes- total number of WRS routes.
* nstrata- total number of strata where this species was detected in WRS data.
* ntime- total number of winters surveys were conducted for WRS data.
Note that the CBC data are contained within the same list and are similar but have a '2' as a suffix.
* count2- CBC counts of that species.
* dist2- total distance driven during a CBC survey.
* time2- The winter when the CBC survey was conducted. E.g. 2004-2005=1, ..., 2020-2021=X
* rt2- route number of CBC survey.
* frst2- first time a CBC route was surveyed.
* strat2- strata index of CBC routes
* ncounts2- total number of CBC counts.
* nroutes2- total number of CBC routes.
* nstrata2- total number of strata where this species was detected in CBC data.
* ntime2- total number of winters surveys were conducted for CBC data.
