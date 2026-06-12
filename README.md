# SDM-Guide-thesis-scripts
R projects for the data processing, modelling, and analyses conducted as part of a master's thesis. These scripts are to be used for learning species distribution modelling (with MaxEnt), although improving copies of the scripts for other projects is encouraged. Larger files will be hosted on Zenodo after the thesis is approved.

# Guide to RProjects linked to StoryMap

Project: **Is Nordic cloudberry moving with climate change? A transmedia guide to species distribution modelling with MaxEnt**

Author: **Nicole K. J. Drewitz** Project Affiliation: University Centre of the Westfjords, Master of Coastal and Marine Management

AI Use: Scripts were drafted with Claude, and sometimes debugged with Claude or Perplexity

------------------------------------------------------------------------

## Overview

*R projects for the data processing, modelling, and analyses conducted as part of a master's thesis. These scripts are to be used for learning species distribution modelling (with MaxEnt), although improving copies of the scripts for other projects is encouraged.*

**SDM Guide**: <https://arcg.is/aP8eq0>

Further info found in the associated **thesis**: Skemman.is

Data Sources: See SDM guide, or thesis

Attribution: CC-BY preferred, some of the data used to create models are CC-BY-NC

Available from GitHub: <https://github.com/NicoleDrewitz/SDM-Guide-thesis-scripts/>

## R Projects
On macOS, the Terminal command to keep the computer on for 2 hours while a long script is running is `caffeinate -i -t 3600`

### A. Data Preparation

> 1.  **Presence/absence points CSV creation**
> 2.  **Select modelling area**
> 3.  **Formatting environmental variables**
> 4.  **Reduce spatial bias in location points CSV**

### B. Build Models

> 1.  **Remove highly correlated variables**
> 2.  **Reduce multicollinearity**
> 3.  **Tune parameter settings**
> 4.  **Stepwise variable removal**
> 5.  **Rerun the final model manually**

### C. Analysis

> 1. **Important environmental variables**
>
> 2.  **Model performance**
>     1.  *value extraction from model*
>     2.  *ROC*
>     3.  *point density*
> 3.  **Moving habitat suitability**
>     1.  *area change*
>     2.  *environmental variable change at presence locations*
> 4.  **(Eco)regional comparisons**
>     1.  *make sampling locations*
>     2.  *extract values at sampling locations*
>     3.  *proportion of highly suitable habitat area*
>     4.  *change in mean suitability*
>     5.  *change of area*
>     6.  *summary tables*

renv lock preserves software version compatability

*Last updated: 2026-05*
