
# B-Build Models: SDM Tutorial *(R Project)*

Project: Is Nordic cloudberry moving with climate change? A transmedia guide to species distribution modelling with MaxEnt

Author: Nicole K. J. Drewitz Project Affiliation: University Centre of the Westfjords, Master of Coastal and Marine Management

**The scripts in this project help:** a. limit the number of environmental variables in a MaxEnt model of habitat suitability, b. iteratively run MaxEnt models, c. select the best model from many using AICc prediction performance comparisons based on habitat suitability at presence vs absence locations.

AI Use: Scripts were drafted with Claude, and sometimes debugged with Claude or Perplexity

SDM Guide: https://arcg.is/aP8eq0

Thesis: Skemman.is

## Prerequisites

Ensure you have R and MaxEnt installed and the following libraries:
* `terra`
* `sf`
* `dismo` to access MaxEnt
* `rJava` to run MaxEnt
* `ggplot2`
* `dplyr`
* `xml2`
* `rvest`
* `raster`
* `gridExtra`
* `corrplot`
* `here` *(optional) to set file paths inside R project*
* `tictoc` *(optional
to track script runtime)*
* check for other libraries before running
each script

## Usage

1.  (Re)place any source files in the `/Input_files` folder, or set new file paths. **Required:** Formatted environmental variables (current/recent climate conditions), modelling species presence locations (CSV), presence/absence locations reserved for model selection (CSVs).

2. Set the file path to your maxent.jar application in each script.

3.  Run scripts matching your data or desired bulk environmental
    variable download from your scripts folder:

    ``` r
    Modelling_Scripts/
    ```

## R Scripts Summary

-   Found in `Modelling_Scripts` folder
-   0-**`MaxEnt untuned settings RUN`** produces a single model run to trial finding variable importance rankings with MaxEnt via permuation importance.
-   Other scripts named with `#-` prefix matching `Modelling Process` steps.

## Outputs
>- Reduced environmental variable dataset _(3-10 most important variables)_

> - Final MaxEnt habitat suitability model _(from MaxEnt application)_
>
> a. Current/recent climate model (ASC)
b. Future climate models (ASCs)
c. binary threshold for suitable habitat (i.e. likely habitat range)
d. Modelling uncertainty rasters
e. html file (and attached CSV data) with permutation importance of variables, variable response curves, AUC, ROC curve, and record of model settings used.

## B. Modelling Process

*See StoryMap timeline and thesis for more details: <https://arcg.is/aP8eq0> | <skemman.is>

>### 1. Remove highly correlated variables 
>
> **Required:** Untuned settings MaxEnt run script.
>
> a. Create a MaxEnt trial run
>
> b. Select a method for removing highly correlated variables, then run the script to choose variables to remove
>
> c. Manually return any variables expected to be ecologically significant, and review correlation

>### 2. Reduce multicollinearity
>
> **Required:** Untuned settings MaxEnt run script.
>
> a. Create a MaxEnt trial run with remaining variables
>
> b. Create a PCA report
>
> c. Manually add to the comparison table, assessing the ecological significance of response curves to select environmental variables to remove.
>
> Can be repeated to reduce noise in response curves. (e.g. 1st scripts run and review = remove 10 variables, then rerun scripts and comparison with fewer variables, and select more to remove.)

>### 3. Tune parameter settings
>
> **Required:** 10 _(recommended)_ environment variable rasters, presence location for modelling, presence/absence locations reserved for model selection.
>
> a. Iterate through fine-tuning parameter settings
>
> b. select best tuning settings with (AICc)
>
> c. optional: view additional plots on AICc change between parameter settings

>### 4. Stepwise variable removal
>
> **Required:** 10 _(recommended)_ environment variable rasters, parameter setting selected from step 3, and presence/absence locations reserved for model selection
>
> a. Iteratively remove the least important variable between model runs.
>
> b. Select the best final model settings from the AICc comparison (sampled from model selection presence/absence locations).

>### 5. Rerun the final model manually
>
> **Required:** Parameter settings selected from step 3, environmental variables for the final model from step 4, and the option to use more time-consuming settings _(see StoryMap)_.
>
> a. Create a bias layer for MaxEnt (raster file) to influence background points selection using the same format as these R scripts.
>
> b. Rerun the final model manually in the MaxEnt application (.jar file) to access the latest version.

*Last updated: 2026-05*