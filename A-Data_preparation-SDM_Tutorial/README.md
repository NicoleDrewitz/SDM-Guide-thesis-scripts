# A-Data Preparation: SDM Tutorial *(R Project)*

Project: Is Nordic cloudberry moving with climate change? A transmedia guide to species distribution modelling with MaxEnt

Author: Nicole K. J. Drewitz Project Affiliation: University Centre of the Westfjords, Master of Coastal and Marine Management

AI Use: Scripts were drafted with Claude, and sometimes debugged with Claude or Perplexity

The scripts in this project help prepare data for species distribution
modelling with MaxEnt or other biostatistical modelling methods.

SDM Guide: https://arcg.is/aP8eq0

Thesis: Skemman.is

## Prerequisites

Ensure you have R installed and the following libraries:
* `terra`
* `sf`
* `readr`
* `httr` *(for bulk downloads)*
* `here` *(optional) to set file paths inside R project*
* `tictoc` *(optional
to track script runtime)*
* check for other libraries before running
each script

## Usage

1.  Replace any source files in the `/Input files` folder.

2.  Run scripts matching your data or desired bulk environmental
    variable download from your scripts folder:

    ``` r
    Data preparation scripts/
    ```

3.  Reset variable extent and cartographic projection as needed.

## R Scripts Summary

-   Found in `Data preparation scripts` folder
-   **`Subsampling CSV files`**
-   Scripts starting in:
    -   **`a - Automated Bulk Download`** for CHELSA bioclim, soilgrids,
        and WEkEO variables.
    -   **`b -` merge** spatially or temporally separated variables.
    -   **`c -`** is for individual and **bulk processing of
        environmental variables** for unified formatting.
    -   **`d -`** environmental variable checks for **missing cell
        values** and possible fix.
    -   **`e -`** After bioclim variables are moved to a folder with
        subfolders for each climate scenario, this script removes file
        name differences between scenarios *(needed for MaxEnt to
        **match variables between climate scenarios)***.
    -   **`O -`** other optional raster file checks.

## Formatted outputs List

- Presence locations (.csv)
- Absence locations (.csv) _(if possible for model selection and analysis)_
- Folder of recent/current environmental variables (.asc)
- Folder(s) of future environmental variables (.asc)
- Optional: vector file for regional analysis (i.e. ecoregions.shp)

Location files can be randomly separated into files for presence-only modelling, model selection, further model performance analysis.

## A. Data Preparation Process

See StoryMap flowcharts for more details: <https://arcg.is/aP8eq0>

> ### 1. Presence/absence points CSV creation
>* (flowchart in StoryMap)
>
>**Required:** Species locations (e.g. GBIF download)\_
>
> Completed manually with QGIS to clean formatting.

> ### 2. Select modelling area
>
> *(flowchart in StoryMap)*
>
> ***Required:** CSV presence locations*
>
> a.  Study area manually created with QGIS.
>
> b.  Modelling area selected from polygons (e.g. ecoregions or
>     watersheds) with study species observations.

> ### 3. Formatting environmental variables
>
> ***Required:** Modelling area vector file*
>
> a.  Bulk download scripts for global bioclimatic and soil variables.
>
> b.  Scripts for unifying environmental variables rasters for
>     modelling.
>
> c.  Script to check for missing values within the modelling area.

> ### 4. Reduce spatial bias in location points CSV
>
> ***Required:** CSV locations files with cleaned formatting, and an
> environmental variable raster (formatted)*
>
> a.  R scripts created to randomly subset CSV files of presence and
>     absence locations by percentages,
>
> b.  Then remove locations if there are more than 1 point within an
>     environmental variable pixel (in GIS).

> ### Other
>
> a.  Read start of ASC files *(to check consistent formatting)*
>
> b.  Unify variable naming format
>
> c.  Check CRS warning logic
>
> *(warning will often appear if projection name is listed slightly
> different, even if both projections are actually the same. Can be
> checked by opening in GIS software)*

*Last updated: 2026-05*
