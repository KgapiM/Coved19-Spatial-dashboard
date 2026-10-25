# COVID-19 Spatial Statistics Dashboard

An interactive R Shiny dashboard for exploring provincial COVID-19 indicators and Community Health Centre locations in South Africa.

**Module:** STAT7006A Spatial Statistics  
**Institution:** University of the Witwatersrand  
**Presentation date:** 26 October 2026

## Project links

- [Open the live dashboard](https://khathucry.shinyapps.io/Covid19/)
- [View the GitHub repository](https://github.com/KgapiM/Coved19-Spatial-dashboard)

## Group members

| Group member | Student number |
| --- | --- |
| Khutso Makobela | 1370561 |
| K Lydia Moukangwe | 554423 |
| Khathutshelo Maumela | 871907 |
| Veronica Sikhosana | 0105163x |

## Research questions

1. Does reported COVID-19 burden show overall spatial association across South Africa’s nine provinces?
2. Which provinces form significant local clusters or spatial outliers?
3. Do Community Health Centre locations depart from Complete Spatial Randomness?

The provincial and facility analyses address complementary questions. They do not establish a causal relationship between facility provision and COVID-19 outcomes.

## Dashboard sections

| Section | Purpose |
| --- | --- |
| Home | Study introduction, headline counts and navigation |
| Map Explorer | Provincial indicators, facility locations and trends |
| Global Moran’s I | Overall provincial spatial association, permutations, neighbour connections and weights |
| Local Moran’s I | Local clusters, spatial outliers and p-value adjustments |
| Point-pattern analysis | Kernel density, quadrat testing and distance-function simulations |
| About | Data information and limitations |

## Methods and justification

### Global Moran’s I

Tests overall spatial autocorrelation in a provincial indicator. Users can choose queen contiguity, rook contiguity or nearest neighbours, with binary or row-standardised weights. Analytical and permutation results support interpretation.

### Local Moran’s I

Examines individual provinces in relation to their neighbours. The dashboard distinguishes High–High and Low–Low patterns from High–Low and Low–High outliers, subject to significance. Unadjusted, Bonferroni and Benjamini–Hochberg false-discovery-rate p-values are available.

### Complete Spatial Randomness tests

The quadrat chi-square test compares facility counts with area-based expectations under a homogeneous Poisson process. K, L, G and F functions provide supporting distance-based analysis. Kernel density estimation describes facility intensity; it is not a significance test.

## Data and preparation

The supplied Excel workbook contains provincial cumulative cases, deaths and recoveries, plus Community Health Centre coordinates. The shared analysis period described in the project is 17 April–4 October 2020.

The preparation script:

- Reads and reshapes the workbook sheets.
- Matches province names and combines the disease series.
- Carries the latest cumulative observation forward over reporting gaps.
- Creates a manually entered 2020 provincial population table.
- Downloads GADM level-one province boundaries using `geodata`.
- Converts valid facility coordinates into spatial points and reports geographic mismatches.
- Saves the prepared objects used by the app.

Population values are described in the preparation script as Statistics South Africa 2020 mid-year estimates and require verification against the original release. The original provenance and redistribution permissions of the disease and facility workbook should also be confirmed. Boundaries are supplied by [GADM](https://gadm.org/); consult its usage conditions.

### Indicators

| Indicator | Calculation |
| --- | --- |
| Cases per 100,000 | Cumulative cases / population × 100,000 |
| Deaths per 100,000 | Cumulative deaths / population × 100,000 |
| Case fatality ratio | Deaths / reported cases × 100 |
| Recovery ratio | Recoveries / reported cases × 100 |
| Centres per 100,000 | Facility count / population × 100,000 |
| Raw cases | Cumulative reported case count |

Facility availability is static in this app and does not measure travel accessibility or capacity. Case fatality is not infection fatality.

## Expected project structure

```text
app.R
00_prepare_data.R
README.md
data/
    ts.rds
    pop.rds
    prov.rds
    chc.rds
```

The original Excel workbook is needed to regenerate prepared objects. Keep any additional assets required by the deployed app in their expected relative folders.

## Run locally

### 1. Open the project

Download the repository or clone it, then open the project folder in RStudio. Use this folder as the working directory so relative `data/` paths resolve correctly.

### 2. Install packages

Run once in the R Console:

```r
packages <- c(
  "shiny", "bslib", "dplyr", "tidyr", "readxl", "ggplot2",
  "sf", "leaflet", "spdep", "spatstat", "geodata"
)
missing <- setdiff(packages, rownames(installed.packages()))
if (length(missing)) install.packages(missing)
```

### 3. Prepare data if required

If all four prepared RDS files are present, proceed to launch. Otherwise, obtain the permitted source workbook. In `00_prepare_data.R`, set the workbook path, including its extension, or use:

```r
xl <- file.choose()
```

Save the script, then run:

```r
source("00_prepare_data.R")
```

Select the workbook if prompted. The boundary download requires internet access. Wait for `Done. Files written to ./data`.

### 4. Launch

```r
shiny::runApp(".")
```

Alternatively, open `app.R` and click **Run App**. Internet access is needed for online basemap tiles.

## Example results

The supplied dashboard screenshots show the following for cumulative cases per 100,000 on 4 October 2020, using queen contiguity, row-standardised weights and a two-sided test:

- Global Moran’s I: **0.144**.
- Expected value: **−0.125**.
- Standardised test score: **1.43**.
- Analytical p-value: **0.154**.
- Permutation p-value: **0.194**, using 999 permutations.
- No significant local clusters or outliers at the 5% level, even before p-value adjustment.

These outputs provide insufficient evidence of provincial spatial association under those settings. Failure to reject the null does not prove randomness.

The dashboard reports 286 facilities with coordinates, with 284 retained inside the simplified point-pattern window. The density map shows uneven intensity, especially around Gauteng. A formal CSR conclusion requires the relevant test output and cannot be inferred from the density map alone.

## Limitations

- Nine provinces restrict statistical power and local stability.
- Neighbour definitions and geographic aggregation can affect results.
- Provincial relationships cannot automatically be attributed to individuals.
- Reporting differences, filled reporting gaps and population denominators affect interpretation.
- The facility inventory may be incomplete, and the observation window excludes some points.
- Homogeneous CSR does not account for population concentration.
- Repeated daily tests and pointwise simulation envelopes require cautious interpretation.

## Troubleshooting

**Workbook path does not exist:** use `file.choose()` or provide the full filename and extension in the preparation script.

**Missing RDS files:** run preparation from the app folder, or place the prepared files under `data/`.

**Package not installed:** install the missing package, then restart the app.

**Simulation not updating:** select the desired function and simulation count, then click **Run simulation** again.

## Acknowledgements and reuse

This repository documents a group coursework project. Retain original author credits and disclose outside assistance according to module requirements. No new licence grant is made by this README; reuse of code and data remains subject to the applicable permissions and source licences.
