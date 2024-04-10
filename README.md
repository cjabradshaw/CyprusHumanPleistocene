# Modelling human arrival and spread in Cyprus 🇨🇾 at the end of the Pleistocene
<a href="https://zenodo.org/doi/10.5281/zenodo.10951687"><img src="https://zenodo.org/badge/DOI/10.5281/zenodo.10951688.svg" alt="DOI"></a>
<a href="https://www.ucy.ac.cy/migrate/"><img align="right" src="www/MIGRATElogoShad.png" width="200" style="margin-top: 20px"></a>

Part of the <a href="https://www.ucy.ac.cy/migrate/">MIGRATE</a> (<strong>M</strong>odell<strong>i</strong>ng Demo<strong>gr</strong>aphy and <strong>A</strong>daptation in the Initial Peopling of the Eastern M<strong>e</strong>diterranean Islandscape) project, under the auspices of the European Union Research and Innovation Foundation for Research, Technological Development and Innovation "Restart 2016-2020".
<br>
<br>
<strong>lead investigator</strong>: Dr <a href="https://www.ucy.ac.cy/directory/en/profile/tmouts01">Theodora Moutsiou</a><br>
<strong>key personnel</strong>: Dr <a href="https://scholar.google.com.au/citations?user=BU25ogMAAAAJ&hl=en">Christian Reepmeyer</a>, Associate Professor <a href="https://www.ucy.ac.cy/directory/en/profile/demest">Stella Demesticha</a>, Dr <a href="https://www.ucy.ac.cy/directory/en/profile/arkasian">Vasiliki Kassianidou</a>, Dr <a href="https://www.cut.ac.cy/faculties/fet/ceg/staff/athos.agapiou/?languageId=1">Athos Agapiou</a>, Dr <a href="https://www.researchgate.net/profile/Zomenia-Zomeni">Zomenia Zomeni</a>, Professor <a href="https://globalecologyflinders.com/people/#DIRECTOR">Corey Bradshaw</a><br>
<strong>collaborators</strong>: Dr <a href="https://globalecologyflinders.com/people/#COORDINATOR">Frédérik Saltré</a>
<br>
## Project overview
Project <a href="https://www.ucy.ac.cy/migrate/">MIGRATE</a> seeks to offer novel insights into population dynamics and range shifts that resulted in dispersals from the Eastern Mediterranean mainland to the island of <a href="https://www.google.com/maps/place/Cyprus/@35.1670135,32.765821,9z/">Cyprus</a> at a critical period (Late Pleistocene, 45-12 ka) through stochastic spatial modelling. This advanced modelling will  enhance our understanding of timing, and climatic and social factors important in the initial colonisation of Cyprus. The proposed project aims to establish new research domains in the field of Cypriot archaeology extending traditional chronological frontiers beyond the Holocene (current Warm Period), encompassing innovative and interdisciplinary methodologies at the forefront of archaeological research.

## Manuscript
Bradshaw, CJA, C Reepmeyer, F Saltré, A Agapiou, V Kassinadiou, S Demesticha, Z Zomeni, M Polidorou, T Moutsiou. 2024. Demographic models predict end-Pleistocene arrival and rapid expansion of pre-agropastoralist humans in Cyprus. <strong><em>Proceedings of the National Academy of Sciences of the USA</em></strong> in press<br>
<br>
and (out-of-date) pre-print:<br>
<br>
Bradshaw, CJA, C Reepmeyer, F Saltré, A Agapiou, V Kassinadiou, S Demesticha, Z Zomeni, M Polidorou, T Moutsiou. <a href="http://doi.org/10.21203/rs.3.rs-3468157/v1">Demographic models predict end-Pleistocene arrival and rapid expansion of pre-agropastoralist humans in Cyprus</a>. <em>Research Square</em> (preprint) doi:10.21203/rs.3.rs-3468157/v1

### Abstract
The antiquity of human dispersal into Mediterranean islands and ensuing coastal adaptation have remained largely unexplored due to the prevailing assumption that the sea was a barrier to movement, and that islands were hostile environments to early hunter-gatherers. Using the latest archaeological data, hindcasted climate projections, and age-structured demographic models, we demonstrate clear evidence for early arrival (14,257 to 13,182 calendar years ago) to Cyprus, and predicted that large groups of people (~ 1,000 to 1,375) arrived in 2–3 main events occurring within < 100 years to ensure low extinction risk. These results indicate that the post-glacial settlement of Cyprus involved only a few large-scale, organised events requiring advanced watercraft technology. Our spatially debiased and Signor-Lipps-corrected estimates indicate rapid settlement of the island within < 200 years, and expansion to a median of 4,000–5,000 people (0.36–0.46 km<sup>−2</sup>) in < 11 human generations (< 300 years). Our results do not support the hypothesis of inaccessible and inhospitable islands in the Mediterranean by pre-agropastoralists, agreeing with analogous conclusions for other parts of the world such as the Indo-Pacific. Our results also highlight the need to revisit these questions in the Mediterranean and test their validity in light of new technologies, field methods, and data. By applying stochastic models based on both temporally and spatially explicit data for the first time to the Mediterranean region, we are able to place Cyprus and large islands in general as attractive and favourable destinations for palaeolithic peoples.

## Scripts
R code by Corey Bradshaw (<a href="http://github.com/cjabradshaw">@cjabradshaw</a>) and Frédérik Saltré (<a href="http://github.com/fredsaltre">@fredsaltre</a>)

### Climate hindcasts
- <code>HadCM3outputs.R</code>: net primary production (kg C m<sup>-2</sup> year<sup>-1</sup>), temperature (°C), and precipitation (mm year<sup>-1</sup>) data (raw and anomalies relative to the present) hindcasted from the <a href="https://www.metoffice.gov.uk/research/approach/modelling-systems/unified-model/climate-models/hadcm3">Hadley Centre Coupled Model version 3</a> (HadCM3) climate model for the 0.5° × 0.5° lat/lon cells covering Cyprus from 20 ka to the present

### Chronology
- <code>calibdate.R</code>: calibrate oldest human archaeological evidence for Cyprus, and use Signor-Lipps correction to estimate window of first entry
- <code>Step-by-Step(CRIWM).R</code>: detailed, step-by-step guide and code to using CRIWM if installing package from Github problematic

### Demography
- <code>CypAspatDemModel.R</code>: stochastic, age-structured human demographic projection model
- <code>CypAspatDemModelMVPsim.R</code>: based on model above, vary founding population size to estimate change in probability of quasi-extinction
- <code>CypAspatDemModelStaggeredEntry.R</code>: vary frequency and intervals of arrival of different group sizes (up to minimum viable population size calculated above) to determine minimum size of entering groups and frequency of immigration required to minimise quasi-extinction
- <code>densities.R</code>: bootstrapped estimates of human densities across Europe/Near East during the Late Pleistocene

### Spatial arrival
- <code>SpatialArrivalEswtimates.R</code>: spatial de-biaising functions to estimate pattern of initial arrival from distribution of archaeological dates across sites

### Source functions
- <code>matrixOperators.r</code>: functions for manipulating matrices for population projections

## Data
### Climate hindcasts
- <em>CyprusRegion(20ka)_NPP(absolutevalues).csv</em>: hindcasted HadCM3 net primary production (absolute values) from 20 ka to present for the general eastern Mediterranean region
- <em>CyprusRegion(20ka)_NPP(anomaliesvalues).csv</em>: hindcasted HadCM3 net primary production (anomalies relative to the present) from 20 ka to present for the general eastern Mediterranean region
- <em>CyprusRegion(20ka)_AnnualMeanTemperature(absolutevalues).csv</em>: hindcasted HadCM3 temperature (absolute values) from 20 ka to present for the general eastern Mediterranean region
- <em>CyprusRegion(20ka)_AnnualMeanTemperature(anomaliesvalues).csv</em>: hindcasted HadCM3 temperature (anomalies relative to the present) from 20 ka to present for the general eastern Mediterranean region
- <em>CyprusRegion(20ka)_AnnualPrecipitation(absolutevalues).csv</em>: hindcasted HadCM3 annual precipitation (absolute values) from 20 ka to present for the general eastern Mediterranean region
- <em>CyprusRegion(20ka)_AnnualPrecipitation(anomaliesvalues).csv</em>: hindcasted HadCM3 annual precipitation (anomalies relative to the present) from 20 ka to present for the general eastern Mediterranean region
  
### Sea level
- <em>cypareaSL.csv</em>: relationship between relative sea level (from the present) and area of the island of Cyprus derived from <a href="http://www.gebco.net">GEBCO</a>
- <em>lambeckESL.csv</em>: global sea level curve from <a href="https://doi.org/10.1073/pnas.1411762111">Lambeck et al. (2014)</a>

### Archaeological dates
- <em>archdates.csv</em>: oldest radiocarbon dates from 10 sites around Cyprus (data from <a href="https://doi.org/10.5281/zenodo.5767862">Near East Radiocarbon Dates (NERD)</a>, <a href="https://doi.org/10.5334/joad.90">Palmisano et al. 2022</a>)

### Demography
- <em>world2013lifetable.csv</em>: human demographic data for age-structure determination (see <a href="http://doi.org/10.1038/s41559-019-0902-6">Bradshaw et al. 2019</a>)
- <em>densities.csv</em>: estimates of human densities standardised to people km<sup>-2</sup> across Europe/Near East during the Late Pleistocene

<p><a href="https://www.ucy.ac.cy"><img align="bottom-left" src="www/UCypruslogo.png" alt="UCyprus logo" height="45" style="margin-top: 20px"></a> &nbsp; <a href="http://www.dainst.org"><img align="bottom-left" src="www/DAIlogo.png" alt="DAI logo" height="55" style="margin-top: 20px"></a> &nbsp; &nbsp; <a href="https://www.flinders.edu.au"><img align="bottom-left" src="www/Flinders_University_Logo_Horizontal_RGB_Master.png" alt="Flinders University logo" height="35" style="margin-top: 20px"></a> &nbsp; <a href="https://globalecologyflinders.com"><img align="bottom-left" src="www/GEL Logo Kaurna New Transp.png" alt="GEL logo" height="55" style="margin-top: 20px"></a> &nbsp; &nbsp; <a href="https://EpicAustralia.org.au"><img align="bottom-left" src="www/CabahFCL.jpg" alt="CABAH logo" height="45" style="margin-top: 20px"></a> &nbsp; <a href="https://www.cut.ac.cy"><img align="bottom-left" src="www/CUTlogoblack.png" alt="CUT logo" height="50" style="margin-top: 20px"></a> &nbsp; <a href="https://www.moa.gov.cy/moa/gsd/gsd.nsf/dmlIndex_en/dmlIndex_en"><img align="bottom-left" src="www/CGSlogo.png" alt="CGS logo" height="45" style="margin-top: 20px"></a></p>
