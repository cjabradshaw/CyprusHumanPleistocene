# Modelling the arrival (and consequences) of humans in Cyprus at the end of the Pleistocene
<img align="right" src="www/MIGRATElogo.jpg" width="200" style="margin-top: 20px">

Part of the MIGRATE (<strong>M</strong>odell<strong>i</strong>ng Demo<strong>gr</strong>aphy and <strong>A</strong>daptation in the Initial Peopling of the Eastern M<strong>e</strong>diterranean Islandscape) project, under the auspices of the European Union Research and Innovation Foundation for Research, Technological Development and Innovation "Restart 2016-2020".
<br>
<br>
<strong>lead investigator</strong>: Dr <a href="https://ucy.academia.edu/TheodoraMoutsiou">Theodora Moutsiou</a><br>
<strong>key personnel</strong>: Dr <a href="https://scholar.google.com.au/citations?user=BU25ogMAAAAJ&hl=en">Christian Reepmeyer</a>, Associate Professor <a href="https://www.ucy.ac.cy/directory/en/profile/demest">Stella Demesticha</a>, Dr <a href="https://www.ucy.ac.cy/directory/en/profile/arkasian">Vasiliki Kassianidou</a>, Dr <a href="https://www.cut.ac.cy/faculties/fet/ceg/staff/athos.agapiou/?languageId=1">Athos Agapiou</a>, Dr <a href="https://www.researchgate.net/profile/Zomenia-Zomeni">Zomenia Zomeni</a>, Professor <a href="https://globalecologyflinders.com/people/#DIRECTOR">Corey Bradshaw</a><br>
<strong>collaborators</strong>: Dr <a href="https://globalecologyflinders.com/people/#COORDINATOR">Frédérik Saltré</a>
<br>
## Project overview
Project MIGRATE seeks to offer novel insights into population dynamics and range shifts that resulted in dispersals from the Eastern Mediterranean mainland to the island of Cyprus at a critical period (Late Pleistocene, 45-12 ka) through stochastic spatial modelling. This advanced modelling will  enhance our understanding of timing, and climatic and social factors important in the initial colonisation of Cyprus. The proposed project aims to establish new research domains in the field of Cypriot archaeology extending traditional chronological frontiers beyond the Holocene (current Warm Period), encompassing innovative and interdisciplinary methodologies at the forefront of archaeological research.

## Scripts
R code by Corey Bradshaw (<a href="http://github.com/cjabradshaw">@cjabradshaw</a>) and Frédérik Saltré (<a href="http://github.com/fredsaltre">@fredsaltre</a>)

### Climate hindcasts
- <code>HadCM3outputs.R</code>: net primary production (kg C m<sup>-2</sup> year<sup>-1</sup>), temperature (° C), and precipitation (mm year<sup>-1) data (raw and anomalies relative to the present) hindcasted from the <a href="https://www.metoffice.gov.uk/research/approach/modelling-systems/unified-model/climate-models/hadcm3">Hadley Centre Coupled Model version 3</a> (HadCM3) climate model for the 0.5° × 0.5° lat/lon cells covering Cyprus from 20 ka to the present

### Chronology
- <code>calibdate.R</code>: calibrate oldest human archaeological evidence for Cyprus, and use Signor-Lipps correction to estimate window of first entry

### Demography
- <code>CypAspatDemModel.R</code>: stochastic, age-structured human demographic projection model
- <code>CypAspatDemModelMVPsim.R</code>: based on model above, vary founding population size to estimate change in probability of quasi-extinction
- <code>CypAspatDemModelStaggeredEntry.R</code>: vary frequency and intervals of arrival of different group sizes (up to minimum viable population size calculated above) to determine minimum size of entering groups and frequency of immigration required to minimise quasi-extinction

### Palaeontology
- <code>dwarf hippo & elephant.R</code>: using dated fossil specimens of dwarf elephant <em>Elaphus cypriotes</em> and dwarf hippopotamus <em>Phanourios minor</em> from Cyprus, estimate Signor-Lipps corrected window of extinction

## Data
### Climate hindcasts
- CyprusRegion(20ka)_NPP(absolutevalues).csv: hindcasted HadCM3 net primary production (absolute values) from 20 ka to present for the general eastern Mediterranean region
- CyprusRegion(20ka)_NPP(anomaliesvalues).csv: hindcasted HadCM3 net primary production (anomalies relative to the present) from 20 ka to present for the general eastern Mediterranean region
- CyprusRegion(20ka)_AnnualMeanTemperature(absolutevalues).csv: hindcasted HadCM3 temperature (absolute values) from 20 ka to present for the general eastern Mediterranean region
- CyprusRegion(20ka)_AnnualMeanTemperature(anomaliesvalues).csv: hindcasted HadCM3 temperature (anomalies relative to the present) from 20 ka to present for the general eastern Mediterranean region
- CyprusRegion(20ka)_AnnualPrecipitation(absolutevalues).csv: hindcasted HadCM3 annual precipitation (absolute values) from 20 ka to present for the general eastern Mediterranean region
- CyprusRegion(20ka)_AnnualPrecipitation(anomaliesvalues).csv: hindcasted HadCM3 annual precipitation (anomalies relative to the present) from 20 ka to present for the general eastern Mediterranean region
  
### Sea level
- cypareaSL.csv: relationship between relative sea level (from the present) and area of the island of Cyprus derived from <a href="http://www.gebco.net">GEBCO</a>
- lambeckESL.csv: global sea level curve from <a href="https://doi.org/10.1073/pnas.1411762111">Lambeck et al. (2014)</a>

### Demography
- world2013lifetable.csv: human demographic data for age-structure determination (see <a href="http://doi.org/10.1038/s41559-019-0902-6">Bradshaw et al. 2019</a>)

### Palaeontology
- elaphus.txt: radiocarbon dates for the dwarf elephant <em>Elaphus cypriotes</em> from Wigand & Simmons (1999) In: Simmons (ed.) <a href="https://link.springer.com/book/10.1007/b109876"><em>Faunal Extinction in an Island Society</em></a>
- phanourios.txt: radiocarbon dates for the dwarf hippopotamus <em>Phanourios minor</em> from <a href="http://doi.org/10.1371/journal.pone.0134429">Zazzo et al. (2015)</a>

<p><a href="https://www.ucy.ac.cy"><img align="bottom-left" src="www/UCypruslogo.png" alt="UCyprus logo" height="50" style="margin-top: 20px"></a> &nbsp; <a href="http://www.dainst.org"><img align="bottom-left" src="www/DAIlogo.png" alt="DAI logo" height="60" style="margin-top: 20px"></a> &nbsp; &nbsp; <a href="https://www.flinders.edu.au"><img align="bottom-left" src="www/Flinders_University_Logo_Horizontal_RGB_Master.png" alt="Flinders University logo" height="40" style="margin-top: 20px"></a> &nbsp; <a href="https://globalecologyflinders.com"><img align="bottom-left" src="www/GEL Logo Kaurna New Transp.png" alt="GEL logo" width="60" style="margin-top: 20px"></a> &nbsp; &nbsp; <a href="https://EpicAustralia.org.au"><img align="bottom-left" src="www/CabahFCL.jpg" alt="CABAH logo" height="50" style="margin-top: 20px"></a> &nbsp; <a href="https://www.cut.ac.cy"><img align="bottom-left" src="www/CUTlogoblack.png" alt="CUT logo" height="50" style="margin-top: 20px"></a></p>
