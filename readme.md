# RobotTeam Papers 
## Robot Team II: Electric Boogaloo (title w.i.p.)
- Discuss real time georectification, generation of reflectance data, etc... 
- Combine Multiple days of observations
- Discuss need for both viewing geometry and solar geometry 
  - make reference to highly nonuniform reflectance as a function of incident angle 
  - make reference to Beer's law as justification for direct determination of concentration of concentration from spectra 
  - in depth discussion of fluorometers (maybe save this for the dissertation)

## Unsupervised Classification of Hyperspectral Imagery for Rapid Characterization of Novel Environments with Autonomous Robotic Teams 
### K-means 
### Self Organizing Maps
- fit an SOM model to the data 
- for each pixel in entire map, assign best matching unit (use distinguishable colors for a $10\times 10$ or $25\times 25$ SOM grid)
- For class, investigate the learned "spectrum" representation and compare against known chemical spectra. 
  - is there a database we could try to use to look up possible species in the reflectance spectra? 
### Generative Topographic Mapping 

## Spectral Indices for Rapid HSI Surveys: Unsupervised and Supervised Methods via SciML
- Apply to PROSPECT database as a simple test case 
- Apply to our own Data 
  - [Mutliple Endmember Spectral Mixture Analysis](https://aviris.jpl.nasa.gov/proceedings/workshops/99_docs/46.pdf)
  - generalize **Spectral Unmixing Models** unmixing models... with gaussian process we could think of an infinite basis of gaussians describing "peaks" in the spectrum". Can we try to kernelize this procedure?

## Synthetic Data Generation for Hyperspectral Imaging with Autonomous Robotic Teams
- Variational Autoencoders 
- Group transformations, e.g. rotations, reflections, translations, cropping, etc... (do these make sense if boat data is point observation) 
- Advanced sampling methods for regions with 

## Uncertainty Quantification via $\partial P$.
- Categorize Methods into two categories: 
  - quantifying uncertainty in collected data 
  - quantifying model uncertainty
- Conformal Prediction (we have a NN code for doing this in flux. Just need to apply it)
- Representativeness Uncertainty i.e. when georectifying HSI images and reducing spatial extend via `ilat` and `ilon` settings, also compute the stdev for each grainy pixel
- Measurements.jl *forward mode* once we have the representativeness uncertainty. 
- Need a way to quantify uncertainty from Boat sensors... 
- Sensativity Analysis with w/ automatic differentiation


# Super Resolution 
## Cloud & Shadow Mask for Sentinel-2
### ML Type 
- Supervised Classification 
### ML Methods 
- Single Pixel w/ Tree Based Methods 
- Deep NN with Convolutional Layers
### Features 
- Sentinel 2 Multi-band Imagery 
- Land Type 
- Viewing Geometry 
- Solar Geometry
### Targets
- Sentinel 2 Cloud Mask + Cloud Shadow Mask

## Cloud & Shadow Fill 
### ML Type 
- supervised regression
### ML Methods 
- pixel based (Would it make sense to do something else here?)
### Features 
- Sentinel 2 Multi-band Imagery 
- Sentinel 2 Cloud Mask & Cloud Shadow Mask 
- Sentinel 1 SAR Variables (GRD or SLC or both?)
- 10 m Digital Elevation Map
- Viewing Geometry 
- Solar Geometry
- Land Type
### Targets
- *Cloudless & Shadowless* Sentinel 2 Multi-band imagery 

## Sentinel 2 RGB *Spatial* Super Resolution 
### ML Type 
- Supervised Regression
### ML Methods 
- This has to be a Deep NN method using convolution to get the upsampling. I don't think we can do this with pixel based models (using tree methods)
- Probably should use a GAN
### Features
- High (spatial) Resolution NAIP RGB Image
- Sentinel 2 Multi-band Imagery 
- Sentinel 1 SAR Variables (GRD or SLC or both?)
- 10 m Digital Elevation Map
- Viewing Geometry 
- Solar Geometry
- Land Type
### Targets
- Sentinel RGB Bands @ NAIP Resolution
### Loss Function Terms
### Notes 
We could make a model that uses all 3 bands (RGB) simultaneously, or we can make a model for a single band that we validate against R, G, and B bands individually. This has the added perc of increasing the training samples. This will be much easier to then apply to *all* bands of the sentinel imagery (and perhaps Sentinel 1, etc...)  independently. We could try: 
- Red, Green, Blue bands separately 
- Black and White converted RGB image
- Data Augmentation via Scaling / Rotation / Reflection

## Sentinel 2 Multiband *Spatial* Super Resolution 
### ML Type 
- Supervised Regression
### ML Methods 
- This has to be a Deep NN method using convolution to get the upsampling. I don't think we can do this with pixel based models (using tree methods)
- Probably should use a GAN
### Features
- High (spatial) Resolution NAIP RGB Image
- Sentinel 2 Multi-band Imagery 
- Sentinel 1 SAR Variables (GRD or SLC or both?)
- 10 m Digital Elevation Map
- Viewing Geometry 
- Solar Geometry
- Land Type
### Targets
- Sentinel RGB Bands @ NAIP Resolution
### Loss Function Terms


## Sentinel 2 Multiband *Spectral* Super Resolution
### ML Type 
- Supervised Regression
### ML Methods 
- This can be pixel based
### Features
- Sentinel 2 Multi-band Imagery 
- Sentinel 1 SAR Variables (GRD or SLC or both?)
- 10 m Digital Elevation Map
- Viewing Geometry 
- Solar Geometry
- UAV Hyperspectral Image
### Targets
- Sentinel *Hyperspectral*  Imagery (i.e. Sentinel at all HSI Bands)
### Loss Function Terms


# Chemical Data Assimilation & ActivePure Work
## ActivePure Research Lab
- Overview of all sensor in sensor matrix 
- Overview of measurement capabilities (list of species, uncertainty levels, etc...)
- Overview of containerized data acquisition pipeline 
  - NodeRed 
  - InfluxDB 
  - Grafana 
  - Quarto
  - Automatic Alerts 
  - Automatic Reports

## Kinetics and Chemical Data Assimilation 
- MCM Implementation in Julia 
- Direct computation of Photolysis rates 
- Combination with Dr. Lary's AutoChem 
- Addition of Ion Chemistry from MIT Lightning disseration
- Visualization of chemical cycles
- SciML methods to infer below detection limits

# Technical Notes 
## Real Time Georectification of Drone Based Imagery
- Georectification of pushbroom HSI 
- Georectifcation of  square visible + thermal FLIR imagery
## Self Organizing Maps 
## Bayesian Optimization with Gaussian Process Regression
## Gaussian Process Regression / Classification in MLJ
## Solar Geometry? (probably not necessary)
