# RobotTeam Papers 
## Robot Team II: Electric Boogaloo (title w.i.p.)
- Discuss real time georectification, generation of reflectance data, etc... 
- Combine Multiple days of observations
- Discuss need for both viewing geometry and solar geometry 
  - make reference to highly nonuniform reflectance as a function of incident angle 
  - make reference to Beer's law as justification for direct determination of concentration of concentration from spectra 
  - in depth discussion of fluorometers (maybe save this for the dissertation)

### Previous outline: 
1. Introduction
   - Human Impact / Importance 
   - Hyperspectral Imaging and Remote sensing
     - Satellites with HSI 
     - Drones with HSI 
     - applications 
       - Spectral Indices + Vegetation Studies
       - Machine Learning with HSI
   - Remote Sensing applications to characterization of bodies of water
     - meteorology 
     - atmospheric sensing 
     - need for in-istu sensing
   - Autonomous Robotic Teams 
     - coordination
     - autonomous data collection 
2. Materials & Methods 
   - Robotic Team 
   - Fluorometers (technology and our included sensors)
   - Data collection procedure
     - AV -> Boat -> AV -> Maps -> Verification w/ Boat
   - Software 
     - Julia & MLJ
     - Docker + Influxdb + NodeRed + Grafana for live visualizations
   - Processing of HSIs 
     - radiance + downwelling irradiance -> reflectances 
       - assumption of *Lambertian* surface
     - georectification 
       - georectification = orthorectify + georeference
     - Other important variables 
       - Viewing geometry (roll, pitch, yaw)
       - Solar geometry i.e. incident geometry (solar azimuth, solar elevation)
         - should visualize this... see the [wikipedia image](https://en.wikipedia.org/wiki/Solar_zenith_angle#/media/File:Solar_Zenith_Angle_min.png)
       - mention grazing incidence reflection and ideal time for data collection
     - Data Collocation procedure 
       - `ilat`, `ilon` and representativeness uncertainty
     - Machine Learning Methods 
       - Individual learners 
       - Feature importance ranking 
         - discuss pros/cons and mentions that Shapley values take too long to compute to be helpful here 
       - Hyperparameter optimization
       - Super learner / Model stacking 
         - discuss training on complimentary folds
       - MLJ model search and data pipeline (i.e. learning networks and DAGs)
         - we should make a graph to visualize the data pipeline for the superlearner stack 
3. Results 
   - Exploratory Data Analysis results 
   - Combining multiple days of data (i.e. the difference between this and the previous paper)
   - Comparison of vanilla models, hpo models and super learner
   - Updated maps of Scotty's ranch 
4. Discussion

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
    

# Technical Notes 
## Real Time Georectification of Drone Based Imagery
    - Georectification of pushbroom HSI 
    - Georectifcation of  square visible + thermal FLIR imagery

## Self Organizing Maps 
## Bayesian Optimization with Gaussian Process Regression
## Solar Geometry? (probably not necessary)
