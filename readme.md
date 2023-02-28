# RobotTeam Papers 
## Robot Team II: Electric Boogaloo (title w.i.p.)
- Discuss real time georectification, generation of reflectance data, etc... 
- Combine Multiple days of observations
- Discuss need for both viewing geometry and solar geometry 
  - make reference to highly nonuniform reflectance as a function of incident angle 
  - make reference to Beer's law as justification for direct determination of concentration of concentration from spectra 
  - in depth discussion of fluorometers (maybe save this for the dissertation)

## Unsupervised Classification of Hyperspectral Imagery for Rapid Characterization of Novel Environments with Autonomous Robotic Teams 
### K-means / Fuzzy K-means
### Self Organizing Maps
- fit an SOM model to the data 
- for each pixel in entire map, assign best matching unit (use distinguishable colors for a $10\times 10$ or $25\times 25$ SOM grid)
- For class, investigate the learned "spectrum" representation and compare against known chemical spectra. 
  - is there a database we could try to use to look up possible species in the reflectance spectra? 
### Generative Topographic Mapping 
### analysis ideas
- The nice feature of both the SOM and the GTM is we can reinterpret them to be spectral-unmixing models. Each SOM node has a feature vector of identical length to the input vector. Similarly, the mean projection $y(x;W)$ in GTM represents the center of a gaussian in *data space*. Thus, we can reinterpret these feature vectors and gaussian centers as representative endmember spectra for the dataset. Further, the topographic properties of these methods ensure we have similarity between classes (at least in the latent space). For the GTM we are guarenteed that the data space projections of our GTM nodes will be similar. **We should see if this holds for SOM**. We can also interpret the cluster means from K-means as our endmembers. 
- For SOM and GTM once we have fit the models, we can look at the map of "winning nodes" (BMU for SOM and Mean for GTM) and perform a secondary clustering (via K-means or DBSCAN). These new clusters can define our endmember-bundles for a spectral-unmixing model.
- Once we have these maps, we should color the point by which day the data came from to see if there are any interesting groups or if the data are distributd across observation days 
- We should make sure we apply the method for unsupervised classification to the dye-released images (see if the maps mirror the diffusion of dye). Pick a day where we have >1 dye flight. Make a map for both cases to illustrate how rapidly fitting an unsupervised model on the fly can enable near real time tracking of plume evolution (good for defense, oil spills, etc...). --> can we then construct a vector field / flow field from the difference and predict the plume evolution? Maybe this would be a good excuse to try lagrangian particle tracking... 
- Provided our GTM/SOM fits, we should try a secondary 
- From nodes/node clusters, can we identify spectral endmembers that represent chemicals we measure (e.g. chlorophyll)? 
- From node activations (SOM) or responsabilities (GTM) can we fit a good model that competes with predictions from full spectra? (dimensionality reduction demonstration)
- cluster viewing geometry separate from refelctances and use resulting viewing geo classes to color GTM/SOM map of reflectance data

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


# Sensor Network + SciML
## Evaluation of local chaos in SharedAirDFWNetwork 
- This gives me an excuse to work with the sensor data 
- Train GTM, SOM, and Variational Autoencoder to produce lower dimensional representation of all data from a central node, e.g. in $\mathbb{R}^2$. 
  - For VAE, test a range of dimensions from the number of sensors down to 2 (better for visualization)
- Analyze the variety of methods from [DataDrivenDiffEq.jl](https://docs.sciml.ai/DataDrivenDiffEq/stable/) to infer dynamics in the low dimensional space 
- Can we infer some kind of Hamiltonian from the data and do a HamiltonianNN approach? 
  - Start of with a standard kinetic-energy style Hamiltonian e.g. $\sum_i \frac{1}{2} \dot{x}_i^2$ where $x_i$ is the 
  - use DataDrivenDiffEq approach to learn the associated potential energy term 
  - alternatively, attempt to capture diurnal cycle (or other relevant time scales) by *learning* coordinate representation that forces dynamics to be uncoupled harmonic oscillators a la Hamilton-Jacobi theory. 
  - Test if this hamiltonian NN model can then be transfered to another central node with an appropriate shift in the "total energy" 
- Attempt to analyze the 2d data to infer Koopman operator. We should treat the original sensor values as observables on which the learned koopman operator acts. This should be doable if the NN is just a function. 
- use DMD appraoch to identify a "forcing" coordinate that can identify when we switch nodes as in [Chaos as an intermittently forced linear system
](https://www.semanticscholar.org/paper/Chaos-as-an-intermittently-forced-linear-system-Brunton-Brunton/2efad3ddade8be144ec43d22a9f2992ab036a923)
- [video on Physics Informed DMD](https://www.youtube.com/watch?v=lx-msllg1kU&ab_channel=SteveBrunton)
- [Deep Learning to Discover Coordinates for Dynamics: Autoencoders & Physics Informed Machine Learning](https://www.youtube.com/watch?v=KmQkDgu-Qp0&list=PLMrJAkhIeNNQ0BaKuBKY43k4xMo6NSbBa&ab_channel=SteveBrunton)

- NOTE: we may need to impute missing values. We shoud do so with either my GPR code or with other ML methods + ConformalPrediction. Provided uncertainty estimates, we should then think about how to propagate errors through our analysis via [Measurements.jl](https://github.com/JuliaPhysics/Measurements.jl), [IntervalArithmetic](https://juliaintervals.github.io/pages/packages/intervalarithmetic/), 

# Other things to think about
- [Sparse Nonlinear Models for Fluid Dynamics with Machine Learning and Optimization](https://www.youtube.com/watch?v=z_CZ_VyMDXE&ab_channel=SteveBrunton)
- [Residual Dynamic Mode Decomposition: A very easy way to get error bounds for your DMD computations](https://www.youtube.com/watch?v=Dc2OYDVp43I&ab_channel=SteveBrunton)
- [Deep Learning of Dynamics and Coordinates with SINDy Autoencoders](https://www.youtube.com/watch?v=WHhDgxkiR9c&t=255s&ab_channel=SteveBrunton)
- [Identifying Dominant Balance Physics from Data](https://www.youtube.com/watch?v=U6eZOzHLSM0&list=PLMrJAkhIeNNRK8orLUJL86k7BcuVX11dl&index=4&ab_channel=SteveBrunton)
