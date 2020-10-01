#fdEIT
The GaussNewton method is one way to map surface measurements to a 
conductivity image. Usually, the minimization function contains only weighted 
differential measurement data and a regularization. This traditional method is 
extended by absolute measurement data to improve fdEIT reconstruction results. 
The key challenge of unknown torso geometries and electrode displacement has 
been addressed for the reconstruction of different lung pathologies.

Folders with the name results are always used for computational results and 
won't be stored in the git repo.

The repository is structured as followed:

- fdEIT_overview.pdf	- conceptional slides about the algorithm

- 0_tissue_properties
	- Models and simulation uses reference conductivites of tissues/organs
	  from gabriel database. Functions to generate frequency specific
	  conductivities are given here.
	  
- 1_models
	- FEM-models can be generated with eidors and netgen. Functions are
	  defined here and should be used for the reconstruction process.

- 3_reconstructions
	- 1_case1_healthyThorax
	  Performs basic reconstruction of simulated data from a healthy thorax
	  model. Hyperparameter and betaparameter sweep is conducted to 
	  optimze via residual/regularization norm.
	- 1_all_pathology_models
		- This folder includes reconstructions of a healthy thorax as well as 9
		  models of different lung pathologies.
		- The premade models and reconstructions are stored in different
		  folders. Reconstructions were made using the skip-4-4 and skip-3-3
		  patterns.
		- Functions are included for creating the different pathology models
		  from scratch, creating a suitable thorax fem, reconstructing the 
		  modelled data and plotting all reconstructions.
		  
- 4_evaluation
	- 5_Noise_evaluation
	- 6_FOM_evaluation
		- Compares the fdEIT algorithm with BP, GN and eidors against the
		  figures of merit. Optimal parameters for fdEIT reconstruction are
		  derived from regularization/residual norm criteria.
		  -> code is under development
	- 3_measurement_evaluation
		- This folder contains tank data that was measured with the AixTOM and
		  SwissTOM measurement devices. Also the conductivity data of potato
		  and pumpking is included. It was measured using the agilent E4980A
		  precision lcr meter.
		- The measurement data is stored in raw and indexed form. For further
		  information on the indexing consult the readme included in the folder.
		- The folder 'scripts' includes all necessary scripts for analysing
		  the different measurement data sets.
	  
- 5_classification:
    - The algorithm itself with different functions including the main
      classification function and functions for image segmentation by 
      magnitude and phase. There's also a function for plotting the result
      as well as a calculated safety factor. A stored colormap defines
      the color coding for the safety factor.
    - Some basic evaluation code that analyses existing reconstructions to
      tell how close they are to the reconstructed model in shape and 
      magnitude values.
    - A test script that is used as a main function to classify different
      reconstructions.

-6_hpselection:
   - Specific files or functions for three hyperparemeter selection approaches 

