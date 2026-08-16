Manuscript: [Insert Paper Title Here]
Authors:    Mariaelena Bottazzi Schenone, Ilaria Mozzetta, 
            Tiziano Iannaccio, and Maurizio Vichi

-------------------------------------------------------------------------------
1. SCOPE AND CORE CONTRIBUTION
-------------------------------------------------------------------------------
This repository provides the computational routines and benchmark datasets associated 
with the manuscript cited above. 

The primary theoretical and methodological contribution of this work is EXCLUSIVELY 
the Interpretability Index (I) for Exploratory Factor Analysis (EFA), implemented 
in the core evaluation function:

    - interp.m : Computes the Interpretability Index (I), along with its 
                 constituents (ECV and weighted-mean Factor Complexity/GS_w), 
                 given a factor pattern matrix and a factor correlation matrix.

The function `interp.m` is the sole routine directly tied to the theoretical 
formulations, empirical claims, and performance guarantees presented in the paper.

-------------------------------------------------------------------------------
2. ADDITIONAL UTILITIES AND REPLICATION DRIVERS
-------------------------------------------------------------------------------
To minimize the user learning curve and enable out-of-the-box replication of 
our empirical results, this repository also includes a comprehensive model-evaluation 
pipeline (e.g., `interpEval.m`, `processData.m`, rotation suites, and MATLAB–R 
bridge scripts), along with 5 empirical benchmark datasets.

Please note the following regarding these supplementary utilities:

  * Scope Boundaries: Methodological issues preceding the post-rotation 
    evaluation step—such as data transformation/cleaning, matrix smoothing, 
    alternative correlation estimations, or standalone EFA extraction logic—are 
    outside the scope of the manuscript.
  
  * EFA & Rotation Engines: Helper routines such as `paf.m` (Principal Axis 
    Factoring), `oblimin.m`, and `srot.m` (sparse rotations called via `rotate.m`) 
    are provided strictly as a convenience for the community to facilitate 
    comparative analysis.
  
  * Pipeline & Data Assumptions: The master driver (`interpEval.m`) and its 
    data-preprocessing workflows (including the MATLAB–R bridge for polychoric 
    correlations) are engineered specifically for the structure of the included 
    datasets. They are not intended as a universal data-processing framework.

For custom data formats, non-standard preprocessing, or alternative extraction 
methods, users are encouraged to interface their own pre-extracted pattern and 
factor correlation matrices directly with `interp.m`.

-------------------------------------------------------------------------------
3. SETUP AND DOCUMENTATION
-------------------------------------------------------------------------------
For step-by-step instructions on setting up the local environment and configuring 
the MATLAB–R bridge required to execute the pipeline on ordinal/Likert data and 
sparse rotations, please consult the accompanying manual:


    README.pdf
