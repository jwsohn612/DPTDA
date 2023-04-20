
# Differentially Private Topological Data Analysis

This is the first work that suggests the differential privacy framework on persistence diagrams.  

- Title: "Differentially Private Topological Data Analysis"
- Authors: Taegyu Kang, Sehwan Kim, Jinwon Sohn, and Jordan Awan

In this work, we introduce privatizing the persistence diagrams on the basis of the exponential mechanism. Because of outliers that could come into play in the mathematical derivation of the sensitivity, we propose using the distance-to-measure (DTM) that is known for robustness in the presence of such outliers. We characterize theoretical properties that DTM brings into the context of the differential privacy. The simulation studies illustrate the sanitizing process via Metropolis-Hasting algorithm coincide with the established theories. Lastly, we apply the algorithm to the real-world data set where the analysis is made on the private topological statistics. 

## Simulation

Refer to 'example.pdf' for testing the proposed method. 

## Required Configuration

- R >= 4.2.1

- Packages : TDA, dplyr, purrr 



