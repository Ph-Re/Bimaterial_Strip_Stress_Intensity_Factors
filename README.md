# Stress Intensity Factor Evaluation for Bi-Material Beams subjected to Three Point Bending and Thermal Loads
Please cite the Publication "Influence of microstructural architecture on the fracture resistance of W-Re focal tracks", if you use this code. See the upcoming doi for further details.

<!-- This code is citable via Zenodo:

[![DOI](https://zenodo.org/badge/DOI/**doi**.svg)](https://doi.org/**doi**)   -->


The simulations were tested on Windows 11 using Abaqus 2019. 

All variations of the given simulation scripts, treat the same basic setup, and will be displayed for one case of interest in a following figure after paper publication:

![Image of the geometry](./setup.jpg)

For the individual simulations, this setup is altered as necessary. The relevant changes are given in the following sections. Common to all simulations are the variables for geometry control:

*L* ... `beam_length`

*S* ... `support_span`

*h<sub>1</sub>* ... `wre_height`

*h<sub>2</sub>* ... `tzm_height`

*a* ... `crack_length`

*F* ... `force_magnitude`

*µ<sub>1</sub>* ... Converted to elastic modulus (`wre_elastic_modulus`) via 2 *µ<sub>1</sub>* (1+$\nu_1$)

$\mu_2$ ... Converted to elastic modulus (`tzm_elastic_modulus`) via 2 *µ<sub>2</sub>* (1+$\nu_2$)

$\nu_1$ ... `wre_poisson_ratio`

$\nu_2$ ... `tzm_poisson_ratio`

Furthermore, not shown in the figure, the coefficients of thermal expansion are defined:

$\alpha_1$ ... `wre_expansion_coeff`

$\alpha_2$ ... `tzm_expansion_coeff`

## 1. Stress Intensity Factor Calculation for Three Point Bending Loads using Singular Elements
A parametric calculation of stress intensity factors of a bi-material beam subjected to three point bending (3PB) load is performed. 

For this the file *bimaterial-beam-thermal-3PB-SIFs.py* is used. To switch to 3PB loads the variable `apply3PBending` is set to **True**, while `applyThermalLoad` is set **False**. Geometry, defined loads and mesh settings (including element types) are adaptable at the beginning of the script. In the script, their individual function is also explained. No settings for plasticity are availabe, as the task is to compare against linear elastic results from the distributed dislocation method, given in the paper.

A for-loop is used to alter the crack length. The range of crack lengths should be adapted directly in the loop at the end of the script, see Line **407** - `crack_length`. 

The `saveAll` variable at the beginning of the document enables the options to save the full model for each crack length, or to delete them after each iteration.

The crack tip is defined using the ABAQUS Engineering Feature. Singular elements are employed at the crack tip, producing the correct square root singularity. The resulting K-Factors are evaluated along 10 contours. The mesh size near the crack tip is adaptive to the distance from other surfaces to ensure enough valid contours.

The results for all contours and additional info is saved in form of pickle file and JSON. The specifics are given at the end of the script.

## 2. Stress Intensity Factor Calculation for Thermal Stresses using Singular Elements
This parametric study follows directly from **1.**

The file *bimaterial-beam-thermal-3PB-SIFs.py* is used with `apply3PBending` is set to **False**, while `applyThermalLoad` is set to **True**.

A special note on the definition of the thermal expansion coefficients is given in the publication and also the code. When assuming isotropic thermal expansion and plane strain condition, then large out-of-plane stresses ($\sigma_{33}$) are generated, as the sample deformation is fully restricted in this direction (think, infinitely thick sample). This is physically unrealistic for typical sample sizes. Nontheless, the stress state at the crack tip can ressemble a plain strain state locally. To consider this case, the thermal expansion coefficients are defined as orthotropic, with vanishing component in the out-of-plane direction.

The same results as in **1.** are generated.

## 3. Stress Intensity Factor Calculation for Combined Loads (Thermal + 3PB) considering Plasticity and Residual Stresses without use of Singular Elements
This simulation considers Plasticity in form of perfectly plastic yielding and Residual Stresses by performing thermal steps.
For this the file *bimaterial-beam-thermal-residualStress-plasticity.py* is used.

The temperature-dependent yield stresses are given for the respective materials via the variables `wre_plasticity` and `tzm_plasticity`.

The simulation performs the following steps:
- CoolDown - The sample is cooled from `T_initial` to `T_coolDown`, representing the initial annealing treatment. No crack is present yet, as the crack is introduced after sample production. To prevent the need of remeshing, the crack is already introduced, but kept closed using an appropriate Contact Interaction.
- HeatUp - First, the Contact Interaction along the crack front is deactivated, then the sample is heated from `T_coolDown` to the testing temperature `T_testing`.
- Establish Contact - Contact Interactions between Plunger / Bending Table and the sample are activated, then a load, typically 1/100 - 1/10 of the total force `force_magnitude`, is applied. This step is performed for better convergence of initial contact.
- Loading - Total force `force_magnitude` is applied.

Thermal expansion definition follows the same considerations as detailed in **2.**

Some care needs to be taken in step control to ensure convergence.

In this simulation the $K_I$ factors are **not** automatically saved to a file.
$K_I$ factors are evaluated. ABAQUS gives a warning for this, which is of no further importance, as they are calculated via J-Integral anyways. 

### Unit System
LENGTH ... mm

FORCE ... N

TEMPERATURE ... K (equiv. °C)

STRESS ... MPa

STRESS INTENSITY FACTOR ... MPa $\sqrt{mm}$


**Codes by Philipp Reindl**
