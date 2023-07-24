# Purpose, Description

This is the flux calculator component of the IOW ESM.
All coupling fields that are communicated between atmospheric and bottom models are passed through or calculated by this component.
The commnuicated information is mapped from the models' grid to the exchange grid and vice versa. 


# Authors

* SK      (sven.karsten@io-warnemuende.de)
* HR      (hagen.radtke@io-warnemuende.de)


# Versions

## 1.03.00 (latest release)

| date        | author(s)   | link      |
|---          |---          |---        |
| 2023-07-24  | SK          | [1.03.00](https://git.io-warnemuende.de/iow_esm/components.flux_calculator/src/branch/1.03.00)     | 

<details>

### changes
* enable bias corrections for exchanged fluxes
  * currently only the evaporation flux can be corrected in an additive fashion
  * however, module `bias_corrections` can be easily extended
  * correction is applied as anual cycle, i.e. for each month there should be a correction file (NetCDF) that contains
    a field that is added to the actual flux
  * these files must be placed in folder called `corrections` that is placed in the input folder called `flux_calculator`
  * the files must be named `mass_evap-MM.nc` for the evaporation correction
  * the `MM` indicates the month and goes from `01` to `12`
* the Fortran code can now call python routines via the Fortrasn module `call_python`
  * currently this used to calculate the current month from the initial start time and the current time step
*build script templates for new target machine have been added
    
### dependencies
* OASIS3-MCT libraries
* see build scripts for more dependencies  
  
### known issues
* none so far

### tested with
* intensively tested on both HLRN machines
  * using example setups available under:
    (coupled) /scratch/usr/mviowmod/IOW_ESM/setups/
              MOM5_Baltic-CCLM_Eurocordex/example_8nm_0.22deg/1.00.00
         and  https://zenodo.org/record/8167743/files/1.00.00.tar.gz (https://doi.org/10.5281/zenodo.8167743)              
* can be built and run on Haumea but output is not intensively tested
  
</details>

<details>
<summary><b><i>older versions</i></b></summary>

## 1.02.00

| date        | author(s)   | link      |
|---          |---          |---        |
| 2022-12-22  | SK          | [1.02.00](https://git.io-warnemuende.de/iow_esm/components.flux_calculator/src/branch/1.02.00)     | 

<details>

### changes
* fluxes can now be calculated according to
  * the MOM5 ocean model (formulars the same as in CCLM but with different transfer coefficients)
  * the RCO ocean model, Meier et al., SMHI REPORTS OCEANOGRAPHY No. 26, August 1999
  * caution this is still a bit experimental
    
### dependencies
* OASIS3-MCT libraries
* see build scripts for more dependencies  
  
### known issues
* none so far

### tested with
* intensively tested on both HLRN machines
  * using example setups available under:
    (coupled) /scratch/usr/mviowmod/IOW_ESM/setups/
              MOM5_Baltic-CCLM_Eurocordex/example_8nm_0.22deg/1.00.00
* can be built and run on Haumea but output is not intensively tested
  
</details>

## 1.01.00

| date        | author(s)   | link      |
|---          |---          |---        |
| 2022-05-31  | SK          | [1.01.00](https://git.io-warnemuende.de/iow_esm/components.flux_calculator/src/branch/1.01.00)     | 

<details>

### changes
* flux calculator can now run in parallel
* each flux calculator instance is responsible only for a small part of the exchange grid
    
### dependencies
* OASIS3-MCT libraries
* see build scripts for more dependencies  
  
### known issues
* none so far

### tested with
* intensively tested on both HLRN machines
  * using example setups available under:
    (coupled) /scratch/usr/mviowmod/IOW_ESM/setups/
              MOM5_Baltic-CCLM_Eurocordex/example_8nm_0.22deg/1.00.00
* can be built and run on Haumea but output is not intensively tested
  
</details>

## 1.00.00 

| date        | author(s)   | link      |
|---          |---          |---        |
| 2022-01-28  | SK, HR      | [1.00.00](https://git.io-warnemuende.de/iow_esm/components.flux_calculator/src/branch/1.00.00)       | 

<details>

### changes
* initial release
* flux calculator can couple MOM5 ocean model and CCLM atmospheric model
* usinf the input file `flux_calculator.nml` it can create a namcouple file that is used by the OASIS3 coupler
    
### dependencies
* OASIS3-MCT libraries
* see build scripts for more dependencies  
  
### known issues
* none so far

### tested with
* intensively tested on both HLRN machines
  * using example setups available under:
    (coupled) /scratch/usr/mviowmod/IOW_ESM/setups/
              MOM5_Baltic-CCLM_Eurocordex/example_8nm_0.22deg/1.00.00
* can be built and run on Haumea but output is not intensively tested
  
</details>

</details>