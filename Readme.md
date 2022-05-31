# Purpose, Description

This is the flux calculator component of the IOW ESM.
All coupling fields that are communicated between atmospheric and bottom models are passed through or calculated by this component.
The commnuicated information is mapped from the models' grid to the exchange grid and vice versa. 


# Authors

* SK      (sven.karsten@io-warnemuende.de)
* HR      (hagen.radtke@io-warnemuende.de)


# Versions

## 1.01.00 (latest release)

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