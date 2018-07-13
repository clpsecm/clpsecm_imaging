# Imaging algorithm of CLP-SECM data
---
## 1. Overview
This code package is dedicated to image reconstruction from scanning chemical microscopic imaging (SECM) data using continuous line probe (CLP). A CLP-SECM scan starts from putting the sample on a rotational stage and locate a line probe
The reconstruction algorithm first formulate the task of locating the reactive species as a variation of LASSO problem, a class of optimization problem, using the 

## 2. Data Management
Create a folder clpsecm_data in the same path parallel to folder clpsecm_algorithm/. Put the data in the folder with the naming convention
```
clpsecm_(MMDDYY)_S(sample_number)_L(number_of_lines)_(data_number).xlsx
```
For example, a CLP-SECM file is generated in date July 04, 2076, with scans the ten dots sample with sample number 3-2. The scan choose 10 lines with different angles, and is the second experiment in that same day under same setting, should has the file name defined as
```
clpsecm_077476_S032_L10_02.xlsx
```


## 3. Overview of Code Package 
      DataSpec


