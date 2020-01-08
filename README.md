# Mhealth
This repository contains the code for the paper "Causal Effect Estimation and Optimal Dose Suggestions in Mobile Health" by Liangyu Zhu, Wenbin Lu and Rui Song (2019)

For questions about the code, please contact Liangyu Zhu: lzhu12@ncsu.edu

## Explanation for the files


- simulation-runmodel.R contains the code for running simulations presented in the paper. It depends on the functions defined in simulation-functions.R

- realdata-preprocess.R contains the code for preprocessing Ohio type 1 diabetes dataset. It depends on the functions defined in realdata-preprocess-functions.R

- realdata-runmodel.R contains the code for running the proposed the method with Ohio type 1 diabetes dataset. It depends on the functions defined in realdata-runmodel-functions.R

Information of Ohio type 1 diabetes dataset can be found in http://smarthealth.cs.ohio.edu/OhioT1DM-dataset.html and in the paper Marling, C. and Bunescu, R.C., 2018, July. The OhioT1DM Dataset For Blood Glucose Level Prediction. In KHD@ IJCAI (pp. 60-63).