# PCF_RNMF
Readme: This fold includes the codes (Matlab and R programing) used for simulated datasets and data analysis  of A novel computational framework for integrating multidimensional data to enhance the accuracy in predicting prognosis of colorectal cancer.

1. Core code
NMF_re.m non-negative matrix factorization
RNMF_re.m random non-negative matrix factorization

2. Main function: The following four files are the main code of 4 simulated datasets, which can be run directly.
simulated dataset 1: main_sim_no_noise_sam45.m
simulated dataset 2: main_sim_no_noise.m
simulated dataset 3: main_sim_with_noise_sam45.m
simulated dataset 4: main_sim_with_noise.m


3. The OS prognosis of colorectal cancer patients collected by TCGA. 
In the "CORE-OS" folder, the data analysis code "OS-code.R" can be run directly, "CORE-OS.RData" is the file that saves the data used by RNMF, and "SC_RSF-OS.R" is the code about the random survival forest screens variables.

4. Recurrence prediction of stage II/III CRC colorectal cancer
In the "CORE-DFI" folder, the data analysis code "DFI-code.R" can be run directly, and "SC_RSF-DFI.R" is the code about the random survival forest screens variables.


