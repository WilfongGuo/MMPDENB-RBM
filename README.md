# Instruction of MMPDENB-RBM package

This package includes Matlab scripts and several datasets for demo of MMPDENB-RBM approach:

‘main_ MMPDENB-RBM.m’ is a Matlab function for the routine of experimental analysis. MMPDENB-RBM aims to identify personalized edge-network biomarkers (multiple PDENBs contains multi-modal PDENBs and PDENBs) for detecting early warning signal of individual patients in cancer and find the drug target genes which can provide effective information for the early treatment of cancer by analyzing the biomarkers.

The input (case: BRCA) include:
(1)Path: The path of the user where the 'Main,m' is located.

(2)Cancer Data: tumor and normal sample data of BRCA.

The output results:
(1)PGIN_BRCA: BRCA patients' personalized gene interaction network construct by SSN method.
i.‘BRCA_i_PGIN.mat’ indicates that the PGIN of the i-th BRCA patient which contains the subnetwork adjacency matrix and the name of the gene in the subnetwork.

(2)PEN_BRCA: BRCA patients' personalized edge network construct by iENA method.
i.‘BRCA_i_PEN.mat’ indicates that the PEN of the i-th BRCA patient which contains the edge-network, gene pair, and the name of the gene in the edge-network.

(3)BRCA_result: Non-dominated solutions of patient samples obtained by MMPDENB-RBM.

i.‘BRCA_sample_i_MMPDENB-RBM_boxchart.mat’ stores non-dominated solutions for the i-th patient by running the evolutionary algorithm 30 times each time.
ii.‘BRCA_sample_i_ MMPDENB-RBM _PF.mat’ indicates that a group Pareto front solutions of i-th patient obtained by performing non-dominated sort on boxchart.
iii.‘BRCA_sample_i_ MMPDENB-RBM _PS.mat’ indicates the PDENB corresponding to Pareto front solutions.
iv.BRCA_PDENB_name.txt: Gene pairs’ name of multi-modal PDENB or PDENB with the highest score  of BRCA patients. For example：i-th patient: patient’s ID. PDENB or multi-modal PDENB: the name of genes in PDENB or multi-modal PDENB.

Suggestions
(1)Hardware suggestions for running this package: Window 10 or above; Matlab 2021 or above; RAM 32G or above.

(2)When users analyzed running this package, please note that:
i.Users should set the path in the program, firstly.
ii.Parameter setting of Popnum, Max_CalNum, and Experiment_num will affect the running time. With default parameters, multi-modal EA in MMPDENB-RBM takes about 40 minutes to identify multiple PDENBs for a BRCA patient. Users can decrease running time by modifying above parameter.
iii.If users want to run their own experimental data, users should add the data of normal samples, tumor samples and mutation samples for the corresponding disease in the 'Code_construct_personalized_network' file.

%    $Id: Main.m Created at 2023-4-27$ 
%   $Copyright (c) 2021 by School of Electrical Engineering, Zhengzhou University, Zhengzhou 450001, China$;
%    $If any problem, please contact zero999396@outlook.com for help. $
