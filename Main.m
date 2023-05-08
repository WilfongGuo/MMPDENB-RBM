clc
clear
%%               MMPDENB-RBM
%                  Multi-modal optimization to identify personalized edge-network biomarkers for detecting early warning signal of individual patients in cancer
%
%                  Attention : Please raad 'ReadMe' file before running MMPDENB-RBM_code and add  folder and subfolders to the path     
%
%                  input  :
%                         cancer_data     :    BRCA or LUNG
%                         path            :    The path of the user where the 'Main,m' is located
%
%                   output :         
%                          PGIN_cancer             : the personalized gene interaction network of cancer patient
%                          PEN_cancer              : the personalized edge network of cancer patient
%                          cancer__result          : non dominated solution of patient samples obtained by MMPDENB-RBM
%                          cancer_PDENB_name       : gene name of MMPDENB or PDENB of cancer patients
%
%                   Proposed by J.J.Liang, Wei feng, Guo, and Zongwei Li in Zhengzhou University on April 27, 2023.
%                   If any problem,pleasse contact zero999396@outlook.com for help.
                  
%%                                                  !!!!!!!!!Attention!!!!!!!!!!!!

%% !!!!!!!!!Extract the compressed file (BRCA_normal_gene_data.7z and BRCA_tunor_gene_data.7z)to the current folder before running!!!!!!!!!!!

%% ***************************Input*************************************


path='F:\MATLAB\MMPDENB-RBM_code\';        % The path of 'Main,m' on the user's computer and '\' need reserve.

unzip('BRCA_tumor.zip');              % Take BRCA as an example
unzip('BRCA_normal.zip');

expression_tumor_fileName = strcat('BRCA_tumor.txt'); 
expression_normal_fileName = strcat('BRCA_normal.txt');

%% *******************The function of MMPDNB****************************

%% default parameters
popnum=300;
Max_CalNum=30000;
Experiment_num=30;

[gen_name,path3]=MMPDENB_RBM(expression_tumor_fileName,expression_normal_fileName,path,popnum,Max_CalNum,Experiment_num);

%% *********Output**********

data2=load([path,expression_tumor_fileName(1:4),'_clinical_stage_information.mat']);
output_file_name=strcat(path3,expression_tumor_fileName(1:4),'_PDENB_name.txt');
fid = fopen(output_file_name,'w');
        for j = 1:length(gen_name)
            k=1;
            fprintf(fid ,num2str(j));
            fprintf(fid ,'%s', ':');
            patient_name=cell2mat(data2.Final_Sample_name(j,1));
            fprintf(fid ,patient_name);
            fprintf(fid ,'%s', ' ');
            while k<=length(gen_name{j,1})
                if length(gen_name{j,1})>1
                    MM=strcat('multi-modal PDENB',num2str(k));
                else
                     MM=strcat('PDENB');
                end
                fprintf(fid ,MM); 
                fprintf(fid ,'%s', ': ');
                l=1;
                while l<=length(gen_name{j,1}{k,1})
                n1=cell2mat(gen_name{j,1}{k,1}(l,1));
                fprintf(fid ,n1);
                fprintf(fid ,'%s', '-');
                n2=cell2mat(gen_name{j,1}{k,1}(l,1));
                fprintf(fid ,n2);
                fprintf(fid ,'%s', '\');
                l=l+1;
                end
            fprintf(fid ,'%s', ';');
            k=k+1;
            end
            fprintf(fid ,'\n');
        end
    fclose(fid);