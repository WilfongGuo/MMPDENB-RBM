function [gen_name,path3]=MMPDENB_RBM(expression_tumor_fileName,expression_normal_fileName,path,popnum,Max_CalNum,experiment_num)
%% *************************tumor data****************************

[tumor,~,~]=importdata(expression_tumor_fileName);
gene_list=tumor.textdata(2:end,1);
Tumor=tumor.data;

%% *************************normal data***************************

[normal,~,~]=importdata(expression_normal_fileName);
Normal=normal.data;

%% *************************network data**************************

if isequal(expression_tumor_fileName(1:4),'BRCA')
    Z=load([path,'Construct_personalized_edge_network\',expression_tumor_fileName(1:4),'_mutation_data.mat']);
    Z=Z.Z;
elseif isequal(expression_tumor_fileName(1:4),'LUNG')
    Z=load([path,'Construct_personalized_edge_network\',expression_tumor_fileName(1:4),'_mutation_data.mat']);
    Z=Z.Z;
end

N1=length(gene_list);
[N2,~]=size(Z);
Net=zeros(N1);
for i=1:N2
    
    Net(Z(i,2),Z(i,1))=1;
    Net(Z(i,1),Z(i,2))=1;
end

%calculate the adjacency matrix of PPI
N1=length(gene_list);
[N2,~]=size(Z);
Net=zeros(N1);
for i=1:N2
    
    Net(Z(i,2),Z(i,1))=1;
    Net(Z(i,1),Z(i,2))=1;
end

%% Construct PGIN and save
% Construct PE and save
path1=strcat(path,'PGIN_',expression_tumor_fileName(1:4),'\');
mkdir(path1);
Ref=Normal;
% for i=2:size(Tumor,2) %all patient samples
%     [subnetwork_genes,subnetwork_adjacency] = construct_PGIN(i,Normal, Tumor,gene_list,Net,Ref);
% 
%     filename=strcat(path1,expression_tumor_fileName(1:4),'_',num2str(i),'_PGIN.mat');
%     save(filename,'subnetwork_genes','subnetwork_adjacency');
% end

% Data processing
reference=strcat(path,expression_tumor_fileName(1:4),'_normal_gene_data.xlsx');
[~,~,~] = delete_non_expression(path1,path1,expression_tumor_fileName(1:4),reference);

% Construct PEN and save
path2=strcat(path,'PEN_',expression_tumor_fileName(1:4),'\');
mkdir(path2);

referenceT=strcat(path,expression_tumor_fileName(1:4),'_tumor_gene_data.xlsx');
referenceN=strcat(path,expression_tumor_fileName(1:4),'_normal_gene_data.xlsx');
for j=2:size(Tumor,2) %all patient samples
    data=load([path1,expression_tumor_fileName(1:4),'_',num2str(j),'_PGIN' '']);
    [subnetwork_genes,edge_net,gg] = construct_edge_network(data,referenceT,referenceN,j);
    filename=strcat(path2,expression_tumor_fileName(1:4),'_',num2str(j),'_PEN.mat');
    save(filename,'subnetwork_genes','edge_net','gg');
end
%% identification PDENB by multi_modal EA in MMPDENB_RBM

path3=strcat(path,expression_tumor_fileName(1:4),'_result\');
mkdir(path3); % output path
%
numf = dir(path2);
for fnum = 1:length(numf)-2 % all patient samples
    data=load([path2,expression_tumor_fileName(1:4),'_',num2str(fnum) '_PEN']);
    
    [PF,PS,Non_dominated_sol]=multi_modal_EA_in_MMPDENB_RBM(data,popnum,Max_CalNum,experiment_num);
    
    filename=strcat(path3,expression_tumor_fileName(1:4),'_sample_',num2str(fnum),'_MMPDENB-RBM_PF.mat');
    save(filename,'PF');
    filename2=strcat(path3,expression_tumor_fileName(1:4),'_sample_',num2str(fnum),'_MMPDENB-RBM_PS.mat');
    save(filename2,'PS');
    filename3=strcat(path3,expression_tumor_fileName(1:4),'_sample_',num2str(fnum),'_MMPDENB-RBM_boxchart.mat');
    save(filename3,'Non_dominated_sol');
end
%% the gene pairs of mutil-modal PDENB or PDENB of patient samples

for fnum=1:length(numf)-2 % all patient samples
    data=load([path2,expression_tumor_fileName(1:4),'_',num2str(fnum) '_PEN']);
    Sol=load([path3,expression_tumor_fileName(1:4),'_sample_',num2str(fnum),'_MMPDENB-RBM_PF.mat']);
    Sol=Sol.PF;
    POP=load([path3,expression_tumor_fileName(1:4),'_sample_',num2str(fnum),'_MMPDENB-RBM_PS.mat']);
    POP=POP.PS;
    maxf=max(Sol(:,2));
    index=Sol(:,2)==maxf;
    individuals=POP(index,:);
    
    MMPDNB=cell(size(individuals,1),1);
    for i=1:size(individuals,1)
        P=find(individuals(i,:)~=0);
        for j=1:length(P)
            node_gene=data.gg(P(j),:);
            MM{j,1}=data.subnetwork_genes{node_gene(1),1};
            MM{j,2}=data.subnetwork_genes{node_gene(2),1};
        end
        MMPDNB{i,1}=MM;
        clear MM
    end
    gen_name{fnum,1}=MMPDNB;
end
end