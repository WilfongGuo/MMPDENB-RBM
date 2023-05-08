function [filename,subnetwork_adjacency,subnetwork_genes] = delete_non_expression(path_in,path_out,cancer,reference)
%% Delete genes whose gene expression is 0 in BRCA/LUNG reference samples

%% Find genes that are not expressed in all reference samples
[BB,~,B]=xlsread(reference);
gene_expression=sum(BB,2);
no_expression=find(gene_expression==0);
no_expression_name=cell(length(no_expression),1);
for i=1:length(no_expression)
    no_expression_name{i,1}=B{no_expression(i)+1,1};
end
B(no_expression+1,:)=[];
filename=strcat(cancer,'_normal');
xlswrite(filename, B);

%% delete

Data=strcat(path_in,'*.mat');

numf = dir(Data);
for fnum=1:length(numf)
    
    data=load([path_in,cancer, '_',num2str(fnum),'_PGIN' '']);
    sample_gene_name=data.subnetwork_genes;
    finger = contains(sample_gene_name,no_expression_name);
    p=find(finger~=0);
    
    
    data.subnetwork_genes(p)=[];
    data.subnetwork_adjacency(p,:)=[];
    data.subnetwork_adjacency(:,p)=[];
    
    
    subnetwork_adjacency=data.subnetwork_adjacency;
    subnetwork_genes=data.subnetwork_genes;
    
    filename=strcat(path_out,cancer, '_',num2str(fnum),'_PGIN.mat');
    save(filename,'subnetwork_adjacency','subnetwork_genes');
end

end

