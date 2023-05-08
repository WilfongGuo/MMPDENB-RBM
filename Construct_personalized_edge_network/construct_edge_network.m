function [subnetwork_genes,edge_net,gg] = construct_edge_network(data,referenceT,referenceN,fnum)
%% 根据参考文献Individual-specific edge-network analysis for disease prediction
%   先选出关联性较强的边（pcc>average pcc）计算具有包含相同节点的边之间的phpcc，从而构建边网络
%   感觉是强强相关

node_net=data.subnetwork_adjacency;

%% Select candidate gene pairs
node_netabs=abs(node_net);
avpcc=mean(node_netabs(node_netabs~=0));
for i=1:length(node_netabs)%% pcc of edge >average pcc was remained
    node_netabs(i,find(node_netabs(i,:)<avpcc))=0;
    node_netabs(i,find(node_netabs(i,:)>avpcc))=1;
end

%% paired gene :gg
Cons=CONS(node_netabs);
cons=Cons;
Cnum=size(cons,1);
g=cell(Cnum,1);
for i=1:Cnum
    g{i}=find(cons(i,:)==1);
end
gg=cell2mat(g);

%% load reference sample
[BB_t,~,B_t]=xlsread(referenceT);
[BB_n,~,B_n]=xlsread(referenceN);

%% the construction of edge_net
edge_net=zeros(Cnum);
for i=1:size(edge_net,1)
    g1p=gg(i,1);
    g2p=gg(i,2);
    p=find(gg(:,1)==g1p); 
    p=setdiff(p,i);
    p1=[];
    for k=1:length(p)
        if edge_net(i,p(k,:))~=0
            p1=[p1;k];
        end
    end
    p(p1)=[];
    if ~isempty(p)
        for j=1:length(p)
            g3p=gg(p(j,:),1);
            g4p=gg(p(j,:),2);
            g1n=data.subnetwork_genes{g1p,1};
            g2n=data.subnetwork_genes{g2p,1};
            g3n=data.subnetwork_genes{g3p,1};
            g4n=data.subnetwork_genes{g4p,1};
            PCC_T=shPCC(g1n,g2n,g3n,g4n,B_t,BB_t,fnum);
            PCC_N=shPCC(g1n,g2n,g3n,g4n,B_n,BB_n,fnum);
            if isnan(log2(abs(PCC_T/PCC_N)))
                edge_net(i,p(j,:))=0;
            else
                edge_net(i,p(j,:))=log2(abs(PCC_T/PCC_N));
            end
            edge_net(p(j,:),i)=edge_net(i,p(j,:));
            
        end
    end
    p=find(gg(:,2)==g1p); 
    p=setdiff(p,i);
    p1=[];
    for k=1:length(p)
        if edge_net(i,p(k,:))~=0
            p1=[p1;k];
        end
    end
    p(p1)=[];
    if ~isempty(p)
        for j=1:length(p)
            g3p=gg(p(j,:),1);
            g4p=gg(p(j,:),2);
            g1n=data.subnetwork_genes{g1p,1};
            g2n=data.subnetwork_genes{g2p,1};
            g3n=data.subnetwork_genes{g3p,1};
            g4n=data.subnetwork_genes{g4p,1};
            PCC_T=shPCC(g1n,g2n,g3n,g4n,B_t,BB_t,fnum);
            PCC_N=shPCC(g1n,g2n,g3n,g4n,B_n,BB_n,fnum);
            if isnan(log2(abs(PCC_T/PCC_N)))
                edge_net(i,p(j,:))=0;
            else
                edge_net(i,p(j,:))=log2(abs(PCC_T/PCC_N));
            end
            edge_net(p(j,:),i)=edge_net(i,p(j,:));
            
        end
    end
    p=find(gg(:,1)==g2p); 
    p=setdiff(p,i);
    p1=[];
    for k=1:length(p)
        if edge_net(i,p(k,:))~=0
            p1=[p1;k];
        end
    end
    p(p1)=[];
    if ~isempty(p)
        for j=1:length(p)
            g3p=gg(p(j,:),1);
            g4p=gg(p(j,:),2);
            g1n=data.subnetwork_genes{g1p,1};
            g2n=data.subnetwork_genes{g2p,1};
            g3n=data.subnetwork_genes{g3p,1};
            g4n=data.subnetwork_genes{g4p,1};
            PCC_T=shPCC(g1n,g2n,g3n,g4n,B_t,BB_t,fnum);
            PCC_N=shPCC(g1n,g2n,g3n,g4n,B_n,BB_n,fnum);
            if isnan(log2(abs(PCC_T/PCC_N)))
                edge_net(i,p(j,:))=0;
            else
                edge_net(i,p(j,:))=log2(abs(PCC_T/PCC_N));
            end
            edge_net(p(j,:),i)=edge_net(i,p(j,:));
            
        end
    end
    p=find(gg(:,2)==g2p); 
    p=setdiff(p,i);
    p1=[];
    for k=1:length(p)
        if edge_net(i,p(k,:))~=0
            p1=[p1;k];
        end
    end
    p(p1)=[];
    if ~isempty(p)
        for j=1:length(p)
            g3p=gg(p(j,:),1);
            g4p=gg(p(j,:),2);
            g1n=data.subnetwork_genes{g1p,1};
            g2n=data.subnetwork_genes{g2p,1};
            g3n=data.subnetwork_genes{g3p,1};
            g4n=data.subnetwork_genes{g4p,1};
            PCC_T=shPCC(g1n,g2n,g3n,g4n,B_t,BB_t,fnum);
            PCC_N=shPCC(g1n,g2n,g3n,g4n,B_n,BB_n,fnum);
            if isnan(log2(abs(PCC_T/PCC_N)))
                edge_net(i,p(j,:))=0;
            else
                edge_net(i,p(j,:))=log2(abs(PCC_T/PCC_N));
            end
            edge_net(p(j,:),i)=edge_net(i,p(j,:));
            
        end
    end
end
subnetwork_genes=data.subnetwork_genes;

end



function A_adjacent=CONS(test_Net)
[z1,z2]=find(triu(test_Net)~=0);
z=[z1,z2];
NNN=length(test_Net);

N1=NNN;
[N2,~]=size(z);
%calculate the adjacency matrix of bipartite graph
A_adjacent=zeros(N2,N1);
for i=1:N2
    
    A_adjacent(i,z(i,1))=1;
    A_adjacent(i,z(i,2))=1;
    
end
end

function value=shPCC(g1n,g2n,g3n,g4n,B,BB,fnum)

index1=find(strcmp(g1n,B(2:end,1)))+1;

index2=find(strcmp(g2n,B(2:end,1)))+1;

index3=find(strcmp(g3n,B(2:end,1)))+1;

index4=find(strcmp(g4n,B(2:end,1)))+1;

C=(BB(index1-1,fnum)-mean(BB(index1-1,:)))*(BB(index2-1,fnum)-mean(BB(index2-1,:)))...
    *(BB(index3-1,fnum)-mean(BB(index3-1,:)))*(BB(index4-1,fnum)-mean(BB(index4-1,:)));
V=(var(BB(index1-1,:))*var(BB(index2-1,:))*var(BB(index3-1,:))*var(BB(index4-1,:))).^0.5;
value=C/V;
end