function [PF,PS,Non_dominated_sol]=multi_modal_EA_in_MMPDENB_RBM(data,popnum,Max_CalNum,experiment_num)
%%     input:
%                   data            :    patient sample's PEN
%                   popnum          :    population size
%                   Max_CalNum      :    the maximum number of function evaluation
%                   experiment_num  :    the number of algorithm rus
%       output:
%                   PF   :    the value of objective function of PDENB obtained by MMPDENB-RBM
%                   PS   :    non dominated PDENB obtained by MMPDENB-RBM
%                   NDS  :    the results of 30 runs of MMPDENB-RBM

test_adjacency=data.edge_net;
Problem.D=size(test_adjacency,1); %dimension
Problem.N=popnum;
Problem.maxFE=Max_CalNum;
Problem.lower=0;
Problem.upper=1;
test_net=zeros(Problem.D,Problem.D);
test_net(test_adjacency~=0)=1;
%% paired edge :g
Cons=CONS(test_net);
Cons(all(Cons==0,2),:)=[];
Cnum=size(Cons,1);
g=cell(Cnum,1);
for i=1:Cnum
    g{i}=find(Cons(i,:)==1);
end
clear Cons
%% 
Dimension=eye(Problem.D,Problem.D);
CV_num=Calcons(Problem.D,Cnum,g,Dimension);
D_score=abs(CV_num-Cnum);
Non_dominated_sol=cell(30,1);
R=cell(30,1);
for EXP_NUM=1:experiment_num
    CalNum=0;
    
    %% Population initialization
    REAL = false;
    
    Dec = ones(Problem.N,Problem.D);
    Mask=creatpop(Problem.N,Problem.D,D_score,g);
    [Mask,~] = unique(Mask,'rows');
    if size(Mask,1)<Problem.N
        newsoultion=Problem.N-size(Mask,1);
        ns=creatpop(newsoultion,Problem.D,D_score,g);
        Mask=[Mask;ns];
    end

    Population = Dec.*Mask;
    [Population,calnum]  = SOLUTION(Population,test_adjacency,g,D_score);
    CalNum=CalNum+calnum;
    
    [Population,Dec,Mask,FrontNo,CrowdDis] = EnvironmentalSelection(Population,Dec,Mask,Problem.N,0,0);
    
    %% Optimization
    rho = 0.5;
    RHO=0.5;
    while(CalNum<Problem.maxFE)
        if CalNum<ceil(Problem.maxFE/2) 
            Site=false(1,ceil(Problem.N));
            if any(Site)
                [rbm,~,allZero,allOne] = ModelTraining(Mask(index,:),Dec(index,:),REAL);
            else
                [rbm,~,allZero,allOne] = deal([]);
            end
            
            Popmarks   = extract_d(Population);
            MatingPool = TournamentSelection_hamming(Problem.N,Problem.N,Popmarks,FrontNo,CrowdDis);
            [OffMask,poss_s_num] = Operator1(MatingPool,rbm,Site,allZero,allOne,D_score,test_adjacency);
            OffMask =logical(OffMask);
            OffDec  =ones(size(OffMask,1),size(OffMask,2));
            Offspring = OffDec.*OffMask;
            [Offspring,calnum]  = SOLUTION(Offspring,test_adjacency,g,D_score);
            CalNum=CalNum+calnum;
            
            [Population,Dec,Mask,FrontNo,CrowdDis,sRatio] = EnvironmentalSelection([Population,Offspring],[Dec;OffDec],[Mask;OffMask],Problem.N,length(Population),poss_s_num);
            rho = (rho+sRatio)/2;
            RHO=[RHO;rho];
        else %% 后半段在充分了解pf的情况下，用RBM去学习pareto最优子空间
            if rho<0.5

                index=FrontNo<ceil(max(FrontNo)/2);
            else
                index=FrontNo==1;
            end
            Site = rho > rand(1,ceil(Problem.N));  
            
            if any(Site)
                [rbm,~,allZero,allOne] = ModelTraining(Mask(index,:),Dec(index,:),REAL);
                
            else
                [rbm,~,allZero,allOne] = deal([]);
            end
            
            Popmarks   = extract_d(Population);
            MatingPool = TournamentSelection_hamming(Problem.N,Problem.N,Popmarks,FrontNo,CrowdDis);
            [OffMask,poss_s_num] = Operator1(MatingPool,rbm,Site,allZero,allOne,D_score,test_adjacency);
            OffMask =logical(OffMask);
            OffDec  =ones(size(OffMask,1),size(OffMask,2));
            Offspring = OffDec.*OffMask;
            [Offspring,calnum]  = SOLUTION(Offspring,test_adjacency,g,D_score);
            CalNum=CalNum+calnum;
            
            [Population,Dec,Mask,FrontNo,CrowdDis,sRatio] = EnvironmentalSelection([Population,Offspring],[Dec;OffDec],[Mask;OffMask],Problem.N,length(Population),poss_s_num);
            rho = (rho+sRatio)/2;
            RHO=[RHO;rho];
        end
    end
    PopObj      = extract(Population);
    Pop         = extract_d(Population);
    [FrontNo,~] = NDSort(PopObj,size(PopObj,1));
    outputpop=Pop(FrontNo==1,:);
    Non_dominated_sol{EXP_NUM}=outputpop;
    R{EXP_NUM}=RHO;
end

Non_dominated_sol=cell2mat(Non_dominated_sol);

[pop,~,~]=unique(Non_dominated_sol,'rows');% De-duplication

functionvalue = Calfunctionvalue(pop,test_adjacency);

[FrontNo,~] = M_non_domination_scd_sort(pop,functionvalue);

POP=pop(FrontNo==1,:);

FV=functionvalue(FrontNo==1,:);

[PF,mod_position]=sortrows(FV);

PS=POP(mod_position,:);

PF(:,2)=-PF(:,2);


end



function pop=creatpop(popnum,D,D_score,g)

pop=zeros(popnum,D);
gg=cell2mat(g);

for i=1:popnum
    mustselectnum_position=randperm(size(gg,1),1);  
    mustselectnum=gg(mustselectnum_position,:);
    pop(i,mustselectnum)=1;
    D_score(mustselectnum)=0;
    canditateD=find(D_score~=0);
    
    min_D_score=D_score;
    gennum=round(0.9*rand*size(canditateD,1));
    for j=1:gennum
        canditateD=find(min_D_score~=0);   
        variables=randperm(size(canditateD,1),2);
        if D_score(canditateD(variables(1)))>D_score(canditateD(variables(2)))
            pop(i,canditateD(variables(1)))=1;
            
            min_D_score(canditateD(variables(1)))=0;
        else
            pop(i,canditateD(variables(2)))=1;
            
            min_D_score(canditateD(variables(2)))=0;
        end
    end
end
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
function CVvalue = Calcons(popnum,Cnum,g,pop)
cv=zeros(Cnum,1);
CVvalue=zeros(popnum,1);
for i=1:popnum
    ind=pop(i,:);
    for j=1:Cnum
        cv(j)=max(1-sum(ind(g{j})),0);
    end
    CVvalue(i)=sum(cv);
end
end
function [Population,calnum]  = SOLUTION(Mask,test_adjacency,g,D_score)
[functionvalue,calnum1] = Calfunctionvalue(Mask,test_adjacency);
[Mask,tt]=FIX2(Mask,functionvalue,g,size(test_adjacency,2),D_score,test_adjacency);
[functionvalue,calnum2]=Calfunctionvalue_afterfix(Mask,test_adjacency,tt,functionvalue);
calnum=calnum1+calnum2;
Population=cell(1,size(Mask,1));
for i=1:size(Population,2)
    Population{1,i}.dec=Mask(i,:);
    Population{1,i}.obj=functionvalue(i,:);
    Population{1,i}.con=0;
    Population{1,i}.add=[];
end
end
function [Functionvalue, calnum]= Calfunctionvalue(pop,test_adjacency)
Functionvalue=zeros(size(pop,1),2);
functionvalue=zeros(size(pop,1),4);
for i=1:size(pop,1)
    Functionvalue(i,1)=sum(pop(i,:));
    
    a=find(pop(i,:)==1);
    matrix=test_adjacency(a,:);
    inmatrix=matrix(:,a);
    genin=nonzeros(inmatrix);
    
    functionvalue(i,1)=abs(mean(genin));
    functionvalue(i,2)=std(genin);
    
    
    gen=nonzeros(matrix);
    genout=setdiff(gen,genin);
    
    functionvalue(i,3)=abs(mean(genout));
    functionvalue(i,4)=functionvalue(i,1)*functionvalue(i,2)/functionvalue(i,3);
    
end
calnum=i;
Functionvalue(:,2)= -round(functionvalue(:,4)*100)/100;
tt= isnan(Functionvalue(:,2));
Functionvalue(tt,1)=size(test_adjacency,2);
Functionvalue(tt,2)=0;
end
function [Functionvalue,calnum] = Calfunctionvalue_afterfix(pop,test_adjacency,tt,functionvalue)

for i=1:size(tt,1)
    functionvalue(tt(i),:)=Calfunctionvalue(pop(tt(i),:),test_adjacency);
end
calnum=i;
Functionvalue=functionvalue;
tt= isnan(Functionvalue(:,2));
Functionvalue(tt,2)=0;
ttt= Functionvalue(:,2)==0;
Functionvalue(ttt,1)=size(test_adjacency,2);
end
function Parents = TournamentSelection_hamming(K,N,Population,FrontNo,SpCrowdDis)
index=zeros(K,1);

hamming_dist=pdist2(Population,Population,'hamming');
[~,site2]=sort(hamming_dist,2);
K1=randperm(N,K);  %The first parent number
K2=site2(K1,2);    %The second parent number

for i=1:K
    if FrontNo(K1(i))<FrontNo(K2(i))
        index(i)=K1(i);
    elseif FrontNo(K1(i))==FrontNo(K2(i))
        if SpCrowdDis(K1(i))>=SpCrowdDis(K2(i))
            index(i)=K1(i);
        else
            index(i)=K2(i);
        end
    else
        index(i)=K2(i);
    end
    
end
Parents=Population(index,:);
end
