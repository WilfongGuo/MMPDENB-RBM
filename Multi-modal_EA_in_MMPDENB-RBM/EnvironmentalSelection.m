function [Population,Dec,Mask,FrontNo,SpCrowdDis,sRatio] = EnvironmentalSelection(Population,Dec,Mask,N,len,num)
% The environmental selection of MOEA/PSL

%------------------------------- Copyright --------------------------------
% Copyright (c) 2022 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

%% Delete duplicated solutions
success = false(1,length(Population));
PopObjs        = extract(Population);
PopDecs        = extract_d(Population);
[~,uni] = unique(PopDecs,'rows');
if length(uni) == 1
    [~,uni] = unique(PopObjs,'rows');
end
Population = Population(uni);
PopObjs    = PopObjs(uni,:);
PopDecs    = PopDecs(uni,:);
Dec        = Dec(uni,:);
Mask       = Mask(uni,:);

% N          = min(N,length(Population));

%% Non-dominated sorting
[FrontNo,SpCrowdDis] = M_non_domination_scd_sort(Mask,PopObjs);
fnum=0;                                                                 
while numel(FrontNo,FrontNo<=fnum+1)<=N                      
    if numel(FrontNo,FrontNo<=fnum+1)==numel(FrontNo,FrontNo<=fnum)
        fnum=fnum-1;
        break
    end
    fnum=fnum+1;
end
MaxFNo=fnum+1;

Next = FrontNo < MaxFNo;


%% Select the solutions in the last front based on their crowding distances
Last     = find(FrontNo==MaxFNo);
[~,Rank] = sort(SpCrowdDis(Last),'descend');
Next(Last(Rank(1:N-sum(Next)))) = true;

%% Calculate the ratio of successful offsprings
success(uni(Next)) = true;
s1     = sum(success(len+1:len+num));
s2     = sum(success(len+num+1:end));
sRatio = (s1+1e-6)./(s1+s2+1e-6);
sRatio = min(max(sRatio,0.1),0.9);

%% Population for next generation
Population = Population(Next);
FrontNo    = FrontNo(Next);
SpCrowdDis   = SpCrowdDis(Next);
Dec        = Dec(Next,:);
Mask       = Mask(Next,:);
end