function [OptMatrix, OptOutVec] = Optimized_Groups(NumFact,N,p,r,GroupMat,Diagnostic)
%
% [OptMatrix, OptOutVec] = Optimized_Groups(NumFact,N,p,r,GroupMat,Diagnostic)
%
% Optimization in the choice of trajectories for the Morris experiment
% clear all

% Inputs
% N:= [1,1]         Total number of trajectories
% p:= [1,1]         Number of levels
% r:= [1,1]         Final number of optimal trjectories
% NumFact:= [1,1]   Number of factors
% LB:= [NumFact,1]  Lower bound of the uniform distribution for each factor
% UB:= [NumFact,1]  Upper bound of the uniform distribution for each factor
% GroupMat:=[NumFact,NumGroups] Matrix describing the groups. Each column represents a group and its elements 
%                               are set to 1 in correspondence of the factors that belong to the fixed group. All
%                               the other elements are zero.
% Diagnostic:= [1,1]            Boolean 1=plot the histograms and compute the
%                               efficiency of the samplign or not, 0
%                               otherwise

if nargin<2 | isempty(N),
  N=100;
end
if nargin<3 | isempty(p),
  p = 4;        % Number of levels
end
if nargin<4 | isempty(r),
  r = 8;        % Number of replicas
end
%GroupMat = [1 0; 1 0; 1 0; 0 1; 0 1; 1 0];
if nargin<1,
  disp('[OutMatrix, OutFact] = Optimized_Groups(NumFact,N,p,r,GroupMat,Diagnostic)')
  return
end
if nargin<5,
  GroupMat = [];
end
if nargin<6|isempty(Diagnostic),
  Diagnostic = 0;
end

LB = zeros(NumFact,1);
UB = ones(NumFact,1);

%[OutMatrix, OutFact] = Sampling_Function(p, NumFact, N, UB, LB);        % Version without Groups
[OutMatrix, OutFact] = Sampling_Function_2(p, NumFact, N, UB, LB, GroupMat);   % Version with Groups

GroupNumber = size(GroupMat,2);
if GroupNumber ~ 0,
    sizeb = GroupNumber + 1;
else
    sizeb = NumFact + 1;
end    

Dist = zeros(N,N);
Diff_Traj = [1:1:N];
% Compute the distance between all pair of trajectories (sum of the distances between points)
% The distance matrix is a matrix N*N
% The distance is defined as the sum of the distances between all pairs of points
% if the two trajectories differ, 0 otherwise
for j =1:N
    for z = j+1:N
    
        MyDist = zeros(sizeb, sizeb);
        for i = 1:sizeb,
            for k = 1:sizeb,    
                MyDist(i,k) = (sum((OutMatrix((j-1)*(sizeb) + i,:) - OutMatrix((z-1)*(sizeb) + k,:)).^2))^0.5;
            end
        end
        
        if size(find(MyDist==0),1) == sizeb,
            % Same trajectory. If the number of zeros in Dist matrix is equal to 
            % (NumFact+1) then the trajectory is a replica. In fact (NumFact+1) is the maximum numebr of 
            % points that two trajectories can have in common
            Dist(j,z) = 0;     
            Dist(z,j) = 0;  
            
            % Memorise the replicated trajectory
            Diff_Traj(1,z) = 0; 
        else
            % Define the distance between two trajectories as 
            % the minimum distance among their points
            Dist(j,z) = sum(sum(MyDist));     
            Dist(z,j) = sum(sum(MyDist));
        end        
    end
end

New_OutMatrix = [];
New_OutFact = [];
% Eliminate replicated trajectories in the sampled matrix
for i = 1:N
    if Diff_Traj(1,i)~=0
        New_OutMatrix = [New_OutMatrix; OutMatrix((i-1)*(sizeb) + 1: (i-1)*(sizeb) + sizeb,:)]; 
        New_OutFact = [New_OutFact; OutFact((i-1)*(sizeb) + 1: (i-1)*(sizeb) + sizeb,:)];
    end
end

% Select in the distance matrix only the rows and columns of different trajectories
Dist_Diff = Dist(find(Diff_Traj),find(Diff_Traj));
New_N = size(find(Diff_Traj), 2);

% Select the optimal set of trajectories
Traj_Vec = zeros(New_N, r);
OptDist = zeros(New_N, r);
for m = 1:New_N 
    
    Traj_Vec(m, 1) = m;

    for z = 2:r
        Max_New_Dist_Diff = 0; 
    
        for j = 1:New_N   
        
            % Check that trajectory j is not already in
            Is_done = 0;
            for h = 1:z
                if j == Traj_Vec(m,h) 
                    Is_done = 1;
                end
            end
        
            if Is_done==0
                New_Dist_Diff = 0;    
            
                % Compute the distance 
                for k = 1:z-1           
                    New_Dist_Diff = New_Dist_Diff + (Dist_Diff(Traj_Vec(m, k),j))^2;  
                end
        
                % Check if the distance is greater than the old one
                if New_Dist_Diff^0.5 > Max_New_Dist_Diff
                    Max_New_Dist_Diff = New_Dist_Diff^0.5;
                    Pippo = j;
                end
            end
        end
    
        % Set the new trajectory
        Traj_Vec(m,z) = Pippo;
        OptDist(m,z) = Max_New_Dist_Diff;
    end
end

% Construct optimal matrix
SumOptDist = sum(OptDist,2);
% Find the maximum distance
Pluto = find(SumOptDist == max(SumOptDist));
Opt_Traj_Vec = Traj_Vec(Pluto(1,1),:);

OptMatrix = [];
OptOutVec = [];  %
for k =1:r
    OptMatrix = [OptMatrix; New_OutMatrix((sizeb)*(Opt_Traj_Vec(1,k)-1) + 1:(sizeb)*(Opt_Traj_Vec(1,k)-1) + sizeb,:)];
    OptOutVec = [OptOutVec; New_OutFact((sizeb)*(Opt_Traj_Vec(1,k)-1) + 1:(sizeb)*(Opt_Traj_Vec(1,k)-1)+ sizeb,:)];
end

if Diagnostic == 1,
    % Clean the trajectories from repetitions and plot the histograms
    HistPlot = zeros(2*r,NumFact);
    for i = 1:NumFact,
        for j = 1:r,
            kk = 1;
        
            % select the first value of the factor
            HistPlot((j-1)*2+kk,i) = OptMatrix((j-1)*(sizeb)+1,i);
            
            % search the second value 
            for ii = 2:sizeb
                if OptMatrix((j-1)*(sizeb)+ii ,i) ~= OptMatrix((j-1)*(sizeb)+1,i),
                    kk = 2;        
                    HistPlot((j-1)*2+kk,i) = OptMatrix((j-1)*(sizeb)+ii ,i);
                end
            end    
        end
    end

	figure('name', 'New Strategy')
	hold on;
	DimPlots = round(NumFact/2);
	for i = 1:NumFact
        subplot(DimPlots,2,i);
        hist(HistPlot(:,i),p);
	end
        
	
	% Plot the histogram for the original samplng strategy
	% Select the matrix
	OrigSample = OutMatrix(1:r*(sizeb),:);
	OriHistPlot = zeros(2*r,NumFact);
	for i = 1:NumFact,
        for j = 1:r,
            kk = 1;
            
            % select the first value of the factor
            OriHistPlot((j-1)*2+kk,i) = OrigSample((j-1)*(sizeb)+1,i);
                
            % search the second value 
            for ii = 2:sizeb
                if OrigSample((j-1)*(sizeb)+ii ,i) ~= OrigSample((j-1)*(sizeb)+1,i),
                    kk = 2;        
                    OriHistPlot((j-1)*2+kk,i) = OrigSample((j-1)*(sizeb)+ii ,i);
                end
            end    
        end
	end
	
	figure('name', 'Old Strategy')
	hold on;
	for i = 1:NumFact
        subplot(DimPlots,2,i);
        hist(OriHistPlot(:,i),p);
	end
	
	% Measure the quality of the sampling strategy
	Levels = [0:(1/(p-1)):1];
	for i = 1:NumFact,
        for j = 1:p,
            % For each facrot and each level count the number of times the factor is on the level
            % This for the new and original sampling
            NumSPoint(i,j) = size(find(abs(HistPlot(:,i)-repmat(Levels(1,j), size(HistPlot,1),1))<1e-5),1);
            NumSOrigPoint(i,j) = size(find(abs(OriHistPlot(:,i)-repmat(Levels(1,j), size(OriHistPlot,1),1))<1e-5),1);
        end
	end
	
	% The optimal sampling has values uniformly distributed across the levels
	OptSampl = 2*r/p;
	QualMeasure = 0;
	QualOriMeasure = 0;
	for i = 1: NumFact,
        for j = 1:p,
            QualMeasure = QualMeasure + abs(NumSPoint(i,j)-OptSampl);
            QualOriMeasure = QualOriMeasure + abs(NumSOrigPoint(i,j)-OptSampl);
        end
	end
	
	QualMeasure = 1 - QualMeasure/(OptSampl*p*NumFact)
	QualOriMeasure = 1 - QualOriMeasure/(OptSampl*p*NumFact)
end	
	
	
	
            
            
            
            
        