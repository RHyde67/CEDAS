function [ ClustersOut ] = CEDAS_Funct( varargin )
%CEDAS_FUNCT
% CEDAS: Clustering of Evolving Data-streams into Arbitrary Shapes
% © R Hyde 2016
% Released under the GNU GPLver3.0
% This version runs on Matlab versions prior to 2015b as it does not utilise Matlab's
% 'graph' functions.
% This is an implementation of the CEDAS clustering algorithm based on the paper:
% Fully Online Clustering of Evolving Data Streams into Arbitrarily Shaped Clusters
% Hyde R, Angelov P, MacKenzie, AR
% submitted to Information Sciences, July 2016
% The function applies the CEDAS clustering algorithm to a single data sample.
% For a data stream, call the function with each new data sample.
% Note this function has been written for clarity and ease of use to
% demonstrate the algorithm and has not been optimized for speed.
% For an example of use, see accompanying wrapper.
% Inputs:
%   Rad: micro-cluster radius
%   MinClusSize: minimum number of sample in a micro-cluster to not be noise
%   Decay: rate at which to decay micro-clusters in number of data samples
% Outputs:
% ClustersOut: micro-cluster information

persistent clusters % keep all cluster information
if ~isstruct(clusters) % create empty structure if none exists
    clusters=struct();
end
%% Read algorithm parameters
Rad=varargin{1}(1,1);
MinClusSize=varargin{1}(1,2);
Decay=varargin{1}(1,3);
%% Read new data sample
NewSample=varargin{2};

ClusterChanged=0; % set flag that no micro-clusters have changed

%% Main Function
if ~isempty(fieldnames(clusters)) % if micro-clusters exist...
    CEDAS_Update_Cluster(); % assign data to micro-cluster and adjust micro-cluster info if reqd
    CEDAS_Kill_Clusters(); % decay and remove old micro-clusters
    if ClusterChanged>size(clusters.Centre,1) % if a cluster was changed, but then removed due to decay, then no clusters have changed
      ClusterChanged=0;
    end
else % create first micro-cluster if none exists
    clusters.Centre(1,:)=NewSample;
    clusters.Count(1,:)=1;
    clusters.Macro(1,:)=1;
    clusters.Life(1,:)=1; % Give cluster life to first cluster
    clusters.Edge(1,:)={1}; % create first global cluster (global clusters are -ve numbers)
end

% if a micro-cluster has changed and has enough data to not be noise then
% update graph
if ClusterChanged~=0 && clusters.Count(ClusterChanged)>MinClusSize
    CEDAS_Update_Graph();
end

ClustersOut=clusters; % return the clusters

function CEDAS_Update_Cluster()
    NumClusts=size(clusters.Centre,1);
    sqDistToAll=sum((repmat(NewSample,NumClusts,1)-clusters.Centre).^2,2); % find square distance to all centres
    [MinDist,MinDistIdx]=min(sqDistToAll); % find minimum sq distance and index of nearest centre
    MinDist=sqrt(MinDist); % find distance
    if MinDist<Rad % if in cluster add to cluster
        ClusterChanged=MinDistIdx;
        clusters.Count(MinDistIdx)=clusters.Count(MinDistIdx)+1; % update Count of samples assigned to cluster
        clusters.Life(MinDistIdx,:)=1; % Renew cluster life to full value

        if MinDist>Rad*0.5 % if outside cluster core adjust cluster info
            clusters.Centre(MinDistIdx,:)=((clusters.Count(MinDistIdx,:)-1)*clusters.Centre(MinDistIdx,:)+NewSample)./clusters.Count(MinDistIdx,:); % update cluster centre to mean of samples
        end

    else % create new cluster
        NumClusts=NumClusts+1; % add new cluster
        clusters.Centre(NumClusts,:)=NewSample;
        clusters.Count(NumClusts,:)=1;
        clusters.Macro(NumClusts,:)=NumClusts;
        clusters.Life(NumClusts,:)=1; % give new cluster some life
        clusters.Edge{NumClusts,:}=max([clusters.Edge{:}])+1;
    end
end % End CEDAS Update Cluster

function CEDAS_Kill_Clusters()
if isfield(clusters,'Life')
    clusters.Life=clusters.Life-(1/Decay); % Decay each cluster
    Dead=find(clusters.Life<0);

    if size(Dead,1)>0;
        for idx1=1:size(Dead,1)
            
            %% remove all references to Dead
            if any([clusters.Edge{:}]==Dead(idx1))
            for idx2=1:size(clusters.Edge,1)
                if any(clusters.Edge{idx2}==Dead(idx1,1))
                    tmp=clusters.Edge{idx2};
                    tmp(tmp==Dead(1,1))=[];
                    clusters.Edge(idx2)={tmp};
                end
            end
            end
            %% Update graph for each Edge micro-cluster
            ClusterChanged=Dead;
            CEDAS_Update_Graph;
            %% delete Dead cluster
            for f=fieldnames(clusters)'
                clusters.(f{1})(Dead(idx1),:)=[];
            end
            %% decrement all cluster references higher than Dead
            for idx2=1:size(clusters.Edge,1)
                if any(clusters.Edge{idx2}>Dead(idx1))
                    tmp=clusters.Edge{idx2};
                    tmp(tmp>Dead(idx1))=tmp(tmp>Dead(idx1))-1;
                    clusters.Edge(idx2)={tmp};

                end
            end
        end
        %% Reassign macro-cluster Number   
        clusters.Macro(:)=0;
        GN=0;
        while any(clusters.Macro==0)
            GN=GN+1;
            tmp1=0;
            tmp2=find(clusters.Macro==0,1,'first');
            while ~isequal(tmp1,tmp2) | isempty(tmp2)
                tmp1=tmp2;
                tmp2=unique([clusters.Edge{tmp1,:}]);
                if size(tmp2,1)==1 & tmp2(1,1)<0
                    tmp2=tmp1;
                else
                tmp2=tmp2(tmp2>0);
                end
            end
            clusters.Macro(tmp2)=GN;
        end

    end
end
end % End CEDAS Kill Cluster

function CEDAS_Update_Graph()
    NumClusts=size(clusters.Centre,1);
    DistToAll=sqrt(sum((repmat(clusters.Centre(ClusterChanged,:),NumClusts,1)-clusters.Centre).^2,2)); % find distance to all Centres
    DistToAll(clusters.Count<MinClusSize)=99; % set all small clusters to far away so not joined
    Rads=Rad*1.5; % sum the radii of changed cluster and all other clusters
    [Intersect,~]=find(DistToAll<Rads); % find where cluster centres are closer than sum of radii
    clusters.Edge(ClusterChanged)={unique([clusters.Edge{ClusterChanged},Intersect'])};
    GN=min(clusters.Macro(clusters.Edge{ClusterChanged}));
    clusters.Macro(ClusterChanged)=GN;
    clusters.Macro(clusters.Edge{ClusterChanged})=GN;

end % End CEDAS Update Graph

end % End Main Function

