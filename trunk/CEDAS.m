function [ C,O,G ] = CEDAS( varargin )
% CEDAS: Clustering of Evolving Data-streams into Arbitrary Shapes
% © R Hyde 2016
% Released under the GNU GPLver3.0
% Requires Matlab 2015b or later as it uses Matlab's 'graph' functions.
% This is an implementation of the CEDAS clustering algorithm based on the paper:
% Fully Online Clustering of Evolving Data Streams into Arbitrarily Shaped Clusters
% Hyde R, Angelov P, MacKenzie, AR
% submitted to Information Sciences, July 2016
% The function applies the CEDAS clustering algorithm to a single data sample.
% For a data stream, call the function with each new data sample.
% Note this function has been written for clarity and ease of use to
% demonstrate the algorithm and has not been optimized for speed.
% For an example of use, see accompanying wrapper.
% Inputs: as described
% Outputs:
% C: cluster information
% O: outliers data not yet clustered
% G: graph structure of micro-clusters

D=varargin{1}; % Data sample
R=varargin{2}; % Micro-cluster radius
C=varargin{3}; % Clusters
F=varargin{4}; % Rate of micro-cluster fade
M=varargin{5}; % Minimum cluster threshold
O=varargin{6}; % Outliers - unclustered data samples
G=varargin{7}; % Graph structure


if isempty(C) % if no microC exist, create first
    O=[O;D];
    [C,O,G]=StartCluster(C,O,R,M,G);
else % else assign and decay
    [C,O,G]=Assign(D,C,R,M,O,G);
    [C,G]=Kill(F,C,R,G);
end
end

function [C,O,G]=StartCluster(varargin)
C=varargin{1}; % Clusters
O=varargin{2}; % Outliers
R=varargin{3}; % Radius
M=varargin{4}; % Minimum threshold
G=varargin{5}; % Graph structure
[d]=pdist2(O,O); % distance between all outliers
[Nin,NC]=max(sum(d<R)); % sum of number of samples within R for all data
    if sum(d(:,NC)<R)>M
    N=1; % new microC number
    [In,~]=find(d(:,NC)<R); % find data within R
    C=struct('C',mean(O(In,:),1)); % make centre mean of clustered data
    C.L(N,:)=1; % Set life
    C.T(N,:)=Nin; % Count number assigned
    C.K(N,:)=sum(d(:,NC)<0.5*R); % Count number in kernel
    G=addnode(G,1); % Add graph structure node
    O(In,:)=[]; % Remove outliers that have been assigned to new microC
    end

end % end function

function [C,O,G]=Assign(varargin)
D=varargin{1}; % Data
C=varargin{2}; % Clusters
R=varargin{3}; % Radius
M=varargin{4}; % Minimum threshold
O=varargin{5}; % Outliers
G=varargin{6}; % Graph structure
[d,I]=pdist2(C.C,D,'euclidean','Smallest',1);
if d<R % if within a microC
    C.L(I,:)=1; % reset life
    C.T(I,:)=C.T(I,:)+1; % increment total membership
    if d<R/2 % if in kernel
        C.K(I,:)=C.K(I,:)+1; % increment kernel membership
        C.C(I,:)=((C.C(I,:)*(C.K(I,:)-1))+D) / C.K(I,:); % recursiveley update centre
        G=Graph(C,G,I,R); % Update Graph
    end
else % Create new microC
    O=[O;D]; % add to outliers (unclustered)
    d=pdist2(O,O); % distance between all outliers
    [~,NC]=max(sum(d<R)); % sum of number of samples within R for all data
    if sum(d(:,NC)<R)>M % if enough data lie within radius to create a new microC
        N=size(C.C,1)+1; % new microC number
        [In,~]=find(d(:,NC)<R); % find data within R
        C.C(N,:)=mean(O(In,:),1); % make centre mean of clustered data
        C.L(N,:)=1; % set life
        C.T(N,:)=size(In,1); % count number assigned
        C.K(N,:)=sum(d(:,NC)<0.5*R); % count number in kernel
        O(In,:)=[]; % Remove outliers that have been assigned to new microC
        G=addnode(G,1); % Add node to graph structure
        G=Graph(C,G,N,R); % Update Graph
    end
    
end
end % End Assign function  

function [C,G]=Kill(varargin)
F=varargin{1}; % Fade (Decay) value
C=varargin{2}; % Clusters
R=varargin{3}; % Radius
G=varargin{4}; % Graph structure
C.L=C.L-F; % Decay all microC
I=find(C.L<0); % Find dead microC
if ~isempty(I)
    C.C(I,:)=[]; % delete microC information from clusters
    C.L(I,:)=[];
    C.T(I,:)=[];
    C.K(I,:)=[];
    G=rmnode(G,I);% delete references to dead microC
end
end % end Kill function   

function [G]=Graph(varargin)
C=varargin{1}; % Clusters
G=varargin{2}; % Graph structure
I=varargin{3}; % List of modified microC
R=varargin{4}; % Radius
d=pdist2(C.C(I,:),C.C); % distances between modified microC and all others
E=find(d<1.5*R); % list of microC that are edges
E(E==I)=[]; % remove self reference
if ~isempty(E) % if any Edges
    E=sort([ones(1,size(E,2))*I;E]',2); % find all edges
    E=setdiff(E,G.Edges.EndNodes,'rows'); % list of edges not present
    if any(E) % add missing edges
        W=C.T(E(:,1))+C.T(E(:,2))-1; % Edge weights, -1 as one will be added
        % in next step, weight is total number of data samples in nodes
        % (not required, but interesting information)
        G=addedge(G,E(:,1),E(:,2),W); % add missing edge
    end
    G.Edges.Weight(findedge(G,I,E(E~=I)))=G.Edges.Weight(findedge(G,I,E(E~=I)))+1; % increment edge weights
end
end % End graph function