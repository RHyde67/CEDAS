function []=CEDAS_IS_Demo()
% Wrapper function to demonstrate CEDAS.
% Requires Matlab 2015b, or later, to run the CEDAS function.
% The CEDAS clustering algorithm is based on the paper:
% Fully Online Clustering of Evolving Data Streams into Arbitrarily Shaped Clusters
% R Hyde, P Angelov, A R MacKenzie
% submitted to Information Sciences, July 2016
% A random data stream is created flowing in a cross formation where two
% groups of data merge and separate in the centre
% The section 'CEDAS Parameters' allows for different algorithm settings
% and the comments describe the effects of different values for this test
% data.
% Set 'Display' flag to 1 to show the data, micro-clusters and graph
% Set 'DisplayRate' to increase speed, e.g. DisplayRate=50 displays only
% after every 50 data samples
% the number displayed on each graph edge is the sum of the number of data
% samples in each node. This is additional to the base algorithm but may be
% of interest.

clear all

%% Read data from a separate function
DataIn=CrossData;

%% Set up display
Display=1; % flag for plotting data, micro-clusters and graph
DisplayRate=50; % number of samples to analyse before displaying results, reduces dealy due to plotting
if Display==1
    NumCols=20; % set the number of available colours for display of Macro-clusters, >20 is hard to distinuish
    Clrs=distinguishable_colors(NumCols);
    figure(1)
    clf
    axis([0 1 0 1])
    hData=scatter(1,1,2,'go'); % create handles to scatter graph
    hData.XData=[]; % empty x data from scatter
    hData.YData=[]; % empty y data from scatter
    hC=[]; % handle for micro-cluster patches
    hT=[]; % handle for micro-cluster text
    Buffer=[]; % buffer is used to hold data for plotting only and is not algorithm related
end

%% Initialise
Outliers=[]; % no outliers exist
Graph=graph(); % initialise graph variable to pass to CEDAS
idx1=0; % data index
Clusters=[];

%% CEDAS parameters
% Radius: 0.01 mostly separate, 0.03 outliers microC appear
% 0.04 good results small microC, 0.1 microC matches data width
Radius=0.05; % CEDAS microC radius
Decay=500; % number of data samples to consider 'recent data'
MinThreshold=1; % Min microC threshold

%% Run data stream
while idx1<=size(DataIn,1) % count through data and loop if desired
    %% Read data sample and loop
    idx1=idx1+1; % increment data sample number
    if idx1==size(DataIn,1)-1 % set loop if desired
        idx1=1;
    end
    
    %% Run CEDAS
    Data=DataIn(idx1,:);
    [Clusters,Outliers, Graph]=CEDAS(Data,Radius,Clusters,1/Decay,MinThreshold,Outliers,Graph);
    
    %% Display data, micro-cluster and graph if required
    if Display==1 & idx1/DisplayRate==floor(idx1/DisplayRate)
        PlotGraph() % display
    elseif Display==1
        Buffer=[Buffer;Data]; % buffer data for displaying at intervals
    end
end

function []=PlotGraph()
hData.XData=[hData.XData,Buffer(:,1)']; % append new data
hData.XData(1:max(0,end-Decay))=[]; % delete expired data
hData.YData=[hData.YData,Buffer(:,2)'];
hData.YData(1:max(0,end-Decay))=[]; % delete expired data
axis([0 1 0 1])
if size(hData.XData,1)>Decay
    hData.XData(1)=[];
    hData.YData(1)=[];
end
 if ~isempty(Clusters)
    figure(1)
    hold on
    axis([0 1 0 1]);
    CP=Clusters.C(~any(isnan(Clusters.C),2),:);
    CC=conncomp(Graph);
    CC=rem(CC,20)+1; % limit number of colours
    NoP=50; % number of points on circle
    t = 2*pi/NoP*(1:NoP);
    delete(hC(size(CP,1)+1:size(hC,2))); % delete unused patches
    hC(size(CP,1)+1:size(hC,2))=[]; % delete patch references
    delete(hT(size(CP,1)+1:size(hT,2))); % delete unused patch text
    hT(size(CP,1)+1:size(hT,2))=[]; % delete patch text references
    for n=1:size(CP,1) % display used patches
        if size(hC,2)>=n
            delete(hC(n));
            delete(hT(n))
        end
       hC(n) = fill(Clusters.C(n,1)+Radius.*cos(t), Clusters.C(n,2)+Radius.*sin(t),'',...
            'FaceColor',Clrs(CC(n),:),'FaceAlpha',0.2,'EdgeColor','none');
    end
    for n=1:size(CP,1)
    hT(n)= text(Clusters.C(n,1),Clusters.C(n,2),num2str(n),'FontSize',15);
    end
    title('Data and Micro-Clusters');
    xlabel('X Data');
    ylabel('Y Data');
    hold off
    figure(2)
    if ~isempty(Graph.Edges)
        hG=plot(Graph,'EdgeLabel',Graph.Edges.Weight);
    else
        hG=plot(Graph);
    end
    hG.NodeColor=Clrs(CC,:);
    hG.Parent.XTickLabel=[];
    hG.Parent.YTickLabel=[];
    hG.LineWidth=2;
    hG.MarkerSize=6;
    title('Micro-Cluster Graph Connections')
    drawnow
 end
end % end function

end % end main function

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
end

function [D1]=CrossData()
% function to create data stream
D1=[];
for idx1=0.1:0.01:0.9
    for idx2=1:5
    D1=[D1;idx1+rand(2,1)/30,0.1+randn(2,1)/30;...
        idx1+rand(2,1)/30,0.9+randn(2,1)/30];
    end
end

for idx1=0.9:-0.01:0.1
    for idx2=1:5
    D1=[D1;idx1+rand(2,1)/30,1-idx1+randn(2,1)/30;...
        idx1+rand(2,1)/30,idx1+randn(2,1)/30];
    end
    
end
end % end function

function colors = distinguishable_colors(n_colors,bg,func)
% DISTINGUISHABLE_COLORS: pick colors that are maximally perceptually distinct
%
% When plotting a set of lines, you may want to distinguish them by color.
% By default, Matlab chooses a small set of colors and cycles among them,
% and so if you have more than a few lines there will be confusion about
% which line is which. To fix this problem, one would want to be able to
% pick a much larger set of distinct colors, where the number of colors
% equals or exceeds the number of lines you want to plot. Because our
% ability to distinguish among colors has limits, one should choose these
% colors to be "maximally perceptually distinguishable."
%
% This function generates a set of colors which are distinguishable
% by reference to the "Lab" color space, which more closely matches
% human color perception than RGB. Given an initial large list of possible
% colors, it iteratively chooses the entry in the list that is farthest (in
% Lab space) from all previously-chosen entries. While this "greedy"
% algorithm does not yield a global maximum, it is simple and efficient.
% Moreover, the sequence of colors is consistent no matter how many you
% request, which facilitates the users' ability to learn the color order
% and avoids major changes in the appearance of plots when adding or
% removing lines.
%
% Syntax:
%   colors = distinguishable_colors(n_colors)
% Specify the number of colors you want as a scalar, n_colors. This will
% generate an n_colors-by-3 matrix, each row representing an RGB
% color triple. If you don't precisely know how many you will need in
% advance, there is no harm (other than execution time) in specifying
% slightly more than you think you will need.
%
%   colors = distinguishable_colors(n_colors,bg)
% This syntax allows you to specify the background color, to make sure that
% your colors are also distinguishable from the background. Default value
% is white. bg may be specified as an RGB triple or as one of the standard
% "ColorSpec" strings. You can even specify multiple colors:
%     bg = {'w','k'}
% or
%     bg = [1 1 1; 0 0 0]
% will only produce colors that are distinguishable from both white and
% black.
%
%   colors = distinguishable_colors(n_colors,bg,rgb2labfunc)
% By default, distinguishable_colors uses the image processing toolbox's
% color conversion functions makecform and applycform. Alternatively, you
% can supply your own color conversion function.
%
% Example:
%   c = distinguishable_colors(25);
%   figure
%   image(reshape(c,[1 size(c)]))
%
% Example using the file exchange's 'colorspace':
%   func = @(x) colorspace('RGB->Lab',x);
%   c = distinguishable_colors(25,'w',func);

% Copyright 2010-2011 by Timothy E. Holy

  % Parse the inputs
  if (nargin < 2)
    bg = [1 1 1];  % default white background
  else
    if iscell(bg)
      % User specified a list of colors as a cell aray
      bgc = bg;
      for i = 1:length(bgc)
	bgc{i} = parsecolor(bgc{i});
      end
      bg = cat(1,bgc{:});
    else
      % User specified a numeric array of colors (n-by-3)
      bg = parsecolor(bg);
    end
  end
  
  % Generate a sizable number of RGB triples. This represents our space of
  % possible choices. By starting in RGB space, we ensure that all of the
  % colors can be generated by the monitor.
  n_grid = 30;  % number of grid divisions along each axis in RGB space
  x = linspace(0,1,n_grid);
  [R,G,B] = ndgrid(x,x,x);
  rgb = [R(:) G(:) B(:)];
  if (n_colors > size(rgb,1)/3)
    error('You can''t readily distinguish that many colors');
  end
  
  % Convert to Lab color space, which more closely represents human
  % perception
  if (nargin > 2)
    lab = func(rgb);
    bglab = func(bg);
  else
    C = makecform('srgb2lab');
    lab = applycform(rgb,C);
    bglab = applycform(bg,C);
  end

  % If the user specified multiple background colors, compute distances
  % from the candidate colors to the background colors
  mindist2 = inf(size(rgb,1),1);
  for i = 1:size(bglab,1)-1
    dX = bsxfun(@minus,lab,bglab(i,:)); % displacement all colors from bg
    dist2 = sum(dX.^2,2);  % square distance
    mindist2 = min(dist2,mindist2);  % dist2 to closest previously-chosen color
  end
  
  % Iteratively pick the color that maximizes the distance to the nearest
  % already-picked color
  colors = zeros(n_colors,3);
  lastlab = bglab(end,:);   % initialize by making the "previous" color equal to background
  for i = 1:n_colors
    dX = bsxfun(@minus,lab,lastlab); % displacement of last from all colors on list
    dist2 = sum(dX.^2,2);  % square distance
    mindist2 = min(dist2,mindist2);  % dist2 to closest previously-chosen color
    [~,index] = max(mindist2);  % find the entry farthest from all previously-chosen colors
    colors(i,:) = rgb(index,:);  % save for output
    lastlab = lab(index,:);  % prepare for next iteration
  end
  
  function c = parsecolor(s)
  if ischar(s)
    c = colorstr2rgb(s);
  elseif isnumeric(s) && size(s,2) == 3
    c = s;
  else
    error('MATLAB:InvalidColorSpec','Color specification cannot be parsed.');
  end
  end

function c = colorstr2rgb(c)
  % Convert a color string to an RGB value.
  % This is cribbed from Matlab's whitebg function.
  % Why don't they make this a stand-alone function?
  rgbspec = [1 0 0;0 1 0;0 0 1;1 1 1;0 1 1;1 0 1;1 1 0;0 0 0];
  cspec = 'rgbwcmyk';
  k = find(cspec==c(1));
  if isempty(k)
    error('MATLAB:InvalidColorString','Unknown color string.');
  end
  if k~=3 || length(c)==1,
    c = rgbspec(k,:);
  elseif length(c)>2,
    if strcmpi(c(1:3),'bla')
      c = [0 0 0];
    elseif strcmpi(c(1:3),'blu')
      c = [0 0 1];
    else
      error('MATLAB:UnknownColorString', 'Unknown color string.');
    end
  end
end

end





