1%% Wrapper to demo CEDAS algorithm on Mackey Glass data streams.

clear
clear functions

%% Read M-G Data
DataIn=csvread('M-G_3D_x2.csv');
DataMG=DataIn(:,1:3);Mins=min(DataMG,[],1);
Maxs=max(DataMG,[],1);
for idx=1:size(DataMG,2)
    DataMG(:,idx) = ( (DataMG(:,idx)-Mins(idx) ) ) / (Maxs(idx) - Mins(idx));
end
ViewAngle=[-22,52];
PA=[min(DataMG(:,1)),max(DataMG(:,1)),min(DataMG(:,2)),max(DataMG(:,2)),min(DataMG(:,3)),max(DataMG(:,3))];
Clrs=distinguishable_colors(20);
PlotSpeed=100;

%% CEDAS parameters
Radius=0.06;
Fade=1000;
MinThreshold=5;
Outliers=[];
Clusters=[];

%% create base micro-cluster sphere data
[x,y,z]=sphere;
x=x*Radius;y=y*Radius;z=z*Radius;

%% Run data stream
tic
for idx1=1001:size(DataIn,1)
    if idx1/1000==floor(idx1/1000)
        sprintf('Sample %i',idx1)
        sprintf('Mean Time %4f',toc/idx1)
    end
    Data=DataMG(idx1,:);
    [Clusters]=CEDAS_Pre2015b([Radius,MinThreshold, Fade],Data);

    %% Plot data and micro-clusters
    if idx1/PlotSpeed==floor(idx1/PlotSpeed) & ~isempty(Clusters)
        figure(1)
        clf
        view(ViewAngle);
        scatter3(DataMG(idx1-Fade:idx1,1),DataMG(idx1-Fade:idx1,2),DataMG(idx1-Fade:idx1,3),1,'ob')
        axis([PA]);
        CP=Clusters.Centre(~any(isnan(Clusters.Centre),2),:);
        CC=Clusters.Macro(~any(isnan(Clusters.Centre),2));
        [~,~,CC]=unique(CC);
        CC=rem(CC,20)+1;
        for idx2=1:size(CP,1)
            Paint=Clrs(CC(idx2),:);
            x1=x+CP(idx2,1);
            y1=y+CP(idx2,2);
            z1=z+CP(idx2,3);
            surface(x1,y1,z1,'FaceColor', Paint,'EdgeColor','none','FaceAlpha',0.6)
            hold on
        end
        drawnow
    end
end
%% Display timings, not relevant if plots displayed
t2=toc;
sprintf('Total Time: %3f',t2)
sprintf('Time per Sample: %3f',t2/idx1)
