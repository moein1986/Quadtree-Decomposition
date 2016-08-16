
clear all
close all
clc

%% Mesh Generation

load('circle1','circle1');
load('circle2','circle2');
load('circle3','circle3');
load('circle4','circle4');

xmin=0; xmax=1; ymin=0; ymax=1; % computational domain

pts=[xmin ymin  ; xmin ymax  ; xmax ymin  ; xmax ymax ];

QuadtreeMesh = QuadtreeDecompFEM(pts,'maxDepth',7,'Gridsize',4);
Nn=max(QuadtreeMesh.ConnectivityArrayEdge(:));  %Number of edges

%% Extract Leaf Elements

FinalElements=find(QuadtreeMesh.BinIsLeaf==0);
QuadtreeMesh.ConnectivityArray(FinalElements,:)=[];
QuadtreeMesh.ConnectivityArrayEdge(FinalElements,:)=[];

QuadtreeMesh.BinEdgeSize(FinalElements,:)=[];
QuadtreeMesh.alpha(FinalElements,:)=[];
QuadtreeMesh.BinBoundaries(FinalElements,:)=[];

 %% Visualizing
 
figure
hold on
axis equal

fill (circle1(:,1),circle1(:,2),'r')
fill (circle2(:,1),circle2(:,2),'r')
fill (circle3(:,1),circle3(:,2),'r')
fill (circle4(:,1),circle4(:,2),'r')

%
pos=QuadtreeMesh.BinBoundaries;
for e=1:length(pos)
    line([pos(e,1) pos(e,3)],[pos(e,2) pos(e,2)],'Color',[0 0 0])
    line([pos(e,3) pos(e,3)],[pos(e,2) pos(e,4)],'Color',[0 0 0])
    line([pos(e,3) pos(e,1)],[pos(e,4) pos(e,4)],'Color',[0 0 0])
    line([pos(e,1) pos(e,1)],[pos(e,4) pos(e,2)],'Color',[0 0 0])


    %                 text(0.5*(pos(e,1)+pos(e,3)),0.5*(pos(e,2)+pos(e,4)),['',num2str(e)],'FontSize',6)
end
% 
line([pts(1,1) pts(3,1)],[pts(1,2) pts(3,2)],'LineWidth',3,'Color',[0 0 0])
line([pts(1,1) pts(2,1)],[pts(1,2) pts(2,2)],'LineWidth',3,'Color',[0 0 0])
line([pts(2,1) pts(4,1)],[pts(2,2) pts(4,2)],'LineWidth',3,'Color',[0 0 0])
line([pts(3,1) pts(3,1)],[pts(1,2) pts(2,2)],'LineWidth',3,'Color',[0 0 0])
 
