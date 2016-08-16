classdef QuadtreeDecompFEM< handle
    
    
    %   Created by Moein Nazari.
    %   1.0     - 2013-06 Initial release
    %
    % Quadtree decomposition algorithm for finite element method
    % application. This class will recursively decompose the computational domain
    % into square elements.
    
    
    % QuadtreeMesh = QuadtreeDecompFEM(pts,maxDepth,Gridsize);
    % input arguments are maxDepth and Gridsize. maxDepth should be set to the
    % maximum level of refinement and Gridsize indicates the initial regular
    % grid's level of refinement.
    
    properties
        Points;
        
        BinCount;
        PointCount=4;
        EdgeCount=4;
        BinBoundaries;
        BinDepths;
        BinGridLevel=zeros(2*10^6,1);
        BinEdgeSize=zeros(2*10^6,2);
        BinIsLeaf=ones(2*10^6,1);
        BinParents = zeros(2*10^6,1);
        HangingNodes=zeros(3,2*10^6); % [  master points   HN number  ]
        HangingEdge=zeros(2*10^6,2); % [slave edge  master edge]
        BinChildren=zeros(2*10^6,4);
        ConnectivityArray=zeros(2*10^6,4);
        ConnectivityArrayEdge=zeros(2*10^6,4);
        alpha=zeros(2*10^6,1);
        BinNeighbor=zeros(2*10^6,4); % Face neighbors [up Down left right]
        Properties;
    end
    
    methods
        
        function this = QuadtreeDecompFEM(pts,varargin)
            
            % Define Boundary Box
            numPts = size(pts,1);
            this.BinBoundaries = [min(pts,[],1) max(pts,[],1)]; %[xmin ymin  xmax ymax ]
            this.Points = pts;
            this.BinDepths = 0;
            this.BinGridLevel=0;
            this.BinParents(1) = 0;
            this.BinCount = 1;
            this.ConnectivityArray(1,:)=[1 2 3 4];
            this.ConnectivityArrayEdge(1,:)=[1 2 3 4];
            
            % Allow custom setting of Properties
            IP = inputParser;
            IP.addParamValue('maxDepth',inf);
            IP.addParamValue('GridSize',1);
            
            
            IP.parse(varargin{:});
            this.Properties = IP.Results;
            
            
            % Return on empty or trivial bins
            if numPts<4, Disp('Boundary box is not defined properly'), return; end
            
            % Start dividing!
            
            this.preallocateSpace;
            this.divideGrid(1);
            SB= find(this.BinIsLeaf(1:this.BinCount));
            for i=1:length(SB)
                UniformFlag(this,SB(i),9.6,9.6,9.6,9.6,1) % this function should be defined based on the input geometry
            end
            this.divide(SB);
            
            %%
            this.deallocateSpace;
            
        end
        
        
        
        function preallocateSpace(this)
            
            numBins = 2*10e6;
            
            this.BinDepths(numBins) = 0;
            this.BinParents(numBins) = 0;
            this.BinBoundaries(numBins,1) = 0;
        end
        function deallocateSpace(this)
            this.BinDepths(this.BinCount+1:end) = [];
            this.BinParents(this.BinCount+1:end) = [];
            this.BinChildren(this.BinCount+1:end,:) = [];
            this.BinBoundaries(this.BinCount+1:end,:) = [];
            this.BinNeighbor(this.BinCount+1:end,:) = [];
            this.BinIsLeaf(this.BinCount+1:end) = [];
            this.BinEdgeSize(this.BinCount+1:end,:) = [];
            this.ConnectivityArray(this.BinCount+1:end,:)=[];
            this.ConnectivityArrayEdge(this.BinCount+1:end,:)=[];
            
            this.HangingNodes(:,this.PointCount+1:end)=[];
            this.HangingNodes(:, find(sum(this.HangingNodes) == 0)) = [];
            
            this.HangingEdge(this.EdgeCount+1:end,:)=[];
            this.HangingEdge(find(sum(this.HangingEdge') == 0),:) = [];
            
            this.alpha(this.BinCount+1:end)=[];
            
        end
        
        function divideGrid(this, startingBins)
            
            % Loop over each bin we will consider for division
            for i = 1:length(startingBins)
                binNo = startingBins(i);
                
                % Prevent dividing beyond the maximum depth
                if this.BinDepths(binNo)+1 >= this.Properties.GridSize
                    
                    continue;
                end
                
                this.BinEdgeSize(binNo,1) = this.BinBoundaries(binNo,3)-this.BinBoundaries(binNo,1); %xmax-xmin  hx
                this.BinEdgeSize(binNo,2) = this.BinBoundaries(binNo,4)-this.BinBoundaries(binNo,2); %ymax-ymin  hy
                
                
                % There is one conditions under which we should divide
                % this bin. 1: Grid Level is not reached to the maximum
                % grid level
                
                if this.BinGridLevel(binNo) < this.Properties.GridSize
                    % leaf elements
                    if this.BinIsLeaf(binNo)
                        this.divideBin(binNo);
                    end
                    this.divideGrid(this.BinChildren(binNo,:));
                    continue;
                    
                    
                    
                    
                end
                
            end
        end
        
        function divide(this, startingBins)
            % Loop over each bin we will consider for division
            for i = 1:length(startingBins)
                binNo = startingBins(i);
                
                % checkMaterial
                UniformFlag(this,binNo,9.6,9.6,9.6,9.6,1);
                %
                this.BinEdgeSize(binNo,1) = this.BinBoundaries(binNo,3)-this.BinBoundaries(binNo,1); %xmax-xmin
                this.BinEdgeSize(binNo,2) = this.BinBoundaries(binNo,4)-this.BinBoundaries(binNo,2); %xmax-xmin
                
                
                % Prevent dividing beyond the maximum depth
                if this.BinDepths(binNo)+1 >= this.Properties.maxDepth
                    % leaf elements
                    
                    continue;
                end
                
                
                % There are two conditions under which we should divide
                % this bin. 1: It's bigger than maxSize. 2: It is
                % non-uniform
                
                if ~this.alpha(binNo)
                    if this.BinIsLeaf(binNo)
                        this.divideBin(binNo);
                    end
                    this.divide(this.BinChildren(binNo,:));
                    continue;
                end
                
                
            end
        end
        
        
        
        
        
        function divideBin(this,binNo)
            % Gather the new points (a bit more efficient to copy once)
            this.BinCount = this.BinCount + 4;
            
            
            this.BinIsLeaf(binNo)=0;
            
            oldMin = this.BinBoundaries(binNo,1:2);
            oldMax = this.BinBoundaries(binNo,3:4);
            
            newDiv = mean([oldMin; oldMax], 1);
            
            % Build the new boundaries of our 4 subdivisions
            minMidMax = [oldMin newDiv oldMax];
            newBounds = minMidMax([...
                1 2 3 4;
                3 2 5 4;
                1 4 3 6;
                3 4 5 6]);
            
            
            newBinInds = this.BinCount-3:this.BinCount;
            this.BinEdgeSize(newBinInds,1)=0.5*this.BinEdgeSize(binNo,1);
            this.BinEdgeSize(newBinInds,2)=0.5*this.BinEdgeSize(binNo,2);
            this.BinChildren(binNo,:)=newBinInds;
            this.BinBoundaries(newBinInds,:) = newBounds;
            this.BinDepths(newBinInds) = this.BinDepths(binNo)+1;
            this.BinGridLevel(newBinInds) = this.BinGridLevel(binNo)+1;
            this.BinParents(newBinInds) = binNo;
            
            UniformFlag(this,newBinInds(1),9.6,9.6,9.6,9.6,1);
            UniformFlag(this,newBinInds(2),9.6,9.6,9.6,9.6,1);
            UniformFlag(this,newBinInds(3),9.6,9.6,9.6,9.6,1);
            UniformFlag(this,newBinInds(4),9.6,9.6,9.6,9.6,1);
            
            
            %% update Face and Edge Neighbor Matrix
            
            
            update_neighbor(this,binNo);
            
            %% plot
            %                         pts=[0 0 0 ; 0 1 0 ; 1 0 0 ; 0 0 1 ];
            %             figure
            %                         boxH = this.plot;
            %                    cols = lines(this.BinCount);
            %                    doplot3 = @(p,varargin)plot3(p(:,1),p(:,2),p(:,3),varargin{:});
            %                    for i = 1:this.BinCount
            % %                         Min = this.BinBoundaries(i,1:3);
            % %                         Max = this.BinBoundaries(i,4:6);
            % %                        A=mean([Min; Max], 1);
            % %                        text(A(1),A(2),A(3),['',num2str(i)],'FontSize',6)
            %
            %                        set(boxH(i),'Color',cols(i,:),'LineWidth', 1+this.BinDepths(i))
            %                        set(boxH(i),'Color',cols(i,:),'LineWidth', 1)
            %
            %                        set(boxH(i),'Color',cols(i,:))
            %
            %                        doplot3(pts(this.PointBins==i,:),'.','Color',cols(i,:))
            %                    end
            %                   axis image, view(3)
            %
            % %
            
            %% Update Connectivity array
            % add CenterPoint
            this.PointCount=this.PointCount+1;
            this.ConnectivityArray(newBinInds(1),3)=this.PointCount;
            this.ConnectivityArray(newBinInds(2),4)=this.PointCount;
            this.ConnectivityArray(newBinInds(3),2)=this.PointCount;
            this.ConnectivityArray(newBinInds(4),1)=this.PointCount;
            
            
            % update corner nodes
            this.ConnectivityArray(newBinInds(1),1)=this.ConnectivityArray(this.BinParents(newBinInds(1)),1);
            this.ConnectivityArray(newBinInds(2),2)=this.ConnectivityArray(this.BinParents(newBinInds(2)),2);
            this.ConnectivityArray(newBinInds(3),4)=this.ConnectivityArray(this.BinParents(newBinInds(3)),4);
            this.ConnectivityArray(newBinInds(4),3)=this.ConnectivityArray(this.BinParents(newBinInds(4)),3);
            
            % update ConnectivityArray for edges inside the divided element
            
            this.ConnectivityArrayEdge(newBinInds(1),4)=this.EdgeCount+1;
            this.ConnectivityArrayEdge(newBinInds(2),3)=this.EdgeCount+1;
            
            this.ConnectivityArrayEdge(newBinInds(1),2)=this.EdgeCount+2;
            this.ConnectivityArrayEdge(newBinInds(3),1)=this.EdgeCount+2;
            
            this.ConnectivityArrayEdge(newBinInds(2),2)=this.EdgeCount+3;
            this.ConnectivityArrayEdge(newBinInds(4),1)=this.EdgeCount+3;
            
            this.ConnectivityArrayEdge(newBinInds(3),4)=this.EdgeCount+4;
            this.ConnectivityArrayEdge(newBinInds(4),3)=this.EdgeCount+4;
            
            this.EdgeCount=this.EdgeCount+4;
            
            
            % update center nodes
            
            updateConArrayCenterNodes(this,newBinInds,binNo);
            
            %             arrayfun(@(direction) updateConArrayCenterNodes(this,newBinInds,binNo,direction) , 1:18);
            
            
            
        end
        
        
        
        function updateConArrayCenterNodes(this,newBinInds,binNo)
            
            for direction=1:4
                switch direction
                    case 1 % up
                        if this.BinNeighbor(binNo,direction)
                            if this.BinIsLeaf(this.BinNeighbor(binNo,direction))
                                this.PointCount=this.PointCount+1;
                                this.ConnectivityArray(newBinInds(3),3)=this.PointCount;
                                this.ConnectivityArray(newBinInds(4),4)=this.PointCount;
                                
                                this.ConnectivityArrayEdge(newBinInds(3),2)=this.EdgeCount+1;
                                this.ConnectivityArrayEdge(newBinInds(4),2)=this.EdgeCount+2;
                                this.EdgeCount=this.EdgeCount+2;
                                
                                this.HangingNodes(:,this.PointCount)=[this.PointCount; this.ConnectivityArray(binNo,3) ; this.ConnectivityArray(binNo,4)  ];
                                this.HangingEdge(this.EdgeCount,:)=[ this.EdgeCount  this.ConnectivityArrayEdge(this.BinNeighbor(binNo,direction),1)];
                                this.HangingEdge(this.EdgeCount-1,:)=[this.EdgeCount-1 this.ConnectivityArrayEdge(this.BinNeighbor(binNo,direction),1)];
                                
                            else
                                this.ConnectivityArray(newBinInds(3),3)=this.ConnectivityArray(this.BinChildren(this.BinNeighbor(binNo,direction),1),2);
                                this.ConnectivityArray(newBinInds(4),4)=this.ConnectivityArray(this.BinChildren(this.BinNeighbor(binNo,direction),1),2);
                                this.HangingNodes(:,this.ConnectivityArray(this.BinChildren(this.BinNeighbor(binNo,direction),1),2))=0;
                                
                                this.ConnectivityArrayEdge(newBinInds(3),2)=this.ConnectivityArrayEdge(this.BinChildren(this.BinNeighbor(binNo,direction),1),1);
                                this.ConnectivityArrayEdge(newBinInds(4),2)=this.ConnectivityArrayEdge(this.BinChildren(this.BinNeighbor(binNo,direction),2),1);
                                
                                this.HangingEdge(this.ConnectivityArrayEdge(this.BinChildren(this.BinNeighbor(binNo,direction),1),1),:)=0;
                                this.HangingEdge(this.ConnectivityArrayEdge(this.BinChildren(this.BinNeighbor(binNo,direction),2),1),:)=0;
                                
                            end
                        else
                            this.PointCount=this.PointCount+1;
                            this.ConnectivityArray(newBinInds(3),3)=this.PointCount;
                            this.ConnectivityArray(newBinInds(4),4)=this.PointCount;
                            
                            this.ConnectivityArrayEdge(newBinInds(3),2)=this.ConnectivityArrayEdge(binNo,2);
                            this.ConnectivityArrayEdge(newBinInds(4),2)=this.EdgeCount+1;
                            this.EdgeCount=this.EdgeCount+1;
                            
                        end
                        
                    case 2 % down
                        if this.BinNeighbor(binNo,direction)
                            if this.BinIsLeaf(this.BinNeighbor(binNo,direction))
                                
                                this.PointCount=this.PointCount+1;
                                this.ConnectivityArray(newBinInds(1),2)=this.PointCount;
                                this.ConnectivityArray(newBinInds(2),1)=this.PointCount;
                                this.HangingNodes(:,this.PointCount)=[this.PointCount; this.ConnectivityArray(binNo,1) ; this.ConnectivityArray(binNo,2)  ];
                                
                                this.ConnectivityArrayEdge(newBinInds(1),1)=this.EdgeCount+1;
                                this.ConnectivityArrayEdge(newBinInds(2),1)=this.EdgeCount+2;
                                this.EdgeCount=this.EdgeCount+2;
                                this.HangingEdge(this.EdgeCount,:)=[this.EdgeCount this.ConnectivityArrayEdge(this.BinNeighbor(binNo,direction),2)];
                                this.HangingEdge(this.EdgeCount-1,:)=[this.EdgeCount-1 this.ConnectivityArrayEdge(this.BinNeighbor(binNo,direction),2)];
                                
                                
                            else
                                this.ConnectivityArray(newBinInds(1),2)=this.ConnectivityArray(this.BinChildren(this.BinNeighbor(binNo,direction),3),3);
                                this.ConnectivityArray(newBinInds(2),1)=this.ConnectivityArray(this.BinChildren(this.BinNeighbor(binNo,direction),3),3);
                                this.HangingNodes(:,this.ConnectivityArray(this.BinChildren(this.BinNeighbor(binNo,direction),3),3))=0;
                                
                                this.ConnectivityArrayEdge(newBinInds(1),1)=this.ConnectivityArrayEdge(this.BinChildren(this.BinNeighbor(binNo,direction),3),2);
                                this.ConnectivityArrayEdge(newBinInds(2),1)=this.ConnectivityArrayEdge(this.BinChildren(this.BinNeighbor(binNo,direction),4),2);
                                
                                this.HangingEdge(this.ConnectivityArrayEdge(this.BinChildren(this.BinNeighbor(binNo,direction),3),2),:)=0;
                                this.HangingEdge(this.ConnectivityArrayEdge(this.BinChildren(this.BinNeighbor(binNo,direction),4),2),:)=0;
                                
                            end
                        else
                            this.PointCount=this.PointCount+1;
                            this.ConnectivityArray(newBinInds(1),2)=this.PointCount;
                            this.ConnectivityArray(newBinInds(2),1)=this.PointCount;
                            
                            this.ConnectivityArrayEdge(newBinInds(1),1)=this.ConnectivityArrayEdge(binNo,1);
                            this.ConnectivityArrayEdge(newBinInds(2),1)=this.EdgeCount+1;
                            this.EdgeCount=this.EdgeCount+1;
                        end
                        
                        
                    case 3 % left
                        if this.BinNeighbor(binNo,direction)
                            if this.BinIsLeaf(this.BinNeighbor(binNo,direction))
                                this.PointCount=this.PointCount+1;
                                this.ConnectivityArray(newBinInds(1),4)=this.PointCount;
                                this.ConnectivityArray(newBinInds(3),1)=this.PointCount;
                                this.HangingNodes(:,this.PointCount)=[this.PointCount; this.ConnectivityArray(binNo,1) ; this.ConnectivityArray(binNo,4)  ];
                                
                                this.ConnectivityArrayEdge(newBinInds(1),3)=this.EdgeCount+1;
                                this.ConnectivityArrayEdge(newBinInds(3),3)=this.EdgeCount+2;
                                this.EdgeCount=this.EdgeCount+2;
                                this.HangingEdge(this.EdgeCount,:)=[this.EdgeCount this.ConnectivityArrayEdge(this.BinNeighbor(binNo,direction),4)];
                                this.HangingEdge(this.EdgeCount-1,:)=[this.EdgeCount-1 this.ConnectivityArrayEdge(this.BinNeighbor(binNo,direction),4)];
                                
                                
                            else
                                this.ConnectivityArray(newBinInds(1),4)=this.ConnectivityArray(this.BinChildren(this.BinNeighbor(binNo,direction),2),3);
                                this.ConnectivityArray(newBinInds(3),1)=this.ConnectivityArray(this.BinChildren(this.BinNeighbor(binNo,direction),2),3);
                                this.HangingNodes(:,this.ConnectivityArray(this.BinChildren(this.BinNeighbor(binNo,direction),2),3))=0;
                                
                                this.ConnectivityArrayEdge(newBinInds(1),3)=this.ConnectivityArrayEdge(this.BinChildren(this.BinNeighbor(binNo,direction),2),4);
                                this.ConnectivityArrayEdge(newBinInds(3),3)=this.ConnectivityArrayEdge(this.BinChildren(this.BinNeighbor(binNo,direction),4),4);
                                this.HangingEdge(this.ConnectivityArrayEdge(this.BinChildren(this.BinNeighbor(binNo,direction),2),4),:)=0;
                                this.HangingEdge(this.ConnectivityArrayEdge(this.BinChildren(this.BinNeighbor(binNo,direction),4),4),:)=0;
                                
                            end
                        else
                            this.PointCount=this.PointCount+1;
                            this.ConnectivityArray(newBinInds(1),4)=this.PointCount;
                            this.ConnectivityArray(newBinInds(3),1)=this.PointCount;
                            
                            this.ConnectivityArrayEdge(newBinInds(1),3)=this.ConnectivityArrayEdge(binNo,3);
                            this.ConnectivityArrayEdge(newBinInds(3),3)=this.EdgeCount+1;
                            this.EdgeCount=this.EdgeCount+1;
                        end
                        
                    case 4 % Right
                        if this.BinNeighbor(binNo,direction)
                            if this.BinIsLeaf(this.BinNeighbor(binNo,direction))
                                this.PointCount=this.PointCount+1;
                                this.ConnectivityArray(newBinInds(2),3)=this.PointCount;
                                this.ConnectivityArray(newBinInds(4),2)=this.PointCount;
                                this.HangingNodes(:,this.PointCount)=[this.PointCount; this.ConnectivityArray(binNo,2) ; this.ConnectivityArray(binNo,3)  ];
                                
                                this.ConnectivityArrayEdge(newBinInds(2),4)=this.EdgeCount+1;
                                this.ConnectivityArrayEdge(newBinInds(4),4)=this.EdgeCount+2;
                                this.EdgeCount=this.EdgeCount+2;
                                this.HangingEdge(this.EdgeCount,:)=[this.EdgeCount this.ConnectivityArrayEdge(this.BinNeighbor(binNo,direction),3)];
                                this.HangingEdge(this.EdgeCount-1,:)=[this.EdgeCount-1 this.ConnectivityArrayEdge(this.BinNeighbor(binNo,direction),3)];
                                
                                
                            else
                                this.ConnectivityArray(newBinInds(2),3)=this.ConnectivityArray(this.BinChildren(this.BinNeighbor(binNo,direction),1),4);
                                this.ConnectivityArray(newBinInds(4),2)=this.ConnectivityArray(this.BinChildren(this.BinNeighbor(binNo,direction),1),4);
                                this.HangingNodes(:,this.ConnectivityArray(this.BinChildren(this.BinNeighbor(binNo,direction),1),4))=0;
                                
                                this.ConnectivityArrayEdge(newBinInds(2),4)=this.ConnectivityArrayEdge(this.BinChildren(this.BinNeighbor(binNo,direction),1),3);
                                this.ConnectivityArrayEdge(newBinInds(4),4)=this.ConnectivityArrayEdge(this.BinChildren(this.BinNeighbor(binNo,direction),3),3);
                                this.HangingEdge(this.ConnectivityArrayEdge(this.BinChildren(this.BinNeighbor(binNo,direction),1),3),:)=0;
                                this.HangingEdge(this.ConnectivityArrayEdge(this.BinChildren(this.BinNeighbor(binNo,direction),3),3),:)=0;
                                
                            end
                        else
                            this.PointCount=this.PointCount+1;
                            this.ConnectivityArray(newBinInds(2),3)=this.PointCount;
                            this.ConnectivityArray(newBinInds(4),2)=this.PointCount;
                            
                            this.ConnectivityArrayEdge(newBinInds(2),4)=this.ConnectivityArrayEdge(binNo,4);
                            this.ConnectivityArrayEdge(newBinInds(4),4)=this.EdgeCount+1;
                            this.EdgeCount=this.EdgeCount+1;
                        end
                        
                end
            end
        end
        
        
        function update_neighbor(this,binNo)
            
            
            
            for direction=1:4 % [up down left right]
                if   (this.BinNeighbor(binNo,direction))
                    Level_diff=this.BinDepths(binNo)-this.BinDepths(this.BinNeighbor(binNo,direction));
                    
                    
                    if Level_diff<1
                        
                        if this.BinIsLeaf(this.BinNeighbor(binNo,direction))
                            updateNeighborMatrixMinus(this,binNo,direction);
                        else
                            updateNeighborMatrixMinusUpdate(this,binNo,direction);
                        end
                        
                        
                    else  % level diff >=2
                        
                        this.divideBin(this.BinNeighbor(binNo,direction))
                        
                        updateNeighborMatrixMinus(this,binNo,direction);
                        %
                        
                        
                        
                    end
                    
                    
                else % this.BinNeighbor(binNo,direction)=0
                    updateNeighborMatrixRegular(this,binNo,direction);
                    
                end
                
            end
        end
        
        
        function updateNeighborMatrixMinusUpdate(this,binNo,direction)
            if direction==1
                this.BinNeighbor(this.BinChildren(binNo,:),direction)=[this.BinChildren(binNo,3) this.BinChildren(binNo,4) this.BinChildren(this.BinNeighbor(binNo,1),1) this.BinChildren(this.BinNeighbor(binNo,1),2) ]';
                
                
                this.BinNeighbor(this.BinChildren(this.BinNeighbor(binNo,direction),1),2)=this.BinChildren(binNo,3);
                this.BinNeighbor(this.BinChildren(this.BinNeighbor(binNo,direction),2),2)=this.BinChildren(binNo,4);
                
                
            elseif direction==2
                this.BinNeighbor(this.BinChildren(binNo,:),direction)=[this.BinChildren(this.BinNeighbor(binNo,2),3) this.BinChildren(this.BinNeighbor(binNo,2),4) this.BinChildren(binNo,1) this.BinChildren(binNo,2)  ]';
                
                
                
                this.BinNeighbor(this.BinChildren(this.BinNeighbor(binNo,direction),3),1)=this.BinChildren(binNo,1);
                this.BinNeighbor(this.BinChildren(this.BinNeighbor(binNo,direction),4),1)=this.BinChildren(binNo,2);
                
                
            elseif direction==3
                this.BinNeighbor(this.BinChildren(binNo,:),direction)=[this.BinChildren(this.BinNeighbor(binNo,3),2)  this.BinChildren(binNo,1) this.BinChildren(this.BinNeighbor(binNo,3),4) this.BinChildren(binNo,3)  ]';
                
                
                this.BinNeighbor(this.BinChildren(this.BinNeighbor(binNo,direction),2),4)=this.BinChildren(binNo,1);
                this.BinNeighbor(this.BinChildren(this.BinNeighbor(binNo,direction),4),4)=this.BinChildren(binNo,3);
                
                
                
            elseif direction==4
                this.BinNeighbor(this.BinChildren(binNo,:),direction)=[this.BinChildren(binNo,2)  this.BinChildren(this.BinNeighbor(binNo,4),1) this.BinChildren(binNo,4) this.BinChildren(this.BinNeighbor(binNo,4),3) ]'  ;
                
                
                this.BinNeighbor(this.BinChildren(this.BinNeighbor(binNo,direction),1),3)=this.BinChildren(binNo,2);
                this.BinNeighbor(this.BinChildren(this.BinNeighbor(binNo,direction),3),3)=this.BinChildren(binNo,4);
                
                
                
            end
            
            
        end
        
        
        
        
        function updateNeighborMatrixMinus(this,binNo,direction)
            
            if direction==1
                this.BinNeighbor(this.BinChildren(binNo,:),direction)=[this.BinChildren(binNo,3) this.BinChildren(binNo,4) this.BinNeighbor(binNo,1) this.BinNeighbor(binNo,1)   ]';
                
                
            elseif direction==2
                this.BinNeighbor(this.BinChildren(binNo,:),direction)=[this.BinNeighbor(binNo,2) this.BinNeighbor(binNo,2) this.BinChildren(binNo,1) this.BinChildren(binNo,2) ]';
                
                
                
            elseif direction==3
                this.BinNeighbor(this.BinChildren(binNo,:),direction)=[this.BinNeighbor(binNo,3)  this.BinChildren(binNo,1) this.BinNeighbor(binNo,3) this.BinChildren(binNo,3)  ]';
                
                
                
            elseif direction==4
                this.BinNeighbor(this.BinChildren(binNo,:),direction)=[this.BinChildren(binNo,2)  this.BinNeighbor(binNo,4) this.BinChildren(binNo,4) this.BinNeighbor(binNo,4) ]'  ;
                
                
                
                
                
                
            end
        end
        
        
        
        
        
        
        function updateNeighborMatrixRegular(this,binNo,direction)
            
            if direction==1
                this.BinNeighbor(this.BinChildren(binNo,:),direction)=[this.BinChildren(binNo,3)  this.BinChildren(binNo,4) this.BinNeighbor(binNo,1) this.BinNeighbor(binNo,1)   ]';
                
            elseif direction==2
                this.BinNeighbor(this.BinChildren(binNo,:),direction)=[this.BinNeighbor(binNo,2) this.BinNeighbor(binNo,2)  this.BinChildren(binNo,1) this.BinChildren(binNo,2)  ]';
                
            elseif direction==3
                this.BinNeighbor(this.BinChildren(binNo,:),direction)=[this.BinNeighbor(binNo,3) this.BinChildren(binNo,1) this.BinNeighbor(binNo,3) this.BinChildren(binNo,3) ]';
                
            elseif direction==4
                this.BinNeighbor(this.BinChildren(binNo,:),direction)=[this.BinChildren(binNo,2) this.BinNeighbor(binNo,4) this.BinChildren(binNo,4) this.BinNeighbor(binNo,4)  ]'  ;
                
                
                
            end
        end
        %
        
        
        
        
        
        function h = plot(this,varargin)
            % OcTree.plot plots bin bounding boxes of an OcTree object
            %
            % H = OT.plot('name',value,...) allows you to specify any
            % properties of the bounding box lines that you would normally
            % supply to a plot(...,'name',value) command, and returns plot
            % object handles (one per bin) to H.
            hold on;
            h = zeros(this.BinCount,1);
            for i = 1:this.BinCount
                binMinMax = this.BinBoundaries(i,:);
                pts = cat(1, binMinMax([...
                    1 2 3; 4 2 3; 4 5 3; 1 5 3; 1 2 3;...
                    1 2 6; 4 2 6; 4 5 6; 1 5 6; 1 2 6; 1 2 3]),...
                    nan(1,3), binMinMax([4 2 3; 4 2 6]),...
                    nan(1,3), binMinMax([4 5 3; 4 5 6]),...
                    nan(1,3), binMinMax([1 5 3; 1 5 6]));
                h(i) = plot3(pts(:,1),pts(:,2),pts(:,3),varargin{:});
            end
        end
        function h = plot3(this,varargin)
            % OcTree.plot plots bin bounding boxes of an OcTree
            %
            % See also OcTree.plot
            h = this.plot(varargin{:});
        end
        
        
        function UniformFlag(this,binNo,alpha1,alpha2,alpha3,alpha4,alpha5)
            
            
            % a 4 circle testcase is implemented here, for other geometris
            % you have to modified this function.
            
            A=[this.BinBoundaries(binNo,1) this.BinBoundaries(binNo,2) ;
                this.BinBoundaries(binNo,1) this.BinBoundaries(binNo,4);
                this.BinBoundaries(binNo,3) this.BinBoundaries(binNo,2);
                this.BinBoundaries(binNo,3) this.BinBoundaries(binNo,4)];
            A=[A;A;A;A];
            
            
            xCenter=[0.8 0.8 0.8 0.8 0.48 0.48 0.48 0.48 0.3 0.3 0.3 0.3 0.7 0.7 0.7 0.7 ]';
            yCenter=[0.8 0.8 0.8 0.8 0.6 0.6 0.6 0.6 0.7 0.7 0.7 0.7 0.2 0.2 0.2 0.2 ]';
            
            dist=sqrt((A(:,1)-xCenter).^2+(A(:,2)-yCenter).^2);
            c1=[dist(1)<.1 ; dist(2)<.1 ; dist(3)<.1 ; dist(4)<.1];
            c2=[dist(5)<.1 ; dist(6)<.1 ; dist(7)<.1 ; dist(8)<.1];
            c3=[dist(9)<.1 ; dist(10)<.1 ; dist(11)<.1 ; dist(12)<.1];
            c4=[dist(13)<.1 ; dist(14)<.1 ; dist(15)<.1 ; dist(16)<.1];
            
            
            C1=sum(c1); C2=sum(c2) ; C3=sum(c3) ; C4=sum(c4);
            
            if C1==0 && C2==0 && C3==0 && C4==0
                this.alpha(binNo)=alpha5;
            elseif C1==4
                this.alpha(binNo)=alpha1;
            elseif C2==4
                this.alpha(binNo)=alpha2;
            elseif C3==4
                this.alpha(binNo)=alpha3;
            elseif C4==4
                this.alpha(binNo)=alpha4;
                
            else
                this.alpha(binNo)=0;
            end
            
        end
        
    end
    
    
end