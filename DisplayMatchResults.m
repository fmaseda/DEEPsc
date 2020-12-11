function DisplayMatchResults(system,Corr,cellnum,varargin)
%% Displays heatmap of results from mapping method in Corr on cell cellnum
%
% Inputs
%   system:         string
%                   -'follicle','mouse','hair': Follicle system
%                   -'zebrafish','satija','seurat': Zebrafish system
%                   -'drosophila','fly','distmap': Drosophila system
%   Corr:           CxP array of correspondence scores from RunMatchingMethods()
%   cellnum:        integer, which cell to display heatmap for (row of Corr)
%
% Required Name-Value pairs:
%   If system='follicle' or 'zebrafish', require Patches and Vertices
%   If system='drosophila', require Centers
%
% Optional inputs
%   figNum:         integer, which figure number to display heatmap in (default=new)
%   Title:          string, text to display in title (default=empty)
%   knownOrigin:    boolean, if origin of cell is known, outline that
%                   location in black (default=true)
%
% ----------
% Example usage
%   DisplayMatchResults('follicle',MyCorr,1)

%% parse input arguments
for k = 1:2:length(varargin)
    switch lower(varargin{k})
        case 'patches'
            Patches = varargin{k+1};
        case 'vertices'
            Vertices = varargin{k+1};
        case 'centers'
            Centers = varargin{k+1};
        case 'fignum'
            figNum = varargin{k+1};
        case 'title'
            titletext = varargin{k+1};
        case 'knownorigin'
            known = varargin{k+1};
    end
end

%% set defaults
if ~exist('titletext','var') || isempty(titletext)
    titletext='';
end
if ~exist('known','var') || isempty(known)
    known=true;
end

%% display heatmap
if ~exist('figNum','var') || isempty(figNum)
    figure()
else
    figure(figNum); clf
end
switch lower(system)
    case {'follicle','mouse','hair'}
        axis([0 1100 -1100 0])
        is3d=false;
    case {'zebrafish','satija','seurat'}
        axis([-1.2 1.2 -.5 1.5])
        is3d=false;
    case {'drosophila','fly','distmap'}
        is3d=true;
    otherwise
        error('First argument is not a valid system name.')
end
if ~is3d
    for i=1:size(Patches,1)
        hold on
        patch(Vertices(Patches{i},1),Vertices(Patches{i},2),1*Corr(cellnum,i))
    end
else
    % Drosophila is 3d, so use scatter3 instead of patches
    % Centers should hold (x,y,z) coords
    scatter3(Centers(:,1),Centers(:,2),Centers(:,3), ...
        20*ones(size(Centers,1),1),Corr(cellnum,:),'o','filled')
    view(0,0)
    axis square
    axis equal
end
map=[linspace(0,1,10)' linspace(0,1,10)' ones(10,1); ...
ones(10,1) linspace(1,0,10)' linspace(1,0,10)'];
colormap(map)
colorbar
title(titletext);
if known
    if is3d
        hold on
        scatter3(Centers(cellnum,1),Centers(cellnum,2),Centers(cellnum,3), ...
            100,Corr(cellnum,cellnum),'o','filled', ...
            'LineWidth',10,'MarkerEdgeColor','black')
    else
        patch(Vertices(Patches{cellnum},1),Vertices(Patches{cellnum},2), ...
            Corr(cellnum,cellnum),'LineWidth',3,'EdgeColor','black')        
    end
end
caxis([0 1])
xticks([]); yticks([]);
set(gca, 'visible', 'off')
set(gca, 'fontsize', 14)
set(findall(gca, 'type', 'text'), 'visible', 'on')
hold off
end