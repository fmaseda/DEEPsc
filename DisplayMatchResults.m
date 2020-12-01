function DisplayMatchResults(system,Corr,figNum,cellnum,titletext,known,Patches,Vertices)
% defaults
if nargin<5
    titletext='';
end
if nargin<6
    known=false;
end

figure(figNum); clf
switch lower(system)
    case {'follicle','mouse','hair'}
        axis([0 1100 -1100 0])
        is3d=false;
    case {'zebrafish','satija','seurat'}
        axis([-1.2 1.2 -.5 1.5])
        is3d=false;
    case {'drosophila','fly'}
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
    % Patches should hold (x,y,z) coords
    scatter3(Patches(:,1),Patches(:,2),Patches(:,3), ...
        20*ones(size(Patches,1),1),Corr(cellnum,:),'o','filled')
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
        scatter3(Patches(cellnum,1),Patches(cellnum,2),Patches(cellnum,3), ...
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