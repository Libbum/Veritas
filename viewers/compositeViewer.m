function compositeViewer(map, tidx, particleType, rundir, logyn, rectsyn)
%Some good default values:
%map = flipud(hot(256)); tidx = 16; particleType = 0; rundir = '../bin/'; logyn = true; rectsyn = false;

%-------------
fid = fopen([rundir 'run.log'],'r');
cfg = fscanf(fid, '%c');
fclose(fid);

res = regexp(cfg,'nx:\s+([0-9\.]+)','tokens');
nx = str2double(res{:});

res = regexp(cfg,'dx:\s+([0-9\.eE-\+]+)','tokens');
dx = str2double(res{:});

res = regexp(cfg,'ratio:\s+([0-9\.]+)','tokens');
r = str2double(res{:});

res = regexp(cfg,'np:\s+([0-9,]+)','tokens');
np = str2double(horzcat(res{:}));
numParticles = numel(np);

res = regexp(cfg,'dp:\s+([0-9eE-\+\.]+)','tokens');
dp = str2double(horzcat(res{:}));
clear cfg
%-------------
rundir = [rundir 'output/rectangleData/'];
files = dir([rundir 'rectangleData_*.txt']);
nfiles = length(files);
times = cell(nfiles,1);
for ii = 1:nfiles
    times(ii) = regexp(files(ii).name,'(?<=_t)[0-9\.eE-\+]+(?=.txt)','match'); %Note: Can't cast this to double here, must keep the string value for file reads. Cast later if required.
end
times = unique(times);
[~,sidx] = sort(str2num(char(times)));
times = times(sidx);

if (tidx > length(times))
    error('tidx not in data');
end

levels = length(dir([rundir 'rectangleData_p' num2str(particleType) '_l*_t' char(times(tidx)) '.txt']));
maxDepth = levels-1;
rmap = lines(levels);
% force level colours
% clrs = lines(levels);
% rmap = zeros(size(clrs));
% rmap(5,:) = clrs(5,:);
% rmap(4,:) = clrs(1,:);
% rmap(3,:) = clrs(4,:);
% rmap(2,:) = clrs(2,:);
% rmap(1,:) = clrs(3,:);

rh = cell(1,levels);
z = cell(1,levels);

pos = get(0,'MonitorPositions');
%pos = pos(end,:); %Gets value of right monitor if 2 or more monitors, or falls back to value of single monitor.
hFig = figure('Visible','off', 'Name','Level Viewer', 'Resize','off', 'Position',[pos(1)+100 100 1000 850]);
movegui(hFig,'center')
ha = axes;
colormap(map);
hold on

ypos = linspace(0.5,np(particleType+1)-0.5,levels);

for ll=0:maxDepth
    xmax = nx.*r.^(maxDepth-ll);
    pmax = np(particleType+1).*r.^(maxDepth-ll);
    
    z{ll+1} = nan(pmax,xmax);
    if exist([rundir 'rectangleData_p' num2str(particleType) '_l' num2str(ll) '_t' char(times(tidx)) '.txt'],'file') == 2
        rects = [];
        ii = 1;
        fid = fopen([rundir 'rectangleData_p' num2str(particleType) '_l' num2str(ll) '_t' char(times(tidx)) '.txt'],'r');
        while ~feof(fid)
            rects(ii,:) = double(cell2mat(textscan(fid, '%*s %u %u %u %u',1)));
            data = textscan(fid, '%u %u %f');
            
            idx = (data{2} + data{1}*pmax)+1;
            for kk=1:length(idx)
                z{ll+1}(idx(kk)) = data{3}(kk);
            end
            ii = ii+1;
        end
        fclose(fid);
        if logyn
            hi{ll+1} = pcolor((0:xmax-1)./r.^(maxDepth-ll),(0:pmax-1)./r.^(maxDepth-ll),real(log10(z{ll+1})));
        else
            hi{ll+1} = pcolor((0:xmax-1)./r.^(maxDepth-ll),(0:pmax-1)./r.^(maxDepth-ll),z{ll+1});
        end
        shading flat
        
        hold on
        if rectsyn
            rh{ll+1} = zeros(size(rects,1),1);
            for ii=1:size(rects,1)
               rh{ll+1}(ii) = rectangle('position',rects(ii,:)./r.^(maxDepth-ll),'EdgeColor',rmap(ll+1,:),'LineWidth',2);
            end
        end
    end
    if rectsyn
        text(-nx/10,ypos(maxDepth-ll+1), ['Level ' num2str(maxDepth-ll)], 'color',rmap(ll+1,:));
    end
end


ch = colorbar;
if logyn
    ch.Label.String = 'log(data)';
    maxv = real(log10(nanmax(cellfun(@nanmax, cellfun(@nanmax, z, 'UniformOutput', 0)))));
    if ~isnan(maxv) && ~isinf(maxv)
        caxis([maxv-10 maxv]);
    end
else
    ch.Label.String = 'data';
    maxv = nanmax(cellfun(@nanmax, cellfun(@nanmax, z, 'UniformOutput', 0)));
end

if exist('hi','var')
    for ll=1:length(hi)
        uistack(hi{ll},'bottom');
    end
end
if exist('rl','var')
    for ll=fliplr(1:length(rh))
        uistack(rh{ll},'top');
    end
end
box on
set(ha,'layer','top')
title(['t = ' char(times(tidx))]);
axis tight
set(hFig, 'Visible','on')
end

