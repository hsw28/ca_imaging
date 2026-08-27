function OUT = same_location_tebc_feasibility(ratNames, varargin)
% SAME_LOCATION_TEBC_FEASIBILITY  Reviewer 3 same-location control diagnostic.
%
% This is a FEASIBILITY analysis only. It does not compute neural activity,
% alter rat structures, or modify any manuscript analysis.
%
% The analysis uses calcium/imaging frame timestamps. Position and speed are
% interpolated to those timestamps. Non-tEBC samples are always restricted to
% speed >= 4 cm/s. Task samples are evaluated both with and without that gate.
%
% Examples
%   OUT = same_location_tebc_feasibility;
%   OUT = same_location_tebc_feasibility({'rat0314','rat0816'});
%   OUT = same_location_tebc_feasibility([], 'Days','all', ...
%       'OutputDir','same_location_feasibility_all_days');
%
% Required fields in each base-workspace rat structure:
%   rat.pos.pos_DATE        [time x y]
%   rat.Ca_ts.CA_ts_DATE    imaging timestamps
%   rat.CS_times.CS_DATE    CS onset times
% Optional:
%   rat.CSUS_id.CSUS_id_DATE [labels; timestamps] or [labels; timestamps]
%   The optional CSUS_id is used to enforce the existing intertrial mask.
%
% Outputs are written to OutputDir:
%   same_location_feasibility_detailed.csv
%   same_location_feasibility_summary.csv
%   same_location_feasibility_matched_bins.csv
%   same_location_feasibility.mat
%   diagnostic PNG figures

if nargin < 1 || isempty(ratNames)
    ratNames = {'rat0222','rat0307','rat0313','rat0314','rat0816'};
end
ratNames = cellstr(ratNames);

p = inputParser;
addParameter(p,'Days','last3_to_An',@(x)ischar(x)||isstring(x)||iscell(x));
addParameter(p,'OutputDir',fullfile(pwd,'same_location_feasibility'),@(x)ischar(x)||isstring(x));
addParameter(p,'SpeedThreshold',4,@(x)isnumeric(x)&&isscalar(x)&&x>=0);
addParameter(p,'FineBinCm',2.5,@(x)isnumeric(x)&&isscalar(x)&&x>0);
addParameter(p,'SpeedBinSize',2,@(x)isempty(x)||(isnumeric(x)&&isscalar(x)&&x>0));
addParameter(p,'RepresentativeSessions',2,@(x)isnumeric(x)&&isscalar(x)&&x>=0);
addParameter(p,'MakePlots',true,@(x)islogical(x)&&isscalar(x));
parse(p,varargin{:});
opt = p.Results;
opt.OutputDir = char(opt.OutputDir);
if ~exist(opt.OutputDir,'dir'), mkdir(opt.OutputDir); end

thresholds = [0.133 0.5 1 2];
epochs = struct('name',{'full_task','trace_only'},'win',{[0 2],[0.25 0.75]});
speedConds = {'all_tEBC','running_tEBC'};
schemes = {'fine_2p5cm','coarse_2x2'};

rows = struct([]);
binRows = struct([]);
heatData = struct([]);
qualityRows = struct([]);

fprintf('\n=== Reviewer 3 same-location non-tEBC FEASIBILITY ===\n');
fprintf('Imaging-frame basis; non-tEBC always speed >= %.1f cm/s.\n',opt.SpeedThreshold);

for ri = 1:numel(ratNames)
    ratName = ratNames{ri};
    if ~evalin('base',sprintf('exist(''%s'',''var'')',ratName))
        warning('%s is not loaded in the base workspace; skipping.',ratName);
        continue
    end
    rat = evalin('base',ratName);
    days = select_days(rat,opt.Days);
    for di = 1:numel(days)
        day = days{di};
        try
            D = load_session(rat,day,opt);
        catch ME
            warning('[%s %s] %s',ratName,day,ME.message);
            continue
        end
        if numel(D.t) < 2 || ~any(D.valid)
            warning('[%s %s] No valid aligned imaging frames.',ratName,day);
            continue
        end

        q = struct('animal',ratName,'session_date',day, ...
            'n_imaging_frames',numel(D.t),'n_valid_xy_speed_frames',nnz(D.valid), ...
            'median_frame_duration_s',median(D.dt,'omitnan'), ...
            'max_frame_gap_s',max(D.rawGap,[],'omitnan'), ...
            'n_frame_gaps_capped',D.nGapCapped, ...
            'intertrial_mask_source',D.intertrialSource, ...
            'fine_edge_source',D.fineEdgeSource);
        qualityRows = append_struct(qualityRows,q);

        for ei = 1:numel(epochs)
            taskAll = in_windows(D.t,D.cs,epochs(ei).win) & D.valid;
            taskRunning = taskAll & D.speed >= opt.SpeedThreshold;
            nonRunning = D.intertrial & D.valid & D.speed >= opt.SpeedThreshold;

            for si = 1:numel(schemes)
                if si == 1, edges = D.fineEdges; else, edges = D.coarseEdges; end
                [bin,nY,nX] = position_bins(D.x,D.y,edges);
                for ci = 1:numel(speedConds)
                    if ci == 1, taskMask = taskAll; else, taskMask = taskRunning; end
                    [r,b] = diagnose_combination(ratName,day,epochs(ei).name, ...
                        schemes{si},speedConds{ci},taskMask,nonRunning,bin, ...
                        D.dt,D.speed,nY*nX,thresholds,opt);
                    rows = append_struct(rows,r);
                    binRows = append_struct_array(binRows,b);
                end

                % Heatmaps use all-tEBC frames. Store full-task and trace-only;
                % representative plotting later chooses sessions by valid frames.
                H = struct('animal',ratName,'session_date',day, ...
                    'task_epoch',epochs(ei).name,'binning_scheme',schemes{si}, ...
                    'nY',nY,'nX',nX,'edges',edges, ...
                    'taskOcc',occupancy_by_bin(bin,taskAll,D.dt,nY*nX), ...
                    'nonOcc',occupancy_by_bin(bin,nonRunning,D.dt,nY*nX), ...
                    'nValid',nnz(D.valid));
                heatData = append_struct(heatData,H);
            end
        end
    end
end

if isempty(rows)
    error(['No sessions were analyzed. Load the rat structures into the MATLAB ' ...
        'base workspace, then rerun this function.']);
end

detailed = struct2table(rows);
if isempty(binRows)
    matchedBins = table(cell(0,1),cell(0,1),cell(0,1),cell(0,1),cell(0,1), ...
        zeros(0,1),zeros(0,1),zeros(0,1),zeros(0,1), ...
        'VariableNames',{'animal','session_date','task_epoch','binning_scheme', ...
        'tEBC_speed_condition','spatial_bin','tEBC_occupancy_s', ...
        'non_tEBC_running_occupancy_s','occ_ratio'});
else
    matchedBins = struct2table(binRows);
end
summary = summarize_table(detailed);
quality = struct2table(qualityRows);

writetable(detailed,fullfile(opt.OutputDir,'same_location_feasibility_detailed.csv'));
writetable(summary,fullfile(opt.OutputDir,'same_location_feasibility_summary.csv'));
writetable(matchedBins,fullfile(opt.OutputDir,'same_location_feasibility_matched_bins.csv'));
writetable(quality,fullfile(opt.OutputDir,'same_location_feasibility_alignment_qc.csv'));

if opt.MakePlots
    make_bar_plot(detailed,'pct_tEBC_bins_with_any_non_tEBC_running', ...
        'tEBC bins with any running non-tEBC occupancy (%)', ...
        fullfile(opt.OutputDir,'pct_bins_with_any_overlap.png'));
    make_bar_plot(detailed,'pct_tEBC_frames_retained_ge_1s_non_tEBC_running', ...
        'tEBC frames retained with >=1 s running non-tEBC occupancy (%)', ...
        fullfile(opt.OutputDir,'pct_frames_retained_ge_1s.png'));
    if ~isempty(matchedBins) && width(matchedBins)>0
        make_ratio_histogram(matchedBins,fullfile(opt.OutputDir,'occupancy_ratio_histograms.png'));
    else
        warning('No matched spatial bins; occupancy-ratio histogram was skipped.');
    end
    make_heatmaps(heatData,opt.RepresentativeSessions,opt.OutputDir);
end

OUT = struct('detailed',detailed,'summary',summary,'matchedBins',matchedBins, ...
    'alignmentQC',quality,'options',opt);
save(fullfile(opt.OutputDir,'same_location_feasibility.mat'),'OUT');

print_recommendations(detailed);
fprintf('\nSaved feasibility outputs to %s\n',opt.OutputDir);
end

function [r,bRows] = diagnose_combination(animal,day,epoch,scheme,speedCond,task,non,bin,dt,speed,nBins,thresholds,opt)
taskOcc = occupancy_by_bin(bin,task,dt,nBins);
nonOcc = occupancy_by_bin(bin,non,dt,nBins);
taskBin = bin(task & isfinite(bin));
taskFrames = accumarray(taskBin,ones(numel(taskBin),1),[nBins 1],@sum,0);
taskOccupied = taskOcc > 0;
matched = taskOccupied & nonOcc > 0;
ratio = nan(nBins,1);
ratio(matched) = min(taskOcc(matched),nonOcc(matched)) ./ max(taskOcc(matched),nonOcc(matched));

nTaskBins = nnz(taskOccupied);
r = struct();
r.animal = animal;
r.session_date = day;
r.task_epoch = epoch;
r.binning_scheme = scheme;
r.tEBC_speed_condition = speedCond;
r.non_tEBC_speed_filter = 'speed_ge_4cm_s';
r.n_tEBC_bins = nTaskBins;
r.n_tEBC_bins_with_any_non_tEBC_running = nnz(matched);
r.pct_tEBC_bins_with_any_non_tEBC_running = pct(nnz(matched),nTaskBins);
for ti = 1:numel(thresholds)
    th = thresholds(ti); tag = threshold_tag(th);
    r.(['n_tEBC_bins_ge_' tag]) = nnz(taskOcc >= th);
    r.(['n_non_tEBC_running_bins_ge_' tag '_among_tEBC_bins']) = nnz(taskOccupied & nonOcc >= th);
    r.(['n_bins_both_ge_' tag]) = nnz(taskOcc >= th & nonOcc >= th);
    r.(['pct_bins_both_ge_' tag]) = pct(nnz(taskOcc >= th & nonOcc >= th),nTaskBins);
end

rr = ratio(isfinite(ratio));
if isempty(rr)
    r.median_occ_ratio=NaN; r.iqr_occ_ratio_low=NaN; r.iqr_occ_ratio_high=NaN;
    r.pct_bins_ratio_ge_0p25=NaN; r.pct_bins_ratio_ge_0p5=NaN; r.pct_bins_ratio_ge_0p75=NaN;
else
    qq=prctile(rr,[25 75]);
    r.median_occ_ratio=median(rr); r.iqr_occ_ratio_low=qq(1); r.iqr_occ_ratio_high=qq(2);
    r.pct_bins_ratio_ge_0p25=100*mean(rr>=.25);
    r.pct_bins_ratio_ge_0p5=100*mean(rr>=.5);
    r.pct_bins_ratio_ge_0p75=100*mean(rr>=.75);
end

totalFrames = nnz(task);
r.total_tEBC_frames = totalFrames;
r.pct_tEBC_frames_retained_any_non_tEBC_running = pct(sum(taskFrames(nonOcc>0)),totalFrames);
r.pct_tEBC_frames_retained_ge_0p5s_non_tEBC_running = pct(sum(taskFrames(nonOcc>=.5)),totalFrames);
r.pct_tEBC_frames_retained_ge_1s_non_tEBC_running = pct(sum(taskFrames(nonOcc>=1)),totalFrames);
r.pct_tEBC_frames_retained_ge_2s_non_tEBC_running = pct(sum(taskFrames(nonOcc>=2)),totalFrames);

% Optional joint spatial + fixed 2 cm/s speed-bin matching. This definition
% already exists in MI_control_rollingPost.m (SpeedBinSize = 2 cm/s).
r.speed_matching_performed = ~isempty(opt.SpeedBinSize);
r.n_tEBC_spatial_speed_bins = NaN;
r.n_tEBC_spatial_speed_bins_with_any_non_tEBC_running = NaN;
r.pct_tEBC_spatial_speed_bins_with_any_non_tEBC_running = NaN;
r.pct_tEBC_frames_retained_joint_spatial_speed = NaN;
if ~isempty(opt.SpeedBinSize) && any(task) && any(non)
    vmax=max(speed(task|non),[],'omitnan');
    lo=floor(min(speed(task|non),[],'omitnan')/opt.SpeedBinSize)*opt.SpeedBinSize;
    hi=ceil(vmax/opt.SpeedBinSize)*opt.SpeedBinSize + opt.SpeedBinSize;
    sEdges=lo:opt.SpeedBinSize:hi;
    sb=discretize(speed,sEdges); nS=numel(sEdges)-1;
    joint=nan(size(bin)); good=isfinite(bin)&isfinite(sb);
    joint(good)=sub2ind([nBins nS],bin(good),sb(good));
    taskJoint=unique(joint(task & isfinite(joint)));
    nonJoint=unique(joint(non & isfinite(joint)));
    isJointMatched=ismember(joint,nonJoint);
    r.n_tEBC_spatial_speed_bins=numel(taskJoint);
    r.n_tEBC_spatial_speed_bins_with_any_non_tEBC_running=nnz(ismember(taskJoint,nonJoint));
    r.pct_tEBC_spatial_speed_bins_with_any_non_tEBC_running= ...
        pct(r.n_tEBC_spatial_speed_bins_with_any_non_tEBC_running,numel(taskJoint));
    r.pct_tEBC_frames_retained_joint_spatial_speed=pct(nnz(task & isJointMatched),totalFrames);
end

ids=find(matched);
bRows=struct([]);
for k=1:numel(ids)
    j=ids(k);
    br=struct('animal',animal,'session_date',day,'task_epoch',epoch, ...
        'binning_scheme',scheme,'tEBC_speed_condition',speedCond, ...
        'spatial_bin',j,'tEBC_occupancy_s',taskOcc(j), ...
        'non_tEBC_running_occupancy_s',nonOcc(j),'occ_ratio',ratio(j));
    bRows=append_struct(bRows,br);
end
end

function D = load_session(rat,day,opt)
pos=get_day_field(rat,'pos',day,{'pos_'});
ts=get_day_field(rat,'Ca_ts',day,{'CA_ts_','Ca_ts_','ts_'});
cs=get_day_field(rat,'CS_times',day,{'CS_'});
pos=double(pos); cs=double(cs(:));
if size(pos,2)<3 && size(pos,1)>=3, pos=pos'; end
if size(pos,2)<3, error('Position must have [time x y] columns.'); end
pt=pos(:,1); px=pos(:,2); py=pos(:,3);
[pt,ia]=unique(pt(:),'stable'); px=px(ia); py=py(ia);
[pt,ord]=sort(pt); px=px(ord); py=py(ord);
t=coerce_timestamps(ts);
t=t(isfinite(t)); t=unique(t(:),'sorted');
if numel(t)<2, error('Fewer than two valid imaging timestamps.'); end

% Existing ca_velocity convention: speed is in row 1, timestamps in row 2.
P=[pt px py];
if exist('ca_velocity','file')==2
    V=ca_velocity(P); vt=double(V(2,:)'); vm=double(V(1,:)');
else
    vt=pt; vm=central_speed(pt,px,py);
end
[vt,iv]=unique(vt,'stable'); vm=vm(iv);
x=interp1(pt,px,t,'linear',NaN); y=interp1(pt,py,t,'linear',NaN);
speed=interp1(vt,vm,t,'linear',NaN);

[dt,rawGap,nCap]=frame_durations(t);
valid=isfinite(x)&isfinite(y)&isfinite(speed);
cs=cs(isfinite(cs));

% Prefer the project's explicit intertrial labels when their timestamp row
% can be recognized; otherwise derive intertrial as outside [CS,CS+2).
intertrial=~in_windows(t,cs,[0 2]); source='derived_outside_CS_0_2s';
if isfield(rat,'CSUS_id') && isstruct(rat.CSUS_id)
    try
        z=get_day_field(rat,'CSUS_id',day,{'CSUS_id_','CSUS_'});
        [lab,lt]=coerce_labels(z);
        if numel(lab)>1 && numel(lt)==numel(lab)
            ilab=interp1(lt,double(lab),t,'nearest',NaN);
            intertrial=(ilab==0); source='existing_CSUS_id_equals_0';
        end
    catch
    end
end

xv=px(isfinite(px)); yv=py(isfinite(py));
if isempty(xv), error('No finite interpolated x/y/speed samples.'); end
xmin=min(xv); xmax=max(xv); ymin=min(yv); ymax=max(yv);
if xmax<=xmin || ymax<=ymin, error('Degenerate spatial bounds.'); end
nXBins=ceil((xmax-xmin)/opt.FineBinCm);
nYBins=ceil((ymax-ymin)/opt.FineBinCm);
fine.x=linspace(xmin,xmax,nXBins+1);
fine.y=linspace(ymin,ymax,nYBins+1);
coarse.x=[xmin (xmin+xmax)/2 xmax]; coarse.y=[ymin (ymin+ymax)/2 ymax];

D=struct('t',t,'x',x,'y',y,'speed',speed,'dt',dt,'rawGap',rawGap, ...
    'nGapCapped',nCap,'valid',valid,'cs',cs,'intertrial',intertrial, ...
    'intertrialSource',source,'fineEdges',fine,'coarseEdges',coarse, ...
    'fineEdgeSource','CA_normalizePosData_session_bounds_ceil_linspace');
end

function t=coerce_timestamps(z)
z=double(z);
if isvector(z), t=z(:); else
    % Existing Ca_ts matrices generally store time in column 2. If only one
    % column is present, use it. Rows are handled by transposing first.
    if size(z,1)<=3 && size(z,2)>size(z,1), z=z'; end
    t=z(:,min(2,size(z,2)));
end
if max(t,[],'omitnan')>1e4, t=t/1000; end
end

function [lab,t]=coerce_labels(z)
z=double(z);
if size(z,1)==2
    a=z(1,:)'; b=z(2,:)';
elseif size(z,2)==2
    a=z(:,1); b=z(:,2);
else
    lab=[]; t=[]; return
end
% Label row has few unique integer values; timestamp row is monotonic.
if numel(unique(a(isfinite(a))))<numel(unique(b(isfinite(b))))
    lab=a; t=b;
else
    lab=b; t=a;
end
[t,ord]=sort(t); lab=lab(ord);
end

function [dt,gaps,nCap]=frame_durations(t)
gaps=diff(t); med=median(gaps(gaps>0),'omitnan');
if ~isfinite(med)||med<=0, error('Invalid imaging timestamp spacing.'); end
dt=[gaps; med];
% A dropped timestamp should not be counted as continuous sampled occupancy.
bad=~isfinite(dt)|dt<=0|dt>3*med; nCap=nnz(bad); dt(bad)=med;
end

function v=central_speed(t,x,y)
v=nan(size(t));
if numel(t)>2
    v(2:end-1)=hypot(x(3:end)-x(1:end-2),y(3:end)-y(1:end-2))./(t(3:end)-t(1:end-2));
    v(1)=v(2); v(end)=v(end-1);
end
end

function x=get_day_field(rat,parent,day,prefixes)
if ~isfield(rat,parent), error('Missing rat.%s.',parent); end
S=rat.(parent); fn=fieldnames(S);
candidates={day};
for i=1:numel(prefixes), candidates{end+1}=[prefixes{i} day]; end %#ok<AGROW>
idx=[];
for i=1:numel(candidates), idx=find(strcmp(fn,candidates{i}),1); if ~isempty(idx),break,end,end
if isempty(idx), idx=find(contains(fn,day),1); end
if isempty(idx), error('Missing %s session field for %s.',parent,day); end
x=S.(fn{idx});
end

function days=select_days(rat,mode)
if exist('autoDateList','file')==2
    try, allDays=cellstr(autoDateList(rat)); catch, allDays={}; end
else, allDays={}; end
if isempty(allDays)
    fn=fieldnames(rat.pos); allDays=regexprep(fn,'^pos_',''); allDays=allDays(:)';
end
if iscell(mode), days=cellstr(mode); return, end
mode=char(mode);
if strcmpi(mode,'all'), days=allDays; return, end
if strcmpi(mode,'last3_to_An') && isfield(rat,'An')
    k=find(strcmp(allDays,rat.An),1);
    if ~isempty(k), days=allDays(max(1,k-2):k); return, end
end
days=allDays(max(1,numel(allDays)-2):end);
end

function tf=in_windows(t,anchors,win)
tf=false(size(t));
for i=1:numel(anchors), tf=tf|(t>=anchors(i)+win(1)&t<anchors(i)+win(2)); end
end

function [bin,nY,nX]=position_bins(x,y,edges)
bx=discretize(x,edges.x); by=discretize(y,edges.y);
nX=numel(edges.x)-1; nY=numel(edges.y)-1;
bin=nan(size(x)); good=isfinite(bx)&isfinite(by);
bin(good)=sub2ind([nY nX],by(good),bx(good));
end

function occ=occupancy_by_bin(bin,mask,dt,nBins)
good=mask&isfinite(bin)&isfinite(dt);
occ=accumarray(bin(good),dt(good),[nBins 1],@sum,0);
end

function tag=threshold_tag(x)
tag=strrep(sprintf('%gs',x),'.','p');
end

function y=pct(a,b)
if b>0, y=100*a/b; else, y=NaN; end
end

function S=append_struct(S,x)
if isempty(S), S=x; else, S(end+1)=x; end
end

function S=append_struct_array(S,x)
if isempty(x), return, end
if isempty(S), S=x; else, S(end+1:end+numel(x))=x; end
end

function summary=summarize_table(T)
isNum=varfun(@isnumeric,T,'OutputFormat','uniform');
metrics=T.Properties.VariableNames(isNum);
metrics=setdiff(metrics,{'n_tEBC_bins'},'stable'); % n bins still useful; restore below
metrics=[{'n_tEBC_bins'} metrics];
[G,epoch,scheme,speedCond]=findgroups(T.task_epoch,T.binning_scheme,T.tEBC_speed_condition);
nG=max(G); out=struct([]);
for g=1:nG
    s=struct('task_epoch',epoch{g},'binning_scheme',scheme{g}, ...
        'tEBC_speed_condition',speedCond{g},'n_sessions',nnz(G==g));
    for m=1:numel(metrics)
        name=metrics{m}; v=T.(name)(G==g); v=v(isfinite(v));
        if isempty(v), vals=nan(1,5); else
            vals=[mean(v),std(v)/sqrt(numel(v)),median(v),min(v),max(v)];
        end
        s.([name '_mean'])=vals(1); s.([name '_sem'])=vals(2);
        s.([name '_median'])=vals(3); s.([name '_min'])=vals(4); s.([name '_max'])=vals(5);
    end
    out=append_struct(out,s);
end
summary=struct2table(out);
end

function make_bar_plot(T,metric,ylabelText,outfile)
sessions=unique(strcat(T.animal," | ",T.session_date),'stable');
conds={'full_task|fine_2p5cm','full_task|coarse_2x2','trace_only|fine_2p5cm','trace_only|coarse_2x2'};
speedNames={'all_tEBC','running_tEBC'};
for sc=1:2
    figure('Color','w','Position',[80 80 1500 650]);
    for c=1:4
        parts=strsplit(conds{c},'|');
        keep=strcmp(T.task_epoch,parts{1})&strcmp(T.binning_scheme,parts{2})&strcmp(T.tEBC_speed_condition,speedNames{sc});
        for i=1:numel(sessions)
            key=strcat(T.animal," | ",T.session_date);
            j=find(keep&key==sessions(i),1); if isempty(j), Y(i,c)=NaN; else, Y(i,c)=T.(metric)(j); end %#ok<AGROW>
        end
    end
    bar(Y); ylim([0 100]); grid on; ylabel(ylabelText);
    xticks(1:numel(sessions)); xticklabels(sessions); xtickangle(45);
    legend({'full fine','full coarse','trace fine','trace coarse'},'Location','bestoutside');
    title(strrep(speedNames{sc},'_',' '));
    [pth,nm,ext]=fileparts(outfile); exportgraphics(gcf,fullfile(pth,[nm '_' speedNames{sc} ext]),'Resolution',200); close(gcf);
    clear Y
end
end

function make_ratio_histogram(B,outfile)
figure('Color','w','Position',[100 100 1100 450]);
schemeNames={'fine_2p5cm','coarse_2x2'};
for i=1:2
    subplot(1,2,i); keep=strcmp(B.binning_scheme,schemeNames{i});
    histogram(B.occ_ratio(keep),0:.05:1); xlim([0 1]); grid on;
    xlabel('min occupancy / max occupancy'); ylabel('Matched spatial bins');
    title(strrep(schemeNames{i},'_',' '));
end
exportgraphics(gcf,outfile,'Resolution',200); close(gcf);
end

function make_heatmaps(H,nRep,outdir)
if nRep<1||isempty(H), return, end
keys=unique(strcat({H.animal}'," | ",{H.session_date}'),'stable');
scores=zeros(size(keys));
for i=1:numel(keys)
    q=strcmp(strcat({H.animal}'," | ",{H.session_date}'),keys{i}); scores(i)=max([H(q).nValid]);
end
[~,ord]=sort(scores,'descend'); keys=keys(ord(1:min(nRep,numel(ord))));
for i=1:numel(keys)
    q=find(strcmp(strcat({H.animal}'," | ",{H.session_date}'),keys{i}));
    figure('Color','w','Position',[50 50 1300 850]);
    for j=1:numel(q)
        h=H(q(j)); base=(j-1)*3;
        A=reshape(h.taskOcc,[h.nY h.nX]); N=reshape(h.nonOcc,[h.nY h.nX]); O=A>0&N>0;
        subplot(numel(q),3,base+1); imagesc(A); axis image; colorbar; title([strrep(h.task_epoch,'_',' ') ' ' strrep(h.binning_scheme,'_',' ') ' tEBC']);
        subplot(numel(q),3,base+2); imagesc(N); axis image; colorbar; title('running non-tEBC');
        subplot(numel(q),3,base+3); imagesc(O); axis image; colorbar; caxis([0 1]); title('overlap mask');
    end
    sgtitle(keys{i}); exportgraphics(gcf,fullfile(outdir,sprintf('occupancy_heatmaps_rep%d.png',i)),'Resolution',200); close(gcf);
end
end

function print_recommendations(T)
fprintf('\n=== Feasibility recommendations (session means) ===\n');
[G,e,b,s]=findgroups(T.task_epoch,T.binning_scheme,T.tEBC_speed_condition);
for g=1:max(G)
    a=mean(T.pct_tEBC_frames_retained_ge_1s_non_tEBC_running(G==g),'omitnan');
    c=mean(T.pct_bins_both_ge_1s(G==g),'omitnan');
    if a>=80 && c>=70, msg='Same-location control appears feasible.';
    else, msg='Same-location control is likely underpowered or biased because spatial occupancy overlap is insufficient.'; end
    fprintf('%s | %s | %s: frames>=1s %.1f%%; bins both>=1s %.1f%%. %s\n',e{g},b{g},s{g},a,c,msg);
end

for sc={'all_tEBC','running_tEBC'}
    ff=mean(T.pct_tEBC_frames_retained_ge_1s_non_tEBC_running(strcmp(T.task_epoch,'full_task')&strcmp(T.binning_scheme,'fine_2p5cm')&strcmp(T.tEBC_speed_condition,sc{1})),'omitnan');
    fc=mean(T.pct_tEBC_frames_retained_ge_1s_non_tEBC_running(strcmp(T.task_epoch,'full_task')&strcmp(T.binning_scheme,'coarse_2x2')&strcmp(T.tEBC_speed_condition,sc{1})),'omitnan');
    if ff<80 && fc>=80
        fprintf('%s: Exact fine-scale spatial matching is underpowered in this open-field dataset, while coarse spatial matching is feasible.\n',sc{1});
    end
end
full=mean(T.pct_tEBC_frames_retained_ge_1s_non_tEBC_running(strcmp(T.task_epoch,'full_task')),'omitnan');
trace=mean(T.pct_tEBC_frames_retained_ge_1s_non_tEBC_running(strcmp(T.task_epoch,'trace_only')),'omitnan');
if trace+10<full
    fprintf('Trace-only retention is substantially lower: the short trace window has too few frames for reliable same-location matching.\n');
end
fprintf('Do not proceed to the final neural analysis unless feasibility is clearly good.\n');
end
