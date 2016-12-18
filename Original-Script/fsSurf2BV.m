function fsSurf2BV(subjName, varargin)
%function fsSurf2BV(subjName, varargin)
%
%Convert surfaces and atlases from freesurfer to BrainVoyager
%surfaces (srf) and patches of interest (poi)
%
%jens.schwarzbach@unitn.it
%
%With default options creates BrainVoyager files for a given participant in
%subdirectory <subjName>:
%brain.vmr
%<subjName>_lh_inflated.srf
%<subjName>_lh_pial.srf
%<subjName>_lh_smoothwm.srf
%<subjName>_lh_sphere.srf
%<subjName>_rh_inflated.srf
%<subjName>_rh_pial.srf
%<subjName>_rh_smoothwm.srf
%<subjName>_rh_sphere.srf
%<subjName>_lh_aparc.a2009s.annot.poi
%<subjName>_lh_aparc.DKTatlas40.annot.poi
%<subjName>_lh_Yeo2011_7Networks_N1000.annot.poi
%<subjName>_rh_aparc.a2009s.annot.poi
%<subjName>_rh_aparc.DKTatlas40.annot.poi
%<subjName>_rh_Yeo2011_7Networks_N1000.annot.poi
%if not yet existing, creates freesurfer annotation files in subject space:
%SUBJECT_DIR/<subjName>/label/lh.Yeo2011_17Networks_N1000.annot
%SUBJECT_DIR/<subjName>/label/rh.Yeo2011_17Networks_N1000.annot
%
%prerequisites:
%%freesurfer
%%surfing toolbox
%%neuroelf
%
%Configurable:
% FREESURFER_HOME
% SUBJECTS_DIR
% projectDir
% hemis
% surfaceTypes
% atlas
% nSubClusters
%
%Example calls
%%for creating all default conversions of freesurfer's subject bert
%fsSurf2BV('bert') %requires write permission to $FREESURFER_HOME/subjects/bert/label
%
%%for creating just the inflated surfaces and pois for bert and the DK40 atlas 
% Cfg.surfaceTypes = {'inflated'};
% Cfg.atlas = {'aparc.DKTatlas40.annot'};
% fsSurf2BV('bert', Cfg)
%
%%for creating just the inflated surfaces and pois for a list of subjects
%%in a given SUBJECTS_DIR whose IDs all start with 'subj_', 
%%and restricting poi-creation to the DK40 atlas 
% Cfg.SUBJECTS_DIR = '~/MRI/myProject/subjects'; %directory and file conventions according to freesurfer
% Cfg.surfaceTypes = {'inflated'};
% Cfg.atlas = {'aparc.DKTatlas40.annot'};
% fsSurf2BV('subj_', Cfg)

if nargin > 1
    Cfg = varargin{1};
else
    Cfg = [];
end

%ENVIRONMENT VARIABLES
if ~isfield(Cfg, 'FREESURFER_HOME'), Cfg.FREESURFER_HOME = '/Applications/freesurfer'; else end;
if ~isfield(Cfg,'SUBJECTS_DIR'), Cfg.SUBJECTS_DIR = '/Applications/freesurfer/subjects'; else end;
%CONFIGURATION
if ~isfield(Cfg, 'projectDir'), Cfg.projectDir = ''; else end;
if ~isfield(Cfg, 'hemis'), Cfg.hemis = {'lh', 'rh'}; else end
if ~isfield(Cfg, 'surfaceTypes'), Cfg.surfaceTypes = {'inflated', 'pial', 'smoothwm', 'sphere'}; else end;
if ~isfield(Cfg, 'atlas'), Cfg.atlas = {'BA.annot', 'aparc.a2009s.annot', 'aparc.DKTatlas40.annot', 'Yeo2011_7Networks_N1000.annot', 'Yeo2011_17Networks_N1000.annot'}; else end;
%; 'BA.annot', Note: 'PALS_B12_Brodmann.annot' fails
if ~isfield(Cfg, 'nSubClusters'), Cfg.nSubClusters = 0; else end;%for segmenting each labeled area into spatially distinct clusters

display(Cfg)

setenv( 'SUBJECTS_DIR', Cfg.SUBJECTS_DIR);
setenv( 'FREESURFER_HOME', Cfg.FREESURFER_HOME); %such as /Applications/freesurfer/bin

if isdir(fullfile(Cfg.SUBJECTS_DIR, subjName))
    %if provided with a path to a directory, interpret this as the
    %subject-directory
    Cfg.currentSubjectName = subjName;
    processSubject(subjName, Cfg)
else
    %interpret this as a search string
    d=dir(fullfile(Cfg.projectDir, subjName));
    numSubs = numel(d);
    fprintf('Found %d subjects matching %s\n', numSubs, subjName);
    h = waitbar(0,'Creating vmrs, srfs, and pois ...');
    for iSub = 1:numel(d)
        waitbar(iSub/numel(d), sprintf('Creating vmrs, srfs, and pois %d/%d', iSub, numSubs));
        Cfg.currentSubjectName = d(iSub).name;
        processSubject(d(iSub).name, Cfg)
    end
    close(h);
    
end

function processSubject(subjName, Cfg)
fprintf(1, '--------------------------------\n');
fprintf(1, 'Processing %s\n', subjName);
fprintf(1, '--------------------------------\n');

Cfg.T1_DIR = fullfile(Cfg.projectDir, subjName, 'T1', '001');
Cfg.bvDir = fullfile(Cfg.projectDir, subjName, 'bv');

%--------------------------------------------------------------------------
%MAKE VMR
%--------------------------------------------------------------------------
fprintf(1, 'MAKE VMR:\n');
T1_fn = fullfile(Cfg.SUBJECTS_DIR, subjName, 'mri', 'brain.mgz');
fprintf(1, 'LOADING %s\n', T1_fn);
T1 = load_mgh(T1_fn);

T1p = permute(T1, [3, 2, 1]);
T1f = flip(T1p, 1);

vmr = xff('new:vmr');
vmr.VMRData = uint8(T1f);
if ~exist(Cfg.bvDir, 'dir')
    mkdir(Cfg.bvDir)
end
vmr_fn = fullfile(Cfg.bvDir, 'brain.vmr');
fprintf(1, 'SAVING %s\n', vmr_fn);
vmr.SaveAs(vmr_fn);
vmr.ClearObject;
fprintf(1, '\n');

%--------------------------------------------------------------------------
%CONVERT SURFACES
%--------------------------------------------------------------------------
fprintf('CONVERT SURFACES:\n')
fprintf(1, '%s\n', Cfg.surfaceTypes{:});
for iSrf = 1:numel(Cfg.surfaceTypes)
    for iHemi = 1:numel(Cfg.hemis)
        strHemi = Cfg.hemis{iHemi};
        %SUBJECT SPECIFIC SURFACE
        fnSurf = fullfile(Cfg.SUBJECTS_DIR, subjName, 'surf', [strHemi, '.', Cfg.surfaceTypes{iSrf}]);
        %READ SURFACE
        [v,f]=freesurfer_read_surf(fnSurf);
        %WRITE BV SURFACE
        outName = fullfile(Cfg.bvDir, [subjName, '_', strHemi, '_', Cfg.surfaceTypes{iSrf}, '.srf']);
        fprintf('WRITING %s (%d vertices)\n\n\n', outName, size(v, 1));
        surfing_write(outName, v, f);
        
    end
end
fprintf(1, '\n');

%--------------------------------------------------------------------------
%CREATE POIS
%--------------------------------------------------------------------------
fprintf(1, 'CREATING POI-FILES FROM ATLASES:\n');
if ~iscell(Cfg.atlas)
    Cfg.atlas = {Cfg.atlas};
end

fprintf(1, '%s\n', Cfg.atlas{:});

for iHemi = 1:numel(Cfg.hemis)
    strHemi = Cfg.hemis{iHemi};
    
    %SUBJECT SPECIFIC SURFACE (inflated may work best)
    fnSurf = fullfile(Cfg.SUBJECTS_DIR, subjName, 'surf', [strHemi, '.', Cfg.surfaceTypes{1}]);
    %READ SURFACE
    [v,~]=freesurfer_read_surf(fnSurf);
    
    for iAtlas=1:numel(Cfg.atlas)
        %e.g. 7 NETWORK ATLAS
        fnAnnot = fullfile(Cfg.SUBJECTS_DIR, subjName, 'label', [strHemi, '.', Cfg.atlas{iAtlas}]);
        fnPoi = fullfile(Cfg.bvDir, [subjName, '_', strHemi, '_', Cfg.atlas{iAtlas}, '.poi']);
        makePoiFromAnnot(fnAnnot, fnPoi, v, Cfg)
    end
end


function makePoiFromAnnot(inAnnotFileName, outPoiFileName, v, Cfg)
%READ ANNOTION (e.g. RS-ATLAS)
if ~exist(inAnnotFileName, 'file')
    fprintf('Mapping atlas to subject space.\n');
    %which hemisphere are we talking about?
    [thisSubjectDir, areaID, ~] = fileparts(inAnnotFileName);
    [thisHemi, rmdr] = strtok(areaID, '.');
    %which atlas are we talking about?
    thisAtlas = rmdr(2:end);
    outAnnotFileName = fullfile(thisSubjectDir, sprintf('%s.%s.annot', thisHemi, thisAtlas));

    %mri_surf2surf --srcsubject fsaverage --trgsubject $1 --hemi lh --sval-annot $SUBJECTS_DIR/fsaverage/label/lh.Yeo2011_7Networks_N1000.annot --tval $SUBJECTS_DIR/$1/label/lh.Yeo2011_7Networks_N1000.annot
    strCmd = sprintf('! /Applications/freesurfer/bin/mri_surf2surf --srcsubject fsaverage --trgsubject %s --hemi %s --sval-annot %s/fsaverage/label/%s.%s --tval %s',...
        Cfg.currentSubjectName, thisHemi, Cfg.SUBJECTS_DIR, thisHemi, thisAtlas, outAnnotFileName);
    fprintf(1, '**********************************************************\n');
    fprintf(1, 'EXECUTING SHELL COMMAND: %s\n', strCmd);
    fprintf(1, '**********************************************************\n');
    eval(strCmd)
end

fprintf(1, '-------------------------------------------------------------------------------\n');
fprintf(1, 'Current atlas (annotation file)\n');
fprintf(1, '%s\n', inAnnotFileName);
fprintf(1, '-------------------------------------------------------------------------------\n');

[vertices, label, colortable] = read_annotation(inAnnotFileName);


%MAKE A POI FILE FROM ANNOTATION
poi = xff('new:poi');
labels = unique(label);
nLabels = numel(labels);

for i = 1:nLabels
    idx = find(label == labels(i));
    idxColorTable = find(labels(i) == colortable.table(:, end));
    if isempty(idxColorTable)
        thisStructName = 'NA';
        thisColor = [0 0 0];
    else
        thisStructName = colortable.struct_names{idxColorTable};
        thisColor = colortable.table(idxColorTable, 1:3);
    end
    if i > 1
        %MAKE A COPY OF FIRST ENTRY
        poi.POI(i) = poi.POI(1);
    end
    %UPDATE WITH CURRENT DATA
    poi.POI(i).Name = thisStructName;
    poi.POI(i).Color = thisColor;
    poi.POI(i).NrOfVertices = length(idx);
    poi.POI(i).Vertices = idx;
    %choose a vertex in the center of an are as reference (this may not be the best one)
    centroid = mean(v(idx, :));
    delta = [(v(idx, 1) - centroid(1)).^2, (v(idx, 2) - centroid(2)).^2, (v(idx, 3) - centroid(3)).^2];
    [~, idxminval] = min(sum(delta, 2));
    poi.POI(i).LabelVertex = idx(idxminval); 
end
poi.NrOfMeshVertices = size(vertices, 1);
poi.NrOfPOIs = nLabels;
fprintf(1, 'SAVING %s\n\n', outPoiFileName);
poi.SaveAs(outPoiFileName);
poi.ClearObject;


if Cfg.nSubClusters > 0
    %MAKE A CLUSTERED POI
    poiClustered = xff('new:poi');
    labelVec = colortable.table(2:end, end); %WE ARE NOT INTERESTED IN THE FIRST LABEL
    structNames = colortable.struct_names(2:end);
    colors = colortable.table(2:end, 1:3);
    nLabels = numel(structNames);
    nClusters = Cfg.nSubClusters;
    poiIdx = 0;
    for i = 1:nLabels
        cases = find(label == labelVec(i));
        if ~isempty(cases)
            fprintf(1, 'Clustering network %d: %s ... \n', i, structNames{i});
            poiIdx = poiIdx + 1;
            if poiIdx > 1
                %MAKE A COPY OF FIRST ENTRY
                poiClustered.POI(poiIdx) = poiClustered.POI(1);
            end
            
            vertexIdx = vertices(cases) + 1;
            poiClustered.POI(poiIdx).Name = structNames{i};
            poiClustered.POI(poiIdx).Color = floor(colors(i, :));
            poiClustered.POI(poiIdx).NrOfVertices = length(cases);
            poiClustered.POI(poiIdx).Vertices = vertexIdx;
            
            %T = clusterdata(v(vertexIdx, :),nClusters);
            T = kmeans(v(vertexIdx, :), nClusters);
            clustervals = unique(T);
            for c = 1:numel(clustervals)
                poiIdx = poiIdx + 1;
                idxThisCluster = find(T == clustervals(c));
                
                if poiIdx > 1
                    %MAKE A COPY OF FIRST ENTRY
                    poiClustered.POI(poiIdx) = poiClustered.POI(1);
                end
                %UPDATE WITH CURRENT DATA
                poiClustered.POI(poiIdx).Name = sprintf('%s_c%d', structNames{i}, c);
                poiClustered.POI(poiIdx).Color = floor(colors(i, :)/c); %EACH CLUSTER BECOMES DARKER
                poiClustered.POI(poiIdx).NrOfVertices = length(idxThisCluster);
                poiClustered.POI(poiIdx).Vertices = vertexIdx(idxThisCluster);
                poiClustered.POI(i).LabelVertex = poiClustered.POI(poiIdx).Vertices(1); %lazy solution
            end
        end
    end
    poiClustered.NrOfMeshVertices = size(vertices, 1);
    poiClustered.NrOfPOIs = poiIdx;
    [strPath, strPre, strPost] = fileparts(outPoiFileName);
    outPoiFileNameClustered = fullfile(strPath, [strPre, '_clustered', strPost]);
    fprintf(1, 'SAVING %s\n', outPoiFileNameClustered);
    
    poiClustered.SaveAs(outPoiFileNameClustered);
    poiClustered.ClearObject;
end

