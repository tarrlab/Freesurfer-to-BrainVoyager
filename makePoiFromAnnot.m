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
