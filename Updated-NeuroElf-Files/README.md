To add label functionality as a feature into the "Load VOI file" function in the NeuroElf GUI when loading volumetric NIfTI Freesurfer atlases created using mri_convert.

1. Replace the ne_loadcluster.m file in the @neuroelf/private folder (make a backup of the old one, just in case)
2. Replace the voi_ImportClusters.m file in @xff/private folder (same)
3. Place the FreeSurferColorLUT.txt file (which is part of FreeSurfer, but just in case) into the folder with the segmentation.
4. Restart NeuroElf (or, to make absolutely sure, first MATLAB and then NeuroElf)
5. Use the "Load VOI file..." function (button), select the segmentation file, threshold 1, say "Yes" (to separate clusters), THEN select the FreeSurferColorLUT.txt file
