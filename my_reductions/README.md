README for `my_reductions/` Folder

Detailed function descriptions.

Data usage types: trial-flattened, time-flattened, neuron-flattened.

============================================


### baselineCheck.m

**Description:** Checks baseline activity across trials/days (quality control).



### bayesDecoder.m

**Description:** Bayesian classifier → decodes trial outcome from neural data.

**Data used:** trial-flattened


### bayesDecoder_REFgrouped.m

**Description:** Bayesian decoder trained on reference day → tested on grouped days.

**Data used:** trial-flattened


### f1score.m

**Description:** Computes F1 score for binary classifier predictions.

**Data used:** n/a (metric function)


### getDayMatrixFromStruct.m

**Description:** Extracts day-specific neural data matrix from master struct.

**Data used:** n/a (data extraction)


### isBadDay.m

**Description:** Identifies bad days (insufficient data or quality).

**Data used:** n/a (QC)


### kernelDecoder.m

**Description:** Kernel SVM → decodes trial outcome from neural data.

**Data used:** trial-flattened


### kernelDecoder2D.m

**Description:** Kernel SVM using 2D PCA projection of neural data.

**Data used:** trial-flattened (PCA 2D)


### kernelDecoder_REF.m

**Description:** Kernel decoder trained on reference day, tested on other days.

**Data used:** trial-flattened


### kernelDecoder_REFgrouped.m

**Description:** Kernel decoder reference day → tested on grouped test days.

**Data used:** trial-flattened


### kernelDecoder_REFgrouped_binned.m

**Description:** Same as above, but with time-binned data.

**Data used:** time-flattened


### linearDecoder.m

**Description:** Linear SVM → decodes trial outcome from neural data.

**Data used:** trial-flattened


### linearDecoder_REF.m

**Description:** Linear decoder trained on reference day → tested on other days.

**Data used:** trial-flattened


### plotAUROC_byPC.m

**Description:** Computes AUROC per PCA dimension → visualizes discriminability in PCA space.

**Data used:** trial-flattened (PCA space)


### plotAUROC_byUMAP.m

**Description:** Computes AUROC per UMAP dimension → visualizes discriminability in UMAP space.

**Data used:** trial-flattened (UMAP space)


### quantifyUMAPManifolds.m

**Description:** Quantifies properties of UMAP manifolds (e.g. separation, overlap, stability).

**Data used:** depends (often time-flattened or trial-flattened)


### runPCA_fromStruct.m

**Description:** Runs PCA on trial-level neural data → returns PCA space.

**Data used:** trial-flattened


### runUMAP_fromStruct.m

**Description:** Runs UMAP on trial-level neural data → returns UMAP space + AUROC.

**Data used:** trial-flattened


### sortgrid.m

**Description:** Sorting utility for grid search results.

**Data used:** n/a (utility)


### testTrainingDayCombos.m

**Description:** Tests multiple training/testing day combinations for decoders.

**Data used:** trial-flattened


### UMAP_cellsByOutcome.m

**Description:** Runs UMAP on neurons — each neuron is a point. Visualizes neuron-level representations of trial outcomes.

**Data used:** neuron-flattened


### UMAP_gridSearch.m

**Description:** Performs parameter grid search for UMAP hyperparameters (n_neighbors, min_dist).

**Data used:** trial-flattened


### UMAP_trajectory_byOutcome.m

**Description:** Plots UMAP trajectories over time by trial outcome (correct/incorrect).

**Data used:** time-flattened


### UMAP_trajectory_significance.m

**Description:** Tests significance of UMAP trajectory separation using permutation testing.

**Data used:** time-flattened


### UMAP_trajectory_variability.m

**Description:** Quantifies within-group (trial type) variability of UMAP trajectories.

**Data used:** time-flattened


### UMAP_trajectory_variability_evolution.m

**Description:** Tracks how UMAP trajectory variability evolves across days.

**Data used:** time-flattened


### UMAPcells_gridSearch.m

**Description:** Grid search for UMAP applied to neuron space (similar to UMAP_cellsByOutcome).

**Data used:** neuron-flattened


### UMAPtrial_crossDayDecoder

**Description:** Runs UMAP per day → visualizes day-to-day evolution of trial-level UMAP space using a reference day.

**Data used:** trial-flattened

### CrossDayDecoding_UMAP_gridSearch

**Description:** Grid search version of UMAPtrial_crossDayDecoder: Runs UMAP per day → visualizes day-to-day evolution of trial-level UMAP space using a reference day.

**Data used:** trial-flattened

### UMAPtrial_dayByDay.m

**Description:** Runs UMAP per day → visualizes day-to-day evolution of trial-level UMAP space, using cross validated own day.

**Data used:** trial-flattened


### UMAPtrial_dayByDay_gridSearch.m

**Description:** UMAP grid search for optimal parameters on trial-level UMAP (day-by-day).

**Data used:** trial-flattened


### UMAPtrial_outcome.m

**Description:** UMAP on trial data → trains SVM to decode trial outcome → reports classifier accuracy.

**Data used:** trial-flattened


# Function Data Type Chart

| Function                                | Trial-flattened   | Time-flattened   | Neuron-flattened   |
|:----------------------------------------|:------------------|:-----------------|:-------------------|
| UMAP_cellsByOutcome.m                   |                   |                  | ✅                 |
| UMAP_gridSearch.m                       | ✅                |                  |                    |
| UMAP_trajectory_byOutcome.m             |                   | ✅               |                    |
| UMAP_trajectory_significance.m          |                   | ✅               |                    |
| UMAP_trajectory_variability.m           |                   | ✅               |                    |
| UMAP_trajectory_variability_evolution.m |                   | ✅               |                    |
| UMAPcells_gridSearch.m                  |                   |                  | ✅                 |
| UMAPtrial_crossDayDecoder               | ✅                |                  |                    |
| UMAPtrial_dayByDay.m                    | ✅                |                  |                    |
| UMAPtrial_dayByDay_gridSearch.m         | ✅                |                  |                    |
| UMAPtrial_outcome.m                     | ✅                |                  |                    |
| baselineCheck.m                         |                   |                  |                    |
| bayesDecoder.m                          | ✅                |                  |                    |
| bayesDecoder_REFgrouped.m               | ✅                |                  |                    |
| CrossDayDecoding_UMAP_gridSearch        | ✅                |                  |                    |
| f1score.m                               |                   |                  |                    |
| getDayMatrixFromStruct.m                |                   |                  |                    |
| isBadDay.m                              |                   |                  |                    |
| kernelDecoder.m                         | ✅                |                  |                    |
| kernelDecoder2D.m                       | ✅                |                  |                    |
| kernelDecoder_REF.m                     | ✅                |                  |                    |
| kernelDecoder_REFgrouped.m              | ✅                |                  |                    |
| kernelDecoder_REFgrouped_binned.m       |                   | ✅               |                    |
| linearDecoder.m                         | ✅                |                  |                    |
| linearDecoder_REF.m                     | ✅                |                  |                    |
| plotAUROC_byPC.m                        | ✅                |                  |                    |
| plotAUROC_byUMAP.m                      | ✅                |                  |                    |
| quantifyUMAPManifolds.m                 | ✅                | ✅               |                    |
| runPCA_fromStruct.m                     | ✅                |                  |                    |
| runUMAP_fromStruct.m                    | ✅                |                  |                    |
| sortgrid.m                              |                   |                  |                    |
| testTrainingDayCombos.m                 | ✅                |                  |                    |
