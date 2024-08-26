# LaneEffects

Here is supplementary information and dataset information for the Babaei et al. manuscript entitled '**Evaluating Ordination and Outlier Methods for Lane Effect Mitigation in Reduced Representation Sequencing Data**'.
Authors: Sarah Babaei, Nathalie M. LeBlanc, Scott A. Pavey, Cassidy C. D'Aloia

See **SupplementaryFiltering/txt** for more information regarding filtering for the original and lane-effect free datasets.

See **LaneEffectMethods.txt** for further information on each method and where to find its corresponding dataset.
Please note that all GENEPOP files have been compressed in .zip folders with matching names.

BayeScan run information can be found in **BayeScan_CombinedDataset** and **BayeScan_GregOnly**. Combined Dataset is the combined dataset, HiSeq+NovaSeq, and the GregOnly is only HiSeq individuals. Parameters for each BayeScan run can be found in the .png file in each folder.

**NonNuclearSNPs.txt:** contains a list of Non-Nuclear SNPs that were removed from all datasets. 

**PCAloadingsSQUARED_neutral.txt**: PCA loadings from PCA
**DAPCLOADINGS_NEUTRAL_byBATCHno3Ps_testmissingindv_0.2_library1rem_NonNucrem.txt**: DAPC loadings from DAPC that ONLY included individuals that were present on both batches. This means this does not include any 3Ps individuals. 
