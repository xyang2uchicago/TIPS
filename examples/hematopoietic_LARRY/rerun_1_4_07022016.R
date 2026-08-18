source('F:/projects/TIPS/source/GSE140802_lineage_tracking/0_Experiment1_invitro_creatSeurat.R')
source('F:/projects/TIPS/source/GSE140802_lineage_tracking/1_Experiment1_invitro.umap_v2.R')
range(larry[['RNA']]$data)

source('F:/projects/TIPS/source/GSE140802_lineage_tracking/2_invitro_NMtrajectory_creatSeurat_v2.R')
range(larry_nm[['RNA']]$data) #[1] 0.00000 6.73168

source('F:/projects/TIPS/source/GSE140802_lineage_tracking/2.2_Experiment1.NMtrajectoryCell_invitro.umap_v2.R')
range(larry_nm[['RNA']]$data) #[1] 0.00000 6.73168


source('F:/projects/TIPS/source/GSE140802_lineage_tracking/3_cluster_annotate_NMtrajectoryCell.R')
range(larry_nm[['RNA']]$data) #[1] 0.00000 6.73168

## on midway3 manually: 
source('F:/projects/TIPS/source/GSE140802_lineage_tracking/4_invitro_NMtrajectory_trajectory_v2.R')
## on midway3 manually:
source('F:/projects/TIPS/source/GSE140802_lineage_tracking/2.1_normcount_2_scanpy.R')

mv timeGuided_NMT_NMbranch*.pdf three_cluster_only/.
mv MST_trajectory_aggregatedBy*.pdf three_cluster_only/.