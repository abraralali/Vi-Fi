# Vi-Fi
## IPSN'22 submission

### Directory over view:

#### "data_collection_preprocessing/" contains:
  * scripts of extracting metadata from ZED .svo files
  * scripts of ground truth labeling
  * scripts of data pre-processing and organizing
  
  
#### "data_prepare_for_train/" contains:
  * scripts of constructing data samples for deep affinity matrix training
  * scripts of spliting training and testing set

#### "aff_mat_train_test/" contains:
  * scripts of model architecture
  * scripts of training and testing

#### "Bipartite_Association/" contains:
 * MATLAB scripts to prepare and process data of the visual tracker, IMU readings and FTM
 * To run the code first run DataPreparation.m then run ZedMain.m
 * You need to change sequences_path and src_path to read data and scripts from your machine path

### TODO:
  * Execution instructions
  * Dataset Link
