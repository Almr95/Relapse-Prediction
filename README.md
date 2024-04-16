# Relapse-Prediction

From PAPER LINK

This set of tools were produced for a study on relapse prediction in childhood leukemia with flow cytometry data. The code is structured according to a specific workflow that can be consulted in the paper linked above, but it is general enough to be applied to any other disease for relapse prediction task (or other outcomes of interest)

The starting point is a database of patients for which we have flow cytometry data at diagnosis. Each flow cytometry experiment is divided into several tubes, each tube with a different set of markers. Each patient is annotated according to the outcome of the disease (relapse VS no relapse). The workflow is the following:

**0 - Manual preprocessing**
    Each flow cytometry file needs to be manually preprocessed to remove debris and dobulets and to explore possible errors in acquisition and other anomalies.

**1 - Preprocessing**
    For each patient, the tubes are imported and preprocessed. This involves compensating, transforming (Logicle), renaming channels and normalizing the data. 

**2 - Merging**
    The set of tubes from each patient is intergrated into a single file by means of file merging algorithms. 

**3 - Patient selection**
    The tradeoff between number of patients and markers is explored and a final selection of both is made.

**4 - Visualization**
    The merged flow cytometry files from the selected patients are imported and visualized together by means of clustering (FlowSOM) and dimensionality reduction (UMAP)

**5 - Feature extraction**
    For each metacluster, the flow cytometry dstributions (marker expression) are summarized by means of mean, median, standard deviation, skewness and kurtosis. Different datasets with this information are created in order to be used for classification-
    Information on cell abundance per cluster is also extracted. Preliminary statistical information is also extracted and plotted.

**6 - Classification**
    Each dataset obtained in the previous step is fed into a classification routine (a nested cross-validation loop that includes four machine learning algorithms: KNN, Naive Bayes, Random Forest, linear SVM) in order to check their predictive value.
    The results are summarized and plotted.


Along with the codes that perform each of this steps (in R markdown), we provide additional R code containing auxiliary code and subroutines. 

The data employed for the study cited at the beginning of this file is available at http://flowrepository.org/id/FR-FCM-Z7A2. The selection of patients used for most of the results are also available here, in the Selection subfolder, along with an excel file containing the outcome annotation. These have already been manually filtered, computationally preprocessed, merged and subsampled (Steps 0-1-2-3). 
