# Evaluation of Fast k-means

## Introduction
This repo holds the source code and scripts for reproduce the key experiments of fast k-means evaluation.
We also upload an exemplar dataset that you can play with in the folder "dataset".

Download our technical report here: https://github.com/tgbnhy/fast-kmeans/unik-vldb.pdf

## Usage


1. If you run in Eclipse, just go to "edu.nyu.unik.expriments.kmeansEfficiency", and click the "run configuration", creat a new java application, and fill the following parameters:

```
./dataset/europediff_169300_2.txt 10 169300 a euro 0 1
```
There are six parameters:
```
arg[0] is the dataset file that you want to cluster
arg[1] is the number of clusters (k)
arg[2] is the number of points in the datafile which will be clustered (|D|)
arg[3] is the edge info file which contains the street name
arg[4] is the name of dataset to distinguish
arg[5] is the start dimension (column) in the dataset file
arg[6] is the end dimension (column) in the dataset file
```
Then, all the result will be recorded into the log file under the "logs" folder.

2. If you want to run from commands (recommended):

```
mvn clean package
```
A file "torch-clus-0.0.1-SNAPSHOT.jar" will be generated under folder "target".
```
 java -Xmx16192M -cp ./target/torch-clus-0.0.1-SNAPSHOT.jar edu.nyu.unik.expriments.kmeansEfficiency ./dataset/europediff_169300_2.txt 10 169300 a euro 0 1
```

## Compared algorithms
| __Algorithm__ | __Paper__ | __Year__ |
|-------------|------------|------------|
| Lloyd         |   Least squares quantization in PCM   | 1987      |
|Ball-tree|The Anchors Hierarchy: Using the Triangle Inequality to Survive High Dimensional Data|2000|
|kd-tree|An efficient k-means clustering algorithm: Analysis and implementation|2002|
| Elkan         | Using the triangle inequality to accelerate k-means | 2003     |
|Hamerly|Making k-means even faster|2010|
|Drake|Accelerated k-means with adaptive distance bounds|2012|
|Annulus|Faster k-means Clustering (Master thesis)|2013|
|Search|Scalable K-Means by ranked retrieval|2014|
|Yinyang|Yinyang K-Means: A Drop-In Replacement of the Classic K-Means with Consistent Speedup|2015|
|Heap|Accelerating Lloyd’s Algorithm for k-Means Clustering|2015|
|Expo|Fast k-means with accurate bounds|2016|
|Drift|Geometric methods to accelerate k -means algorithms|2016|
|Vector|Speeding up k -means by approximating Euclidean distances via block vectors|2016|
|Regroup|Two Modifications of Yinyang K-means Algorithm|2017|
|Cover-tree|A Dual-Tree Algorithm for Fast k-means Clustering With Large k|2017|


## Datasets
Any dataset with csv format can be clustered, just need to specify the columns that have numeric values.
For example, you can download some datasets from UCI (https://archive.ics.uci.edu/ml/datasets.php), and put into the "dataset" folder.

| __Dataset__ | __link__ | __Dimension__ |
|-------------|------------|------------|
| BigCross         | https://s3.amazonaws.com/h2o-training/clustering/BigCross.data.gz     | 57      |
| Conflong         | http://networkrepository.com/ConfLongDemo-JSI.php | 3     |
|Covtype|https://archive.ics.uci.edu/ml/datasets/covertype|55|
|Europe|http://cs.joensuu.fi/sipu/datasets/europediff.txt|2|
|KeggDirect|https://archive.ics.uci.edu/ml/datasets/KEGG+Metabolic+Relation+Network+(Directed) |24|
|KeggUndirect|https://archive.ics.uci.edu/ml/datasets/KEGG+Metabolic+Reaction+Network+(Undirected) |29|
|NYC-Taxi||2|
|Skin|https://archive.ics.uci.edu/ml/datasets/skin+segmentation |4|
|Power|https://archive.ics.uci.edu/ml/datasets/Individual+household+electric+power+consumption |9|
|RoadNetwork|https://archive.ics.uci.edu/ml/datasets/3D+Road+Network+(North+Jutland,+Denmark)|4|
|US-Census|https://archive.ics.uci.edu/ml/machine-learning-databases/census1990-mld/ |68|
|Mnist|http://yann.lecun.com/exdb/mnist/|784|

## Results Interpretation
### Log files
In the log file, according to the name, you can go to check multiple metrics we use.
For example, in "logs/vldb_log1/euro_169300_2_10_BallMetric_30.log", you can see multiple lines.
In each line, we present multiple methods' corresponding performance.
The first line shows the overall running time, for the other lines, please check class "kmeansAlgorihtm" line 1667, function writelogs to see the metric if you are intereted in.
For all the methods we tested, you can check line 2184 function testExisting and line 2272 function testIndex to see the full names.

### Terminal
After you run the command, you will also observe logs from the terminal. It mainly shows the running time of each iteration calling various methods.


## Citation
If you use our code for research work, please cite our paper as below:
```
@article{wang2019unik,
  title={On the Efficiency of {K}-Means Clustering: Evaluation, Optimization, and Algorithm Selection},
  author={Wang, Sheng and Bao, Zhifeng and Sun, Yuan},
  year={2020},
}
```
