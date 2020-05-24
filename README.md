# Evaluation of Fast k-means

## Introduction
This repo holds the source code and scripts for reproduce the key experiments of fast k-means evaluation.
We also upload an exemplar dataset that you can play with in the folder "dataset"

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
| Cat         | Soccer     | Apple      |
| Dog         | Basketball | Orange     |


## Datasets
Any dataset with csv format can be clustered, just need to specify the columns that have numeric values.
For example, you can download some datasets from UCI (https://archive.ics.uci.edu/ml/datasets.php), and put into the "dataset" folder.

| __Dataset__ | __link__ | __Fruits__ |
|-------------|------------|------------|
| Cat         | Soccer     | Apple      |
| Dog         | Basketball | Orange     |


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
