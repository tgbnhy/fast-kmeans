# Evaluation of Fast k-means

## Introduction
This repo holds the source code and scripts for reproduce the key experiments of fast k-means evaluation.

## Usage


1. If you run in Eclipse, just go to "au.edu.rmit.trajectory.expriments.kmeansEfficiency", and click the "run configuration", creat a new java application, and fill the following parameters:

```
./dataset/NYC_TAXI_Location_2013_Oct_15004557_6.csv 10 3500455 a NYC  8 9
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

 #run the tdrive clustering for efficiency comparision.
```
 java -Xmx16192M -cp ./torch-clus-0.0.1-SNAPSHOT.jar ./dataset/NYC_TAXI_Location_2013_Oct_15004557_6.csv 10 3500455 a NYC  8 9
```

## Datasets
Any dataset with csv format can be clustered, just need to specify the columns that have numeric values.
For example, you can download and put into the "dataset folder"

## Results Interpretation
In the log file, according to the name, you can go to check multiple metric we use.


## Citation
If you use our code for research work, please cite our paper as below:
```
@article{wang2019fast,
  title={Fast large-scale trajectory clustering},
  author={Wang, Sheng and Bao, Zhifeng and Culpepper, J Shane and Sellis, Timos and Qin, Xiaolin},
  journal={Proceedings of the VLDB Endowment},
  volume={13},
  number={1},
  pages={29--42},
  year={2019},
  publisher={VLDB Endowment}
}
```
