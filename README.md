# Approximate Cluster Identities (ACI)

A python package to visualise the approximate within and between cluster identities of short sequences as assigned by e.g. mmseqs2, cd-hit or panaroo.

# Installation
```
pip install approximate-cluster-identities
```

# Usage

```
aci -h

Create visualisations of approximate between and within cluster nucleotide identities for short sequences.

positional arguments:
  input_fasta           Input FASTA file of all sequences.
  input_json            Input JSON file with cluster assignments ({<sequence header>: <cluster assignment>}).

optional arguments:
  -h, --help            show this help message and exit
  --clusterGML CLUSTERGML
                        Output path of GML clustering file to view with Cytoscape or similar.
  --distanceTable DISTANCETABLE
                        Output path of CSV of distances (may take a long time).
  --clusterPlot CLUSTERPLOT
                        Output path of jointplot to visualise between and within cluster identities.
  --kmerSize KMERSIZE   Kmer size (default: 9).
  --windowSize WINDOWSIZE
                        Minimiser window size (default: 20).
  --threshold THRESHOLD
                        Jaccard similarity threshold (default: 0.9).
  --threads THREADS     Threads for sketching and jaccard distance calculations (default: 1).
  --shorter             Assess identity relative to the shorter sequence.
```

# Methods

We calculate sequence identities by pairwise calculation of jaccard distances using minimizers of size ```--kmerSize``` where 1 *k*-mer is sampled from a window of a total of ```--windowSize``` *k*-mers that slides across the input sequences. Increasing ```--windowSize``` will decrease the number of minimizers per sequence, decreasing the sensitivity of the identity calculations but increasing the speed of the programme. This tool is designed to give you an idea of how variable a large number of short sequences are within and between clusters to choose an appropriate sequencing clustering tool and its parameters.

# Example output

Example cluster plots for data in ```test/``` using ```--windowSize 1``` and ```--windowSize 100```.

### Window size = 1

#### Mean identities
![Mean identities](images/cluster_distances_window_1.mean.png)
#### Mode identities
![Mode identities](images/cluster_distances_window_1.mode.png)
#### Median identities
![Median identities](images/cluster_distances_window_1.median.png)
#### Range identities
![Range of identities](images/cluster_distances_window_1.range.png)

### Window size = 100

#### Mean identities
![Mean identities](images/cluster_distances_window_100.mean.png)
#### Mode identities
![Mode identities](images/cluster_distances_window_100.mode.png)
#### Median identities
![Median identities](images/cluster_distances_window_100.median.png)
#### Range identities
![Range of identities](images/cluster_distances_window_100.range.png)
