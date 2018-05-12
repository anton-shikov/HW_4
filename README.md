# HW_4. De Brujin Graph visualizer. 
A simple phyton program for genome assemble.

## Getting Started

This tool allows you to assemble genomes and visualize De Bruin Graph in different types and formats.

### Prerequisites

You need to install python3 with Biopython and Graphvis modules to run this tool.

### Installing

To install this tool clone this repository to your PC.

```
~$ git clone https://github.com/anton-shikov/HW_4
```

## Running and using tool

Using this tool is not complicted. After downloading repository launch terminal and enter this repository or run programm using absolute path to .py file.
Use this following command to execute tool:
```
~$ python De_Brujin.py -sq hw_4_5_dataset.fasta -ks 55 -g a -d p -gi pd -o hw_4_5_dataset_graph -O HW4-5_contigs.fasta -i 5 -t 2
```
Information about flags:  
-sq input sequence in fasta format  
-ks the size of k-mer (3 for default)  
-g graph visualization type: f for full; a for abridged (f for default)  
-d DNA-strand for analysis: p for (+)-strand; m for (-)-strand (p for default)  
-gi graph output format: pn for png; pd for pdf; sv for svg (pd for default)  
-o graph output name  
-O contigs output name
-i number of iteration for graph compressing\cutting tails  
-t coverage threshold for cutting tails  


Output format: Graph in dot format and in png/pdf/svg format, contigs in fasta format.

By the way, blasting contigs has shown that reads for homework were from Hepatitis B virus

## Author

* **Anton Shikov** - *Initial work* - [anton-shikov](https://github.com/anton-shikov)


## License

This project is free and available for everyone.

## Acknowledgments

Eugene Bakin for python course.
