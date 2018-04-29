# HW_4. De Brujin Graph visualizer. 
A simple phyton program for genome assemble.

## Getting Started

This tool allows you to assemble genomes and visualize De Bruin Graph in different ways and formats

### Prerequisites

You need to install python3 with Biopython and Graphvis modules to run this tool.

### Installing

To install this tool clone this repository to your PC.

```
~$ git clone https://github.com/anton-shikov/HW_4
```

## Running and using tool

Using this tool is not complicted. After downloading repository launch terminal and enter this repository or run programm using absolute path to .py file.
Use this following to execute tool:
```
~$  python De_Brujin.py -sq test.fasta -ks 15 -g f -d p -gi pd -o testgraph 
```
Information about flags:  
-sq input sequence in fasta format  
-ks the size of k-mer (3 for default)
-g graph visualization type: f for full; a for abridged (f for default) 
-d DNA-strand for analysis: p for (+)-strand; m for (-)-strand (p for default)
-gi graph output format: pn for png; pd for pdf; sv for svg (pd for default)
-o output name

Output format: Graph in dot format amd in png/pdf/svg format.

## Author

* **Anton Shikov** - *Initial work* - [anton-shikov](https://github.com/anton-shikov)


## License

This project is free and available for everyone.

## Acknowledgments

Eugene Bakin for python course.
