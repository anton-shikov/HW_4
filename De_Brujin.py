# -*- coding: utf-8 -*-
"""
Created on Fri Apr 27 16:06:51 2018

@author: anton
"""

from Bio import SeqIO
from collections import defaultdict
from graphviz import Digraph
import argparse

class Vertex:
    
    def __init__(self, seq):
        self.seq = seq
        self.coverage = 1
        self.in_edges = defaultdict()
        self.out_edges = defaultdict()
        
    def increase_coverage(self):
        self.coverage += 1

class Edge:
    
    def __init__(self,k1,k2):
        self.seq = k1 + k2[-1]
        self.coverage = 0
    
    def calc_coverage(self,c1,c2):
        self.coverage = (c1+c2)/2


class Graph:

    def __init__(self,k):
        self.vertices = defaultdict()
        self.k = k
        
    def add_read(self,read):
        read_lng = len(read)
        if read_lng < self.k:
            return
            
        kmer = read[:self.k]
        if kmer in self.vertices:
            self.vertices[kmer].increase_coverage()
        else:
            self.vertices[kmer] = Vertex(kmer)
        
        for next_kmer_indx in range(1,read_lng-self.k+1,1):
            next_kmer = read[next_kmer_indx:(next_kmer_indx+self.k)]
            if next_kmer in self.vertices:
                self.vertices[next_kmer].increase_coverage()
            else:
                self.vertices[next_kmer] = Vertex(next_kmer)
            
            new_edge = Edge(kmer,next_kmer)
            
            self.vertices[next_kmer].in_edges[kmer]  = [new_edge]
            
            self.vertices[kmer].out_edges[next_kmer] = [new_edge]

            kmer = next_kmer
    
    def calc_init_edge_coverage(self):
        
        for current_vertex in self.vertices.keys():
            for next_vertex in self.vertices[current_vertex].out_edges.keys():
                self.vertices[current_vertex].out_edges[next_vertex][0].calc_coverage(self.vertices[current_vertex].coverage,self.vertices[next_vertex].coverage)
    
    def graph_vis(self, gtype, form, output):
        
        dot = Digraph(comment='De_Bruijn_Graph')
        if gtype == 'f':
            for kmer, node_obects in self.vertices.items():
                dot.node(kmer, kmer)
                for out_kmers, out_obects in node_obects.out_edges.items():
                    dot.edge(kmer, out_kmers, '{} coverage={}'.format(out_obects[0].seq, out_obects[0].coverage))
        else:
            for kmer, node_obects in self.vertices.items():
                dot.node(kmer, 'coverage={}'.format(node_obects.coverage))
                for out_kmers, out_obects in node_obects.out_edges.items():
                    dot.edge(kmer, out_kmers, 'coverage={} length={}'.format(out_obects[0].coverage,len(out_obects[0].seq)))
        if form == 'pn':
            dot.format = 'png'
        elif form == 'pd':
            dot.format = 'pdf'
        elif form == 'sv':
            dot.format = 'svg'
        dot.render('{}'.format(output), view=True)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='De_Bruijn_Graph visualization')
    parser.add_argument('-sq', '--sequence', help='Please enter sequence file', type=str)
    parser.add_argument('-ks', '--kmersize', help='Please enter desirable k-mer size', default=3, type=int)
    parser.add_argument('-g', '--graphtype', help='Please enter graph visualization type: f for full; a for abridged', default='f', type=str)
    parser.add_argument('-d', '--strand', help='Please enter DNA strand: p for (+)-strand; m for (-)-strand', default='p', type=str)
    parser.add_argument('-gi', '--graphimage', help='Please enter graph output: pn for png; pd for pdf; sv for svg', default='pd', type=str)
    parser.add_argument('-o', '--graphout', help='Please enter output file name', type=str)
    
    args = parser.parse_args()
    
    sq, ks, g, d, gi, o  = args.sequence, args.kmersize, args.graphtype, args.strand,  args.graphimage, args.graphout

    my_graph = Graph(ks)
    
    if d=='p':
        with open(sq, "r") as handle:
            for record in SeqIO.parse(handle, "fasta"):
                my_graph.add_read(str(record.seq))
    else:
        with open(sq, "r") as handle:
            for record in SeqIO.parse(handle, "fasta"):
                my_graph.add_read(str(record.reverse_complement().seq) )
    my_graph.calc_init_edge_coverage()         
    my_graph.graph_vis(g, gi, o)