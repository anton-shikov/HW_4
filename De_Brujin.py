# -*- coding: utf-8 -*-
"""
Created on Fri Apr 27 16:06:51 2018

@author: anton
"""

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
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
        
    def merge (self, vertex, prev_vertex, after_vertex):
        self.prev_vertex,self.after_vertex =  prev_vertex,after_vertex
        self.vertices[self.prev_vertex].out_edges[self.after_vertex], self.vertices[self.after_vertex].in_edges[self.prev_vertex] = [Edge(self.prev_vertex, self.after_vertex)],[Edge(self.prev_vertex, self.after_vertex)]
        
        self.new_seq = self.vertices[vertex].in_edges[self.prev_vertex][0].seq + self.vertices[self.after_vertex].in_edges[vertex][0].seq[self.k:]
        self.vertices[self.prev_vertex].out_edges[self.after_vertex][0].seq, self.vertices[self.after_vertex].in_edges[self.prev_vertex][0].seq = self.new_seq,self.new_seq
        
        self.new_cov = (self.vertices[vertex].out_edges[self.after_vertex][0].coverage * (len(self.vertices[vertex].out_edges[self.after_vertex][0].seq) - self.k + 1) + self.vertices[vertex].in_edges[self.prev_vertex][0].coverage * (len(self.vertices[vertex].in_edges[self.prev_vertex][0].seq) - self.k + 1) - self.vertices[vertex].coverage) / (len(self.vertices[self.prev_vertex].out_edges[self.after_vertex][0].seq) - self.k + 1)
        self.vertices[self.prev_vertex].out_edges[self.after_vertex][0].coverage, self.vertices[self.after_vertex].in_edges[self.prev_vertex][0].coverage = self.new_cov, self.new_cov 
        
        del self.vertices[vertex]
        del self.vertices[self.prev_vertex].out_edges[vertex]
        del self.vertices[self.after_vertex].in_edges[vertex]
       
    def compress(self):
        self.to_delete = [vertex for vertex in self.vertices.keys() if (len(self.vertices[vertex].out_edges.keys()) == 1) and (len(self.vertices[vertex].in_edges.keys()) == 1)]
        for vertex in self.to_delete: 
            prev_vertex, after_vertex = list(self.vertices[vertex].in_edges.keys())[0],list(self.vertices[vertex].out_edges.keys())[0]
            self.merge(vertex, prev_vertex, after_vertex)
    def cut_tails(self, lim):
        self.to_cut = [vertex for vertex in self.vertices.keys() if ((len(self.vertices[vertex].out_edges.keys()) == 0) and (len(self.vertices[vertex].in_edges.keys()) == 1)) or ((len(self.vertices[vertex].out_edges.keys()) == 1) and (len(self.vertices[vertex].in_edges.keys()) == 0)) or ((len(self.vertices[vertex].out_edges.keys()) == 0) and (len(self.vertices[vertex].in_edges.keys()) == 0))]
        for vertex in self.to_cut: 
            if self.vertices[vertex].coverage <= lim:
                if len(self.vertices[vertex].out_edges.keys()) == 1:
                    self.vertices[list(self.vertices[vertex].out_edges.keys())[0]].in_edges.pop(vertex)
                if len(self.vertices[vertex].in_edges.keys()) == 1:
                    self.vertices[list(self.vertices[vertex].in_edges.keys())[0]].out_edges.pop(vertex)
                del self.vertices[vertex]    
                
    def get_contigs(self, out):
        i=0
        contigs = []
        for vertex in self.vertices.keys():
            for out_vertex in self.vertices[vertex].out_edges.keys():
                contigs.append(SeqRecord(Seq(self.vertices[vertex].out_edges[out_vertex][0].seq), id="contig {}".format(i)))
                i+=1
        SeqIO.write(contigs, out, "fasta")

        
               
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='De_Bruijn_Graph visualization')
    parser.add_argument('-sq', '--sequence', help='Please enter sequence file', type=str)
    parser.add_argument('-ks', '--kmersize', help='Please enter desirable k-mer size', default=3, type=int)
    parser.add_argument('-g', '--graphtype', help='Please enter graph visualization type: f for full; a for abridged', default='f', type=str)
    parser.add_argument('-d', '--strand', help='Please enter DNA strand: p for (+)-strand; m for (-)-strand', default='p', type=str)
    parser.add_argument('-gi', '--graphimage', help='Please enter graph output: pn for png; pd for pdf; sv for svg', default='pd', type=str)
    parser.add_argument('-o', '--graphout', help='Please enter output file name', type=str)
    parser.add_argument('-O', '--contigsout', help='Please enter contigs file name', type=str)
    parser.add_argument('-i', '--iterations', help='Please enter number of iterations for graph compression', type=int)  
    parser.add_argument('-t', '--threshold', help='Please enter coverage theshold for cutting', type=int)  
       
    args = parser.parse_args()
    
    sq, ks, g, d, gi, o, O, t, i  = args.sequence, args.kmersize, args.graphtype, args.strand,  args.graphimage, args.graphout, args.contigsout,  args.threshold, args.iterations

    my_graph = Graph(ks)

    if d=='p':
        with open(sq, "r") as handle:
            for record in SeqIO.parse(handle, "fasta"):
                my_graph.add_read(str(record.seq))
    else:
        with open(sq, "r") as handle:
            for record in SeqIO.parse(handle, "fasta"):
                my_graph.add_read(str(record.reverse_complement().seq) )
    
    if i is not None:
        for k in range (i):
            my_graph.compress()
            my_graph.cut_tails(t)
        
    my_graph.calc_init_edge_coverage()  
    if O is not None:  
        my_graph.get_contigs(O)
    if o is not None and g is not None and gi is not None:
        my_graph.graph_vis(g, gi, o)
