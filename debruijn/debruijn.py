#!/bin/env python3
# -*- coding: utf-8 -*-
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    A copy of the GNU General Public License is available at
#    http://www.gnu.org/licenses/gpl-3.0.html

"""Perform assembly based on debruijn graph."""

from operator import itemgetter
from random import randint
import argparse
import os
import statistics
import sys
import networkx as nx
import matplotlib
import random as rd
rd.seed(9001)


__author__ = "Julien Pichon"
__copyright__ = "Universite Paris Diderot"
__credits__ = ["Julien Pichon"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Julien Pichon"
__email__ = "julien.pichon@@cri-paris.org"
__status__ = "Developpement"

def isfile(path):
    """Check if path is an existing file.
        :Parameters:
        path: Path to the file
    """
    if not os.path.isfile(path):
        if os.path.isdir(path):
            msg = "{0} is a directory".format(path)
        else:
            msg = "{0} does not exist.".format(path)
        raise argparse.ArgumentTypeError(msg)
    return path


def get_arguments():
    """Retrieves the arguments of the program.
        Returns: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage=
                                    "{0} -h"
                                    .format(sys.argv[0]))
    parser.add_argument('-i', dest='fastq_file', type=isfile,
                        required=True, help="Fastq file")
    parser.add_argument('-k', dest='kmer_size', type=int,
                        default=21, help="K-mer size (default 21)")
    parser.add_argument('-o', dest='output_file', type=str,
                        default=os.curdir + os.sep + "contigs.fasta",
                        help="Output contigs in fasta file")
    return parser.parse_args()


def read_fastq(fastq_file):
    """Take a fastq file as argument.
        Return a generator of sequences
    """
    with open(fastq_file, "r") as filin:
        for line in filin:
            yield next(filin).strip()
            next(filin)
            next(filin)


def cut_kmer(seq, k):
    """Take a sequence and an integer as argument.
        Return a generator of sequence chunks of length = k.
    """
    for i in range(len(seq) - k + 1):
        yield seq[i:i+k]


def build_kmer_dict(fastq_file, k):
    """Take a fastq file and an interger as argument.
        Return a dictionnary of each kmer as key and
        its occurence as vallue.
    """
    kmer_dict = {}
    for seq in read_fastq(fastq_file):
        for kmer in cut_kmer(seq, k):
            if kmer not in kmer_dict:
                kmer_dict[kmer] = 1
            else:
                kmer_dict[kmer] += 1
    return kmer_dict


def build_graph(kmer_dict):
    """Take a dictionnary as input with kmer as key and
        occurences as value.
        Return a directionnal graph with kmer as nodes and occurences
        as length
    """
    kmer_tree = nx.DiGraph()
    for kmer in kmer_dict:
        kmer_tree.add_edge(kmer[0:len(kmer) - 1], kmer[1:len(kmer)], weight = kmer_dict[kmer])
    return kmer_tree


def get_starting_nodes(kmer_tree):
    """Take a graph as argument and return the
        starting nodes list.
    """
    starting_nodes = []
    for node in kmer_tree.nodes:
        if len(list(kmer_tree.predecessors(node))) == 0:
            starting_nodes.append(node)
    return starting_nodes


def get_sink_nodes(kmer_tree):
    """Take a graph as argument and return the
        ending nodes list.
    """
    sink_nodes = []
    for node in kmer_tree.nodes:
        if len(list(kmer_tree.successors(node))) == 0:
            sink_nodes.append(node)
    return sink_nodes


def get_contigs(kmer_tree, starting_nodes, sink_nodes):
    """Take as argument a graph, the lists of starting and ending nodes.
        Find the different paths and build a contig, keeping only the last letter
        of each kmer.
        Return a list of tuple with the sequence and its length.
    """
    all_contigs = []
    for input_node in starting_nodes:
        for output_node in sink_nodes:
            for path in nx.all_simple_paths(kmer_tree, input_node, output_node):
                contig = path[0]
                for kmer in path[1:]:
                    contig = contig + kmer[-1]
                all_contigs.append((contig, len(contig)))
    return all_contigs


def fill(text, width=80):
    """Split text with a line return to respect fasta format"""
    return (os.linesep.join(text[i:i+width] for i in range(0, len(text), width)))


def save_contigs(contigs, fasta):
    """Take a list of contigs and a output file name as arguments.
        Create a fasta file in the local depositery.
    """
    with open(fasta, "w") as filout:
        for number, contig in enumerate(contigs):
            filout.write(">contig_" +  str(number + 1) + " len=" + str(contig[1])
            + "\n" + str(fill(contig[0])) + "\n")


def std(value_list):
    """calculate standart deviation of a list."""
    standart_deviation = statistics.stdev(value_list)
    return standart_deviation


def path_average_weight(kmer_tree, path):
    """Take a graph and a path as argument.
    Caculate and return the average weigth of the path.
    """
    total_weight = 0
    for kmer in path:
        total_weight += kmer_tree.out_degree(kmer, weight = "weight")
    average_weight = total_weight / (len(path) - 1)
    return average_weight


def remove_paths(kmer_tree, path_list, delete_entry_node, delete_sink_node):
    """Take a graph, a list of paths, the boleans deletion or not of entry and
    sink nodes as arguments.
    Return a graph without the unwanted paths.
    """
    entry = 1
    sink = -1
    if delete_entry_node == True:
        entry = 0
    if delete_sink_node == True:
        sink = None
    for path in path_list:
        kmer_tree.remove_nodes_from(path[entry:sink])
    return kmer_tree


def select_best_path(kmer_tree, path_list, length_list, weight_list,
    delete_entry_node = False, delete_sink_node = False):
    """Take a graph, a list of paths, a list of length paths, a list of average length
        of each path and the boleans deletion or not of entry and sink nodes as
        arguments.
        Return a graph without the unwanted paths.
    """
    max_weight = 0
    best_path_len = 0
    best_path_index = -1
    for list_index, weight in enumerate(weight_list):
        if weight > max_weight:
            max_weight = weight
            best_path_len = length_list[list_index]
            best_path_index = list_index
        elif weight == max_weight:
            if best_path_len < length_list[list_index]:
                best_path_len = length_list[list_index]
                best_path_index = list_index
            elif best_path_len == length_list[list_index]:
                best_path_index = rd.choice([best_path_index, list_index])
    if best_path_index == -1:
        best_path_index = rd.randint(0, len(path_list))
    kmer_tree = remove_paths(kmer_tree, path_list[:best_path_index]
    + path_list[best_path_index + 1:], delete_entry_node, delete_sink_node)
    return kmer_tree


def solve_bubble(kmer_tree, ancestor_node, descendant_node):
    """Take a graph, and ancestor and descendant node as arguments.
        Return a graph without the bubble between the specified nodes.
    """
    bubble_path = []
    bubble_len_path = []
    bubble_weight = []
    for path in nx.all_simple_paths(kmer_tree, ancestor_node, descendant_node):
        bubble_path.append(path)
        bubble_len_path.append(len(path))
        bubble_weight.append(path_average_weight(kmer_tree, path))
    kmer_tree = select_best_path(kmer_tree, bubble_path, bubble_len_path, bubble_weight)
    return kmer_tree


def simplify_bubbles(kmer_tree):
    """Take a graph as argument and return a graph without bubbles."""
    bubble_nodes = []
    for node in kmer_tree:
        ancestor_node = [i for i in kmer_tree.predecessors(node)]
        if len(ancestor_node) >= 2:
            ancestor = nx.lowest_common_ancestor(kmer_tree, ancestor_node[0], ancestor_node[1])
            bubble_nodes.append([ancestor, node])
    for node_couples in bubble_nodes:
        kmer_tree = solve_bubble(kmer_tree, node_couples[0], node_couples[1])
    return kmer_tree


def solve_entry_tips():
    pass


def solve_out_tips():
    pass


#==============================================================
# Main program
#==============================================================
def main():
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()
    kmer_dict = build_kmer_dict(args.fastq_file, args.kmer_size)
    kmer_tree = build_graph(kmer_dict)
    kmer_tree = simplify_bubbles(kmer_tree)
    starting_nodes = get_starting_nodes(kmer_tree)
    sink_nodes = get_sink_nodes(kmer_tree)
    all_contigs = get_contigs(kmer_tree, starting_nodes, sink_nodes)
    save_contigs(all_contigs, args.output_file)


if __name__ == '__main__':
    main()
