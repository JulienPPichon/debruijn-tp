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
import sys
import networkx as nx
import matplotlib
import random
random.seed(9001)
import statistics

__author__ = "Your Name"
__copyright__ = "Universite Paris Diderot"
__credits__ = ["Your Name"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Your Name"
__email__ = "your@email.fr"
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
	with open(fastq_file, "r") as filin:
		for line in filin:			
			yield next(filin).strip()
			next(filin)
			next(filin)
			

def cut_kmer(seq, k):
	for i in range(len(seq) - k + 1):
		yield seq[i:i+k]


def build_kmer_dict(fastq_file, k):
	kmer_dict = {}
	for seq in read_fastq(fastq_file):
		for kmer in cut_kmer(seq, k):
			if kmer not in kmer_dict:
				kmer_dict[kmer] = 1
			else:
				kmer_dict[kmer] += 1
	return kmer_dict


def build_graph(kmer_dict):
	kmer_tree = nx.DiGraph()
	for kmer in kmer_dict:
		kmer_tree.add_edge(kmer[0:len(kmer) - 1], kmer[1:len(kmer)], weight = kmer_dict[kmer])
	return kmer_tree
		

def get_starting_nodes(kmer_tree):
	starting_nodes = []
	for node in kmer_tree.nodes:
		if len(list(kmer_tree.predecessors(node))) == 0:
			starting_nodes.append(node)
	return starting_nodes


def get_sink_nodes(kmer_tree):
	sink_nodes = []
	for node in kmer_tree.nodes:
		if len(list(kmer_tree.successors(node))) == 0:
			sink_nodes.append(node)
	return sink_nodes
	

def get_contigs(kmer_tree, starting_nodes, sink_nodes):
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
	with open("../results/" + fasta, "w") as filout:
		for number, contig in enumerate(contigs):
			filout.write(">contig_" +  str(number + 1) + " len=" + str(contig[1]) + "\n" + str(fill(contig[0])) + "\n")
			

def std(value_list):
	std = statistics.stdev(value_list)
	return std
	
	
def path_average_weight(kmer_tree, path):
	total_weight = 0
	for kmer in path:
		total_weight += kmer_tree.out_degree(kmer, weight = "weight")
	average_weight = total_weight / (len(path) - 1)
	return average_weight
	

def remove_paths(kmer_tree, path_list, delete_entry_node, delete_sink_node):
	entry = 1
	sink = -1
	if delete_entry_node == True:
		entry = 0
	if delete_sink_node == True:
		sink = None
	for path in path_list:
		kmer_tree.remove_nodes_from(path[entry:sink]) 
	return kmer_tree
    
    
def select_best_path(kmer_tree, path_list, length_list, weight_list, delete_entry_node = False, delete_sink_node = False):
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
	kmer_tree = remove_paths(kmer_tree, path_list[:best_path_index] + list_path[best_path_index + 1:], delete_entry_node, delete_sink_node)
	return kmer_tree


def solve_bubble(kmer_tree, ancestor_node, descendant_node):
	bubble_path = []
	bubble_len_path = []
	bubble_weight = []
	for path in nx.all_simple_paths(kmer_tree, ancestor_node, descendant_node):
		bubble_path.append(path)
		bubble_len_path.append(len(path))
		bubble_weight.append(path_average_weight(kmer_tree, path))
	kmer_tree = select_best_path(kmer_tree, bubble_path, bubble_len_path, bubble_weight)
	return kmer_tree

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
	starting_nodes = get_starting_nodes(kmer_tree)
	sink_nodes = get_sink_nodes(kmer_tree)
	all_contigs = get_contigs(kmer_tree, starting_nodes, sink_nodes)
	save_contigs(all_contigs, args.output_file) 


if __name__ == '__main__':
	main()
