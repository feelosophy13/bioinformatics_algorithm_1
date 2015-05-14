import numpy as np
import collections  # for Counter()
import operator  # for itemgetter()

from module1 import *
from module2 import *


## this function finds DNA patterns that is off by one nucleotide at most
def get_immediate_neighbors(ptrn):
    nucleotides = ['A', 'C', 'G', 'T']
    ptrn_len = len(ptrn)
    neighbors = [ptrn]

    ## for each observed nucleotide in k-mer
    for i in range(0, ptrn_len):
        nucleotide = ptrn[i]
        
        ## for each nucleotide that is NOT the observed nucleotide
        for nuc in list(set(nucleotides) - set(nucleotide)):
            
            ## generate a 1-mismatch pattern by replacing the observed nucleotide with a mismatch
            neighbor = list(ptrn)
            neighbor[i] = nuc
            neighbor = ''.join(neighbor)

            ## add to a list of neighbors
            neighbors.append(neighbor)

    return neighbors


## this function is equivalent to get_aprx_match_ptrns function but possibly more efficient (because it doesn't list out all the permutations)
## this function finds DNA patterns with at most d mismatches
def get_neighbors(ptrn, d):
    if d == 0:
        return [ptrn]
    if len(ptrn) == 1:
        return ['A', 'C', 'G', 'T']
    neighbors = []
    suffix = ptrn[1::]
    suffix_neighbors = get_neighbors(suffix, d)
    for suffix_neighbor in suffix_neighbors:
        Hamming_dist = calc_Hamming_dist(suffix, suffix_neighbor)
        if Hamming_dist < d:
            for nucleotide in ['A', 'C', 'G', 'T']:
                neighbors.append(nucleotide + suffix_neighbor)
        else:
            neighbors.append(ptrn[0] + suffix_neighbor)

    return neighbors

    
## this function returns a list of all k-mers seen in a single DNA string
def find_all_kmers_in_DNA_str(DNA_str, k):
    n_iter = len(DNA_str) - k + 1
    kmers = []
    for i in range(0, n_iter):
        kmer = DNA_str[i:i+k]
        kmers.append(kmer)
    return kmers


## this function returns a list of all k-mers seen in a list of DNA strings
def find_all_kmers_in_DNA_list(DNA_list, k):
    kmers = []
    for DNA_str in DNA_list:
        kmers_per_DNA_str = find_all_kmers_in_DNA_str(DNA_str, k)
        kmers += kmers_per_DNA_str
    return kmers 


## this function returns True if a kmer is in a DNA string with most d mismatches
def check_kmer_exists_in_DNA_str(kmer, DNA_str, d):
    k = len(kmer)
    n_iter = len(DNA_str) - k + 1
    for i in range(0, n_iter):
        kmer_window = DNA_str[i:i+k]
        Hamming_dist = calc_Hamming_dist(kmer, kmer_window)
        if Hamming_dist <= d:
            return True
    return False


## this function returns True if a kmer is in all of the strings in a list of DNA strings, each with most d mismatches
def check_kmer_exists_in_all_DNA_str(kmer, DNA_list, d):
    for DNA_str in DNA_list:
        kmer_exists_in_DNA_str = check_kmer_exists_in_DNA_str(kmer, DNA_str, d)
        if not kmer_exists_in_DNA_str:
            return False
    return True

## this function takes in a list of DNA strings and returns a list of unique motifs that are different from each other with at most d mismatches (brute force approach)
def find_motifs(DNA_list, k, d):
    ptrns = []
    kmers_in_DNA_list = find_all_kmers_in_DNA_list(DNA_list, k)
    for kmer in kmers_in_DNA_list:
        kmer_neighbors = get_neighbors(kmer, d)
        for kmer_neighbor in kmer_neighbors:
            kmer_exists_in_all_DNA = check_kmer_exists_in_all_DNA_str(kmer_neighbor, DNA_list, d)
            if kmer_exists_in_all_DNA:
                ptrns.append(kmer_neighbor)
    ptrns = list(set(ptrns))
    return ptrns


## this function calculates the lowest distance between a pattern and a DNA string
def calc_dist_btwn_ptrn_and_DNA_str(ptrn, DNA_str):
    k = len(ptrn)
    n_iter = len(DNA_str) - k + 1
    dist = float('Inf')  # initialized to inifinity
    for i in range(0, n_iter):
        ptrn_window = DNA_str[i:i+k]
        Hamming_dist = calc_Hamming_dist(ptrn, ptrn_window)
        if Hamming_dist < dist:
            dist = Hamming_dist
    return dist


## this function calculates the cumulative distance between a pattern and a list of DNA strings
def calc_cum_dist_btwn_ptrn_and_all_DNA_str(ptrn, DNA_list):
    cum_dist = 0
    for DNA_str in DNA_list:
        dist = calc_dist_btwn_ptrn_and_DNA_str(ptrn, DNA_str)
        cum_dist += dist
    return cum_dist


## this function returns a median string, a pattern with the smallest cumulative Hamming distance against a list of DNA strings, that has a length of k
## this function is extremely slow for large values for k and is the cause of bottleneck for Greedy Motif Search algorithm
## also, there could be multiple median strings; this function returns only one 
def find_median_str(DNA_list, k):
    min_cum_dist = float('Inf')
    n_iter = 4 ** k
    for i in range(0, n_iter):
        ptrn = conv_base10_num_to_DNA_str(i, k)
        cum_dist = calc_cum_dist_btwn_ptrn_and_all_DNA_str(ptrn, DNA_list)
        if cum_dist < min_cum_dist:
            min_cum_dist = cum_dist
            median_str = ptrn
    return median_str


## this function returns the probability a specific k-mer will form based on k-mer's profile matrix 
def get_kmer_prob_based_on_prof_mtrx(kmer, profile_matrix):
    k = len(kmer)
    prob = 1
    for i in range(0, k):
        nuc = kmer[i]
        if nuc=='A':
            row = profile_matrix[0]
        elif nuc=='C':
            row = profile_matrix[1]
        elif nuc=='G':
            row = profile_matrix[2]
        else:
            row = profile_matrix[3]
        prob *= row[i]
    return prob


## this function returns the most probable k-mer from a DNA string based on a k-mer's profile matrix
# 1. Get a list of all k-mers from given text. (n - k + 1)
# 2. Calculate probabilities of all k-mers in above list. 
# 3. Find a k-mer with maximum probability value.

# The function returns an empty string if all the kmers in the string have zero probability
def find_most_prob_kmer_based_on_prof_mtrx(DNA_str, k, profile_matrix, zero_prob_adj=True):
    k = len(profile_matrix[0])
    n_iter = len(DNA_str) - k + 1
    most_prob_kmer = ''
    max_prob = 0

    if zero_prob_adj:
        profile_matrix = np.array(profile_matrix) + 1
        profile_matrix = profile_matrix.tolist()
    for i in range(0, n_iter):
        kmer = DNA_str[i:i+k]
        kmer_prob = get_kmer_prob_based_on_prof_mtrx(kmer, profile_matrix)
        if kmer_prob > max_prob:
            max_prob = kmer_prob
            most_prob_kmer = kmer
    return most_prob_kmer


## this function converts nucleotide count dictionary to a list
## e.g. {'A': 0, 'C':2, 'G':1, 'T':3} to [0, 2, 1, 3]
## 0th value represents A-count, 1st value C-count, and so forth
def conv_nuc_cnt_dict_to_list(nuc_cnt_dict):
    nuc_to_row_index = {'A':0, 'C':1, 'G':2, 'T':3}
    nuc_cnt_list = [0] * 4
    for nuc in nuc_cnt_dict.keys():
        nuc_cnt = nuc_cnt_dict[nuc]
        ind = nuc_to_row_index[nuc]
        nuc_cnt_list[ind] = nuc_cnt
    return nuc_cnt_list


## this function converts nucleotide count list to dictionary (opposite for conv_nuc_cnt_dict_to_list function)
def conv_nuc_cnt_list_to_dict(nuc_cnt_list):
    nuc_cnt_dict = {}
    nuc_cnt_dict['A'] = nuc_cnt_list[0]
    nuc_cnt_dict['C'] = nuc_cnt_list[1]
    nuc_cnt_dict['G'] = nuc_cnt_list[2]
    nuc_cnt_dict['T'] = nuc_cnt_list[3]
    return nuc_cnt_dict


## this function calculates the sum of Hamming distances of all elements of motifs to another string
## this function assumes that all motifs and the compared string have the same length
def calc_cum_dist_btwn_motifs_and_str(motifs, str):
    cum_dist = 0
    for motif in motifs:
        dist = calc_dist_btwn_ptrn_and_DNA_str(motif, str)
        cum_dist += dist
    return cum_dist


## this function calculates the sum of Hamming distances of all elements of motifs to motifs' median string 
## this function assumes that all motifs have the same length
def calc_cum_dist_btwn_motifs_and_median_str(motifs):
    k = len(motifs[0])
    median_str = find_median_str(motifs, k)
    cum_dist = calc_cum_dist_btwn_motifs_and_str(motifs, median_str)
    return cum_dist


## this function calculates the sum of Hamming distances of all elements of motifs to motifs' consensus string
## this function assumes that all motifs have the same length
def calc_cum_dist_btwn_motifs_and_consensus(motifs):
    consensus_str = find_consensus_str_from_motifs(motifs)
    cum_dist = calc_cum_dist_btwn_motifs_and_str(motifs, consensus_str)
    return cum_dist


## this function takes in a list of motifs and outputs a consensus string
## this function assumes that all the motifs have the same length
def find_consensus_str_from_motifs(motifs):
    nuc_cnt_mtrx = create_nuc_cnt_mtrx_from_motifs(motifs)
    nuc_cnt_mtrx = np.array(nuc_cnt_mtrx)
    consensus_str = ''
    for col in nuc_cnt_mtrx.T:
        nuc_cnt_dict = conv_nuc_cnt_list_to_dict(col)
        dom_nuc = max(nuc_cnt_dict.iteritems(), key=operator.itemgetter(1))[0]  # dominant nucleotide is the one that appears with the highest frequency at the ith position out of all motifs
        consensus_str += dom_nuc    
    return consensus_str


## this function returns a nucleotide count matrix (2-dimensional list) from a list of motifs of the same length
## this function uses numpy
def create_nuc_cnt_mtrx_from_motifs(motifs):
    k = len(motifs[0])
    t = len(motifs)

    ## convert 2-dim list to 2-dim numpy array
    motifs = [list(motif) for motif in motifs]
    motifs = np.array(motifs)

    ## initialize (transpose of) nucleotide count matrix
    nuc_cnt_mtrx_T = np.zeros([k, 4])
    
    ## for each column in motifs matrix, count by nucleotide (A, C, T, G)
    for i in range(0, len(motifs.T)):
        col = motifs.T[i]
        nuc_cnt_dict = dict(collections.Counter(col))
        nuc_cnt_list = conv_nuc_cnt_dict_to_list(nuc_cnt_dict)
        nuc_cnt_mtrx_T[i] = nuc_cnt_list

    ## return 
    nuc_cnt_mtrx = nuc_cnt_mtrx_T.T
    nuc_cnt_mtrx = nuc_cnt_mtrx.tolist()
    return nuc_cnt_mtrx


## this function returns a profile matrix (2-dimensional list) from a list of motifs of the same length
## this function assumes that all motifs have the same length
def create_prof_mtrx_from_motifs(motifs, zero_prob_adj=False):

    ## calculate nucleotide count matrix
    nuc_cnt_mtrx = create_nuc_cnt_mtrx_from_motifs(motifs)
    nuc_cnt_mtrx = np.array(nuc_cnt_mtrx)
    if zero_prob_adj:
        nuc_cnt_mtrx += 1

    ## calculate profile matrix
    prof_mtrx = nuc_cnt_mtrx / nuc_cnt_mtrx.sum(axis=0, keepdims=True)
    prof_mtrx = prof_mtrx.tolist()

    ## return
    return prof_mtrx


## this function read a text file and loads each line into list
def load_file_lines_onto_list(file):
    f = open(file, 'r')
    list = []
    for line in f:
        list.append(line.strip('\n|\r|\r\n'))
    f.close()
    return list


## this function returns a set of motifs that form the lowest sum of Hammings distance through greedy search algorithm
## this function assumes that all DNA strings in the DNA_list have the same length
## http://www.mrgraeme.co.uk/greedy-motif-search/
def find_best_motifs_greedy(DNA_list, k, zero_prob_adj=True):
    t = len(DNA_list)
    best_motifs = []
    for DNA_str in DNA_list:
        best_motifs.append(DNA_str[0:k])
    base_DNA_str = DNA_list[0]
    other_DNA_str = DNA_list[1:t]
    
    for i in range(0, len(base_DNA_str)-k+1):
        kmer_base_DNA_str = base_DNA_str[i:i+k]
        motifs = [kmer_base_DNA_str]
        for DNA_str in other_DNA_str:
            profile_matrix = create_prof_mtrx_from_motifs(motifs, zero_prob_adj)
            most_prob_kmer = find_most_prob_kmer_based_on_prof_mtrx(DNA_str, k, profile_matrix)
            if most_prob_kmer != '':
                motifs.append(most_prob_kmer)
            else:
                motifs.append(DNA_str[0:k])  # if no most probable k-mer is found, just add the first k-mer from the DNA string
        current_motifs_cum_dist = calc_cum_dist_btwn_motifs_and_consensus(motifs)
        best_motifs_cum_dist = calc_cum_dist_btwn_motifs_and_consensus(best_motifs)
        if best_motifs_cum_dist > current_motifs_cum_dist:
            best_motifs = motifs
    return best_motifs



"""
print get_immediate_neighbors('ACGT')

print get_neighbors('ACG', 1)

print find_all_kmers_in_DNA_str('ABCDEFG', 4)

DNA_list = ['ABCDE', 'VWXYZ']
print find_all_kmers_in_DNA_list(DNA_list, 3)

print check_kmer_exists_in_DNA_str('ACA', 'ACCA', 1)
print check_kmer_exists_in_DNA_str('CCC', 'CATT', 1)

kmer = 'AAA'
DNA_list = ['AAAA', 'AACA', 'ACAA']
d = 1
print check_kmer_exists_in_all_DNA_str(kmer, DNA_list, d)

DNA_list = ['ATTTGGC', 'TGCCTTA', 'CGGTATC', 'GAAAATT']
k = 3
d = 1
print find_motifs(DNA_list, k, d)

print calc_dist_btwn_ptrn_and_DNA_str('AAA', 'TTACCTTAAC')

ptrn = 'AAA'
DNA_list = ['TTACCTTAAC', 'GATATCTGTC', 'ACGGCGTTCG', 'CCCTAAAGAG', 'CGTCAGAGGT']
print calc_cum_dist_btwn_ptrn_and_all_DNA_str(ptrn, DNA_list)

DNA_list = ['AAATTGACGCAT', 'GACGACCACGTT', 'CGTCAGCGCCTG', 'GCTGAGCACCGG', 'AGTTCGGGACAG']
k = 3
print find_median_str(DNA_list, k)

kmer = 'ACGGGGATTACC'
profile_matrix = \
[[.2, .2,  0,  0,  0,  0, .9, .1, .1, .1, .3,  0],
 [.1, .6,  0,  0,  0,  0,  0, .4, .1, .2, .4, .6],  
 [0,   0,  1,  1, .9, .9, .1,  0,  0,  0,  0,  0], 
 [.7, .2,  0,  0, .1, .1,  0, .5, .8, .7, .3, .4]]
print get_kmer_prob_based_on_prof_mtrx(kmer, profile_matrix)

DNA_str = 'ACCTGTTTATTGCCTAAGTTCCGAACAAACCCAATATAGCCCGAGGGCCT'
k = 5
profile_matrix = \
[[0.2, 0.2, 0.3, 0.2, 0.3],
 [0.4, 0.3, 0.1, 0.5, 0.1],
 [0.3, 0.3, 0.5, 0.2, 0.4],
 [0.1, 0.2, 0.1, 0.1, 0.2]]
print find_most_prob_kmer_based_on_prof_mtrx(DNA_str, k, profile_matrix)

DNA_str = 'AAGAATCAGTCA'
k = 3
profile_matrix = \
[[0.0, 0.0, 0.0], 
 [0.0, 0.0, 1.0], 
 [1.0, 1.0, 0.0], 
 [0.0, 0.0, 0.0]]
print find_most_prob_kmer_based_on_prof_mtrx(DNA_str, k, profile_matrix)

motifs = ['ACA', 'ATA', 'AAA', 'CCA']
print calc_cum_dist_btwn_motifs_and_median_str(motifs)

nuc_cnt_dict = {'A': 0, 'C':2, 'G':1, 'T':3}
print conv_nuc_cnt_dict_to_list(nuc_cnt_dict)

nuc_cnt_list = [0, 2, 1, 3]
print conv_nuc_cnt_list_to_dict(nuc_cnt_list)

motifs = ['ATC', 'ACC']
print create_nuc_cnt_mtrx_from_motifs(motifs)

motifs = ['ACA', 'ATA', 'AAA', 'CCA']
print create_prof_mtrx_from_motifs(motifs)

motifs = \
['TCGGGGGTTTTT', 
'CCGGTGACTTAC', 
'ACGGGGATTTTC', 
'TTGGGGACTTTT', 
'AAGGGGACTTCC', 
'TTGGGGACTTCC', 
'TCGGGGATTCAT', 
'TCGGGGATTCCT', 
'TAGGGGAACTAC', 
'TCGGGTATAACC']
print find_consensus_str_from_motifs(motifs)
print calc_cum_dist_btwn_motifs_and_consensus(motifs)

k = 3
DNA_list = \
['GGCGTTCAGGCA', 
 'AAGAATCAGTCA',
 'CAAGGAGTTCGC',
 'CACGTCAATCAC',
 'CAATAATATTCG']
print find_best_motifs_greedy(DNA_list, k, zero_prob_adj=False)
print find_best_motifs_greedy(DNA_list, k, zero_prob_adj=True)

DNA_list = load_file_lines_onto_list('dataset.txt')
k = 12
print find_best_motifs_greedy(DNA_list, k, zero_prob_adj=False)
print find_best_motifs_greedy(DNA_list, k, zero_prob_adj=True)

"""
