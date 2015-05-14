import random

from module1 import *
from module2 import *
from module3 import *


## this function create a list of most probable motifs with length k from a list o DNA strings, evaluated for DNA string based on k-mer's profile matrix
def create_motifs_based_on_prof_mtrx(prof_mtrx, DNA_list, k, zero_prob_adj=True):
    motifs = []
    for DNA_str in DNA_list:
        most_prob_kmer = find_most_prob_kmer_based_on_prof_mtrx(DNA_str, k, prof_mtrx, zero_prob_adj)
        motifs.append(most_prob_kmer)
    return motifs


## this function stochastically tries to find a set of motifs that minimizes the sum of Hamming's distance to its consensus
def find_good_motifs_rand_srch(DNA_list, k):
    n = len(DNA_list[0])
    t = len(DNA_list)

    ## randomly select k-mer motifs in each string in DNA_list
    motifs = []
    for DNA_str in DNA_list:
        rand_ind_strt = random.randrange(0, n-k+1)
        rand_ind_end = rand_ind_strt + k
        motif = DNA_str[rand_ind_strt:rand_ind_end]
        motifs.append(motif)
    best_motifs = motifs

    ## loop until cumulative Hamming distance stops minimizing
    while True:
        prof_mtrx = create_prof_mtrx_from_motifs(motifs)
        motifs = create_motifs_based_on_prof_mtrx(prof_mtrx, DNA_list, k)
        cur_motifs_cum_dist = calc_cum_dist_btwn_motifs_and_consensus(motifs)
        best_motifs_cum_dist = calc_cum_dist_btwn_motifs_and_consensus(best_motifs)
        if cur_motifs_cum_dist < best_motifs_cum_dist:
            best_motifs = motifs
        else:
            return best_motifs


## this function runs find_good_motifs_rand_srch() iteratively to find the best set of motifs found in those iterations 
def find_good_motifs_rand_srch_iter(DNA_list, k, n_iter=10):
    lowest_cum_dist = float('Inf')
    for i in range(0, n_iter):
        motifs_rand_srch = find_good_motifs_rand_srch(DNA_list, k)
        cum_dist = calc_cum_dist_btwn_motifs_and_consensus(motifs_rand_srch)
        if cum_dist < lowest_cum_dist:
            best_motifs = motifs_rand_srch
            lowest_cum_dist = cum_dist
    return best_motifs


## this function normalizes the elements in a list of probability distribution so they add up to 1
def normalize_prob_dist(prob_dist):
    prob_sum = sum(prob_dist)
    prob_dist = np.array(prob_dist) / prob_sum
    prob_dist = prob_dist.tolist()
    return prob_dist




## this function stochastically converges to find a set of motifs that minimizes the sum of Hamming's distance to its consensus by utilizing randomized search and Gibbs sampling
def find_good_motifs_rand_srch_Gibbs_sample(DNA_list, k, N):
    n = len(DNA_list[0])
    t = len(DNA_list)

    ## randomly select k-mer motifs in each string in DNA_list
    motifs = []
    for DNA_str in DNA_list:
        rand_ind_strt = random.randrange(0, n-k+1)
        rand_ind_end = rand_ind_strt + k
        motif = DNA_str[rand_ind_strt:rand_ind_end]
        motifs.append(motif)
    best_motifs = motifs

    ## loop until cumulative Hamming distance stops minimizing
    for i in range(0, N):
        j = np.random.choice(t)
        motifs_except_jth = motifs[0:j] + motifs[j+1::]
        prof_mtrx = create_prof_mtrx_from_motifs(motifs_except_jth)

        DNA_str = DNA_list[j]
        new_motif_jth = find_most_prob_kmer_based_on_prof_mtrx(DNA_str, k, prof_mtrx)

        motifs[j] = new_motif_jth
        cur_motifs_cum_dist = calc_cum_dist_btwn_motifs_and_consensus(motifs)
        best_motifs_cum_dist = calc_cum_dist_btwn_motifs_and_consensus(best_motifs)
        if cur_motifs_cum_dist < best_motifs_cum_dist:
            best_motifs = motifs
    return best_motifs


## this function calls find_good_motifs_rand_srch_Gibbs_sample() iteratively and finds a set of motifs that yields the lowest sum of Hamming's distance
def find_good_motifs_rand_srch_Gibbs_sample_iter(DNA_list, k, N, n_iter):
    lowest_cum_dist = float('Inf')
    for i in range(0, n_iter):
        motifs_rand_srch = find_good_motifs_rand_srch_Gibbs_sample(DNA_list, k, N)
        cum_dist = calc_cum_dist_btwn_motifs_and_consensus(motifs_rand_srch)
        if cum_dist < lowest_cum_dist:
            best_motifs = motifs_rand_srch
            lowest_cum_dist = cum_dist
    return best_motifs


"""
##
DNA_list = ['TTACCTTAAC', 'GATGTCTGTC', 'ACGGCGTTAG', 'CCCTAACGAG', 'CGTCAGAGGT']
prof_mtrx = \
[[0.8,   0,   0, 0.2], 
 [  0, 0.6, 0.2,   0],
 [0.2, 0.2, 0.8,   0],
 [  0, 0.2,   0, 0.8]]
print create_motifs_based_on_prof_mtrx(prof_mtrx, DNA_list, 4, zero_prob_adj=True)

##
k = 8
DNA_list = \
['CGCCCCTCTCGGGGGTGTTCAGTAAACGGCCA',
 'GGGCGAGGTATGTGTAAGTGCCAAGGTGCCAG',
 'TAGTACCGAGACCGAAAGAAGTATACAGGCGT',
 'TAGATCAAGTTTCAGGTGCACGTCGGTGAACC',
 'AATCCACCAGCTCCACGTGCAATGTTGGCCTA']
x = find_good_motifs_rand_srch_iter(DNA_list, k, n_iter=1000)
print x
print calc_cum_dist_btwn_motifs_and_consensus(x)

y = ['TCTCGGGG', 'CCAAGGTG','TACAGGCG', 'TTCAGGTG', 'TCCACGTG']
print y
print calc_cum_dist_btwn_motifs_and_consensus(y)

##
k = 15
DNA_list = load_file_lines_onto_list('dataset2.txt')
x = find_good_motifs_rand_srch_iter(DNA_list, k)
print x
print calc_cum_dist_btwn_motifs_and_consensus(x)

## 
prob_dist = [0.1, 0.2, 0.3]
print normalize_prob_dist(prob_dist)

##
k = 8
N = 100
n_iter = 20
DNA_list = \
['CGCCCCTCTCGGGGGTGTTCAGTAAACGGCCA',
 'GGGCGAGGTATGTGTAAGTGCCAAGGTGCCAG',
 'TAGTACCGAGACCGAAAGAAGTATACAGGCGT',
 'TAGATCAAGTTTCAGGTGCACGTCGGTGAACC',
 'AATCCACCAGCTCCACGTGCAATGTTGGCCTA']
x = find_good_motifs_rand_srch_Gibbs_sample_iter(DNA_list, k, N, n_iter)
print x
print calc_cum_dist_btwn_motifs_and_consensus(x)

y = ['TCTCGGGG', 'CCAAGGTG', 'TACAGGCG', 'TTCAGGTG', 'TCCACGTG']
print y
print calc_cum_dist_btwn_motifs_and_consensus(y)

##
k = 15
N = 200
n_iter = 20
DNA_list = load_file_lines_onto_list('dataset3.txt')
x = find_good_motifs_rand_srch_Gibbs_sample_iter(DNA_list, k, N, n_iter)
print x
print calc_cum_dist_btwn_motifs_and_consensus(x)

"""
