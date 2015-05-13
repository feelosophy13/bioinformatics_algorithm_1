from module1 import *
import itertools


## this function outputs a list of skew, where each skew is the difference between the number of guanines and that of cytosines (#G - #C) at (i+1)-th index in the DNA string
def create_GC_skews_list(DNA_str):
    skew = 0
    result = [0]
    for nucleotide in DNA_str:
        if nucleotide == 'G':
            skew += 1
        elif nucleotide == 'C':    
            skew -= 1
        result.append(skew)
    return result


## this function outputs the indices where minimum G-C skews occur in a DNA string
def get_min_GC_skew_indices(DNA_str):
    DNA_str = DNA_str[:-1]
    GC_skews_list = create_GC_skews_list(DNA_str)
    min_skew_val = min(GC_skews_list)
    min_skew_indices = []
    for i in range(0, len(GC_skews_list)):
        if min_skew_val == GC_skews_list[i]:
            min_skew_indices.append(i)
    return min_skew_indices


## this function calculates the number of mis-matches in two given DNA strings
def calc_Hamming_dist(DNA_str1, DNA_str2):
    if len(DNA_str1) == len(DNA_str2):
        n_mismatches = 0
        for i in range(0, len(DNA_str1)):
            if DNA_str1[i] != DNA_str2[i]:
                n_mismatches += 1
        return n_mismatches
    return None

## this function returns start indices of a pattern if it approximately matches a k-mer in a DNA string
def find_ptrn_start_indices_aprx_match(ptrn, DNA_str, d):
    indices = []
    ptrn_len = len(ptrn)
    k_iter = len(DNA_str) - ptrn_len + 1
    for i in range(0, k_iter):
        obs_kmer = DNA_str[i:i+ptrn_len]
        Hamming_distance = calc_Hamming_dist(obs_kmer, ptrn)
        if Hamming_distance <= d:
            indices.append(i)
    return indices


## this function counts the number of matches of a single (specific) k-mer that approximately matches the pattern in a DNA string
def get_kmer_freq_cnt_aprx_match(DNA_str, ptrn, d):
    freq_cnt = 0
    ptrn_len = len(ptrn)
    n_iter = len(DNA_str) - ptrn_len + 1
    for i in range(0, n_iter):
        kmer = DNA_str[i:i+ptrn_len]
        if calc_Hamming_dist(kmer, ptrn) <= d:
            freq_cnt += 1
    return freq_cnt


## this function generates approximate-match patterns whose Hamming distances to the original pattern are less than or equal to d; includes the observed pattern in the list; (THERE MUST EXIST A MORE EFFICIENT METHOD)
def get_aprx_match_ptrns(ptrn, d):
    aprx_match_ptrns = []
    ptrn_len = len(ptrn)
    perms = get_ordered_kmer_permutations(ptrn_len)
    for perm in perms:
        Hamming_distance = calc_Hamming_dist(perm, ptrn)
        if Hamming_distance <= d:
            aprx_match_ptrns.append(perm)
    return aprx_match_ptrns


## this function returns the most frequent k-mers (with mismatches and reverse complements) in a DNA string
## if greedy is set to True, then the function evaluates all possible k-mers
## if greedy is set to False, then the function evaluates only the k-mers observed in the DNA string
## if incl_rev_comp is set to True, then the function also evaluates the reverse complements of k-mers
def get_most_freq_kmers_aprx_match(DNA_str, k, d, greedy=False, incl_rev_comp=False):
    freq_array = get_ordered_kmers_freq_array(DNA_str, k)
    aprx_match_freq_dict = {}

    ## for every possible k-mer
    for i in range(0, len(freq_array)):

        ## if not greedy, don't evaluate k-mer that hasn't been observed in the DNA string
        if not greedy:
            if freq_array[i] == 0: 
                continue

        ## find k-mer and its approximate matches count
        kmer = conv_base10_num_to_DNA_str(i, k)
        aprx_match_ptrns = get_aprx_match_ptrns(kmer, d)  # includes the observed k-mer pattern itself
        aprx_match_ptrns_base10 = [conv_DNA_str_to_base10_num(ptrn) for ptrn in aprx_match_ptrns]
        aprx_match_cnt = sum([freq_array[base10] for base10 in aprx_match_ptrns_base10])
        kmer_list = [kmer]

        ## k-mer reverse complement and its approximate matches
        if incl_rev_comp:
           kmer_rev_comp = get_rev_complement(kmer)
           aprx_match_ptrns_rev_comp = [get_rev_complement(ptrn) for ptrn in aprx_match_ptrns]
           aprx_match_ptrns_base10_rev_comp = [conv_DNA_str_to_base10_num(ptrn) for ptrn in aprx_match_ptrns_rev_comp]
           aprx_match_cnt_rev_comp = sum([freq_array[base10] for base10 in aprx_match_ptrns_base10_rev_comp])
           aprx_match_cnt += aprx_match_cnt_rev_comp
           kmer_list = [kmer, kmer_rev_comp]

        ## add k-mers to the match count dictionary
        if aprx_match_cnt in aprx_match_freq_dict.keys():
            aprx_match_freq_dict[aprx_match_cnt] += kmer_list 
        else:
            aprx_match_freq_dict[aprx_match_cnt] = kmer_list

    max_cnt = max(aprx_match_freq_dict.keys())
    most_freq_kmers = aprx_match_freq_dict[max_cnt]
    return most_freq_kmers



"""
DNA = 'CATGGGCATCGGCCATACGCC'
print create_GC_skews_list(DNA)

DNA = 'TAAAGACTGCCGAGAGGCCAACACGAGTGCTAGAACGAGGGGCGTAAACGCGGGTCCGAT'
print create_GC_skews_list(DNA)
print get_min_GC_skew_indices(DNA)

DNA1 = 'GGGCCGTTGGT'
DNA2 = 'GGACCGTTGAC'
print calc_Hamming_dist(DNA1, DNA2)

ptrn = 'ATTCTGGA'
DNA = 'CGCCCGAATCCAGAACGCATTCCCATATTTCGGGACCACTGGCCTCCACGGTACGGACGTCAATCAAAT'
d = 3
print find_ptrn_start_indices_aprx_match(ptrn, DNA, d)

DNA = 'TTTAGAGCCTTCAGAGG'
ptrn = 'GAGG'
d = 2
print get_kmer_freq_cnt_aprx_match(DNA, ptrn, d)

DNA = 'ACGTTGCATGTCGCATGATGCATGAGAGCT'
k = 4
d = 1
print get_most_freq_kmers_aprx_match(DNA, k, d)

DNA = 'ACGTTGCATGTCGCATGATGCATGAGAGCT'
k = 4
d = 1
print get_most_freq_kmers_aprx_match(DNA, k, d, greedy=False, incl_rev_comp=False)
print get_most_freq_kmers_aprx_match(DNA, k, d, greedy=False, incl_rev_comp=True)

"""
