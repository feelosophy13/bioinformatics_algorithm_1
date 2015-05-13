import operator
import itertools


## function to count the number of a single (specific) k-mer pattern in a DNA string
def get_kmer_freq_cnt(DNA_str, kmer_ptrn):
    kmer_len = len(kmer_ptrn)
    n_iter = len(DNA_str) - kmer_len + 1
    freq_cnt = 0
    for i in range(0, n_iter):
        if DNA_str[i:i+kmer_len] == kmer_ptrn:
            freq_cnt += 1
    return freq_cnt


## function to count all k-mer patterns in a DNA string 
def get_kmers_freq_cnts(DNA_str, k):
    freq_cnts = {}
    DNA_len = len(DNA_str)
    n_iter = DNA_len - k + 1
    for i in range(0, n_iter):
        kmer = DNA_str[i:i+k]
        if kmer in freq_cnts:
            freq_cnts[kmer] += 1
        else:
            freq_cnts[kmer] = 1
    return freq_cnts


## function to determine the most frequent k-mer in a DNA string
## http://stackoverflow.com/questions/268272/getting-key-with-maximum-value-in-dictionary
## http://stackoverflow.com/questions/18272160/access-multiple-elements-of-list-knowing-their-index
## http://stackoverflow.com/questions/6294179/how-to-find-all-occurrences-of-an-element-in-a-list

def get_most_freq_kmers(DNA_str, k):
    freq_cnts = get_kmers_freq_cnts(DNA_str, k)

    kmers = list(freq_cnts.keys())  # list of k-mers found (order paired with cnts list)
    cnts = list(freq_cnts.values())  # list of k-mer counts (order paired with kmers list)
    max_cnt = max(cnts)  # max k-mer count discovered
    max_cnt_kmers_indices = [i for i, x in enumerate(cnts) if x == max_cnt]  # indices of desired k-mers in kmers list
    most_freq_kmers = [kmers[j] for j in max_cnt_kmers_indices]

    return most_freq_kmers


## function that takes in a string of DNA sequence and outputs the reverse complement of pattern (that reads from 5' to 3')
def get_rev_complement(DNA_str):
    complementary = ''
    for nucleotide in DNA_str:
        if nucleotide == 'A':
            complementary = 'T' + complementary
        elif nucleotide == 'T':
            complementary = 'A' + complementary
        elif nucleotide == 'C':
            complementary = 'G' + complementary
        else:
            complementary = 'C' + complementary
    return complementary


## function that takes in a substring and a DNA string and outputs indices of the substring pattern in the DNA string.
def find_ptrn_start_indices(ptrn, DNA_str):
    indices = []
    ptrn_len = len(ptrn)
    k_iter = len(DNA_str) - ptrn_len + 1 
    for i in range(0, k_iter):
        if DNA_str[i:i+ptrn_len] == ptrn:
            indices.append(i)
    return indices


## function to convert a base-10 number to another base
## e.g. binary to decimal or hexadecimal to binary
## http://stackoverflow.com/questions/2267362/convert-integer-to-a-string-in-a-given-numeric-base-in-python
def conv_base(num, b, numerals="0123456789abcdefghijklmnopqrstuvwxyz"):
    return ((num == 0) and numerals[0]) or (conv_base(num // b, b, numerals).lstrip(numerals[0]) + numerals[num % b])


## function to convert base-4 number to base-10 number
def conv_base4_num_to_base10_num(base4):
    base4 = str(base4)
    n_digits_base4 = len(base4)
    final_base10 = 0
    for i in range(0, n_digits_base4):
        base4_digit = int(base4[i])
        exp = n_digits_base4 - i - 1
        base10 = base4_digit * (4 ** exp)
        final_base10 += base10
    return final_base10


## function to convert DNA string to base-4 number representation
## i.e. replace DNA nucleotides with digits 0 through 3 
## e.g. ACGT to 0124 or AAAA to 0000 
def conv_DNA_str_to_base4_num(DNA_str):
    s = {'A':0, 'C':1, 'G':2, 'T':3} 
    base4_num = ''
    for nucleotide in DNA_str:
        digit = s[nucleotide]
        base4_num += str(digit)
    return base4_num


## function to convert DNA string to base-10 number
def conv_DNA_str_to_base10_num(DNA_str):
    DNA_base4_num = conv_DNA_str_to_base4_num(DNA_str)
    DNA_base10_num = conv_base4_num_to_base10_num(DNA_base4_num)
    return DNA_base10_num
   

## function to convert base-4 numeric representation of a DNA string to its original DNA string form
def conv_base4_num_to_DNA_str(DNA_base4_num):
    DNA_base4_num = str(DNA_base4_num)
    s = {0:'A', 1:'C', 2:'G', 3:'T'}
    DNA_str = ''
    for digit in DNA_base4_num:
        nucleotide = s[int(digit)]
        DNA_str += nucleotide
    return DNA_str


## function to convert base-10 numeric representation of DNA string to original DNA string, where k is the length of the original DNA string
def conv_base10_num_to_DNA_str(DNA_base10_num, k):
    DNA_base4_num = conv_base(DNA_base10_num, 4)
    DNA_str = conv_base4_num_to_DNA_str(DNA_base4_num)
    if len(DNA_str) < k:
        prepend_DNA_str = 'A' * (k - len(DNA_str))
        DNA_str = prepend_DNA_str + DNA_str
    return DNA_str


## function to create an ordered list of all possible permutations of k-mers
def get_ordered_kmer_permutations(k):
    nucleotides = ['A', 'C', 'G', 'T']
    perms = [''.join(p) for p in itertools.product(nucleotides, repeat=k)]
    return perms


## function to calculate frequency array of ordered k-mers
def get_ordered_kmers_freq_array(DNA_str, k):
    #ordered_kmer_perms = get_ordered_kmer_permutations(k)
    kmer_freq_cnts = [0] * (4 ** k)
    for i in range(0, len(DNA_str) - k + 1):
        kmer = DNA_str[i:i+k]
        kmer_base10 = conv_DNA_str_to_base10_num(kmer)  # k-mer's numeric representation is also the k-mer index in the frequency array
        kmer_freq_cnts[kmer_base10] += 1
    return kmer_freq_cnts 


## function to find clustered k-mers that appears at least t times in a window frame of length L
def find_clumped_kmers_inefficient(DNA_str, k, L, t):

    ## initialize an empty list to store clumped k-mers
    clumped_kmers = []

    ## number of all possible k-mer permutations 
    n = 4 ** k

    ## initialize a list whose index represents a base-10 numeric representation of a k-mer and whose value represents a condition whether the k-mer is clumped (0 or 1)
    clump_array = [0] * n
    
    ## for each window, find clumped k-mers
    for i in range(0, len(DNA_str) - L + 1):

        ## get a window frame of length L into the DNA string
        DNA_str_window = DNA_str[i:i+L]

        ## generate k-mer frequency array (whose index represents a base-10 representation of a k-mer) per window 
        freq_array = get_ordered_kmers_freq_array(DNA_str_window, k)

        ## for each k-mer, check if it occurs at least t times; if it is, mark it as clumped
        for j in range(0, n):
            if freq_array[j] >= t:
                clump_array[j] = 1

    ## scan clumped k-mer conditional array; for each k-mer, if it's marked as clumped, then add to the list of clumped k-mers
    for i in range(0, n):
        if clump_array[i] == 1:
            kmer = conv_base10_num_to_DNA_str(i, k)
            clumped_kmers.append(kmer)
        
    return clumped_kmers


## function to find clustered k-mers that appears at least t times in a window frame of length L
def find_clumped_kmers(DNA_str, k, L, t):

    ## initialize an empty list to store clumped k-mers
    clumped_kmers = []

    ## number of all possible k-mer permutations 
    n = 4 ** k

    ## initialize a list whose index represents a base-10 numeric representation of a k-mer and whose value represents a condition whether the k-mer is clumped (0 or 1)
    clump_array = [0] * n
    
    ## first window into DNA string
    DNA_str_window = DNA_str[0:L]
    
    ## generate k-mer frequency array (whose index represents a base-10 representation of a k-mer) for the first window
    freq_array = get_ordered_kmers_freq_array(DNA_str_window, k)

    ## check the frequency array for the first window to see if any clumped k-mers exist
    for i in range(0, n):
        if freq_array[i] >= t:
            clump_array[i] = 1

    ## for each moving window frame
    for i in range(1, len(DNA_str) - L):

        ## decrement count for the first k-mer pattern from the previous window
        prev_window_first_ptrn = DNA_str[i-1:i-1+k]
        j = conv_DNA_str_to_base10_num(prev_window_first_ptrn)
        freq_array[j] = freq_array[j] - 1
        
        ## increment count for the last (new) k-mer pattern from the current window
        cur_window_last_ptrn = DNA_str[i+L-k:i+L]
        j = conv_DNA_str_to_base10_num(cur_window_last_ptrn)
        freq_array[j] = freq_array[j] + 1
        
        ## after these updates, if frequency count indicates k-mer count greater than t, mark the k-mer
        if freq_array[j] >= t:
            clump_array[j] = 1

    ## scan the clumped array to find clumped k-mers
    for i in range(0, len(clump_array)):
        if clump_array[i]==1:
            clumped_kmer = conv_base10_num_to_DNA_str(i, k)
            clumped_kmers.append(clumped_kmer)

    ## return data
    return clumped_kmers




"""
## function to calculate the probability that a string pattern (ptrn) appears t or more times in a random string of length N formed from an alphabet of A letters
def get_str_ptrn_prob(N, A, ptrn, t)
    pass
"""



"""
DNA = 'ACAACTATGCATACTATCGGGAACTATCCT'
kmer = 'ACTAT'
print get_kmer_freq_cnt(DNA, kmer)

DNA = 'ACGTTGCATGTCGCATGATGCATGAGAGCT'
k = 4
print get_most_freq_kmers(DNA, k)

DNA = 'AAAACCCGGT'
print get_rev_complement(DNA)

ptrn = 'ATAT'
DNA = 'GATATATGCATATACTT'
print find_ptrn_start_indices(ptrn, DNA)

print conv_DNA_str_to_base4_num('ACGT')

print conv_base4_num_to_base10_num(13)

print conv_DNA_str_to_base10_num('AAAAAA')
print conv_DNA_str_to_base10_num('AAAAAC')
print conv_DNA_str_to_base10_num('ATGCAA')
print conv_DNA_str_to_base10_num('TTTTTT')

print conv_base(10 , 4)

print conv_base4_num_to_DNA_str(123)
print conv_base4_num_to_DNA_str('0123')

print conv_base10_num_to_DNA_str(912, 5)
print conv_base10_num_to_DNA_str(912, 7)

print get_ordered_kmer_permutations(2)

DNA = 'ACGCGGCTCTGAAA'
k = 2
print get_ordered_kmers_freq_array(DNA, k)

DNA = 'ACACGTACGT'
k = 2
L = 4
t = 2
print find_clumped_kmers_inefficient(DNA, k, L, t)
print find_clumped_kmers(DNA, k, L, t)

"""




