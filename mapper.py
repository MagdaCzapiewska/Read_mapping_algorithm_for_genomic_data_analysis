import sys
from Bio import SeqIO
from collections import deque
import math

# Using an edit-distance-like dynamic programming formulation, we can
# look for approximate occurrences of p in t.

import numpy

# Assume x is the string labeling rows of the matrix and y is the
# string labeling the columns

def trace(D, x, y):
    ''' Backtrace edit-distance matrix D for strings x and y '''
    i, j = len(x), len(y)
    xscript = []
    while i > 0:
        diag, vert, horz = sys.maxsize, sys.maxsize, sys.maxsize
        delt = None
        if i > 0 and j > 0:
            delt = 0 if x[i-1] == y[j-1] else 1
            diag = D[i-1, j-1] + delt
        if i > 0:
            vert = D[i-1, j] + 1
        if j > 0:
            horz = D[i, j-1] + 1
        if diag <= vert and diag <= horz:
            # diagonal was best
            xscript.append('R' if delt == 1 else 'M')
            i -= 1; j -= 1
        elif vert <= horz:
            # vertical was best; this is an insertion in x w/r/t y
            xscript.append('I')
            i -= 1
        else:
            # horizontal was best
            xscript.append('D')
            j -= 1
    # j = offset of the first (leftmost) character of t involved in the
    # alignment
    return j, (''.join(xscript))[::-1] # reverse and string-ize

def kEditDp(p, t):
    ''' Find and return the alignment of p to a substring of t with the
        fewest edits.  We return the edit distance, the offset of the
        substring aligned to, and the edit transcript.  If multiple
        alignments tie for best, we report the leftmost. '''
    D = numpy.zeros((len(p)+1, len(t)+1), dtype=int)
    # Note: First row gets zeros.  First column initialized as usual.
    D[1:, 0] = range(1, len(p)+1)
    for i in range(1, len(p)+1):
        for j in range(1, len(t)+1):
            delt = 1 if p[i-1] != t[j-1] else 0
            D[i, j] = min(D[i-1, j-1] + delt, D[i-1, j] + 1, D[i, j-1] + 1)
    # Find minimum edit distance in last row
    mnJ, mn = None, len(p) + len(t)

    for j in range(len(t)+1):
        if D[len(p), j] < mn:
            mnJ, mn = j, D[len(p), j]
    # Backtrace; note: stops as soon as it gets to first row
    off, xcript = trace(D, p, t[:mnJ])
    # Return edit distance, offset into T, edit transcript
    return mn, off, mnJ, xcript, D

def encode_kmer(kmer):
    encoding = {'A': 0b00, 'C': 0b01, 'G': 0b10, 'T': 0b11}
    val = 0
    for base in kmer:
        val = (val << 2) | encoding[base]
    return val

def roll_kmer(prev_val, new_base, k):
    encoding = {'A':0b00, 'C':0b01, 'G':0b10, 'T':0b11}
    mask = (1 << (2*k)) - 1 # keep only k-mer bits
    val = ((prev_val << 2) & mask) | encoding[new_base]
    return val

def get_minimizers_index(seq, k=15, w=10):
    if len(seq) < k:
        return {}

    kmer_val = encode_kmer(seq[:k])
    # COMMENT_ON_CHANGES:
    # In the first version of the program I precomputed a full list of k-mer hashes
    # and then iterated over it, which caused excessive memory usage.
    # In the new approach I compute rolling k-mer hashes on-the-fly.
    
    # OLD_CODE_1_START - precomputing a full list of k-mer hashes
    # kmer_vals = [kmer_val]

    # Compute rolling hash for the remaining k-mers
    # for i in range(1, len(seq)-k+1):
    #     kmer_val = roll_kmer(kmer_val, seq[i+k-1], k)
    #     kmer_vals.append(kmer_val)
    # OLD_CODE_1_END

    dq = deque()  # stores pairs (hash, position)
    index = {}
    last_added = set()
    prev_min_val = None
    
    # OLD_CODE_2_START - iterating over the list of hash values for all k-mers
    # for i, val in enumerate(kmer_vals):
    # OLD_CODE_2_END
    # NEW_CODE_1_START - computing rolling k-mer hashes on-the-fly
    for i in range(len(seq)-k+1):
        if i > 0:
            kmer_val = roll_kmer(kmer_val, seq[i+k-1], k)
        val = kmer_val
        # NEW_CODE_1_END
        # Remove elements outside the window
        while dq and dq[0][1] <= i - w:
            dq.popleft()

        # Remove elements larger than current (cannot be minimizer)
        while dq and dq[-1][0] > val:
            dq.pop()

        dq.append((val, i))

        # First minimizer appears after the first w k-mers
        if i >= w - 1:
            # smallest hash at the front of deque
            min_val = dq[0][0]
            if min_val != prev_min_val:
                last_added = set()
                prev_min_val = min_val

            # Add all elements in deque with the same min_val
            # Iterate over all positions with this min_val
            j = 0
            while j < len(dq) and dq[j][0] == min_val:
                pos = dq[j][1]
                if (min_val, pos) not in last_added:
                    if min_val not in index:
                        index[min_val] = []
                    index[min_val].append(pos)
                    last_added.add((min_val, pos))
                j += 1

    return index

def main():
    if len(sys.argv) != 4:
        print("Usage: python3 mapper.py reference.fasta reads.fasta output.txt")
        sys.exit(1)

    ref_path = sys.argv[1]
    reads_path = sys.argv[2]
    out_path = sys.argv[3]

    reference_records = [rec for rec in SeqIO.parse(ref_path, "fasta")]
    if not reference_records:
        print("Error: reference file is empty!")
        sys.exit(1)
    reference = reference_records[0]

    ref_minimizers = get_minimizers_index(str(reference.seq))

    reads = [r for r in SeqIO.parse(reads_path, "fasta")]
    if not reads:
        print("Warning: reads file is empty.")

    with open(out_path, "w") as fout:
        for r in reads:

            read_min = get_minimizers_index(str(r.seq))

            offsets = []
            for h in read_min:
                if h in ref_minimizers:
                    for r_pos in read_min[h]:
                        ref_positions = ref_minimizers[h]

                        if len(ref_positions) != 1:
                            continue

                        ref_pos = ref_positions[0]
                        diff = ref_pos - r_pos
                        offsets.append(diff)
            if offsets:
                offsets.sort()
                window_size = math.ceil(len(r.seq) * 0.2)

                best_count = 0
                best_clusters = []
                j = 0

                for i in range(len(offsets)):
                    while offsets[i] - offsets[j] > window_size:
                        j += 1

                    count = i - j + 1

                    if count > best_count:
                        best_count = count
                        best_clusters = [offsets[j]]
                    elif count == best_count:
                        best_clusters.append(offsets[j])

                read_seq = str(r.seq)
                best_dp_score = float('inf')
                best_dp_result = None

                for cluster_offset in best_clusters:
                    candidate_start = max(0, cluster_offset - window_size)
                    candidate_end = min(cluster_offset + len(read_seq) + window_size, len(reference.seq))

                    ref_region = str(reference.seq[candidate_start:candidate_end])

                    edit_dist, aln_start, aln_end, transcript, D = kEditDp(read_seq, ref_region)
                    aln_start_global = candidate_start + aln_start
                    aln_end_global = candidate_start + aln_end

                    if edit_dist < best_dp_score:
                        best_dp_score = edit_dist
                        best_dp_result = (aln_start_global, aln_end_global, cluster_offset)

                ed_too_large = math.ceil(len(read_seq) * 0.2)
                if best_dp_result is not None and best_dp_score < ed_too_large:
                    best_start, best_end, used_offset = best_dp_result
                    fout.write(f"{r.id}\t{best_start}\t{best_end}\n")

if __name__ == "__main__":
    main()
