Program mapper.py is an implementation of a read mapping algorithm that:
- is designed to work on reads of length around 1kbp with error rate 5 − 10%,
- may fail to map some reads, but avoids incorrect alignments,
- is efficient (i.e. fast) and has good quality (i.e. high proportion of mapped reads).

Program is executable using syntax:
`python3 mapper.py reference.fasta reads.fasta output.txt`

Input data are in fasta format. Output file consists of one line for each mapped read, consisting of read
identifier, start and end positions of read alignment, separated by tabs.

It is assumed that:
- Majority of reads come from the reference sequence, but contain errors (substitutions, insertions and deletions of single nucleotides) resulting from the sequencing process.
- Errors occur independently at each position at the assumed error rate, so the total number of errors may slightly exceed 10% of the read length.
- Input data may contain unmappable reads that cannot be mapped to the reference sequence with error rate < 20%.

A read included in the output file is considered to be correctly mapped if it is mappable and the reported
mapping coordinates differ from the actual coordinates of the fragment it comes from by <= 20bp.
