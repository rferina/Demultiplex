#!/usr/bin/env python
# Author: Rachel Ferina rferina@uoregon.edu

# Check out some Python module resources:
#   - https://docs.python.org/3/tutorial/modules.html
#   - https://python101.pythonlibrary.org/chapter36_creating_modules_and_packages.html
#   - and many more: https://www.google.com/search?q=how+to+write+a+python+module

'''This module is a collection of useful bioinformatics functions
written during the Bioinformatics and Genomics Program coursework.
These functions will be useful on FASTA and FASTQ files.'''

__version__ = "0.7"         # Read way more about versioning here:
                            # https://en.wikipedia.org/wiki/Software_versioning

from lib2to3.pytree import convert
from xmlrpc.client import boolean

DNA_bases = "AGTCCCATGGCCCNNA"
RNA_bases = "UUACGAGGUUAUAGNU"

def convert_phred(letter: str) -> int:
    """Inputs a single character, and returns it as a phred score."""
    return ord(letter) - 33


def qual_score(phred_score: str) -> float:
    """Takes in an unmodified string of phred quality
    scores. Returns the average quality score of the
    input phred string.
    """
    phred_sum = 0
    for letter in phred_score:
        # convert letter to phred score
        score = convert_phred(letter)
        phred_sum += score
    # calculate average
    avg_phred = phred_sum / len(phred_score)
    return avg_phred


def validate_base_seq(seq: str, RNAflag: bool=False) -> bool:
    '''Takes in a string of DNA or RNA. Returns True if string has
    only As, Gs, Cs and Ts (or Us if RNAflag), otherwise returns False. Case insensitive.'''
    # make sequence uppercase
    seq = seq.upper()
    # if the sequence isn't RNA, count AGTC and see if equal to length of sequence
    if RNAflag == False:
        if seq.count('A') + seq.count('G') + seq.count('T') + seq.count('C') == len(seq):
            return True
    # if the sequence is RNA, count AGUC and see if equal to length of sequence
    else:
        if seq.count('A') + seq.count('G') + seq.count('U') + seq.count('C') == len(seq):
            return True
    return False


def gc_content(seq: str) -> float:
    '''Takes in a string (DNA). Returns GC content of the DNA sequence as a decimal between 0 and 1.'''
    # make sequence uppercase
    seq = seq.upper()
    # count Gs
    Gs = seq.count("G")
    # count Cs
    Cs = seq.count("C")
    return (Gs+Cs) / len(seq)


def oneline_fasta(file):
    '''Makes FASTA sequences on one line. Writes out to the file
    fa_one_line.fa. Returns the number of records so they can be
    manually compared to the number of header lines in the output file,
    to confirm the output file is accurate.'''
    # make dict with headers as keys and sequences as values
    seq_dict = {}
    with open(file, 'r') as fa:
        line_count = 0
        for line in fa:
            line_count +=1
            line = line.strip('\n')
            # only get header lines
            if line[0] == '>':
                header_line = line
            # populate dict with seq lines (non-header lines)
            else:    
                if header_line not in seq_dict:
                    seq_dict[header_line] = line
                else:
                    seq_dict[header_line] += line
    # write out to file
    fa_one_line = open('fa_one_line.fa', 'w')
    for keys,vals in seq_dict.items():
        fa_one_line.write(str(keys) + '\n' + str(vals) + '\n')
    fa_one_line.close()
    return len(seq_dict)


def reverse_complement(str):
    '''
    Takes in a string of a DNA sequence, and returns the reverse
    complement of the sequence in a new string. N's don't have
    a reverse complement, and remain the same.
    '''
    rev_str = ''
    comp_dict = {'G': 'C', 'C': 'G', 'A': 'T', 'T': 'A', 'N': 'N'}
    comp_list = []
    position = len(str)
    for nuc in range(len(str)):
        position -= 1
        rev_str += str[position]
    for base in rev_str:
        comp_list.append(comp_dict[base])
    rev_comp = ''.join(comp_list)
    return rev_comp


if __name__ == "__main__":
    # write tests for functions above
    print(reverse_complement('AGAG'))
    assert convert_phred('A') == 32, "incorrect 'A' convert_phred score"
    assert convert_phred('F') == 37, "incorrect 'F' convert_phred score"
    assert convert_phred("@") == 31, "incorrect '@' convert_phred score"
    print('passed convert_phred tests')

    assert qual_score('HJIC2@JFFH$$') == 30.166666666666668, "qual_score produced incorrect average"
    print('passed qual_score test')

    assert validate_base_seq("ACTCGCCT", False) == True, "Validate base seq does not work on DNA"
    assert validate_base_seq("UACAUG", True) == True, "Validate base seq does not work on RNA"
    assert validate_base_seq("CTGUUA",False) == False, "Validate base seq worked on invalid sequence"
    assert validate_base_seq(DNA_bases, False) == False, 'Validate base seq does not work on DNA_bases'
    assert validate_base_seq(RNA_bases, True) == False, 'Validate base seq does not work on DNA_bases'
    print('validate_base_seq passed DNA and RNA tests')

    assert gc_content("GCGCCCG") == 1
    assert gc_content("TTTATAAA") == 0
    assert gc_content("GGCCTATA") == 0.5
    assert gc_content("GCTATAAT") == 0.25
    print("gc_content passed GC content tests")

    file_fa = oneline_fasta('Danio_rerio.GRCz11.pep.all.fa')
    assert path.exists('fa_one_line.fa') == True, "fa_one_line_.fa file not created"
    print('oneline_fasta passed tests')

    assert reverse_complement('GTAGCGTA') == 'TACGCTAC', 'Reverse complement is incorrect'
    assert reverse_complement('TGTTCCGT') == 'ACGGAACA', 'Reverse complement is incorrect'
    assert reverse_complement('AGAGTCCA') == 'TGGACTCT', 'Reverse complement is incorrect'
    print('reverse_complement passed tests')