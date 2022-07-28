#!/usr/bin/env python
# Author: Rachel Ferina rferina@uoregon.edu

# Check out some Python module resources:
#   - https://docs.python.org/3/tutorial/modules.html
#   - https://python101.pythonlibrary.org/chapter36_creating_modules_and_packages.html
#   - and many more: https://www.google.com/search?q=how+to+write+a+python+module

'''This module is a collection of useful bioinformatics functions
written during the Bioinformatics and Genomics Program coursework.
These functions will be useful on FASTA and FASTQ files.'''

__version__ = "0.6"         # Read way more about versioning here:
                            # https://en.wikipedia.org/wiki/Software_versioning

from lib2to3.pytree import convert
from xmlrpc.client import boolean


DNA_bases = ""
RNA_bases = ""

def convert_phred(letter: str) -> int:
    """Converts a single character into a phred score"""
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


def validate_base_seq(seq: str, RNAflag: bool) -> bool:
    '''Takes in a string of DNA or RNA. Returns True if string is composed
    of only As, Ts (or Us if RNAflag), Gs, Cs; otherwise false. Case insensitive.'''
    # make sequence uppercase
    seq = seq.upper()
    # if the sequence isn't RNA, count AGTC and see if equal to length of sequence
    if RNAflag is False:
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


def init_list(lst: list, value: float=0.0) -> list:
    '''This function takes an empty list and will populate it with
    the value passed in "value". If no value is passed, initializes list
    with 101 values of 0.0.'''
    # account for no value parameter entered
    if len(lst) == 0:
        lst = [0.0] * 101
    else:
        lst = [value] * 101
    return lst


def populate_list(file: str) -> tuple[list, int]:
    """Opens a FASTQ file and decodes Phred quality scores to numbers
    accounting for Phred+33. Sums the quality scores for each position
    and counts the total number of lines in the FASTQ file.
    Returns the array and line count."""

    lst = init_list([])
    with open(file, 'r') as fq:
        line_count = 0
        for line in fq:
            line = line.strip('\n')
            line_count += 1
            # obtain lines with quality scores
            if line_count % 4 == 0:
                # specify for position when converting phred score
                for count, letter in enumerate(line):
                    lst[count] += bioinfo.convert_phred(letter)
    return (lst, line_count)




if __name__ == "__main__":
    # write tests for functions above
    assert convert_phred('A') == 32, "incorrect 'A' convert_phred score"
    assert convert_phred('F') == 37, "incorrect 'F' convert_phred score"
    assert convert_phred("@") == 31, "incorrect '@' convert_phred score"
    print('passed convert_phred tests')

    assert qual_score('HJIC2@JFFH$$') == 30.166666666666668, "qual_score produced incorrect average"
    print('passed qual_score test')

    assert validate_base_seq("ACTCGCCT", False) == True, "Validate base seq does not work on DNA"
    assert validate_base_seq("UACAUG", True) == True, "Validate base seq does not work on RNA"
    assert validate_base_seq("CTGUUA",False) == False, "Validate base seq worked on invalid sequence"
    print('validate_base_seq passed DNA and RNA tests')

    assert gc_content("GCGCCCG") == 1
    assert gc_content("TTTATAAA") == 0
    assert gc_content("GGCCTATA") == 0.5
    assert gc_content("GCTATAAT") == 0.25
    print("gc_content passed GC content tests")

    file_fa = oneline_fasta('Danio_rerio.GRCz11.pep.all.fa')
    assert path.exists('fa_one_line.fa') == True, "fa_one_line_.fa file not created"
    print('oneline_fasta passed tests')