#!/usr/bin/env python

#
#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

# Originally made by Eva Strauch
#   touch up by Brian Coventry


import optparse, time, pdb, sys, time, numpy
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import  MeltingTemp


def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

def main():

    #input
    parser = optparse.OptionParser()
    parser.add_option('--in', action = 'store', type = 'string', dest = 'infile', help = 'input file containing fasta sequences, DNA')
    parser.add_option('--table', action = 'store', type = 'string', dest = 'table', help = 'top frequency codons')
 
    #optional:
    parser.add_option('--start', action = 'store', type = 'int', default = '1', dest = 'start', help = 'start with library')
    parser.add_option('--end', action = 'store', type = 'int', default = '0', dest = 'end', help = 'end position with library')

    (option, args) = parser.parse_args()
    start = option.start - 1

        #open the files for input and output
    input_file = open(option.infile, 'r')
    out = open(( option.infile + '_SSM.tab'), 'w')
    #sum_out = open(( option.infile + '_SSM_primers_sum.tab'), 'w')
    end = option.end    

    #def build_codons_dict(input_file):
    tab  = open(option.table)
    codon_tab = {}
    for line in tab:
            f = line.strip().split()
            codon_tab[f[0]] = format(f[1])

    #out = open('ssm_' + .tab' , 'w')

    #iterating over fasta sequences
    for cur_record in SeqIO.parse(input_file, "fasta") :
        name = cur_record.name
        myseq = cur_record.seq
        seq = myseq #nterm + myseq + cterm
        proteinseq = ''
        
        #if end is not defined default to the end of the sequence
        if option.end == 0:
            end = len(myseq)
        else:
            end = option.end - 3

        for pos in range(start, end, 3):
            codon = ''
            fiveprime = ''
            threeprime = ''

                        #sanity check get the codons and translate them...
            for i in range(0,3):
                codon += seq[pos + i]

            curr_aa = Seq(codon).translate()
            eprint("replacing codon: ", codon , "for" , curr_aa , (pos+3)//3 , "base position", pos+1 )
                                    
            #for frame check later
            proteinseq += Seq(codon).translate()
            
            for i in range(start,pos):
                fiveprime += seq[i]            
            #print fiveprime     
        
            for i in range(pos+3,len(seq)):
                      threeprime += seq[i]

            #now iterate through codons
            for aa in codon_tab:
                if str(aa) == str(curr_aa):
                    eprint("wt sequence")

                else:   
                    print(fiveprime + codon_tab[aa] + threeprime, name + '__' + str(pos//3+1) + '__' + aa ) 

            


        eprint("frame check, this is the protein sequence that was in frame for this:\n" , proteinseq)    
                


if __name__ == '__main__':
    main()    
