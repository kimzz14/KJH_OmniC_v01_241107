import re
from itertools import groupby
from joblib import Parallel, delayed
import gzip
class SamLine:
    def __init__(self, line):
        self.line = line.rstrip()
        self.data_LIST = self.line.split('\t')

        self.QNAME = self.data_LIST[ 0]
        self.FLAG  = self.data_LIST[ 1]
        self.RNAME = self.data_LIST[ 2]
        self.POS   = self.data_LIST[ 3]
        self.MAPQ  = self.data_LIST[ 4]
        self.CIGAR = self.data_LIST[ 5]
        self.RNEXT = self.data_LIST[ 6]
        self.PNEXT = self.data_LIST[ 7]
        self.TLEN  = self.data_LIST[ 8]
        self.SEQ   = self.data_LIST[ 9]
        self.QUAL  = self.data_LIST[10]

        self.FLAG  = int(self.FLAG)
        #Define strand
        if self.FLAG & 16 == 0:
            self.strand = '+'
        else:
            self.strand = '-'

        #Define CIGAR_LIST
        self.CIGAR_LIST = [] 
        for match in re.findall(r'(\d+)([A-Za-z]+)', self.CIGAR):
            self.CIGAR_LIST += [(int(match[0]), match[1])]

        #Read optional fields
        self.option_DICT = {}
        for data in self.data_LIST[11:]:
            TAG, TYPE, VALUE = data.split(':')
            if TYPE == 'i':
                self.option_DICT[TAG] = int(VALUE)
            else:
                self.option_DICT[TAG] = VALUE

        #Define score
        self.score = self.option_DICT['AS']

        #Define sorted_CIGAR_LIST
        if self.strand == '+':
            self.sorted_CIGAR_LIST = self.CIGAR_LIST
        else:
            self.sorted_CIGAR_LIST = self.CIGAR_LIST[::-1]

    def isUsable(self):
        if self.CIGAR == '*': return False
        if len(self.CIGAR_LIST) > 2: return False
        if self.sorted_CIGAR_LIST[0][1] != 'M': return False

        return True
    
    def toPrint(self):
        print('  '.join(map(str, ['   ', self.score, self.strand, self.sorted_CIGAR_LIST, self.line])))
    
    def toString(self, isFirst):
        if isFirst == True:
            self.data_LIST[ 1] = self.FLAG | 64
        else:
            self.data_LIST[ 1] = self.FLAG | 128

        return '\t'.join(map(str, self.data_LIST))

def filter(line_LIST):
    samLine_LIST = []
    for line in line_LIST:
        samLine = SamLine(line)

        if samLine.isUsable() == False:
            continue

        samLine_LIST += [samLine]

        if len(samLine_LIST) == 1:
            if samLine.score < samLine.sorted_CIGAR_LIST[0][0] - 5:
                return None
        elif len(samLine_LIST) == 2:
            break

    if len(samLine_LIST) == 0:
        return None
    elif len(samLine_LIST) == 1:
        return samLine_LIST[0]
    else:
        if samLine_LIST[0].score - samLine_LIST[1].score > 15:
            return samLine_LIST[0]
        else:
            return None

import subprocess
class SAMTOOLS:
    def __init__(self):
        pass

    def merge(self, outfile, prefix):
        command_LIST  = ['samtools']
        command_LIST += ['merge']
        command_LIST += ['-o', outfile]
        command_LIST += ['-o']


        command = 'samtools merge -o {0} {1}*.sam.gz'.format(outfile, prefix)
        #print(command)
        process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        process.wait()
    

from datetime import datetime

def run_batch(batchN, inSAM, prefix, isFirst):
    def run_single(batchIDX):
        fin = gzip.open(inSAM, 'rt')
        fout = gzip.open(prefix + '.' + str(batchIDX).zfill(6) + '.sam.gz', 'wt')
        for line in fin:
            fout.write(line)
            if line.startswith('@PG') == True:
                break
        
        readIDX = -1
        for readID, line_LIST in groupby(fin, lambda line: line.rstrip('\n').split('\t')[0]):
            readIDX += 1
            if readIDX%batchN != batchIDX: continue

            samLine = filter(line_LIST)
            if samLine == None: continue
            fout.write(samLine.toString(isFirst=isFirst) + '\n')
        fout.close()
        fin.close()
    
    Parallel(n_jobs=batchN)(delayed(run_single)(batchIDX) for batchIDX in range(0, batchN))

from optparse import OptionParser
import sys
#option parser
parser = OptionParser(usage="""Run annotation.py \n Usage: %prog [options]""")
parser.add_option("-t","--threadN",action = 'store',type = 'int',dest = 'threadN',help = "")
parser.add_option("-1","--in1",action = 'store',type = 'string',dest = 'inSAM1',help = "")
parser.add_option("-2","--in2",action = 'store',type = 'string',dest = 'inSAM2',help = "")
parser.add_option("-o","--out",action = 'store',type = 'string',dest = 'outSAM',help = "")

(opt, args) = parser.parse_args()
if opt.threadN == None or opt.inSAM1 == None or opt.inSAM2 == None or opt.outSAM == None:
    print('Basic usage')
    print('')
    print('     python filter_sam.py -t 24 -1 test_1.sam.gz -2 test_2.sam.gz -o test.uniq.bam')
    print('')
    sys.exit()

threadN = opt.threadN
inSAM1 = opt.inSAM1
inSAM2 = opt.inSAM2
outSAM = opt.outSAM

tmpPrefix = './tmp/tmp.filter.' + datetime.now().strftime("%Y%m%d%H%M%S")


run_batch(24, inSAM1, tmpPrefix + '_1', isFirst = True)
run_batch(24, inSAM2, tmpPrefix + '_2', isFirst = False)

samtools = SAMTOOLS()
samtools.merge(outSAM, tmpPrefix)