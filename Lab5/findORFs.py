#!/usr/bin/env python3
# Name: Your full name (aspandit)
# Group Members: List full names (CATS usernames) or “None”

from sequenceAnalysis import FastAreader as FastaReader

########################################################################
# CommandLine
########################################################################
class CommandLine() :
    '''
    Handle the command line, usage and help requests.

    CommandLine uses argparse, now standard in 2.7 and beyond. 
    it implements a standard command line argument parser with various argument options,
    a standard usage and help.

    attributes:
    all arguments received from the commandline using .add_argument will be
    avalable within the .args attribute of object instantiated from CommandLine.
    For example, if myCommandLine is an object of the class, and requiredbool was
    set as an option using add_argument, then myCommandLine.args.requiredbool will
    name that option.
 
    '''
    
    def __init__(self, inOpts=None) :
        '''
        Implement a parser to interpret the command line argv string using argparse.
        '''
        
        import argparse
        self.parser = argparse.ArgumentParser(description = 'Program prolog - a brief description of what this thing does', 
                                             epilog = 'Program epilog - some other stuff you feel compelled to say', 
                                             add_help = True, #default is True 
                                             prefix_chars = '-', 
                                             usage = '%(prog)s [options] -option1[default] <input >output'
                                             )
        self.parser.add_argument('-lG', '--longestGene', action = 'store', nargs='?', const=True, default=False, help='longest Gene in an ORF')
        self.parser.add_argument('-mG', '--minGene', type=int, choices= (0,100,200,300,500,1000), default=100, action = 'store', help='minimum Gene length')
        self.parser.add_argument('-s', '--start', action = 'append', default = ['ATG'],nargs='?', 
                                 help='start Codon') #allows multiple list options
        self.parser.add_argument('-t', '--stop', action = 'append', default = ['TAG','TGA','TAA'],nargs='?', help='stop Codon') #allows multiple list options
        self.parser.add_argument('-v', '--version', action='version', version='%(prog)s 0.1')  
        if inOpts is None :
            self.args = self.parser.parse_args()
        else :
            self.args = self.parser.parse_args(inOpts)

########################################################################
# Main
# Here is the main program
# 
#
########################################################################

class ORFfinder:
    """
    Parse the open reading frames (ORFs) found within a given DNA sequence and its reverse complement.
    Following the path of DNA polymerase, the following cases will be returned as ORFs for all frames

    START-STOP
    START-End Of Line
    Beginning Of line-STOP
    No START or STOP found (whole sequence)

    cases where multiple start codons are encountered before a stop are handled based on lG or longestGene:

    if lG is True:
    START1-START2-STOP -> START1-STOP

    if lG is False:
    START1-START2-STOP -> START1-STOP, START2-STOP

    The main method to call is the .parse() method, which takes a list of start codons to consider, a list
    of stop codons to consider, and whether lG is True or False, as well as a doesReturn parameter that dictates
    if the method only modifies the orfs attribute in place or also returns it.
    """

    complement = {
        "A": "T",
        "T": "A",
        "C": "G",
        "G": "C"
    }

    def __init__(self, topStrand):
        """constructor for ORFfinder objects, intializes a top (base) strand of
        dna and a bottom (reverse complement)"""
        self.top = topStrand.replace(' ', '')  # 5 -> 3 (read L -> R)
        self.bottom  = ''.join([self.complement[base] for base in self.top])  # 3 -> 5 (Read R -> L)
        self.orfs = []


    def _parseTop(self, startCodons=['ATG'], stopCodons=['TAA', 'TGA', 'TAG'], lG=False) -> None:
        """generates tuples of format (frame, start pos, end pos, length) for all valid
        ORFs in the base dna strand and appends them to self.orfs"""
        
        # boolean conditions for each frame to know if a start/stop codon has been
        # encountered in that frame
        hitStart1, hitStart2, hitStart3 = False, False, False
        hitStop1, hitStop2, hitStop3 = False, False, False

        # define helper function for use in base iteration
        def appendToORF(stack, frame, hitStart, i) -> None:
            """nested helper function which handles appending orfs when a stop codon is hit"""
            if hitStart:
                if not lG:
                    # handle when not longestGene
                    for start in stack:
                        # print(f'+{frame} {start}..{i+3} {i+4-start}', self.top[i+2:i+3-start])  # debug print
                        self.orfs.append((frame, start, i+3, i+4-start, 'yes stop & start'))
                else:  
                    # handle when only longestGene
                    # print(stack)
                    # print(f'+{frame} {stack}..{i+3} {i+4}')
                    self.orfs.append((frame, min(stack), i+3, i+4-min(stack), 'yes stop & start'))
            else:
                # print(f'+{frame} 1..{i+3} {i+3}')  # debug print
                self.orfs.append((frame, 1, i+3, i+3, 'yes stop no start'))
    
        # define helper function to handle no start/stop codon edge case
        def nothingInGene(frame, hitStop, hitStart) -> None:
            """handles the condition where no start or stop codon is found in a
            frame and the entire frame is then treated as an ORF"""
            if not hitStop and not hitStart:
                # print(f'+{frame} 1..{len(self.top)} {len(self.top)}')  # debug print
                self.orfs.append((frame, 1, len(self.top), len(self.top), 'nothing in gene'))

        # the following stacks hold start codon indices in respective frames until hitting a stop codon
        stack1, stack2, stack3 = [], [], []

        # iterates through every base
        for i in range(len(self.top)-2):
            codon = self.top[i:i+3]
            frame = i % 3 + 1
            # print('frame:', frame, ', codon:', codon, ', pos:', i+1)  # debug print

            # handles building stacks for each frame
            if codon in startCodons:
                # print(codon)
                # print('^ start codon')  # debug print
                match frame:
                    case 1:
                        hitStart1 = True
                        stack1.append(i+1)
                    case 2:
                        hitStart2 = True
                        stack2.append(i+1)
                    case 3:
                        hitStart3 = True
                        stack3.append(i+1)

            # handles building orfs for each frame & stack
            if codon in stopCodons:
                # print(codon)
                # print('^ stop codon')  # debug print
                match frame:
                    case 1:
                        # print('stack for frame', frame, stack1)  # debug print
                        if stack1:
                            hitStop1 = True
                            appendToORF(stack1, 1, hitStart1, i)
                            stack1 = []
                    case 2:
                        # print('stack for frame', frame, stack2)  # debug print
                        if stack2:
                            hitStop2 = True
                            appendToORF(stack2, 2, hitStart2, i)
                            stack2 = []
                    case 3:
                        # print('stack for frame', frame, stack3)  # debug print
                        if stack3:
                            hitStop3 = True
                            appendToORF(stack3, 3, hitStart3, i)
                            stack3 = []
        
        # handle cases for each frame where neither a start nor stop
        # codon are found
        nothingInGene(1, hitStop1, hitStart1)
        nothingInGene(2, hitStop2, hitStart2)
        nothingInGene(3, hitStop3, hitStart3)


    def _parseBottom(self, startCodons=['ATG'], stopCodons=['TAA', 'TGA', 'TAG'], lG=False) -> None:
        """generates tuples of format (frame, start pos, end pos, length) for all valid
        ORFs in the complement dna strand and appends them to self.orfs"""
        
        # since we are iterating in the opposite direction (->) as DNA Pol III (<-)
        # in this case, the start and stop codons are read in reverse and should be
        # counted in reverse
        startCodons = [codon[::-1] for codon in startCodons]
        stopCodons = [codon[::-1] for codon in stopCodons]

        # boolean conditions for each frame to know if a start codon has been
        # encountered in that frame
        hitStart1, hitStart2, hitStart3 = False, False, False

        # stacks for holding start codon positions
        start1, start2, start3 = [], [], []

        # keeping track of the stop codon position for each frame, init at 1 since w/o stop codon
        # the reading would end at the end of the gene or the first base in the r. complement strand
        stopPos1, stopPos2, stopPos3 = 1, 1, 1

        # nested helper func for appending orf
        def appendORF(startStack: list, stopPos: int, frame: int) -> None:
            """given a stack of start positions, a stop position, and a frame number, this helper func
            will load the necessary tuple(s) into self.orfs"""
            # 'if startstack' already guaranteed
            # frame is already matched
            if lG:  # only get the longest orf 
                startCodonPos = max(startStack) # since stop comes before start, longest frag = max start
                self.orfs.append((-frame, stopPos, startCodonPos, startCodonPos-stopPos+1))
            else:
                for start in startStack:
                    self.orfs.append(-frame, stopPos, start, start-stopPos+1)

        # -----------
        # begin iterating through each nucleotide in the sequence
        # -----------
        for i in range(len(self.bottom)-2):
            codonR = self.bottom[i:i+3]
            frame = (len(self.bottom) - i) % 3 + 1

            # q (queueing start codon positions)
            if codonR in startCodons:
                match frame:
                    case 1:
                        hitStart1 = True
                        start1.append(i+3)  # ex. ATG (pos of G) frame 1
                    case 2:
                        hitStart2 = True
                        start2.append(i+3)  # ex. ATG (pos of G) frame 2
                    case 3:
                        hitStart3 = True
                        start3.append(i+3)  # ex. ATG (pos of G) frame 3
            
            # if s: d & s = [], x = i
            if codonR in stopCodons:
                match frame:
                    case 1:
                        if start1:
                            appendORF(start1, stopPos1, 1)  # adds tuples to self.orfs (if s: d)
                        start1 = []  # s = []
                        stopPos1 = i+1 # ex. TAA (pos of T) frame 1 (x = i)
                    case 2:
                        if start2:
                            appendORF(start2, stopPos2, 2)  # adds tuples to self.orfs (if s: d)
                        start2 = []  # s = []
                        stopPos2 = i+1  # ex. TAA (pos of T) frame 2 (x = i)
                    case 3:
                        if start3:
                            appendORF(start3, stopPos3, 3)  # adds tuples to self.orfs (if s: d)
                        start3 = []  # s = []
                        stopPos3 = i+1  #ex. TAA (pos of T) frame 3 (x = i)
        
        # in each frame if there are any orfs that would be appended but haven't bc eol was reached
        # before a stop codon (when orfs are appended) they are appended here
        if start1:
            appendORF(start1, stopPos1, 1)
        if start2:
            appendORF(start2, stopPos2, 2)
        if start3:
            appendORF(start3, stopPos3, 3)

        # handles cases where there are no start codons present at all in the gene, stop codons may
        # or may not be present (if not, stop codon pos is treated as beginning of gene candidate)
        if not hitStart1:
            appendORF([len(self.bottom)], stopPos1, 1)
        if not hitStart2:
            appendORF([len(self.bottom)], stopPos2, 2)
        if not hitStart3:
            appendORF([len(self.bottom)], stopPos3, 3)

    
    def _sortORFs(self, presorted=False):
        """sorts the ORFs by their lengths, they are already presorted by their frames"""
        if not presorted:
            self.orfs.sort(key=lambda x: x[0])
        self.orfs.sort(key=lambda x: x[3], reverse=True)


    def parse(self, startCodons=['ATG'], stopCodons=['TAA', 'TGA', 'TAG'], lG=False, doesReturn=False):
        """Parses the ORFfinder sequence as well as its reverse complement in place within self.orfs"""
        self._parseTop(startCodons, stopCodons, lG)
        self._sortORFs()
        self._parseBottom(startCodons, stopCodons, lG)
        self._sortORFs(presorted=True)

        if doesReturn:
            return self.orfs
   
def main(inFile = None, options = None):
    '''
    Find some genes.  
    '''
    thisCommandLine = CommandLine(options)
    reader = FastaReader(inFile)
    for fasta in reader.readFasta():
        orfFinder = ORFfinder(fasta[1])
        orfFinder.parse(
            startCodons=thisCommandLine.args.start, 
            stopCodons=thisCommandLine.args.stop,
            lG=thisCommandLine.args.longestGene
            )

        print(fasta[0])
        for frame in orfFinder.orfs:
            if frame[3] >= thisCommandLine.args.minGene:
                print(f'{frame[0]:+d} {frame[1]:>5d}..{frame[2]:>5d} {frame[3]:>5d}')
    
    ###### replace the code between comments.
    # thisCommandLine.args.longestGene is True if only the longest Gene is desired
    # thisCommandLine.args.start is a list of start codons
    # thisCommandLine.args.stop is a list of stop codons
    # thisCommandLine.args.minGene is the minimum Gene length to include
    #
    #######
    
if __name__ == "__main__":
    main() # delete this stuff if running from commandline
