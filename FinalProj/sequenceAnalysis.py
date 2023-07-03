import sys
from numpy import Inf


class NucParams:
    rnaCodonTable = {
        # RNA codon table
        # U
        'UUU': 'F', 'UCU': 'S', 'UAU': 'Y', 'UGU': 'C',  # UxU
        'UUC': 'F', 'UCC': 'S', 'UAC': 'Y', 'UGC': 'C',  # UxC
        'UUA': 'L', 'UCA': 'S', 'UAA': '-', 'UGA': '-',  # UxA
        'UUG': 'L', 'UCG': 'S', 'UAG': '-', 'UGG': 'W',  # UxG
        # C
        'CUU': 'L', 'CCU': 'P', 'CAU': 'H', 'CGU': 'R',  # CxU
        'CUC': 'L', 'CCC': 'P', 'CAC': 'H', 'CGC': 'R',  # CxC
        'CUA': 'L', 'CCA': 'P', 'CAA': 'Q', 'CGA': 'R',  # CxA
        'CUG': 'L', 'CCG': 'P', 'CAG': 'Q', 'CGG': 'R',  # CxG
        # A
        'AUU': 'I', 'ACU': 'T', 'AAU': 'N', 'AGU': 'S',  # AxU
        'AUC': 'I', 'ACC': 'T', 'AAC': 'N', 'AGC': 'S',  # AxC
        'AUA': 'I', 'ACA': 'T', 'AAA': 'K', 'AGA': 'R',  # AxA
        'AUG': 'M', 'ACG': 'T', 'AAG': 'K', 'AGG': 'R',  # AxG
        # G
        'GUU': 'V', 'GCU': 'A', 'GAU': 'D', 'GGU': 'G',  # GxU
        'GUC': 'V', 'GCC': 'A', 'GAC': 'D', 'GGC': 'G',  # GxC
        'GUA': 'V', 'GCA': 'A', 'GAA': 'E', 'GGA': 'G',  # GxA
        'GUG': 'V', 'GCG': 'A', 'GAG': 'E', 'GGG': 'G'  # GxG
    }
    dnaCodonTable = {key.replace('U', 'T'): value for key, value in rnaCodonTable.items()}

    def __init__(self, inString='') -> None:
        """constructor, sets up instance vars nucStr (the string nucleic acid) and codonList 
        (nucStr broken into its codons) for use by later methods"""
        self.nucStr = inString.upper()
        self.codonList = [self.nucStr[base:base + 3].replace('T', 'U') for base in range(0, len(self.nucStr), 3)]
        # breaks inString into a list of codons, replacing all Ts with Us as output only wants RNA

    def addSequence(self, inSeq: str) -> None:
        """adds a sequence to both the nucleotide string and codon list by appending
        to the end of said string or list"""
        self.nucStr += inSeq.upper()
        newSeq = [inSeq[base:base+3] for base in range(0, len(inSeq), 3)]  # breaks inSeq into codons
        self.codonList.extend(newSeq)  

    def aaComposition(self) -> dict:
        """Returns a dictionary of amino acids and their counts in the nucleic acid"""
        if self.codonList:
            aaList = [self.rnaCodonTable[codon.replace('T', 'U')] for codon in self.codonList]
        return {aa: aaList.count(aa) for aa in set(aaList)}

    def nucComposition(self) -> dict:
        """returns a dictionary of nucleotides and their counts in the nucleic acid"""
        return {nuc: self.nucStr.count(nuc) for nuc in set(self.nucStr)}

    def codonComposition(self) -> dict:
        """returns a dictionary of rna codons and their counts"""
        codonList = [codon.replace('T', 'U') for codon in self.codonList]
        return {codon: self.codonList.count(codon.replace('U', 'T')) for codon in set(codonList) if self.validCodon(codon)}

    def nucCount(self) -> int:
        """returns the total number of valid nucleotides"""
        return sum(self.nucComposition().values())

    def validCodon(self, codon) -> bool:
        """helper method designed for codonComposition to validate a codon within a list comprehension"""
        for base in codon:
            if base not in 'ACGU':
                return False
        return True


class FastAreader:
    ''' 
    Define objects to read FastA files.
    
    instantiation: 
    thisReader = FastAreader ('testTiny.fa')
    usage:
    for head, seq in thisReader.readFasta():
        print (head,seq)
    '''

    def __init__(self, fname=None):
        '''contructor: saves attribute fname '''
        self.fname = fname

    def doOpen(self):
        ''' Handle file opens, allowing STDIN.'''
        if self.fname is None:
            return sys.stdin
        else:
            return open(self.fname)

    def readFasta(self):
        ''' Read an entire FastA record and return the sequence header/sequence'''
        header = ''
        sequence = ''

        with self.doOpen() as fileH:

            header = ''
            sequence = ''

            # skip to first fasta header
            line = fileH.readline()
            while not line.startswith('>'):
                line = fileH.readline()
            header = line[1:].rstrip()

            for line in fileH:
                if line.startswith('>'):
                    yield header, sequence
                    header = line[1:].rstrip()
                    sequence = ''
                else:
                    sequence += ''.join(line.rstrip().split()).upper()

        yield header, sequence


class ProteinParam:
    """class containing multiple data tables and methods for parsing strings of amino acids in their single character form"""
    # These tables are for calculating:
    #     molecular weight (aa2mw), along with the mol. weight of H2O (mwH2O)
    #     absorbance at 280 nm (aa2abs280)
    #     pKa of positively charged Amino Acids (aa2chargePos)
    #     pKa of negatively charged Amino acids (aa2chargeNeg)
    #     and the constants aaNterm and aaCterm for pKa of the respective termini
    #  Feel free to move these to appropriate methods as you like

    # As written, these are accessed as class attributes, for example:
    # ProteinParam.aa2mw['A'] or ProteinParam.mwH2O

    # defines each amino acid by single character symbol with its molecular weight in daltons or grams per mole
    aa2mw = {
        'A': 89.093, 'G': 75.067, 'M': 149.211, 'S': 105.093, 'C': 121.158,
        'H': 155.155, 'N': 132.118, 'T': 119.119, 'D': 133.103, 'I': 131.173,
        'P': 115.131, 'V': 117.146, 'E': 147.129, 'K': 146.188, 'Q': 146.145,
        'W': 204.225, 'F': 165.189, 'L': 131.173, 'R': 174.201, 'Y': 181.189
    }

    mwH2O = 18.015  # defines the molecular weight of water in daltons or grams per mole
    aa2abs280 = {'Y': 1490, 'W': 5500, 'C': 125}  # defines absorption values for Y, W,and C

    aa2chargePos = {'K': 10.5, 'R': 12.4, 'H': 6}  # defines the pKa of positively charged amino acids
    aa2chargeNeg = {'D': 3.86, 'E': 4.25, 'C': 8.33, 'Y': 10}  # defines the pKa of negatively charged amino acids
    aaNterm = 9.69  # pKa of N-terminus
    aaCterm = 2.34  # pKa of C-terminus

    def __init__(self, protein):
        """constructor for ProteinParam, initializes the input protein and its composition in a dict"""
        self.protein = ''.join([aa for aa in protein if aa in self.aa2mw.keys()]).upper()
        # defines the protein as the input string with all erroneous characters removed
        self.aaComposition_ = {aa: self.protein.count(aa) for aa in self.aa2mw}
        # generates dict of key:value pairs in aa2mw if those keys are in proteins

    def aaCount(self) -> int:
        """gets the length of the purified protein parameter"""
        return len(self.protein)

    def pI(self):
        """finds the isoelectric value for a given pH. """
        return_pH = -1
        isoelectric = Inf
        for pH in range(0, 1401):
            pH /= 100  # brings it down to the scale of 0.00 to 14.00
            newCharge = abs(self._charge_(pH))  # we want closest to 0, so measure distance from 0 with abs()
            if newCharge < isoelectric:
                isoelectric = newCharge
                return_pH = pH
        if return_pH >= 0:
            return return_pH
        raise AttributeError("For some reason, return_pH was never executed. Investigate this.")

    def aaComposition(self) -> dict:
        """returns a dictionary of the amino acid composition of the protein with 1 letter shortened amino acids as keys
        and their molecular weight as their value"""
        return self.aaComposition_

    def _charge_(self, pH):
        """Implements the equation to calculate the net charge on an amino acid"""
        posCharge, negCharge = 0, 0

        for aa, pKa in self.aa2chargePos.items():  # Calculate positive charge
            posCharge += self.aaComposition_[aa] * (10 ** pKa) / (10 ** pKa + 10 ** pH)

        for aa, pKa in self.aa2chargeNeg.items():  # Calculate negative charge
            negCharge += self.aaComposition_[aa] * (10 ** pH) / (10 ** pKa + 10 ** pH)

        posCharge += 10 ** self.aaNterm / (10 ** self.aaNterm + 10 ** pH)  # increments n terminus charge
        negCharge += 10 ** pH / (10 ** self.aaCterm + 10 ** pH)  # increments c terminus charge

        return posCharge - negCharge  # returns net charge

    def molarExtinction(self, oxidizing=True) -> float:
        """Implements the molar extinction equation, with a base True case to account for oxidative conditions and if false, reductive conditions"""
        molarExt = 0
        for aa in self.aa2abs280.keys():
            if oxidizing and aa != 'C':
                molarExt += self.aaComposition()[aa] * self.aa2abs280[aa]
            elif not oxidizing:
                molarExt += self.aaComposition()[aa] * self.aa2abs280[aa]

        return molarExt

    def massExtinction(self, oxidizing=True) -> float:
        """divides the molar extinction by molecular weight to get the mass extinction"""
        myMW = self.molecularWeight()
        return self.molarExtinction(oxidizing=oxidizing) / myMW if myMW else 0.0

    def molecularWeight(self) -> float:
        """
        Formula: MW(H2O) + Sum for all inputted amino acids(N(aa) * (MW(aa) - MW(H2O)))
        """
        molWeight = self.mwH2O  # initialize mol weight w/ 1 water
        for aa in self.protein:
            molWeight += (self.aa2mw[aa] - self.mwH2O)
            # adds the difference of the amino acid mol. weight and the mol. weight of water to molWeight
        return molWeight


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
        "G": "C",
        "N": "N"
    }

    def __init__(self, topStrand):
        """constructor for ORFfinder objects, intializes a top (base) strand of
        dna and a bottom (reverse complement)"""
        topStrand = topStrand.upper().replace('U', 'T')
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
                    self.orfs.append((-frame, stopPos, start, start-stopPos+1))

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
   

def main(fileName=None):
    """Driver file- generates NucParams object, parses the data, and formats it nicely in output file"""
    myReader = FastAreader(fileName)
    myNuc = NucParams()
    for head, seq in myReader.readFasta():
        myNuc.addSequence(seq)

    # sort codons in alpha order, by Amino Acid
    aaComp = myNuc.aaComposition()
    codonComp = myNuc.codonComposition()
    nucs = [(codon, myNuc.rnaCodonTable[codon], codonComp[codon]/aaComp[myNuc.rnaCodonTable[codon]])
            for codon in codonComp.keys()]
    nucs.sort(key=lambda x: x[0])  # first sort by codons
    nucs.sort(key=lambda x: x[1])  # then sort by amino acids
    # amino acid sorting will override codon sorting, but for lines with same aa, codon sorting will remain

    # calculate sequence length and gc usage
    seqLen = myNuc.nucCount() / 1000000  # 1e6 bytes / Mb
    gcUsage = (myNuc.nucComposition()['G'] + myNuc.nucComposition()['C']) / myNuc.nucCount()

    # write everything to output
    with open('output.out', 'w') as out:
        # write header to file
        out.write("sequence length = " + str(round(seqLen, 2)) + " Mb\n\n")
        out.write("GC content = " + str(round(gcUsage*100, 1)) + "%\n\n")
        
        # loop through nucs and output each codon's data to output file
        for nucI in nucs:
            codon, aa, val = nucI[0], nucI[1], nucI[2]
            out.write(f'{codon} : {aa} {val*100:5.1f} ({codonComp[codon]:6d})\n')

    # NOTE: as this is supposed to be incorporated into my "personal package", output only printed to file, 
    # not stdout


if __name__ == "__main__":
    main("testGenome.fa")  # no param = use stdin
