import sys
from numpy import Inf


# STATUS: Not Started
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
    dnaCodonTable = {key.replace('U','T'):value for key, value in rnaCodonTable.items()}

    def __init__ (self, inString=''):
        self.nucString = inString.upper()
        self.codonList = [inString[i:i + 3] for i in range(0, len(inString), 3)]

        if 'T' in inString:
            self.isDNA = True
        else:
            self.isDNA = False
        

    def addSequence (self, inSeq: str) -> None:
        self.codonList.insert(0, inSeq)
        

    def aaComposition(self) -> dict:
        """Returns """
        if self.isDNA:
            aaList = [self.dnaCodonTable[codon] for codon in self.codonList]
        else:
            aaList = [self.rnaCodonTable[codon] for codon in self.codonList]

        return {aa: self.codonList.count(aa) for aa in aaList}


    def nucComposition(self):
        if self.isDNA:
            return {nuc:self.nucString.count(nuc) for nuc in self.nucString if nuc in "ACGTN"}
        return {nuc:self.nucString.count(nuc) for nuc in self.nucString if nuc in "ACGUN"}
    

    def codonComposition(self):
        cleanCodonList = [codon for codon in self.codonList if ]
        return self.codonComp
    

    def nucCount(self):
        return sum(self.nucComp.values())
    

# STATUS: Not Started
class FastAreader :
    ''' 
    Define objects to read FastA files.
    
    instantiation: 
    thisReader = FastAreader ('testTiny.fa')
    usage:
    for head, seq in thisReader.readFasta():
        print (head,seq)
    '''
    def __init__ (self, fname=None):
        '''contructor: saves attribute fname '''
        self.fname = fname
            
    def doOpen (self):
        ''' Handle file opens, allowing STDIN.'''
        if self.fname is None:
            return sys.stdin
        else:
            return open(self.fname)
        
    def readFasta (self):
        ''' Read an entire FastA record and return the sequence header/sequence'''
        header = ''
        sequence = ''
        
        with self.doOpen() as fileH:
            
            header = ''
            sequence = ''
            
            # skip to first fasta header
            line = fileH.readline()
            while not line.startswith('>') :
                line = fileH.readline()
            header = line[1:].rstrip()

            for line in fileH:
                if line.startswith ('>'):
                    yield header,sequence
                    header = line[1:].rstrip()
                    sequence = ''
                else :
                    sequence += ''.join(line.rstrip().split()).upper()

        yield header,sequence


# Done Previously
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

    aa2mw = {  # defines each amino acid by single character symbol with its molecular weight in daltons or grams per mole
        'A': 89.093,  'G': 75.067,  'M': 149.211, 'S': 105.093, 'C': 121.158,
        'H': 155.155, 'N': 132.118, 'T': 119.119, 'D': 133.103, 'I': 131.173,
        'P': 115.131, 'V': 117.146, 'E': 147.129, 'K': 146.188, 'Q': 146.145,
        'W': 204.225,  'F': 165.189, 'L': 131.173, 'R': 174.201, 'Y': 181.189
        }

    mwH2O = 18.015  # defines the molecular weight of water in daltons or grams per mole
    aa2abs280= {'Y':1490, 'W': 5500, 'C': 125}  # defines absorption values for Y, W,and C

    aa2chargePos = {'K': 10.5, 'R':12.4, 'H':6}  # defines the pKa of positively charged amino acids
    aa2chargeNeg = {'D': 3.86, 'E': 4.25, 'C': 8.33, 'Y': 10}  # defines the pKa of negatively charged amino acids
    aaNterm = 9.69  # pKa of N-terminus
    aaCterm = 2.34  # pKa of C-terminus

    def __init__ (self, protein):
        """constructor for ProteinParam, initializes the input protein and its composition in a dict"""
        self.protein = ''.join([aa for aa in protein if aa in self.aa2mw.keys()]).upper()
        # defines the protein as the input string with all erroneous characters removed
        self.aaComposition_ = {aa : self.protein.count(aa) for aa in self.aa2mw}
        # generates dict of key:value pairs in aa2mw if those keys are in proteins

    def aaCount(self) -> int:
        """gets the length of the purified protein parameter"""
        return len(self.protein)

    def pI (self):
        """finds the isoelectric value for a given pH. """
        isoelectric = Inf
        for pH in range(0, 1401):
            pH /= 100  # brings it down to the scale of 0.00 to 14.00
            newCharge = abs(self._charge_(pH))  # we want closest to 0, so measure distance from 0 with abs()
            if newCharge < isoelectric: 
                isoelectric = newCharge
                return_pH = pH
        if return_pH:
            return return_pH
        raise AttributeError("For some reason, return_pH was never executed. Investigate this.")

    def aaComposition(self) -> dict:
        """returns a dictionary of the amino acid composition of the protein with 1 letter shortened amino acids as keys
        and their molecular weight as their value"""
        return self.aaComposition_

    def _charge_ (self, pH):
        """Implements the equation to calculate the net charge on an amino acid"""
        posCharge, negCharge = 0, 0

        for aa,pKa in self.aa2chargePos.items(): # Calculate positive charge
            posCharge += self.aaComposition_[aa] * (10**pKa) / (10**pKa + 10**pH)
            
        for aa,pKa in self.aa2chargeNeg.items(): # Calculate negative charge
            negCharge += self.aaComposition_[aa] * (10**pH) / (10**pKa + 10**pH)
        
        posCharge += 10**self.aaNterm / (10**self.aaNterm + 10**pH) # increments n terminus charge
        negCharge += 10**pH / (10**self.aaCterm + 10**pH) # increments c terminus charge
        
        return posCharge - negCharge # returns net charge


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
    

def main (fileName=None):
    myReader = FastAreader(fileName) 
    myNuc = NucParams()
    for head, seq in myReader.readFasta() :
        myNuc.addSequence(seq)
        
    # sort codons in alpha order, by Amino Acid
    
    # calculate relative codon usage for each codon and print
    for nucI in nucs:
        ...
        print ('{:s} : {:s} {:5.1f} ({:6d})'.format(nuc, aa, val*100, thisCodonComp[nuc]))

if __name__ == "__main__":
    main() # no param = use stdin
    