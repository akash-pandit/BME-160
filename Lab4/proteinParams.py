#!/usr/bin/env python3
# Name: Akash Pandit
# Group Members: Preeti Rajarathinam(prrajara), Delson Hays (dhhays), Natalie Cellucci(ncellucc), Riya Tiloda (rtiloda)

from numpy import inf  # used in ProteinParam.pI()

"""
This program defines a class ProteinParam which is built around a string of amino acids represented as single characters.
Many operations can be performed on these, as listed by the below methods

self.aaCount() - find the number of amino acids in the input string.

self.pI() - finds the isoelectric point, the pH value where the protein is most electrically neutral.

self.aaComposition() - returns a dictionary mapping each unique amino acid in the protein to its count.

self._charge_(pH) - an internal method that calculates the total electric charge of the protein at a specific pH.

self.molarExtinction(oxidizing) - calculates the molar extinction coefficient for the protein. In an oxidizing environment, 
    cysteine is left out of the calculations as its present in cystine residues which do not contribute. In a reducing 
    environment, this does not occur and cysteine remains in the molar extinction coefficient calculations.

self.massExtinction(oxidizing) - divides the molar extinction value from self.molarExtinction by the molar mass to get the
    protein's mass extinction value. The oxidizing parameter value is passed to self.molarExtinction().

self.molecularWeight() - calculates the molecular weight of the protein.

SAMPLE INPUT:
VLSPADKTNVKAAW

SAMPLE OUTPUT:
Number of Amino Acids: 14
Molecular Weight: 1499.7
molar Extinction coefficient: 5500.00
mass Extinction coefficient: 3.67
Theoretical pI: 9.88
Amino acid composition:
A = 21.43%
C = 0.00%
D = 7.14%
E = 0.00%
F = 0.00%
G = 0.00%
H = 0.00%
I = 0.00%
K = 14.29%
L = 7.14%
M = 0.00%
N = 7.14%
P = 7.14%
Q = 0.00%
R = 0.00%
S = 7.14%
T = 7.14%
V = 14.29%
W = 7.14%
Y = 0.00%

This program also includes a driver function main() to build a ProteinParam object and output
the results of the object's methods.
"""

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
        isoelectric = inf
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

# Please do not modify any of the following.  This will produce a standard output that can be parsed
    
import sys
def main():
    inString = input('protein sequence?')
    while inString :
        myParamMaker = ProteinParam(inString)
        myAAnumber = myParamMaker.aaCount()

        print ("Number of Amino Acids: {aaNum}".format(aaNum = myAAnumber))
        print ("Molecular Weight: {:.1f}".format(myParamMaker.molecularWeight()))
        print ("molar Extinction coefficient: {:.2f}".format(myParamMaker.molarExtinction()))
        print ("mass Extinction coefficient: {:.2f}".format(myParamMaker.massExtinction()))
        print ("Theoretical pI: {:.2f}".format(myParamMaker.pI()))
        print ("Amino acid composition:")
        
        if myAAnumber == 0 : myAAnumber = 1  # handles the case where no AA are present 
        
        for aa,n in sorted(myParamMaker.aaComposition().items(), key= lambda item:item[0]):
            print ("\t{} = {:.2%}".format(aa, n/myAAnumber))
    
        inString = input('protein sequence?')

if __name__ == "__main__":
    main()