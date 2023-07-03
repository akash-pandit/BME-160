from sequenceAnalysis import ORFfinder, FastAreader
import os
import matplotlib.pyplot as plt

class CLReader:
    """
    A custom Command Line Reader for this program.

    Arguments:

    -t (--target) the target file of comparison

    -c1 (--compareTo1) one of 2 files the target gets compared with 

    -c2 (--compareTo2) the other of 2 files the target gets compared with

    -g (--graph) optional flag for if you want a visual representation of matching orfs
    """
    def __init__(self):
        import argparse

        self.p = argparse.ArgumentParser()

        self.p.add_argument('-t', '--target', help="Provide the target FASTA file path")
        self.p.add_argument('-c1', '--compareTo1', 
                            help="Provide the first FASTA file for ORF comparison with the target file")
        self.p.add_argument('-c2', '--compareTo2', 
                            help="Provide the second FASTA file for ORF comparison with the target file")
        self.p.add_argument('-g', '--graph', help="Visualize the data that led to the given output", action='store_true')

        self.args = self.p.parse_args()

        # handles missing arg -t / --target
        if not self.args.target:
            print("\nExecution Halted:\nPlease provide a FASTA to compare with. (-t file.fa)")
            exit()

        # handles invalid filepath for arg -t / --target
        if not os.path.exists(self.args.target):
            print("\nExecution Halted:\nPlease provide an existing filepath to argument -t.")
            exit()

        # handles missing arg -c1 / --compareTo1
        if not self.args.compareTo1:
            print("\nExecution Halted:\nPlease provide 2 FASTA files for comparing the target file to. (-c1 file.fa for first file)")
            exit()

        # handles invalid filepath for arg -c1 / --compareTo1
        if not os.path.exists(self.args.compareTo1):
            print("\nExecution Halted:\nPlease provide an existing filepath to argument -c1.")
            exit()

        # handles missing arg -c2 / --compareTo2
        if not self.args.compareTo2:
            print("\nExecution Halted:\nPlease provide 2 FASTA files for comparing the target file to. (-c2 file.fa for second file)")
            exit()

        # handles invalid filepath for arg -c2 / --compareTo2
        if not os.path.exists(self.args.compareTo2):
            print("\nExecution Halted:\nPlease provide an existing filepath to argument -c2.")
            exit()
            

class Comparison:
    def __init__(self, targetObj: ORFfinder, c1_Obj: ORFfinder, c2_Obj: ORFfinder, 
                 startList=['ATG'], stopList=['TAA', 'TAG', 'TGA']) -> None:
        """Create ORFfinder objects based on the input sequences from their respective files"""

        self.target = targetObj
        self.c1, self.c2 = c1_Obj, c2_Obj
        self.c1counter, self.c2counter = 0, 0
        self.startList, self.stopList = startList, stopList

        self.target.parse()
        self.c1.parse()
        self.c2.parse()

    def _getSetLens(self) -> tuple:
        """Helper method that returns sets of lengths of.orfs that could possibly match to another file's lengths"""

        # get orf lens for target file
        lenTarget = set(orf[3] for orf in self.target.orfs)
        # get orf lens for c1 file only if also in lenTarget
        lenC1 = set(orf[3] for orf in self.c1.orfs if orf[3] in lenTarget)
        # get orf lens for c2 file only if also in lenTarget
        lenC2 = set(orf[3] for orf in self.c1.orfs if orf[3] in lenTarget)
        # only keep lenTarget lens if in either lenC1 or lenC2
        lenTarget = set(orfLen for orfLen in lenTarget if orfLen in lenC1.union(lenC2))

        # return a tuple of the 3 sets
        return lenTarget, lenC1, lenC2
    
    def _getORFSeq(self, orfTuple: tuple, fullSeq: str) -> str:
        """Helper method that returns the actual sequence of the ORF"""

        # fullSeq = ORFfinder.top (str of bases)
        # orfTuple -> (frame, startPos, endPos, len)
        frame, startPos, endPos = orfTuple[0], orfTuple[1], orfTuple[2]

        if frame > 0:
            # handles getting the ORF from the base strand
            return fullSeq[startPos-1:endPos]
        if frame < 0:
            # handles getting the ORF from the reverse complement
            returnSeq = fullSeq[startPos-1:endPos:][::-1]
            return ''.join(ORFfinder.complement[i] for i in returnSeq)

    def _genDicts(self, targSet: set, setC1: set, setC2: set) -> tuple:
        """generate and populate dictionaries with keys being unique lengths in a respective
        ORFfinder object and """

        # generate dictionaries with empty sets for all unique lengths in each file's set of lengths
        dictTarget = {length : set() for length in targSet}
        dictC1 = {length : set() for length in setC1}
        dictC2 = {length : set() for length in setC2}

        # define helper func that adds the orflens from a given ORFfinder to each len's respective set
        # in the given dictFile
        def appendORF(orfFinder: ORFfinder, dictFile: dict) -> None:
            """helper function which just shortens repeated code"""

            def notGeneFrag(seq) -> bool:
                """helper function which returns whether or not the orf has a start & stop codon
                makes sure gene fragments aren't included"""
                hasStart, hasStop = False, False

                for codon in self.startList:
                    if codon in seq:
                        hasStart = True
                for codon in self.stopList:
                    if codon in seq:
                        hasStop = True
                return hasStart and hasStop
                
            # go through each orf
            for orf in orfFinder.orfs:
                orfLen = orf[3]
                orfSeq = self._getORFSeq(orf, orfFinder.top)
                if notGeneFrag(orfSeq):
                    
                    if dictFile.get(orfLen) != None:
                        dictFile[orfLen].add(orfSeq)

        appendORF(self.target, dictTarget)
        appendORF(self.c1, dictC1)
        appendORF(self.c2, dictC2)

        # return the 3 populated dicts in a tuple
        return dictTarget, dictC1, dictC2
    
    def compareFiles(self) -> None:
        """compares the ORFs of each comparison file to the target file and updates the object
        field incrementers for each matching orf left"""
        lenTargetSet, lenC1Set, lenC2Set = self._getSetLens()
        dTarget, dC1, dC2 = self._genDicts(lenTargetSet, lenC1Set, lenC2Set)
        self.dC1Graph = {length : len(dC1[length]) for length in dC1.keys()}
        self.dC2Graph = {length: len(dC2[length]) for length in dC2.keys()}

        for length in dTarget.keys():
            # generate sets of sequences for that particular length
            targetSeqSet, c1SeqSet, c2SeqSet = dTarget.get(length), dC1.get(length), dC2.get(length)
            # iterate through the target sequences
            for seq in targetSeqSet:
                # check first if the file sequence set exists then check if target
                # sequence is in that set of sequences (they have it in common)
                if c1SeqSet != None and seq in c1SeqSet:
                    # if so increment the counter
                    self.c1counter += 1
                    self.dC1Graph[length] += 1

                # same as above here but for c2
                if c2SeqSet != None and seq in c2SeqSet:
                    self.c2counter += 1
                    self.dC2Graph[length] += 1

    def genGraph(self) -> None:
        """generate a set of 3 subplots, 1 histogram displaying the different counter values and
        2 scatterplots for looking at distribution of matches"""
        def graphPipeline(sub, vals, freqs, xlabel, ylabel, title, type='scatter'):
            """coagulate the actual generation of a particular graph"""
            if type == 'scatter':
                axs[sub].scatter(vals, freqs)
            elif type == 'hist':
                axs[sub].hist(vals)
            axs[sub].set_title(title)
            axs[sub].set_xlabel(xlabel)
            axs[sub].set_ylabel(ylabel)

        # generate a figure with 1 subplot tall and 3 wide
        axs = plt.subplots(1, 3, figsize=(10, 4))[1]

        values = []
        for _ in range(self.c1counter):
            values.append('C1')
        for _ in range(self.c2counter):
            values.append('C2')
        graphPipeline(0, values, None, "File", "Frequency", "ORF Matches By Comparison (C1/C2) File", type='hist')

        values = self.dC1Graph.keys()
        frequencies = self.dC1Graph.values()
        graphPipeline(1, values, frequencies, 'ORF Lengths', "Num of ORFs w/ Length", "Target - C1 ORF Matches")

        values = self.dC2Graph.keys()
        frequencies = self.dC2Graph.values()
        graphPipeline(2, values, frequencies, 'ORF Lengths', 'Num of ORFs w/ Length', "Target - C2 ORF Matches")

        plt.tight_layout()
        plt.show()
    

def main():
    """Execution of compare.py through command line"""
    cl = CLReader()
    targetName, c1Name, c2Name = cl.args.target, cl.args.compareTo1, cl.args.compareTo2
    targetFile, c1File, c2File = FastAreader(targetName), FastAreader(c1Name), FastAreader(c2Name)

    targetSeq = ''.join([seq[1] for seq in targetFile.readFasta()])
    c1Seq = ''.join([seq[1] for seq in c1File.readFasta()])
    c2Seq = ''.join([seq[1] for seq in c2File.readFasta()])
    # [0] -> gene title sequence, [1] -> the actual sequence

    target, c1, c2 = ORFfinder(targetSeq), ORFfinder(c1Seq), ORFfinder(c2Seq)
    comparer = Comparison(target, c1, c2)
    comparer.compareFiles()
    if comparer.c1counter > comparer.c2counter:
        print(f"{targetName} is likely to be more closely related to {c1Name}")
        print(f"matches found: {comparer.c1counter} ORFs matched")
        print(f"matches for {c2Name}: {comparer.c2counter}")
    if comparer.c2counter > comparer.c1counter:
        print(f"{targetName} is likely to be more closely related to {c2Name}")
        print(f"matches found: {comparer.c2counter} ORFs matched")
        print(f"matches for {c1Name}: {comparer.c2counter}")

    if cl.args.graph:
        comparer.genGraph()


if __name__ == '__main__':
    """File usage:
    py compare.py -t fastaT -c1 fastaA    -c2 fastaB   -g
                  --target  --compareTo1  --compareTo2"""
    main()