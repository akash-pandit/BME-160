#!/usr/bin/env python3
# Name: Akash Pandit (aspandit)
# Group: None


"""
    Input for STDIN (in order to compile):
    python3 findUnique.py <bos-tRNA-7.fa >bos-tRNA-7-output.txt
"""


from sequenceAnalysis import FastAreader

class findUnique:
    def __init__(self):
        """Reads tRNA from the inputted Fasta file. Adds the power set of each sequence to a list"""
        self.pSets, self.uniques, self.headerDict = [], [], {}

    def genPSets(self, sequence: str) -> None:
        """Appends all unique tRNA substrings (a powerset)"""
        powerSets = set()
        for index in range(len(sequence)):
            size = len(sequence)
            while size > index:
                powerSets.add(sequence[index:size])
                size -= 1

        self.pSets.append(powerSets)

    def findUniques(self):
        """
        Looks for duplicate tRNA sequences and removes from all but the appriopiate set to
        obtain a unique set of tRNA subsequences.
        """
        for pSet in self.pSets:
            union = set()
            copyList = self.pSets.copy()
            copySet = pSet.copy()  # Make duplicate of the possible power set
            copyList.remove(copySet)  # Remove power set from list.

            for powerSet in copyList: # removes the duplicate elements between union and copySet
                union = union.union(powerSet)  # The union of all the other tRNA sets.

            copySet.difference_update(union)  # Update the copySet, removing elements found in union.

            newSet = copySet.copy()
            for s1 in copySet:  # checks prior strings to identify more possible substrings
                uniqueSet = copySet.copy()
                uniqueSet.remove(s1)
                for s2 in uniqueSet:
                    if s1 in s2:
                        if len(s1) < len(s2):
                            newSet.discard(s2)  # Gets rid of the larger substrings.

            self.uniques.append(newSet)

    def print(self) -> None:
        """Prints to STDOUT (CL in this case)"""
        for index in range(0,len(self.headerDict)):
            header, sequence = self.headerDict[index]
            print(header)
            print(sequence)
            for i in range(len(sequence)): 
                for substr in self.uniques[index]:
                    if substr == sequence[len(substr)]:
                        print('.'*i + substr)

def main():
    """Runs above code to find and print header, sequence and substrings from the fasta file"""

    # init objects & counter var
    tRNA = findUnique()
    fastaFile = FastAreader()
    count = 0

    # generate psets & increment counter
    for header, sequence in fastaFile.readFasta():
        filteredSequence = sequence.replace('-', '').replace('_', '')  # makes sure only valid characters are being used
        tRNA.headerDict[count] = [header, filteredSequence]
        tRNA.genPSets(filteredSequence)  # generate & append pset to self.psets
        count += 1

    # find & print unique seqs
    tRNA.findUniques()
    tRNA.print()


if __name__ == "__main__":
    main()