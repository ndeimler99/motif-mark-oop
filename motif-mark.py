#!/usr/bin/env python

import argparse
import cairo
import re
import itertools

OFFSET = 50
COLORS = [[255, 0, 0], [0, 0, 255], [0, 255, 0], [213, 0, 255], [0, 255, 230], [255, 255, 0], [129, 95, 55]]

class Gene:
    '''This class is used for storing the information about each read in the fasta file; primarily its length and number.
    Number Corresponds to which number gene it is in the fasta file.'''

    def __init__(self, start, length, number):
        '''initiate gene object. start will always be 0..., but whatever. Length is the length of the sequence in bases (correlates to pixels);
        this function stores necesary info in class variables'''
        self.start = start
        self.length = length
        self.number = number
    
    def draw_gene(self, context):
        '''takes class variables and plots gene. One base correlates to one pixel in length'''
        context.set_source_rgb(0,0,0)
        context.move_to(self.start + OFFSET, self.number * 150 - 10)
        context.line_to(self.length + OFFSET, self.number * 150 - 10)
        context.stroke()
        pass

class Exon:
    '''This class is very similar to Gene class.  However, this takes in the actual sequence and determines where the exon is within that sequence
    It will then plot a rectangle at that position with a width the same length of the exon sequence'''

    def __init__(self, sequence, thickness, number):
        '''Assign and initiate class variables'''
        self.sequence = sequence
        self.thickness = thickness
        self.number = number

    def determine_start_pos(self):
        '''split the sequence by lower and uppercase. We know each gene contains one section of uppercase - this is the exon
        The start position of the exon will be the length of the pre-exon segment.'''
        sequences = re.split(r"([A-Z]+)", self.sequence)
        start = len(sequences[0]) + 1
        length = len(sequences[1])
        return start, length

    def draw_exon(self, context):
        '''using class variables draw the exon'''
        context.set_source_rgb(0,0,0)
        start, length = self.determine_start_pos()
        context.rectangle(start + OFFSET, self.number * 150 - (self.thickness/2 + 10), length, self.thickness)
        context.fill()

class FastaHeader:
    '''This class is similar to the above two classes. It takes in necessary information and draws a header'''

    def __init__(self, text, number):
        '''assign and initate class variables.'''
        self.header = text
        self.number = number

    def draw_header(self, context):
        '''draw header'''
        context.move_to(OFFSET + 10, self.number * 150 - 110)
        context.show_text(self.header)
        context.stroke()

class Motif:
    '''This class is a bit more complicated as it has to calculate the starting position of each motif within the sequence'''

    def __init__(self, sequence, motif, motif_list, number, gene_number):
        '''sequence is the entire sequene of pre, exon, and post-exon regions; motif is the actual motif sequence; 
        motif_list is a list of the possible motifs due to ambiguity; number is the motif number depending on the number of motifs in motif file;
        gene_number is the gene number infasta file we are on'''
    
        self.sequence = sequence
        self.motif = motif
        self.motif_list = motif_list
        self.gene_number = gene_number
        self.number = number
        
        

    def get_locations(self, sequence):
        '''return the start locations for motif object'''
        locations = []
        sequence = sequence.upper()
        for m in self.motif_list:
            for i in range(0, len(sequence)-len(self.motif) + 1):
                #what is the sequence of the sliding window
                check = sequence[i: i + len(self.motif)]
                #is the sliding window sequence the same sequence as the motif we are checking
                if check == m:
                    #if it is append the start position of the motif
                    locations.append(i)
        # #return a list of starting locations
        return locations

    def draw_motif(self, context):
        '''draw all motifs in this object'''
        locations = self.get_locations(self.sequence)
        for local in locations:
            context.rectangle(local + OFFSET, 150 * self.gene_number - 10 - 5 * self.number, len(self.motif), 5)
        context.set_source_rgb(COLORS[self.number][0], COLORS[self.number][1], COLORS[self.number][2])
        context.fill()
        #plot rectangles based off of height number, starting locations, and length of motif
        pass

class Legend:
    '''this class will draw and store legend information. In reality a better design would be storing motifs as dictionary here with the colors as keys and having a getColor method'''
    def __init__(self, motif):
        self.motif = motif

    def draw(self, context):
        context.set_source_rgb(0,0,0)
        context.move_to(20, 20)
        context.show_text("Legend")
        n = 0
        x=20
        for m in self.motif:
            context.set_source_rgb(COLORS[n][0], COLORS[n][1], COLORS[n][2])
            context.rectangle(x, 50, 20, 10)
            context.fill()
            context.set_source_rgb(0, 0, 0)
            context.move_to(x + 25, 60)
            context.show_text(m)
            context.fill()
            x += 80        
            n += 1

class geneGroup:                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      
    '''This is a wrapper class, passes off information to other classes as needed'''
    def __init__(self, sequence, header, number, motif):
        self.sequence = sequence
        self.number = number
        self.gene = Gene(0, len(sequence), number)
        self.exon = Exon(sequence, 15, number)
        self.fastaHeader = FastaHeader(header, number)
        self.motif = motif

    def motif_finder(self, motif, context):
        i=0
        for m in motif:
            new_motif = Motif(self.sequence, m, motif[m], i, self.number)
            new_motif.draw_motif(context)
            i+=1
            
    def draw(self, context):
        gene = self.gene
        exon = self.exon
        fastaHeader = self.fastaHeader
        gene.draw_gene(context)
        exon.draw_exon(context)
        fastaHeader.draw_header(context)
        self.motif_finder(self.motif, context)

        


#gene class - gene object stores start, length, line thickeness
#exon class - exon start, length, and line thickness
#motif class - start, length, width of motif, color of motif
#FastaHeader class - start pos, text, font?
#gene group class composed of other ojbects, one gene object, one exon object, and a group of motif objects and a fasta header object and rank (y pos in final figure)

#argparse function to take user input
def get_args():
    """argparse function that takes user input"""
    parser = argparse.ArgumentParser(description = "motif-mark: for easy visualization of the presence of motifs in your sequences")
    #action = append will take any duplicated -f flag and append it to the end of a list, allows for multiple input fasta files
    parser.add_argument("-f", "--files", help = "what is your input fasta files?", required = True, action="append")
    #a txt file containing motifs, each motif must be on its own line
    parser.add_argument("-m", "--motifs", help="What are your motifs? Each motif should appear on a new line", required=True)
    parser.add_argument("-n", "--number", help="What is the number of genes you would like to plot on each plot?", required=True)
    return parser.parse_args()

#assign user input to global script variables
args = get_args()
fasta_files = args.files
motif_file = args.motifs
number_of_genes = int(args.number)

def parser(fh, context):
    '''docstring goes here'''
    
    #initialize an empty list that will hold the genegroup objects
    #open the file passed to the method
    with open(fh, "r") as file:
        #create empty header and sequence strings
        header = ""
        sequence = ""
        gene_number = 2
        #for each line in the file
        for line in file:
            if gene_number > number_of_genes + 1:
                print("Maximum Number of Genes Exceeded, please fix Command Call")
                return False
            #remove the newline character
            line = line.strip()
            #if line is a header line it would start with a ">"
            if line.startswith(">"):
                #if the len(header) > 0 then we know that we are not on the first entry so a sequence exists
                #the very first line would be a header line, but we do not want to split the sequence, as we have not yet reached the sequence
                if len(header) > 0:
                    if "U" in sequence.upper():
                        gene = geneGroup(sequence, header, int(gene_number), RNA)
                    else:
                        gene = geneGroup(sequence, header, int(gene_number), DNA)
                    gene.draw(context)

                    gene_number+=1
                #assign the header line to header
                header = line
                #reset the sequence
                sequence = ""
            #else the line is not a header
            else:
                #concatenate the new line to the last line
                sequence += line 
        #split and assign the very last sequence in the file
        #this is necessary because the last sequence split will not be triggered by a new header
        if "U" in sequence.upper():
            gene = geneGroup(sequence, header, int(gene_number), RNA)
        else:
            gene = geneGroup(sequence, header, int(gene_number), DNA)
        gene.draw(context)
    return True

#In the case of ambiguous codes the value is a list of possible nucleotides
iupac_conversions = {"A":"A", "T":"T", "U":"T", "C":"C", "G":"G", "R":["A","G"], "Y":["C", "T"], "S":["G", "C"], "W":["A", "T"], "K":["G", "T"], "M":["A", "C"],"B":["C", "G", "T"],"D":["A", "G", "T"], "H":["A", "C", "T"], "V":["A", "C", "G"],"N":["A", "C", "G", "T"]}

#motif_maker takes the input of a .txt file that contains one motif per line (while this program will technically work with as many motifs as you want, the default color scheme allows for ten colors)
def motif_maker(motif_file):
    '''docstring goes here'''
    #create two empty dictionaries, one for DNA and one for RNA motifs
    #the keys will be the actual motif sequence with the value of all possible motif sequences (after adjusting them for ambiguous nucleotides)    
    DNA = {}
    RNA = {}

    #open the motif file
    with open(motif_file, "r") as fh:
        #go through the file line by line
        for line in fh:
            #remove the newline character
            line = line.strip()
            #make the motif sequence uppercase
            line = line.upper()
            #create new dictionary entries for each motif
            DNA[line] = []
            RNA[line] = []

    #for every motif in the dictionary
    for motif in DNA:
        #intialize empty temporary lists
        new_dna_motif=[]
        new_rna_motif=[]
        
        #use some fancy list comprehension
        #for letter in motif, create an itertools.product permutation replacing that letter with each other possible nucleotide
        for i in itertools.product(*[iupac_conversions[j] for j in motif]):
            #by default iupac_conversions dictionary creates DNA motif sequences
            new_dna_motif.append("".join(i))
            #replace every thymine with uracil for the RNA motif sequences
            new_rna_motif.append("".join(i).replace("T", "U"))

        #add all possible motif sequences into respective dictionaries
        DNA[motif] = new_dna_motif
        RNA[motif] = new_rna_motif
    
    #return dictionaries
    return DNA, RNA

DNA, RNA = motif_maker(motif_file)

#for each fasta file run the program on it
for file in fasta_files:
    #extract file name for use in naming plot
    figname = re.findall('(.*)\.fa', file)[0] + ".svg"

    #create a pycairo surface
    surface = cairo.SVGSurface(figname, 1000, 150*(number_of_genes+1))

    #create context on the pycairo surface
    context = cairo.Context(surface)
    geneGroups = parser(file, context)
    legend = Legend(DNA)
    legend.draw(context)
    context.stroke()
    surface.finish()
    
