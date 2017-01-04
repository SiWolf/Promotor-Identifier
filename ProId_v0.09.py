# ----------------------------------------------------------------------------------------------------------------------------------
# Author: Silver A. Wolf
# Last Modified: Mi, 04.01.2017
# Title: In silico Promotor Detektion von Oberflächen-assoziierten Faktoren und ausgewählten Resistenzgenen in Staphylococcus aureus
# Version: 0.0.9a for Python 3.5.2 (Windows)
# ----------------------------------------------------------------------------------------------------------------------------------

#Imports:
from Bio import SeqIO
from Bio.Blast import NCBIXML
from fnmatch import fnmatch, fnmatchcase
import os
import platform
import re

#Functions:
def fetch_p_region(genome, gene_name, gene_begin, gene_end, gene_strand, contig_name):
    #Gibt die Promoter-Region und Gensequenz nach Eingabe der Blast-Informationen zurück
    fasta_file = genome+".fasta"
    fasta_path = os.path.join(genomeDir, fasta_file)

    if p_region_mode == "short":
        p_length = 201
    else:
        p_length = 551

    for contig in SeqIO.parse(fasta_path, "fasta"):
        if contig.id == contig_name:
            if gene_strand == "1":
                gene = contig.seq[gene_begin-1:gene_end]
                if gene_begin < p_length:
                    promoter = contig.seq[0:gene_begin-1]
                    print("Warnung - Genom " + genome + ": Gen " + gene_name + " liegt am Rande des Contigs, die Promoterregion wurde teilweise abgeschnitten. (Fehlercode 001)")
                else:
                    promoter = contig.seq[gene_begin-p_length:gene_begin-1]
            else:
                gene = contig.seq[gene_end-1:gene_begin].reverse_complement()
                if gene_end < p_length:
                    promoter = contig.seq[0:gene_end-1].reverse_complement()
                    print("Warnung - Genom " + genome + ": Gen " + gene_name + " liegt am Rande des Contigs, die Promoterregion wurde teilweise abgeschnitten. (Fehlercode 002)")
                else:
                    promoter = contig.seq[gene_begin:gene_begin+p_length-1].reverse_complement()

    p_and_g_region = [promoter, gene]
    return(p_and_g_region)

def final_results(fastaList):
    #Gibt Ergebnisse der bisherigen Analysen in einer CSV-Datei aus
    final_path = os.path.join(genomeDir, "all_promoters_results.csv")
    final_results = open(final_path, 'w')

    final_results.write("Genom")
    
    for gene in Gene_Names:
        final_results.write("\t" + gene + "\t" + gene + "-Promotor-Region")
    
    final_results.write("\n")

    for fasta in fastaList:
        final_results.write(fasta)
        genome_name = fasta.split('.')[0]
        gene_allels_path = os.path.join(genomeDir, genome_name + "_genes.txt")
        promotor_allels_path = os.path.join(genomeDir, genome_name + "_genes_promoters.txt")
        
        for gene in Gene_Names:
            allel_list = find_allel_types(fastaList, gene)
            gene_allels = open(gene_allels_path, 'r')
            promotor_allels = open(promotor_allels_path, 'r')
            p_allel = g_allel = " "
            
            for line in promotor_allels:
                if line[0] == ">":
                    stripped_line = line.strip()
                    promotor_name = stripped_line.split("Promoter-Region (Gene): ")[1]
                    if promotor_name == gene:
                        sequence = promotor_allels.readline().split("\n")[0]
                        p_allel = allel_list[0].index(sequence)+1

            for line in gene_allels:
                if line[0] == ">":
                    stripped_line = line.strip()
                    gene_name = stripped_line.split("Gene: ")[1]
                    if gene_name == gene:
                        sequence = gene_allels.readline().split("\n")[0]
                        g_allel = allel_list[1].index(sequence)+1

            final_results.write("\t" + str(g_allel) + "\t" + str(p_allel))
            gene_allels.close()
            promotor_allels.close()
            
        final_results.write("\n")

    final_results.close()

def find_allel_types(fastaList, gene):
    #Identifiziert Allel-Varianten von Genen bzw. deren Promotor-Regionen
    gene_allels = []
    promotor_allels = []
    
    for fasta in fastaList:
        genes_file_path = os.path.join(genomeDir, fasta.split('.')[0] + "_genes.txt")
        genes_file = open(genes_file_path, 'r')
        promoters_file_path = os.path.join(genomeDir, fasta.split('.')[0] + "_genes_promoters.txt")
        promoters_file = open(promoters_file_path, 'r')
        
        for line in promoters_file:
            if line[0] == ">":
                stripped_line = line.strip()
                promotor_name = stripped_line.split("Promoter-Region (Gene): ")[1]
                if promotor_name == gene:
                    sequence = promoters_file.readline().split("\n")[0]
                    if promotor_allels.count(sequence) == 0:
                        promotor_allels.append(sequence)

        for line in genes_file:
            if line[0] == ">":
                stripped_line = line.strip()
                gene_name = stripped_line.split("Gene: ")[1]
                if gene_name == gene:
                    sequence = genes_file.readline().split("\n")[0]
                    if gene_allels.count(sequence) == 0:
                        gene_allels.append(sequence)

        genes_file.close()
        promoters_file.close()

    if promotor_allels and gene_allels:
        genes_allel_file_path = os.path.join(genomeDir, "genes_" + gene + "_allels.txt")
        genes_allel_file = open(genes_allel_file_path, 'w')
        promoters_allel_file_path = os.path.join(genomeDir, "genes_" + gene + "_allels_promoters.txt")
        promoters_allel_file = open(promoters_allel_file_path, 'w')

        for gene_entry in gene_allels:
            genes_allel_file.write("> Gene: " + gene + " | Allel-Type: " + str(gene_allels.index(gene_entry) + 1) + "\n" + str(gene_entry) + "\n")

        for promotor_entry in promotor_allels:
            promoters_allel_file.write("> Gene: " + gene + " | Promoter-Type: " + str(promotor_allels.index(promotor_entry) + 1) + "\n" + str(promotor_entry) + "\n")

        genes_allel_file.close()
        promoters_allel_file.close()

    allel_list = [promotor_allels , gene_allels]
    return(allel_list)

def gene_reader():
    #Ließt die vorgegebene Multifasta-Datei genes.txt ein um die Gen-Liste zurückzugeben
    gene_file_path = os.path.join(genomeDir, "genes.txt")
    gene_file = open(gene_file_path, 'r')

    gene_list = []
    
    for line in gene_file:
        if line[0] == ">":
            stripped_line = line.strip()
            gene_list.append(stripped_line.split("Gene: ")[1])
        else:
            if line[0:3] not in ["ATG", "GTG", "TTG"]:
                print("Warnung - Genom Referenzgenom: Gen " + gene_list[-1] + " beginnt mit unbekanntem Startcodon " + line[0:3] + ". (Fehlercode 004)")

    gene_file.close()
    return(gene_list)

def makeblastDB(genomeName, genomeFile):
    #BLAST-Befehl zum Erstellen einer lokalen BLAST-Datenbank
    blastDBFile = genomeName+"_BlastDB"
    blastDBcommand = 'makeblastdb -dbtype nucl -title S_aureus_ST398 -out ' + blastDBFile + ' -in ' + genomeFile
    os.system('cd '+ genomeDir + ' & ' + blastDBcommand)
    #os.system(blastDBcommand+' & timeout 15')
    return(blastDBFile)

def analyseBlastResult(genome, blastfile):
    #Ließt eine BLAST-Ergebnisdatei aus und parst die Ergebnisse in eine CSV und die Promoter.txt
    blast_results_path = os.path.join(genomeDir, genome + "_parse_results.csv")
    blast_results = open(blast_results_path, 'w')

    blast_xml_path = os.path.join(genomeDir, blastfile)
    blast_xml = open(blast_xml_path)

    genes_path = os.path.join(genomeDir, genome + "_genes.txt")
    genes = open(genes_path, 'w')   

    promoters_path = os.path.join(genomeDir, genome + "_genes_promoters.txt")
    promoters = open(promoters_path, 'w')

    blastRecords = NCBIXML.parse(blast_xml)

    blast_results.write("Genom: " + genome + "\t" + "Hits" + "\t" + "Significant Hits" + "\t" + "Beginn" + "\t" + "Ende" + "\t" + "Länge" + "\t" + "Strang" + "\t" + "Contig" + "\n")

    for blastRecord in blastRecords:
        first_hit = hits = significant_hits = 0
        begin = end = length = strand = contig = ""

        gene = blastRecord.query.split("Gene: ")[1]
        blast_results.write(gene + ":" + "\t")
        
        if blastRecord.alignments == []:
            blast_results.write(str(hits + significant_hits) + "\t" + str(significant_hits) + "\t" + str(begin) + "\t" + str(end) + "\t" + str(length) + "\t" + str(strand) + "\t" + contig + "\n")
        else:
            for align in blastRecord.alignments:
                for hsp in align.hsps:
                    hsp_length = len(hsp.match)
                    if significantHit(hsp_length, blastRecord.query_letters, hsp.identities):
                        significant_hits += 1
                        if first_hit == 0:
                            begin = hsp.sbjct_start
                            end = hsp.sbjct_end
                            if begin == 1 or end == 1:
                               print("Warnung - Genom " + genome + ": Gen " + gene + " liegt am Rande des Contigs, Gen wird ignoriert da die Promoterregion vollständig abgeschnitten wurde. (Fehlercode 003)")
                               break
                            length = abs(begin - end)
                            strand = hsp.frame[1]
                            contig = align.hit_def.split(" ")[0]
                            promoter_and_gene_regions = fetch_p_region(genome, gene, begin, end, str(strand), contig)
                            promoters.write("> Genome: " + genome + " | Promoter-Region (Gene): " + gene + "\n" + str(promoter_and_gene_regions[0]) + "\n")
                            genes.write("> Genome: " + genome + " | Gene: " + gene + "\n" + str(promoter_and_gene_regions[1]) + "\n")
                            first_hit += 1
                    else:
                        hits += 1
            if significant_hits > 1:
                print("Warnung - Genom " + genome + ": Gen " + gene + " besitzt mehrere signifikante BLAST-Treffer. Nur der erste Treffer wird gespeichert. (Fehlercode 005)")
            blast_results.write(str(hits + significant_hits) + "\t" + str(significant_hits) + "\t" + str(begin) + "\t" + str(end) + "\t" + str(length) + "\t" + str(strand) + "\t" + contig + "\n")

    blast_results.close()
    blast_xml.close()
    genes.close()
    promoters.close()

def regular_expression_search():
    #Verwendet reguläre Ausdrücke zur Motifsuche
    all_promoters_path = os.path.join(genomeDir, "all_promoters.txt")
    all_promoters = open(all_promoters_path, 'r')
    regular_expression_results_path = os.path.join(genomeDir, "all_promoters_regular_expression_search.csv")
    regular_expression_results = open(regular_expression_results_path, 'w')

    motif_names = ["Discrimination Region", "-10 Element (Proximal Promotor Box/Pribnow Box)", "Extended -10 Element", "-35 Element (Distal Promotor Box)", "UP Element"]
    motif_list = ["GGG", "TATAAT", "TGTG", "TTGACA", "AAA"]
    mutated_motifs = ["(CGG|GCG|GGC)", "(AATAAT|TTTAAT|TAAAAT|TATTAT|TATATT|TATAAA)", "(AGTG|TCTG|TGAG|TGTC)", "(ATGACA|TAGACA|TTCACA|TTGTCA|TTGAGA|TTGAGT)", "(TAA|ATA|AAT)"]
   
    regular_expression_results.write("Genom" + "\t" + "Promoter-Region (Gene)")

    for motif in motif_names:
        regular_expression_results.write("\t" + motif)

    regular_expression_results.write("\n")

    fasta_name = ""

    for line in all_promoters:
        if line[0] == ">":
            stripped_line = line.strip()
            promoter_name = stripped_line.split("Promoter-Region (Gene): ")[1]
            
            if fasta_name == stripped_line.split(" ")[2]:
                regular_expression_results.write("\t" + promoter_name)
            else:
                fasta_name = stripped_line.split(" ")[2]
                regular_expression_results.write(fasta_name + ".fasta \t" + promoter_name)
        else:
            current_rA = 0
            for rA in motif_list:
                first_match = 0
                matches = ""
                r_line = line[::-1]
                rA_check = re.findall(rA, r_line[0:50])
                rA_mutated_check = re.findall(mutated_motifs[current_rA], r_line[0:50])
                if rA_check:
                    rA_search = re.finditer(rA, r_line[0:50])
                    for match in rA_search:
                        if first_match == 0:
                            hit = r_line[match.span()[0]:match.span()[1]]
                            matches = (hit + str(match.span()))
                            first_match += 1
                        else:
                            hit = r_line[match.span()[0]:match.span()[1]]
                            matches = (matches + "; " + hit + str(match.span()))
                    regular_expression_results.write("\t" + matches)
                #Suche nach Mutationen
                elif rA_mutated_check:
                    rA_mutated_search = re.finditer(mutated_motifs[current_rA], r_line[0:50])
                    for match in rA_mutated_search:
                        if first_match == 0:
                            hit = r_line[match.span()[0]:match.span()[1]]
                            matches = (hit + str(match.span()))
                            first_match += 1
                        else:
                            hit = r_line[match.span()[0]:match.span()[1]]
                            matches = (matches + "; " + hit + str(match.span()))
                    regular_expression_results.write("\t" + matches)
                else:
                    regular_expression_results.write("\t")
                current_rA += 1
            regular_expression_results.write("\n")

    all_promoters.close()
    regular_expression_results.close()

def significantHit(matchLength, geneLength, identities):
    #Überprüfung ob ein Hit wirklich siginifikant ist (nach globalem Threshold und Coverage)
    percentage = float(identities) / matchLength
    #Überprüfung der Identität und Coverage der Sequenzen
    return  percentage > threshold and float(geneLength) * mincoverage <= matchLength

#Main:
def main():  
    fileList = os.listdir(genomeDir)
    genomeList = [name for name in fileList if fnmatch(name, '*.fasta')]

    all_promoters_path = os.path.join(genomeDir, "all_promoters.txt")
    all_promoters = open(all_promoters_path, 'w')
   
    for fasta in genomeList:
        genomeName = fasta.split('.')[0]
        
        blastDB = makeblastDB(genomeName, fasta)
        blastResultsFile = genomeName+"_BlastResults.xml"
        os.system('cd ' + genomeDir + ' & blastn -query genes.txt' + ' -db ' + blastDB + ' -out ' + blastResultsFile + ' -outfmt=5')
        analyseBlastResult(genomeName, blastResultsFile)

        genome_promoters_path = os.path.join(genomeDir, genomeName + "_genes_promoters.txt")
        genome_promoters = open(genome_promoters_path, 'r')
        all_promoters.write(genome_promoters.read())

    genome_promoters.close()
    all_promoters.close()

    regular_expression_search()

    #MUSCLE alignments
    for gene in Gene_Names:
        all_promoters = open(all_promoters_path, 'r')
        first_gene = 0
        for line in all_promoters:
            if line[0] == ">":
                stripped_line = line.strip()
                promotor_name = stripped_line.split("Promoter-Region (Gene): ")[1]
                if promotor_name == gene:
                    if first_gene == 0:
                        gene_promoters = open(os.path.join(genomeDir, "all_promoters_"+gene+".txt"), 'w')
                        first_gene += 1
                    gene_promoters.write(line)
                    gene_promoters.write(all_promoters.readline())
        
        gene_promoters.close()
        os.system('cd ' + genomeDir + ' & ' + muscle_client + ' -in all_promoters_' + gene + '.txt -out all_promoters_' + gene + '_alignment.txt')

    all_promoters.close()

    final_results(genomeList)
  
#Global constants:
#script_dir = os.path.dirname(__file__)
genomeDir = "Genomes"
Gene_Names = gene_reader()
p_region_mode = "large"
mincoverage = 0.9
threshold = 0.7
user_system = platform.system()

if user_system == "Windows":
    muscle_client = "muscle3.8.31_i86win32.exe"
else:
    muscle_client = "muscle3.8.31_i86linux64"

main()

# To-Do:
# Liegen die gesuchten Gene in einem Operon/an erster Position in einem Operon?
# Ergebnisse vergleichen mit Online Promotor Tools (z.B. pEPPER)/bekannten Promoterregionen aus Papern
# Programm auf Linux oder Mac testen
# Gefundene/Nicht-gefundene Gene mit RidomSeqSphere abgleichen
# Warnings exportieren
# Alignments in Geneious laden (Bereiche betrachten - Wo liegen diese? Gibt es (neue) Motife?)
# Verkürtzte Promotor-Regionen markieren?
