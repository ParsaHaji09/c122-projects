##########################################################################################
# GOAL: find which genome each read comes from, knowing that the reads come from a handful
# genomes
##########################################################################################
""" imports """
import os
import random
import re
##########################################################################################
""" functions """

def map_genomes(reads, genomes):
    reads_to_genome = {}
    genome_to_reads = {}

    for genome_number, genome in genomes.items():
        kmer_genome_positions = {}
        for i in range(len(genome) - 16 + 1):
            kmer = genome[i:i + 16]
            if kmer not in kmer_genome_positions:
                kmer_genome_positions[kmer] = [i]
            else:
                kmer_genome_positions[kmer].append(i)

        for i, read in enumerate(reads):
            kmers_read = [read[j:j+16] for j in range(0, len(read) - 16 + 1)]
            for kmer in kmers_read:
                if kmer in kmer_genome_positions:
                    if i not in reads_to_genome:
                        reads_to_genome[i] = []
                    reads_to_genome[i].append(genome_number)

                    if genome_number not in genome_to_reads:
                        genome_to_reads[genome_number] = []
                    genome_to_reads[genome_number].append(i)
                    break

    return reads_to_genome, genome_to_reads

def low_genomes(genome_to_reads):
    small_genomes = []
    for genome, reads in genome_to_reads.items():
        if len(reads) < 9998:
            small_genomes.append(genome)
    return small_genomes
##########################################################################################
""" main """

reads = []
file_path = './project1c_reads.fasta'
with open(file_path, 'r') as f:
    for line in f:
        line = line.strip()
        if line and line[0] != '>':
            reads.append(line[:48])


genomes = {}
genome_directory = './project1c'
for file_name in os.listdir(genome_directory):
    file_path = os.path.join(genome_directory, file_name)
    with open(file_path, 'r') as f:
        genome = ''
        genome_number = None
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                genome_number = re.search(r"Genome_Number_(\d+)", line).group(1)
            else:
                genome += line
        genomes[genome_number] = genome

reads_to_genome, genome_to_reads = map_genomes(reads, genomes)
small_genomes = low_genomes(genome_to_reads)
for read_id, genomes in reads_to_genome.items():
    reads_to_genome[read_id] = [genome for genome in genomes if genome not in small_genomes]

sorted_reads = sorted(reads_to_genome.keys())
for index in sorted_reads:
    genomes = reads_to_genome[index]
    if not genomes:
        genome_number = 1
    elif len(genomes) > 1:
        genome_number = random.choice(genomes)
    else:
        genome_number = genomes[0]
    with open('result.txt', 'a') as file: 
        file.write(">read_" + str(index) + " Genome_Number_" + str(genome_number) + f"\n")