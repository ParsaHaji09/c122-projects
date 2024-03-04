##########################################################################################
# GOAL: Determine if there were insertions, deletions or subsitutions in reference genome 
# based on error-reads
##########################################################################################

reads = []
file_path = './project1b-u_with_error_paired_reads.fasta'
with open(file_path, 'r') as f:
    for line in f:
        line = line.strip()
        if line and line[0] != '>':
            reads.append(line[:48])


genome = ''
genome_directory = './project1c'

file_path = './project1b-u_reference_genome.fasta'
with open(file_path, 'r') as f:
    genome_number = None
    for line in f:
        line = line.strip()
        if line.startswith('>'):
            continue
        else:
            genome += line

##########################################################################################

#build hash_table
def create_dictionary(larger_string):
    result_dict = {}
    for i in range(len(larger_string) - 15):
        substring = larger_string[i:i+16]
        if substring in result_dict:
            result_dict[substring].append(i)
        else:
            result_dict[substring] = [i]
    return result_dict

def hamming_distance(p, q):
    """Calculate the Hamming distance between two strings."""
    if len(p) != len(q):
        raise ValueError("Input strings must be of the same length.")
    count = 0
    indexes = []
    changes = []
    for i in range(len(q)):
        if p[i] != q[i]:
            count += 1
            indexes.append(i)
            changes.append((p[i], q[i]))
    return count, indexes, changes

def get_deletion(error, index):
    reference_snippet = genome[index:index+len(error)]
    for i in range(len(error)):
        temp_error = error[:i] + "-" + error[i:]
        temp_error = temp_error[1:]
        count,_,_ = hamming_distance(reference_snippet, temp_error)
        if count <= 2:
            return ("D", index + i, (error[i]))
    return ("N", index, (error[0],0))

def get_insertion(error, index):
    reference_snippet = genome[index:index+len(error)]
    for i in range(len(error)):
        temp_error = error[:i] + error[i+1:]
        temp_reference = reference_snippet [1:]
        count,_,_ = hamming_distance(temp_reference, temp_error)
        if count <= 2:
            return ("I", index + 1 + i, (error[i]))
    return ("N", index, (error[0],0))


def configure_sequence(a, preferred_index):
    # given that it won't be a perfect match 2/3rds
    # check for substitution or indel
    # return tuple of format (s, index)
    count, indexes, changes = hamming_distance(genome[preferred_index: preferred_index + len(a)], a)
    
    result = []
    # probably a subsitution
    if count <= 2:
        for i,j in zip(indexes,changes):
            result.append (("S", i + preferred_index, j))
    # indel case
    else:
        
        insertion = get_insertion(a, preferred_index)
        deletion = get_deletion(a, preferred_index)

        if insertion[0] != "N":
            result.append(insertion)
        else:
            result.append(deletion)

    return result

def check_continuous(a1, a2):
    for i in a1:
        for j in a2:
            if j-i == 16:
                return i,j
    return -1,-1

##########################################################################################

hash_16 = create_dictionary(genome)
results_list = {}

for error_read in reads:
    
    if len(error_read) > 47 and len(error_read) < 51:
        array_one = []
        array_two = []
        array_three = []

        missing_portion = ""
        missing_index = -1

        if error_read[0:16] in hash_16:
            array_one = hash_16[error_read[0:16]]
        if error_read[16:32] in hash_16:
            array_two = hash_16[error_read[16:32]]
        if error_read[32:48] in hash_16:
            array_three = hash_16[error_read[32:48]]

        section3_i1, section3_i2 = check_continuous (array_one, array_two)
        section2_i1, section2_i3 = check_continuous (array_one, array_three)
        section1_i2, section1_i3 = check_continuous (array_two, array_three)

        if array_one and array_two and section3_i1 != -1:
            i3 = section3_i2 + 16
            if i3 in array_three:
                continue
            else:
                missing_portion = error_read[32:48]
                missing_index = i3
                
        elif array_two and array_three and section1_i2 != -1:
            i1 = section1_i2 - 16
            if i1 in array_one:
                continue
            else:
                missing_portion = error_read[0:16]
                missing_index = i1

        elif array_one and array_three and section2_i1 != -1:
            i2 = section2_i1 + 16
            if i2 in array_two:
                continue
            else:
                missing_portion = error_read[16:32]
                missing_index = i2
        
        if missing_portion != "":
            res_list = configure_sequence(missing_portion, i3)
            for res_item in res_list:
                if res_item[0] != "N":
                    if res_item in results_list:
                        results_list[res_item] += 1
                    else:
                        results_list[res_item] = 1


sorted_dict = dict(sorted(results_list.items(), key=lambda item: (0 if item[0][0] == 'S' else (1 if item[0][0] == 'I' else 2), item[0])))

with open('result.txt', 'w') as file:       
    for result,freq in sorted_dict.items():
        if freq > 1:
            change_type = result[0]
            index = result[1]
            changes_str = ""
            if change_type == "S":
                changes_str = str(result[2][0]) + " " + str(result[2][1])
            else:
                changes_str = str(result[2][0])
            string = f">{change_type}{index} {changes_str}\n"

            file.write(string)