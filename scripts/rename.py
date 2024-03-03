def process_fasta(fasta_file, names_file, output_file):
    # Read the names from the txt file
    with open(names_file, 'r') as names_fh:
        new_names = [line.strip() for line in names_fh]

    # Open the input fasta file and the output file
    with open(fasta_file, 'r') as fasta_fh, open(output_file, 'w') as output_fh:
        current_name_index = 0

        for line in fasta_fh:
            if line.startswith('>'):
                # Replace the existing header with the new name
                print(current_name_index)

                # Create new header and add reference sequence chromosome name
                new_header = '>' + new_names[current_name_index] + ' CP004010.2'
                output_fh.write(new_header + '\n')
                current_name_index += 1
            else:
                # Remove asterisks from sequence lines
                cleaned_sequence = line.replace('*', 'N')
                output_fh.write(cleaned_sequence)

if __name__ == "__main__":
    
    # Replace 'merged.fasta', 'strains.csv', and 'output.fasta' with your file paths
    process_fasta('merged.fasta', 'strains.csv', 'output.fasta')