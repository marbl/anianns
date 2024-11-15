import numpy as np
import argparse

# Function to read TSV file and convert it to a NumPy array
def read_tsv_to_numpy(file_path):
    return np.genfromtxt(file_path, delimiter='\t', dtype=object, encoding=None)

# Function to process the data
def process_data(data):
    output = []

    # Initialize variables to keep track of the current group
    current_query_name = None
    current_query_start = None
    current_reference_start = None
    current_reference_end = None
    identities = []
    rows_merged = 0

    def finalize_group():
        if current_query_name is not None:
            avg_identity = sum(identities) / len(identities)
            output.append([current_query_name, current_query_start, current_reference_start, current_reference_end, round(avg_identity, 2), rows_merged])

    for row in data:
        query_name, query_start, query_end, _, reference_start, reference_end, identity = row
        
        query_start = int(query_start)
        query_end = int(query_end)
        reference_start = int(reference_start)
        reference_end = int(reference_end)
        identity = float(identity)

        # If starting a new group
        if current_reference_start is None or query_start != current_query_start or reference_start > current_reference_end + 30001:
            # Finalize the previous group
            finalize_group()
            
            # Start a new group
            current_query_name = query_name
            current_query_start = query_start
            current_reference_start = reference_start
            current_reference_end = reference_end
            identities = [identity]
            rows_merged = 1
        else:
            # Continue the current group
            current_reference_end = reference_end
            identities.append(identity)
            rows_merged += 1

    # Finalize the last group
    finalize_group()

    # Filter out groups with 3 or fewer rows merged
    output = [row for row in output if row[5] > 3]

    # Convert the output to a NumPy array for consistency
    output_array = np.array(output, dtype=object)
    # Sort the output by the reference start column (3rd column) and then by the query start column (2nd column)
    output_array = output_array[np.lexsort((output_array[:, 1], output_array[:, 3]))]
    
    return output_array

# Function to print the output in the desired format
def print_output(output_array):
    for row in output_array:
        print("\t".join(map(str, row)))

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process a TSV file and output the filtered results.')
    parser.add_argument('input_file', type=str, help='Path to the input TSV file')
    parser.add_argument('output_file', type=str, help='Path to the output TSV file')
    args = parser.parse_args()

    data = read_tsv_to_numpy(args.input_file)
    output_array = process_data(data)
    #print_output(output_array)

    # Save filtered output to a file
    np.savetxt(args.output_file, output_array, fmt='%s', delimiter='\t')