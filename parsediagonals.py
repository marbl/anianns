import pandas as pd
import argparse

# Function to parse command line arguments
def parse_args():
    parser = argparse.ArgumentParser(description="Filter and merge intervals from input file")
    parser.add_argument('input_file', help='Input file containing data')
    parser.add_argument('filtered_output_file', help='Output file for filtered data')
    parser.add_argument('merged_output_file', help='Output file for merged data')
    return parser.parse_args()

# Main function to process the data
def process_data(input_file, filtered_output_file, merged_output_file):
    # Read the data from the input file into a DataFrame
    df = pd.read_csv(input_file, sep="\t", header=None)

    # Remove the b'' from the first column
    df[0] = df[0].str.replace("b'", "").str.replace("'", "")

    # Count occurrences of each unique value in column 4
    value_counts = df[3].value_counts()

    # Add a new column with the counts
    df['count_col4'] = df[3].map(value_counts)
    df['diff_diag'] = df[3] - df[1] + 1

    # Filter the DataFrame based on the conditions
    filtered_df = df[(df['count_col4'] >= 0.5 * df[5]) & (df['count_col4'] > 3)]

    # Group by unique values in column 3 and get the row with the largest diff_diag for each group
    result_df = filtered_df.loc[filtered_df.groupby(3)['diff_diag'].idxmax()]

    # Save only columns 1, 2, and 4 to the filtered output file
    result_df[[0, 1, 3]].to_csv(filtered_output_file, sep='\t', header=False, index=False)

    # Read the filtered data from the filtered output file into a DataFrame
    try:
        df_filtered = pd.read_csv(filtered_output_file, sep="\t", header=None)
        df_filtered.columns = ['chromosome', 'start', 'end']

        # Sort the DataFrame by 'start' column
        df_filtered = df_filtered.sort_values(by='start')

        # Initialize variables to store merged intervals
        merged_intervals = []

        # Initialize variables to track current interval
        current_chromosome = df_filtered.iloc[0]['chromosome']
        current_start = df_filtered.iloc[0]['start']
        current_end = df_filtered.iloc[0]['end']

        # Iterate through the sorted DataFrame and merge intervals
        for index, row in df_filtered.iloc[1:].iterrows():
            if row['chromosome'] == current_chromosome and row['start'] <= current_end + 1:
                # Merge intervals
                current_end = max(current_end, row['end'])
            else:
                # Store the merged interval
                merged_intervals.append((current_chromosome, current_start, current_end))
                # Update current interval
                current_chromosome = row['chromosome']
                current_start = row['start']
                current_end = row['end']

        # Append the last merged interval
        merged_intervals.append((current_chromosome, current_start, current_end))

        # Convert merged intervals into a DataFrame
        merged_df = pd.DataFrame(merged_intervals, columns=['chromosome', 'start', 'end'])

        # Save the merged intervals to the merged output file
        merged_df.to_csv(merged_output_file, sep='\t', header=False, index=False)
        print(f"Filtered data saved to {filtered_output_file}")
        print(f"Merged data saved to {merged_output_file}")
    except pd.errors.EmptyDataError as err:
        print(f"{input_file} empty.")


if __name__ == '__main__':
    # Parse command line arguments
    args = parse_args()

    # Process data with provided input and output filenames
    process_data(args.input_file, args.filtered_output_file, args.merged_output_file)