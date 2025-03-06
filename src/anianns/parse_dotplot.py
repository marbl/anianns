import numpy as np
import pysam
from pybedtools import BedTool


def row_non_zero_stats(matrix):
    # Set diagonal elements to zero
    np.fill_diagonal(matrix, 0)

    # Total number of non-zero elements per row
    non_zero_counts = np.count_nonzero(matrix, axis=1)

    # Total sum of non-zero elements per row
    non_zero_sums = np.sum(matrix * (matrix != 0), axis=1)

    return non_zero_counts, non_zero_sums


def save_to_file(
    non_zero_counts,
    non_zero_sums,
    window_size,
    final,
    j,
    band,
    filename="non_zero_stats.csv",
):
    with open(filename, "w") as f:
        print(len(non_zero_counts), len(non_zero_sums))
        for i in range(len(non_zero_counts)):
            start = (i * window_size) + (j * band)
            end = min((start + window_size), final)
            f.write(f"{start}-{end}, {non_zero_counts[i]}, {non_zero_sums[i]} \n")
        f.close()


def find_non_zero_length(matrix, window):
    n = matrix.shape[0]
    lengths = []
    coords = []

    # Process rows from the top to the bottom
    for i in range(n):  # Start from the first row (top) to the last row (bottom)
        count = 0
        zero_streak = 0
        col = i  # Diagonal column for the current row
        start_range = col
        end_range = col

        # Check elements on the left side of the diagonal
        for j in range(col - 1, -1, -1):  # Move left from the diagonal
            if matrix[i, j] != 0:
                count += 1
                start_range = j
                zero_streak = 0
            else:
                zero_streak += 1
            if zero_streak == 3:  # Stop after 3 consecutive zeros
                break

        # Reset zero streak for the right side
        zero_streak = 0

        # Check elements on the right side of the diagonal
        for j in range(col + 1, n):  # Move right from the diagonal
            if matrix[i, j] != 0:
                count += 1
                zero_streak = 0
                end_range = j
            else:
                zero_streak += 1
            if zero_streak == 3:  # Stop after 3 consecutive zeros
                break

        lengths.append(count)
        coords.append((start_range * window, end_range * window))

    return lengths, coords


def merge_coordinates(data):
    """
    Merges coordinates and keeps the entry with the highest count if two entries share the same x or y axis.

    Parameters:
    data (Counter): A Counter object with tuples of (x, y) as keys and counts as values.

    Returns:
    list: A list of tuples with merged coordinates and their counts.
    """
    # Filter out entries with count <= 3
    filtered_data = [(key, count) for key, count in data.items() if count > 3]

    # Sort by y coordinate, then by count (descending)
    filtered_data.sort(key=lambda item: (item[0][1], -item[1]))

    result = []
    seen_y = set()
    for (x, y), count in filtered_data:
        # If y already exists, skip this entry (as the highest count for this y is already added)
        if y not in seen_y:
            result.append(((x, y), count))
            seen_y.add(y)

    # Now process for x axis
    result.sort(
        key=lambda item: (item[0][0], -item[1])
    )  # Sort by x, then by count (descending)
    final_result = []
    seen_x = set()
    for (x, y), count in result:
        # If x already exists, skip this entry (as the highest count for this x is already added)
        if x not in seen_x:
            final_result.append(((x, y), count))
            seen_x.add(x)

    # Sort the final result by x-coordinate for better readability
    final_result.sort(key=lambda item: item[0][0])
    return final_result


def mask_diagonals(matrix, coordinates):
    x = 0
    for i in range(matrix.shape[0]):  # Use an index to iterate over rows
        start = int(coordinates[x][0])
        end = int(coordinates[x][1])
        matrix[i, start : end + 1] = 0  # Mask the specified range in the row
        x += 1
    return matrix


def mask_offdiagonals(matrix, coordinates):
    for i in range(matrix.shape[0]):  # Iterate over indices, not rows directly
        start = int(coordinates[i][0])  # Get start of the range
        end = int(coordinates[i][1])  # Get end of the range

        # Set all elements before `start` to 0
        matrix[i, 0:start] = 0

        # Set all elements after `end` to 0
        matrix[i, end + 1 : matrix.shape[1]] = 0

    return matrix


def create_bed_file(intervals, bed_filename):
    with open(bed_filename, "a") as bed_file:
        for interval in intervals:
            chr, start, end = interval
            bed_file.write(f"{chr}\t{start}\t{end}\n")


def mask_fasta_with_bed(fasta_file, seq_name, bed_file, output_fasta):
    # Load the BED file
    bed = BedTool(bed_file)

    # Open the FASTA file
    fasta = pysam.FastaFile(fasta_file)

    # Prepare the output file
    with open(output_fasta, "w") as out_fasta:
        for seq in seq_name:
            print(seq)
            # Get the original sequence
            sequence = list(fasta.fetch(seq))

            # Filter BED entries for the current sequence
            for interval in bed.filter(lambda x: x.chrom == seq_name):
                start, end = int(interval.start), int(interval.end)
                # Replace the region with 'N'
                for i in range(start, end):
                    sequence[i] = "N"

            # Write the masked sequence to the output file
            masked_sequence = "".join(sequence)
            out_fasta.write(f">{seq_name}\n")
            # Split the sequence into lines of 60 characters for FASTA formatting
            for i in range(0, len(masked_sequence), 60):
                out_fasta.write(masked_sequence[i : i + 60] + "\n")

    print(f"Masked FASTA file written to {output_fasta}")
