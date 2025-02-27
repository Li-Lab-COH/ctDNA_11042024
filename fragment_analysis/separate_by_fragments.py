#!/usr/bin/env python3

import os
import sys
import glob
import pysam
import pandas as pd   # <-- We'll use pandas to build the final matrix.

def split_by_fragment_size(
        anumbers,
        input_dir,
        output_dir,
        bins=[(40, 150), (125, 155), (150, 180), (170, 220), (220, 345)]
):
    """
    For each A-number in anumbers:
      1. Locate the corresponding BAM file(s) in input_dir.
      2. Create subfolders in output_dir for each bin range: e.g. 170_220/bam/.
      3. Split reads by fragment size and write to anumber_binrange.bam in the correct bin folder.
      4. Keep track of counts (# of reads in original and in each bin) for summary.
    """

    # Create a list to store the row of counts for each anumber.
    # Each row will look like:
    # {
    #    "anumber": Axxx,
    #    "original_sequence_count": ...,
    #    "bin_40_150": ...,
    #    "bin_125_155": ...,
    #    ...
    # }
    results = []

    # Create subfolders for each bin range ahead of time
    for (b_start, b_end) in bins:
        bin_label = f"{b_start}_{b_end}"
        bin_dir = os.path.join(output_dir, bin_label)
        bam_dir = os.path.join(bin_dir, "bam")
        os.makedirs(bam_dir, exist_ok=True)

    # Process each A-number
    for anumber in anumbers:
        # Find the matching BAM file
        pattern = os.path.join(input_dir, f"{anumber}*.bam")
        input_bam_path = glob.glob(pattern)

        if len(input_bam_path) != 1:
            print(f"[WARNING] Expected 1 BAM file for {anumber}, but found {len(input_bam_path)} in {input_dir}",
                  file=sys.stderr)
            continue  # Skip this A-number if there are 0 or more than 1 matching files

        input_bam_path = input_bam_path[0]
        print(f"[INFO] Processing {anumber} with input file: {input_bam_path}")

        # Open the input BAM once
        in_bam = pysam.AlignmentFile(input_bam_path, "rb")

        # -------------------------
        # Initialize counters
        # -------------------------
        total_reads = 0
        bin_counts = {(b_start, b_end): 0 for (b_start, b_end) in bins}

        # Create output BAM handles for each bin range
        out_bams = {}
        for (b_start, b_end) in bins:
            bin_label = f"{b_start}_{b_end}"
            bam_dir = os.path.join(output_dir, bin_label, "bam")
            out_name = f"{anumber}_{bin_label}.bam"
            out_path = os.path.join(bam_dir, out_name)
            out_bams[(b_start, b_end)] = pysam.AlignmentFile(out_path, "wb", template=in_bam)

        # Read the input BAM, distribute reads into the appropriate bin
        for read in in_bam.fetch(until_eof=True):
            if read.is_unmapped:
                continue

            total_reads += 1  # increment total read count

            tlen = abs(read.template_length)
            for (b_start, b_end) in bins:
                if b_start <= tlen <= b_end:
                    out_bams[(b_start, b_end)].write(read)
                    bin_counts[(b_start, b_end)] += 1
                    break  # Stop once we've written it to the matching bin

        # Close input and output BAM files
        in_bam.close()
        for bam_obj in out_bams.values():
            bam_obj.close()

        # -------------------------
        # Store row of counts
        # -------------------------
        row_data = {
            "anumber": anumber,
            "original_sequence_count": total_reads
        }

        # For each bin, add a column named bin_<start>_<end>
        for (b_start, b_end) in bins:
            bin_label = f"bin_{b_start}_{b_end}"
            row_data[bin_label] = bin_counts[(b_start, b_end)]

        results.append(row_data)

        print(f"[INFO] Finished splitting {anumber}")

    # -------------------------
    # Create/Print the DataFrame
    # -------------------------
    df = pd.DataFrame(results)

    # Get current timestamp in YYYYMMDD_HHMMSS format
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

    # Create output filename with timestamp
    output_filename = f"fragment_size_summary_{timestamp}.csv"

    # Save DataFrame to output directory
    output_path = os.path.join(output_dir, output_filename)
    df.to_csv(output_path, index=False)

    print(f"\n[INFO] Summary saved to: {output_path}")


def main():
    if len(sys.argv) < 4:
        print("Usage: python separate_by_fragments.py <comma_separated_anumbers> <input_dir> <output_dir>")
        sys.exit(1)

    # 1) Parse command-line args
    anumbers_str = sys.argv[1]  # e.g. "A14891,A14892"
    anumbers = anumbers_str.split(",")
    input_dir = sys.argv[2]
    output_dir = sys.argv[3]

    # 2) Define your bin ranges
    bins = [
        (40, 150),
        (125, 155),
        (150, 180),
        (170, 220),
        (220, 345),
    ]

    # 3) Run
    split_by_fragment_size(anumbers, input_dir, output_dir, bins)


if __name__ == "__main__":
    main()
