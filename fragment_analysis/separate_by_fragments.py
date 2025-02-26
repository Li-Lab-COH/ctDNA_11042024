#!/usr/bin/env python3

import os
import sys
import pysam

# I will need to customize this code to work as a testing and then the full thing
# the way it will work is that the input arguments will be supplied to this file by
# anothe script. The __
def split_by_fragment_size(
        # TODO: add a variable that takes the list of Anumbers for a job array and - For testing
        input_bam_path,
        output_dir, # Should be folder location where the fragment_range/bam and bed files will be located
        bins=[(40, 150), (125, 155), (150, 180), (170, 220), (220, 345)]
):
    """
    Splits a BAM file into multiple BAMs based on fragment-size bins.
    """

    # Open the input BAM
    in_bam = pysam.AlignmentFile(input_bam_path, "rb")

    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    # TODO: have to add a bit that makes makes a parent directory given the current bin range
    # TODO: then make a subfolder called bam. This is the folder where all the output is going to go into for the two Anumbers supplied

    # Prepare output BAM files (one per bin)
    out_bams = {}
    for b_start, b_end in bins:
        bin_label = f"{b_start}_{b_end}"
        out_path = os.path.join(output_dir, f"{os.path.basename(input_bam_path).replace('.bam', '')}_{bin_label}.bam")
        out_bams[(b_start, b_end)] = pysam.AlignmentFile(out_path, "wb", template=in_bam)

    # Iterate through reads, write to correct bin
    for read in in_bam.fetch(until_eof=True):
        if read.is_unmapped:
            continue

        # Template length can be negative depending on read orientation
        # so use abs(tlen)
        tlen = abs(read.template_length)
        for b_start, b_end in bins:
            if b_start <= tlen <= b_end:
                out_bams[(b_start, b_end)].write(read)
                break

    # Close all files
    in_bam.close()
    for bam_obj in out_bams.values():
        bam_obj.close()


def main():
    if len(sys.argv) < 3:
        print("Usage: python separate_by_fragments.py <input.bam> <output_dir>")
        sys.exit(1)
    #  TODO: Uncomment to work on the cluster
    # input_bam_path = sys.argv[1]
    # output_dir = sys.argv[2]

    input_bam_path = "/Users/janzules/Roselab/ctDNA_11042024/data/testing_binned_analysis/mapped_bam_subset"
    output_dir = "/Users/janzules/Roselab/ctDNA_11042024/data/testing_binned_analysis/binned_sequences"
    # Define your custom bins
    bins = [
        (40, 150),
        (125, 155),
        (150, 180),
        (170, 220),
        (220, 345),
    ]

    split_by_fragment_size(input_bam_path, output_dir, bins)


if __name__ == "__main__":
    main()
