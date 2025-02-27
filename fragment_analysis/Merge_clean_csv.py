#!/usr/bin/env python3

import os
import sys
import glob
import pandas as pd
from datetime import datetime

def main(directory):
    """
    Combine all CSV files matching fragment_size_summary*.csv
    in `directory` into one final CSV, remove the partial CSV files,
    and save the final CSV with a timestamp in the filename.
    """

    # 1. Find the relevant CSV files
    pattern = os.path.join(directory, "fragment_size_summary*.csv")
    csv_files = glob.glob(pattern)

    if not csv_files:
        print(f"[WARNING] No CSV files matching {pattern}")
        sys.exit(0)

    print("[INFO] Found the following CSV files to combine:")
    for f in csv_files:
        print(f" - {f}")

    # 2. Read all CSVs into a list of DataFrames
    dfs = []
    for csv_file in csv_files:
        df = pd.read_csv(csv_file)
        dfs.append(df)

    # 3. Concatenate all DataFrames into a single DataFrame
    final_df = pd.concat(dfs, ignore_index=True)

    # 4. Remove the original CSV files
    for csv_file in csv_files:
        os.remove(csv_file)
        print(f"[INFO] Deleted {csv_file}")

    # 5. Save the merged DataFrame with a timestamp
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    output_filename = f"final_fragment_size_summary_{timestamp}.csv"
    output_path = os.path.join(directory, output_filename)
    final_df.to_csv(output_path, index=False)

    print(f"[INFO] Final combined CSV saved to {output_path}")


if __name__ == "__main__":
    """
    Usage: python combine_csvs.py <directory>
    """
    if len(sys.argv) < 2:
        print("Usage: python combine_csvs.py <directory>")
        sys.exit(1)

    outdir = sys.argv[1]
    main(outdir)
