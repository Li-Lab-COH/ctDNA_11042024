#!/bin/bash
#SBATCH --job-name=adapter_clean_test    # Job name
#SBATCH --output=/home/janzules/ctDNA_11042024/code/slurmOutput/adapter_clean_test.log
#SBATCH --error=/home/janzules/ctDNA_11042024/code/slurmOutput/adapter_clean_test.err
#SBATCH --ntasks=1              
#SBATCH --cpus-per-task=16      
#SBATCH --mem=120G              
#SBATCH --time=24:00:00         


#region Initializing
#--------------------------------------------------------------------#
#Module loads
module load fastp
module load FastQC


# Define file locations
# fastqz_list="/home/janzules/ctDNA_11042024/code/addresses/all_files.txt"
# fastp_output="/labs/yunroseli_grp/Jon/32_samples_10182024/data/fast_star_all/fastp"

fastqz_list="/home/janzules/ctDNA_11042024/code/addresses/test_address.txt"
fastp_output="/home/janzules/ctDNA_11042024/data/testing/fastp_output/cleaned"
fasp_reports="/home/janzules/ctDNA_11042024/data/testing/fastp_output/reports"

# # In case I have do rerun a test
# rm -rf "$fastp_output"
# mkdir -p "$fastp_output"
# mkdir -p "$star_temp_dir"
# mkdir -p "$output_dir"
# Load the necessary modules



#--------------------------------------------------------------------#
#endregion

#region saving space for parallel cleaning and fastqc

# BLANK

#endregion



echo "+++++++++++++++++++++++++++++++++++++++++++++++"
echo "total files: $(cat "$fastqz_list" | wc -l)"
echo "+++++++++++++++++++++++++++++++++++++++++++++++"

# Read through the list of fastq files and process them in pairs
while IFS= read -r line; do
    # Check if the file is an R1 file (assuming R1 and R2 convention is used)
    if [[ "$line" == *"_R1_"* ]]; then
        r1_file="$line"
        r2_file="${line/_R1_/_R2_}"  # Infer the corresponding R2 file

        # Check if both R1 and R2 files exist
        if [[ -f "$r1_file" && -f "$r2_file" ]]; then
            echo " "
            echo " "
            echo "-------------------------"
            echo "Processing R1: $r1_file"
            echo "Processing R2: $r2_file"
            echo "-------------------------"
            echo " "
            echo " "

            # Extract the A number from the filename for output file naming
            anumber=$(echo "$r1_file" | sed -E 's/.*_(A[0-9]+)_.*R1_.*/\1/')

            # Define output filenames based on the extracted A number
            output_r1="$fastp_output/${anumber}_R1.fastq.gz"
            output_r2="$fastp_output/${anumber}_R2.fastq.gz"
            # Defining reports file
            json_report="$fasp_reports/${anumber}.json"
            html_report="$fasp_reports/${anumber}.html"

            # FastQC report generation before running fastp
            echo " "
            echo "Generating FastQC report before cleaning for: $anumber"
            echo " "

            fastqc "$r1_file" "$r2_file" --outdir "$fasp_reports" --threads 16
            mv "$fasp_reports/$(basename "$r1_file" .fastq.gz)_fastqc.html" "$fasp_reports/${anumber}_R1_before_fastp.html"
            mv "$fasp_reports/$(basename "$r1_file" .fastq.gz)_fastqc.zip" "$fasp_reports/${anumber}_R1_before_fastp.zip"
            mv "$fasp_reports/$(basename "$r2_file" .fastq.gz)_fastqc.html" "$fasp_reports/${anumber}_R2_before_fastp.html"
            mv "$fasp_reports/$(basename "$r2_file" .fastq.gz)_fastqc.zip" "$fasp_reports/${anumber}_R2_before_fastp.zip"

            echo " "
            echo "FastQC report generated, continuing with cleaning for: $anumber"
            echo " "


            # # Run fastp with default settings for paired-end reads
            # test 1 default with quality filtering disabled, but explicitly saying what I want to
            fastp -i $r1_file -o $output_r1 \
            -I $r2_file -O $output_r2 \
            -w 16 \
            --disable_quality_filtering \
            --trim_poly_g \
            --adapter_sequence AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
            --adapter_sequence_r2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
            --cut_tail \
            --cut_mean_quality 20 \
            --cut_window_size 4 \
            --json $json_report \
            --html $html_report \
            

            if [[ -f "$output_r1" && -f "$output_r2" ]] ; then
                echo "Successfully processed: $r1_file and $r2_file -> $output_r1 and $output_r2"
            fi

            
        else
            echo "Warning: R1 and/or R2 files not found for $line"
        fi

        # Check if fastp was successful
        if [[ -f "$output_r1" && -f "$output_r2" ]]; then
            echo "Successfully cleaned: $r1_file and $r2_file -> $output_r1 and $output_r2"

            # FastQC report generation after running fastp
            echo " "
            echo "Generating FastQC report after cleaning for: $anumber"
            echo "cleaned R1"
            echo "$output_r1"
            echo "cleaned R2"
            echo "$output_r2"
            echo " "

            
            fastqc "$output_r1" "$output_r2" --outdir "$fasp_reports" --threads 16
            mv "$fasp_reports/$(basename "$output_r1" .fastq.gz)_fastqc.html" "$fasp_reports/${anumber}_R1_after_fastp.html"
            mv "$fasp_reports/$(basename "$output_r1" .fastq.gz)_fastqc.zip" "$fasp_reports/${anumber}_R1_after_fastp.zip"
            mv "$fasp_reports/$(basename "$output_r2" .fastq.gz)_fastqc.html" "$fasp_reports/${anumber}_R2_after_fastp.html"
            mv "$fasp_reports/$(basename "$output_r2" .fastq.gz)_fastqc.zip" "$fasp_reports/${anumber}_R2_after_fastp.zip"
        else
            echo "Fastp processing failed for: $r1_file and $r2_file"
        fi
    fi
done < "$fastqz_list"
