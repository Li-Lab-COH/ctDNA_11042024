

name_mapping_file="/home/janzules/ctDNA_11042024/code/addresses/human_samples_name_mapping.csv"


map_qual_dir="$output_dir/map_qual"


declare -A id_mapping
while IFS=',' read -r tgen_id patient_id; do
   # Skip the header
   if [[ "$tgen_id" != "TGen_Sample_ID" ]]; then
       id_mapping["$tgen_id"]="$patient_id"
   fi
done < "$name_mapping_file"

