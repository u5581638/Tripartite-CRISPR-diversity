

module load python3/3.10.4
find split_sequences_0_5kb/ -name "*_pfam_descriptions.txthits.csv" -type "f" | xargs -n 1 -I {} -P 1 python3 sequence_name_adder.py {}
find split_sequences_0_5kb/ -name "*_crisprcas_hits.txt_hmmscan.csv" -type "f" | xargs -n 1 -I {} -P 1 python3 sequence_name_adder.py {}
find split_sequences_5_10kb/ -name "*_crisprcas_hits.txt_hmmscan.csv" -type "f" | xargs -n 1 -I {} -P 1 python3 sequence_name_adder_5_10kb.py {}
find split_sequences_5_10kb/ -name "*_pfam_descriptions.txthits.csv" -type "f" | xargs -n 1 -I {} -P 1 python3 sequence_name_adder_5_10kb.py {}

find split_sequences_5_10kb/ -name "*pfam_descriptions.txthits.csvheader_added.csv" -type "f" | xargs -n 1 -I {} -P 1  cat {} >> pfam_5_10kb_header_hits2.csv
find split_sequences_0_5kb/ -name "*pfam_descriptions.txthits.csvheader_added.csv" -type "f" | xargs -n 1 -I {} -P 1  cat {} >> pfam_0_5kb_header_hits2.csv
find split_sequences_0_5kb/ -name "*_crisprcas_hits.txt_hmmscan.csvheader_added.csv" -type "f" | xargs -n 1 -I {} -P 1  cat {} >> padlocplus_0_5kb_header_hits2.csv
find split_sequences_5_10kb/ -name "*_crisprcas_hits.txt_hmmscan.csvheader_added.csv" -type "f" | xargs -n 1 -I {} -P 1  cat {} >> padlocplus_5_10kb_header_hits2.csv