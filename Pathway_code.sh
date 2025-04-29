############################################## Pathway Analysis #################################################




############################################## Pathway Merging ##################################################
humann_join_tables -i ./ -o merged_pathabundance.tsv --file_name pathabundance
humann_join_tables -i ./ -o merged_genefamilies.tsv --file_name genefamilies
humann_join_tables -i ./ -o merged_pathcoverage.tsv --file_name pathcoverage
humann_renorm_table -i merged_pathabundance.tsv -o merged_pathabundance_cpm.tsv --units cpm
awk -F'\t' '{sum=0; for(i=2; i<=NF; i++) sum+=$i; if(sum>0) print $0}' merged_pathabundance.tsv > merged_pathabundance_filtered.tsv
awk -F'\t' 'NR==1 || ($2+$3+$4+$5+$6+$7+$8 != 0)' merged_pathabundance.tsv > filtered_pathabundance.tsv 
humann_renorm_table -i filtered_pathabundance.tsv -o merged_pathabundance_cpm1.tsv --units cpm
