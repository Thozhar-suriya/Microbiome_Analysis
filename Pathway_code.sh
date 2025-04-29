############################################## Pathway tools Installation #########################################
conda activate biobakery3
wget http://huttenhower.sph.harvard.edu/humann_data/chocophlan/full_chocophlan.v201901_v31.tar.gz
humann3_databases --download chocophlan full /home/vamouda/Suriya/Meta_JIP/humann --update-config yes
pip install metaphlan2
bowtie2 --version
tar -xvzf diamond-linux64.tar.gz
ls
which diamond
conda install -c bioconda diamond=2.0.15
which diamond
mv diamond ~/Paul/miniconda3/envs/biobakery3/bin/
chmod +x ~/Paul/miniconda3/envs/biobakery3/bin/diamond
diamond-version

humann3_databases --download uniref : uniref90_diamond full /home/vamouda/Suriya/Meta_JIP/humann --update-config 
wget http://huttenhower.sph.harvard.edu/humann_data/uniprot/uniref_annotated/uniref50_annotated_v201901b_full.tar.gz
humann3_databases --download uniref uniref50_diamond /home/vamouda/Suriya/Meta_JIP/humann --update-config yes
############################################## Pathway Analysis #################################################
# Create output and merged directories
mkdir -p ../out/merged_humann3

# Process each sample
for sample in ../Complete1/*.fastq; do
    sample_name=$(basename "$sample" .fastq)
    echo "Processing $sample_name..."

    # Run HUMAnN3 with outputs in sample-specific folders
    humann3 -i "$sample" \
            -o ../out/"$sample_name"_humann3 \
            --threads 64 \
            --metaphlan-options "--nproc 64" \
            --diamond-options "--threads 64"

    # Create symlinks for merging in ../out/merged_humann3
    for file_type in genefamilies pathabundance pathcoverage; do
        # Source file (original output)
        src_file="../out/${sample_name}_humann3/${sample_name}_humann3_${file_type}.tsv"
        
        # Symlink name (to avoid conflicts)
        symlink_name="../out/merged_humann3/${sample_name}_${file_type}.tsv"
        
        # Create symlink if source file exists
        if [ -f "$src_file" ]; then
            ln -sf "$src_file" "$symlink_name"
        else
            echo "WARNING: File $src_file not found!"
        fi
    done
done

echo "Done! Merged files are in ../out/merged_humann3"
############################################## Pathway Merging ##################################################
humann_join_tables -i ./ -o merged_pathabundance.tsv --file_name pathabundance
humann_join_tables -i ./ -o merged_genefamilies.tsv --file_name genefamilies
humann_join_tables -i ./ -o merged_pathcoverage.tsv --file_name pathcoverage
humann_renorm_table -i merged_pathabundance.tsv -o merged_pathabundance_cpm.tsv --units cpm
awk -F'\t' '{sum=0; for(i=2; i<=NF; i++) sum+=$i; if(sum>0) print $0}' merged_pathabundance.tsv > merged_pathabundance_filtered.tsv
awk -F'\t' 'NR==1 || ($2+$3+$4+$5+$6+$7+$8 != 0)' merged_pathabundance.tsv > filtered_pathabundance.tsv 
humann_renorm_table -i filtered_pathabundance.tsv -o merged_pathabundance_cpm1.tsv --units cpm
