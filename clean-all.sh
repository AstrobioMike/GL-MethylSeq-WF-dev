rm -rf FastQC_Outputs/ work/ .nextflow* samples.txt Filtered_Sequence_Data/ Bismark* Reference_Genome_Files/ Metadata/ Differential_Methylation_Analysis_Data/ Resource_Usage/ Methylation_Call_Data/
rm -rf .nf-test/ tests/FastQC_Outputs/ tests/work/ tests/.nextflow* tests/samples.txt tests/Filtered_Sequence_Data/ tests/Bismark* tests/Reference_Genome_Files/ tests/Metadata/ tests/Differential_Methylation_Analysis_Data/ tests/Resource_Usage/ tests/Methylation_Call_Data/

if [ ! -z ${1} ]; then

    if [ $1 == 'all' ]; then
    
        rm -rf test-data/

    fi

fi
