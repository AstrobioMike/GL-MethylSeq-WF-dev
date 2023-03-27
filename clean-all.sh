rm -rf MultiQC_Outputs/ work/ .nextflow* samples.txt Trimmed_Sequence_Data/ Bismark* Reference_Genome_Files/ Metadata/ MethylKit_Outputs/ Resource_Usage/
rm -rf .nf-test/ tests/MultiQC_Outputs/ tests/work/ tests/.nextflow* tests/samples.txt tests/Filtered_Sequence_Data/ tests/Bismark* tests/Reference_Genome_Files/ tests/Metadata/ tests/MethylKit_Outputs/ tests/Resource_Usage/

if [ ! -z ${1} ]; then

    if [ $1 == 'all' ]; then
    
        rm -rf test-data/

    fi

fi
