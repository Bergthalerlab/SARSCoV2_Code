#!/bin/bash
#SBATCH --output /scratch/lab_bergthaler/2020_SARS-CoV-2_Evolution/phylogeny/CeMM_build_global_new_v19/slurmoutput/Phylo.%j.%N.out 
#SBATCH --error /scratch/lab_bergthaler/2020_SARS-CoV-2_Evolution/phylogeny/CeMM_build_global_new_v19/slurmoutput/Phylo.%j.%N.err
#SBATCH --job-name=19CoV2
#SBATCH --partition=shortq
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --time=12:00:00
#SBATCH --mem=40000


echo "Enviromental variables"
echo "======================"

echo $SLURM_SUBMIT_DIR
echo $SLURM_JOB_NAME
echo $SLURM_JOB_PARTITION
echo $SLURM_NTASKS
echo $SLURM_NPROCS
echo $SLURM_JOB_ID
echo $SLURM_JOB_NUM_NODES
echo $SLURM_NODELIST
echo $SLURM_CPUS_ON_NODE

echo "======================"

## NEXTSTRAIN COMMAND LIST

# filtering was done manually, subsampling then with subsampling.py down to 8000 sequences
augur filter \
  --sequences /scratch/lab_bergthaler/2020_SARS-CoV-2_Evolution/phylogeny/CeMM_build_global_new_v19/data/sequences.fasta \
  --metadata /scratch/lab_bergthaler/2020_SARS-CoV-2_Evolution/phylogeny/CeMM_build_global_new_v19/data/metadata.tsv \
  --exclude /scratch/lab_bergthaler/2020_SARS-CoV-2_Evolution/phylogeny/CeMM_build_global_new_v19/config/exclude.txt \
  --output /scratch/lab_bergthaler/2020_SARS-CoV-2_Evolution/phylogeny/CeMM_build_global_new_v19/results/filtered.fasta \
  --group-by division year month \
  --min-length 25000

augur align \
  --sequences /scratch/lab_bergthaler/2020_SARS-CoV-2_Evolution/phylogeny/CeMM_build_global_new_v19/data/filteredsubsampled.fasta \
  --reference-sequence /scratch/lab_bergthaler/2020_SARS-CoV-2_Evolution/phylogeny/CeMM_build_global_new_v19/config/reference.fasta \
  --output /scratch/lab_bergthaler/2020_SARS-CoV-2_Evolution/phylogeny/CeMM_build_global_new_v19/results/aligned.fasta \
  --fill-gaps \
  --debug \
  --nthreads 16

  # masking first 130 and last 50 nucleotides of alignment
augur mask \
  --sequences /scratch/lab_bergthaler/2020_SARS-CoV-2_Evolution/phylogeny/CeMM_build_global_new_v19/results/aligned.fasta \
  --mask /scratch/lab_bergthaler/2020_SARS-CoV-2_Evolution/phylogeny/CeMM_build_global_new_v19/config/mask.bed \
  --output /scratch/lab_bergthaler/2020_SARS-CoV-2_Evolution/phylogeny/CeMM_build_global_new_v19/results/aligned_masked.fasta \
  --no-cleanup

# build the tree => MIGHT NEED TO ADJUST NUMBER OF AVAILABLE CPUS UNDER NTHREADS! MUCH FASTER WITH MORE CPUs

augur tree \
  --alignment /scratch/lab_bergthaler/2020_SARS-CoV-2_Evolution/phylogeny/CeMM_build_global_new_v19/results/aligned_masked.fasta \
  --output /scratch/lab_bergthaler/2020_SARS-CoV-2_Evolution/phylogeny/CeMM_build_global_new_v19/results/tree_raw.nwk \
  --vcf-reference /scratch/lab_bergthaler/2020_SARS-CoV-2_Evolution/phylogeny/CeMM_build_global_new_v19/config/reference.fasta \
  --nthreads 16

# Tree refinement => if this gives an error then it's probably about the name of Wuhan-Hu-1/2019 
augur refine \
  --tree /scratch/lab_bergthaler/2020_SARS-CoV-2_Evolution/phylogeny/CeMM_build_global_new_v19/results/tree_raw.nwk \
  --alignment /scratch/lab_bergthaler/2020_SARS-CoV-2_Evolution/phylogeny/CeMM_build_global_new_v19/results/aligned_masked.fasta \
  --metadata /scratch/lab_bergthaler/2020_SARS-CoV-2_Evolution/phylogeny/CeMM_build_global_new_v19/data/metadata.tsv \
  --output-tree /scratch/lab_bergthaler/2020_SARS-CoV-2_Evolution/phylogeny/CeMM_build_global_new_v19/results/tree.nwk \
  --output-node-data /scratch/lab_bergthaler/2020_SARS-CoV-2_Evolution/phylogeny/CeMM_build_global_new_v19/results/branch_lengths.json \
  --timetree \
  --coalescent skyline \
  --date-inference marginal \
  --date-confidence \
  --no-covariance \
  --precision 3 \
  --date-format "%d/%m/%Y" \
  --divergence-unit mutations \
  --clock-filter-iqd 4 \
  --root "Wuhan-Hu-1/2019"

#Reconstruct Ancestral Traits
augur traits \
  --tree /scratch/lab_bergthaler/2020_SARS-CoV-2_Evolution/phylogeny/CeMM_build_global_new_v19/results/tree.nwk \
  --metadata /scratch/lab_bergthaler/2020_SARS-CoV-2_Evolution/phylogeny/CeMM_build_global_new_v19/data/metadata.tsv \
  --output-node-data /scratch/lab_bergthaler/2020_SARS-CoV-2_Evolution/phylogeny/CeMM_build_global_new_v19/results/traits.json \
  --columns country location region division \
  --confidence
  

#Infer Ancestral Sequences
augur ancestral \
  --tree /scratch/lab_bergthaler/2020_SARS-CoV-2_Evolution/phylogeny/CeMM_build_global_new_v19/results/tree.nwk \
  --alignment /scratch/lab_bergthaler/2020_SARS-CoV-2_Evolution/phylogeny/CeMM_build_global_new_v19/results/aligned_masked.fasta \
  --output-node-data /scratch/lab_bergthaler/2020_SARS-CoV-2_Evolution/phylogeny/CeMM_build_global_new_v19/results/nt_muts.json \
  --infer-ambiguous \
  --inference joint 


#Identify Amino-Acid Mutations
augur translate \
  --tree /scratch/lab_bergthaler/2020_SARS-CoV-2_Evolution/phylogeny/CeMM_build_global_new_v19/results/tree.nwk \
  --ancestral-sequences /scratch/lab_bergthaler/2020_SARS-CoV-2_Evolution/phylogeny/CeMM_build_global_new_v19/results/nt_muts.json \
  --vcf-reference /scratch/lab_bergthaler/2020_SARS-CoV-2_Evolution/phylogeny/CeMM_build_global_new_v19/config/reference.fasta \
  --vcf-reference-output /scratch/lab_bergthaler/2020_SARS-CoV-2_Evolution/phylogeny/CeMM_build_global_new_v19/results/reference_translation.fasta \
  --reference-sequence /scratch/lab_bergthaler/2020_SARS-CoV-2_Evolution/phylogeny/CeMM_build_global_new_v19/config/MN908947_annotations.gb \
  --output /scratch/lab_bergthaler/2020_SARS-CoV-2_Evolution/phylogeny/CeMM_build_global_new_v19/results/aa_muts.json \
  --alignment-output /scratch/lab_bergthaler/2020_SARS-CoV-2_Evolution/phylogeny/CeMM_build_global_new_v19/results/aligned_aa_%GENE.fasta


#Clade definitions
augur clades \
  --tree /scratch/lab_bergthaler/2020_SARS-CoV-2_Evolution/phylogeny/CeMM_build_global_new_v19/results/tree.nwk \
  --mutations /scratch/lab_bergthaler/2020_SARS-CoV-2_Evolution/phylogeny/CeMM_build_global_new_v19/results/nt_muts.json /scratch/lab_bergthaler/2020_SARS-CoV-2_Evolution/phylogeny/CeMM_build_global_new_v19/results/aa_muts.json \
  --clades /scratch/lab_bergthaler/2020_SARS-CoV-2_Evolution/phylogeny/CeMM_build_global_new_v19/config/clades.tsv \
  --output-node-data /scratch/lab_bergthaler/2020_SARS-CoV-2_Evolution/phylogeny/CeMM_build_global_new_v19/results/clades.json


#Export the Results

augur export v2\
  --tree /scratch/lab_bergthaler/2020_SARS-CoV-2_Evolution/phylogeny/CeMM_build_global_new_v19/results/tree.nwk \
  --metadata /scratch/lab_bergthaler/2020_SARS-CoV-2_Evolution/phylogeny/CeMM_build_global_new_v19/data/metadata.tsv \
  --node-data /scratch/lab_bergthaler/2020_SARS-CoV-2_Evolution/phylogeny/CeMM_build_global_new_v19/results/branch_lengths.json \
              /scratch/lab_bergthaler/2020_SARS-CoV-2_Evolution/phylogeny/CeMM_build_global_new_v19/results/traits.json \
              /scratch/lab_bergthaler/2020_SARS-CoV-2_Evolution/phylogeny/CeMM_build_global_new_v19/results/nt_muts.json \
              /scratch/lab_bergthaler/2020_SARS-CoV-2_Evolution/phylogeny/CeMM_build_global_new_v19/results/aa_muts.json \
	      /scratch/lab_bergthaler/2020_SARS-CoV-2_Evolution/phylogeny/CeMM_build_global_new_v19/results/clades.json \
  --colors /scratch/lab_bergthaler/2020_SARS-CoV-2_Evolution/phylogeny/CeMM_build_global_new_v19/config/colors.tsv \
  --lat-longs /scratch/lab_bergthaler/2020_SARS-CoV-2_Evolution/phylogeny/CeMM_build_global_new_v19/config/lat_longs.tsv \
  --auspice-config /scratch/lab_bergthaler/2020_SARS-CoV-2_Evolution/phylogeny/CeMM_build_global_new_v19/config/auspice_config.json \
  --output /scratch/lab_bergthaler/2020_SARS-CoV-2_Evolution/phylogeny/CeMM_build_global_new_v19/auspice/SARSCoV2_NextstrainAustria.json