#!/bin/bash
#SBATCH --output ./LoFreq.%j.%N.out 
#SBATCH --error ./LoFreq.%j.%N.err
#SBATCH --job-name=runBSF0757_LoFreq
#SBATCH --partition=longq
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=7
#SBATCH --time=70:00:00
#SBATCH --mem=100000


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

# # create a mpileup of all positions in the genome based on all the bam files
bcftools mpileup --threads 3 -Ov -f ./RefSeq_sequence_Wuhan-Hu1.fa -q 0 -Q 15 -a DP,AD,ADF,ADR,SP -I -d 100000000  -o ./BSF_0766_HLYWGDRXX_MPileup.vcf -b ./BSF_0766_HLYWGDRXX_viralBamListForLoFreq.txt
bgzip ./BSF_0766_HLYWGDRXX_downSample_MPileup.vcf
tabix -p vcf ./BSF_0766_HLYWGDRXX_MPileup.vcf.gz

#LoFreq has been run on each file independently in the VirSeq_VariantCaller.py
#merge multiple vcf files into one
bcftools merge -i "DP:sum,DP4:sum,AF:max,SB:max" --force-samples -m none -O v ./results_pipeline/*/*_samp.lofreq.bed.vcf.gz > ./runBSF0757_merge_samp_lofreq.vcf
bgzip -f ./runBSF0757_merge_samp_lofreq.vcf
tabix -f -p vcf ./runBSF0757_merge_samp_lofreq.vcf.gz

#annotate based on samtools mpileup depth of coverage info
#normalize the indels. Left-align indels. -m is for splitting the multiallelic sites into biallelic records
# have to merge and unmerge entries to get DP adn AD2 for all unannotated samples and loci (else just first occurence)
bcftools annotate -a ./runBSF0757_MPileup_GoodNames.vcf.gz --collapse all -c "+FORMAT/DP,FORMAT/AD2:=FORMAT/AD" ./runBSF0757_merge_samp_lofreq_GoodNames.vcf.gz | bcftools norm -f ./RefSeq_sequence_Wuhan-Hu1.fa -m+any | bcftools norm -f ./RefSeq_sequence_Wuhan-Hu1.fa -m-any > ./runBSF0757_all_samp_norm_lofreq.vcf

bcftools norm -f ./RefSeq_sequence_Wuhan-Hu1.fa -m+any  ./runBSF0757_all_samp_norm_lofreq.vcf >  ./runBSF0757_lofreq2_samp_norm.vcf

#filtering on AF and annotation with snpeff
bcftools view -i "FORMAT/AF>0.001" ./runBSF0757_lofreq2_samp_norm.vcf | bcftools norm -f ./RefSeq_sequence_Wuhan-Hu1.fa -m+any  > ./runBSF0757_lofreq2_samp_norm_AF0001.vcf

#annotate with SNPEff
java -jar /home/apopa/Desktop/work/software/snpEff/snpEff.jar -v sarsCov2_lukas -c /home/apopa/Desktop/work/software/snpEff/snpEff.config -dataDir /home/apopa/Desktop/work/software/snpEff/data -no-intergenic -no "INTRAGENIC" -no-downstream -no-upstream -stats ./runBSF0757_SNPEFF_stats.html ./runBSF0757_lofreq2_samp_norm_AF0001.vcf > ./runBSF0757_lofreq2_samp_norm_AF0001_SNPEff.vcf

bcftools norm  -f ./RefSeq_sequence_Wuhan-Hu1.fa -m-any ./runBSF0757_lofreq2_samp_norm_AF0001_SNPEff.vcf | bcftools norm  -f ./RefSeq_sequence_Wuhan-Hu1.fa -m+snps | bcftools view -i 'TYPE="snp"' >  ./runBSF0757_lofreq2_snpeff_snp_only.vcf

#extract the informations from the vcf file
#cat ./runBSF0757_lofreq2_samp_norm_AF0001_SNPEff.vcf | /home/apopa/Desktop/work/software/snpEff/scripts/vcfEffOnePerLine.pl | java -jar /home/apopa/Desktop/work/software/snpEff/SnpSift.jar extractFields - CHROM POS REF AF "ANN[*].GENE" "ANN[*].ALLELE" "ANN[*].EFFECT" "ANN[*].AA" "ANN[*].IMPACT" > ./runBSF0757_lofreq2_samp_norm_AF0001_OneLine.tab
#java -jar /home/apopa/Desktop/work/software/snpEff/SnpSift.jar extractFields  -s "," -e "." ./runBSF0757_lofreq2_samp_norm_AF0001_SNPEff.vcf CHROM POS REF AF "ANN[*].GENE" "ANN[*].ALLELE" "ANN[*].EFFECT" "ANN[*].AA" > ./runBSF0757_lofreq2_samp_norm_AF0001_snpeff.tab

java -jar /home/apopa/Desktop/work/software/snpEff/SnpSift.jar extractFields ./runBSF0757_lofreq2_samp_norm_AF0001_SNPEff.vcf CHROM POS REF "ANN[0].GENE" "ANN[0].ALLELE" "ANN[0].EFFECT" "ANN[0].AA" "ANN[1].ALLELE" "ANN[1].EFFECT" "ANN[1].AA" "ANN[2].GENE" "ANN[2].ALLELE" "ANN[2].EFFECT" "ANN[2].AA"> ./runBSF0757_lofreq2_samp_norm_AF0001_snpeff.tab
bcftools query -Hf "%CHROM\t%POS\t%REF\t%ALT[\t%AF\t%DP\t%SB\t%PQ]\n"  ./runBSF0757_lofreq2_samp_norm_AF0001.vcf > ./runBSF0757_lofreq2_samp_norm_AF0001.stats.tab
bcftools query -Hf "%CHROM\t%POS\t%REF\t%ALT[\t%AF]\n"  ./runBSF0757_lofreq2_samp_norm_AF0001.vcf > ./runBSF0757_lofreq2_samp_norm_AF0001.vcf.afs.tab
bcftools query -Hf "%CHROM\t%POS\t%REF\t%ALT[\t%DP]\n"  ./runBSF0757_lofreq2_samp_norm_AF0001.vcf > ./runBSF0757_lofreq2_samp_norm_AF0001.vcf.dp.tab
bcftools query -Hf "%CHROM\t%POS\t%REF\t%ALT[\t%PQ]\n"  ./runBSF0757_lofreq2_samp_norm_AF0001.vcf > ./runBSF0757_lofreq2_samp_norm_AF0001.vcf.pq.tab

paste ./runBSF0757_lofreq2_samp_norm_AF0001.vcf.afs.tab <(cut -f5- ./runBSF0757_lofreq2_samp_norm_AF0001.vcf.dp.tab) <(cut -f5- ./runBSF0757_lofreq2_samp_norm_AF0001.vcf.pq.tab) >  ./runBSF0757_lofreq2_samp_norm_AF0001.afs.dp.pq.tab
