#!/usr/bin/env python
""" clean, align, variant call pipeline """
""" started on 24.03.2020 """
""" Alexandra Popa, Lukas Endler enriched from https://covid19.galaxyproject.org/genomics/4-variation/#workflows """
"""  it needs:  module load htslib/1.9
				module load samtools/1.9
				module unload gcc/6.1.0
				module load gcc/7.1.0
				module load fastqc/0.11.8
				module load python/3.7.2
				module load bwa/0.7.17
				module load bcftools/1.9
				module load vcftools/0.1.13
"""
""" 1. BBMerge (do not remove adapters)"""
""" 2. FASTQC"""
""" 2'. Repair pairing if needed"""
""" 3. BWA-MEM pipeline """
""" 4. MultiQC stats """
""" 5. Extract viral sequences """
""" 6. Get Coverages """
""" 7. iVar to correct for the primers """
""" 8. Realignment Viterbi """
""" 9. Mark indel qualities """
""" 10. Call the variants with LoFreq"""

from argparse import ArgumentParser
import os
import sys
import subprocess
import re
import pypiper
#from envbash import load_envbash
#load_envbash('/scratch/lab_bergthaler/pipelines/virSeqVariantCaller/moduleLoad.sh')


########################
### Define software path ###
########################

#for trimming
BBDUK = './software/bbmap/bbduk.sh'
#for read length normalization
BBNORM = './software/bbmap/bbnorm.sh'
#for error correction of overalaping regions
BBMERGE = './software/bbmap/bbmerge.sh'
#for de novo assembly
SPADES = './software/SPAdes-3.14.0-Linux/bin/spades.py'
#for checking the assembly and making stats
QUAST = './software/quast-5.0.2/quast.py'
#bamUtils for clipping overlaping edges
bamUTIL = './software/bamUtil/bin/bam'
#for primer removal from bam
iVar = './software/miniconda3/bin/ivar'
#LoFREQ
LOFREQ = './software/lofreq_star-2.1.2/bin/lofreq'
#SEQTK
SEQTK = './software/seqtk/seqtk'
#this file contains both primers and adapters for different sequencing techniques
ADAPTERS = './software/bbmap/resources/adapters.fa'
#the sequences of the 98 set of primers used to amplify the viral sequences
PRIMERS = './genomes/SARS-CoV-2/artic-ncov2019-master/primer_schemes/nCoV-2019/V2/nCoV-2019.primers.fa'
#the sequences of the 98 set of primers used to amplify the viral sequences
PRIMERSBED = './genomes/SARS-CoV-2/artic-ncov2019-master/primer_schemes/nCoV-2019/V2/nCoV-2019.scheme_modif.bed'
#Genome
refGENOME = "./genomes/hg38_SARSCoV2/indices_for_BWA/hg38_sars_cov2.fa"
refSARSCoV2 = "./genomes/SARS-CoV-2/RefSeq_sequence_Wuhan-Hu1.fa"
#Annotations
refANNOT = "./genomes/hg38_sars_cov2/NC_045512.2.gff3"
#NC_045512.2.rewrite.gtf


########################
### Argument Parsing ###
########################
parser = ArgumentParser(description='Pypiper arguments VirSeq.')
#parser = pypiper.add_pypiper_args(parser, all_args=True)

parser.add_argument('-y', '--sample_yaml', dest='config_file', help='yaml config file with sample attributes')
parser.add_argument('-dp', '--data_path', dest='data_path', help='path to sequencing data file')
parser.add_argument('-n', '--sample_name', dest='sample_name', help='name of the sample')
parser.add_argument('-r', dest='results_folder', help='path to folder to store results in')
parser.add_argument('-fc', dest='flowcell', help='the flow cell id')

args = parser.parse_args()

##################
### Initialize ###
##################

outfolder = os.path.abspath(os.path.join(args.results_folder,args.sample_name))


pm = pypiper.PipelineManager(name = 'VirSeq', outfolder = outfolder, args = args) # initialize pipeline manager instance

################
### To Fastq ###
################

pm.timestamp('### BAM to FASTQ: ')


trimAdapPair_1_fq = args.sample_name + '_npa_1.fastq.gz'
pathAdapTrimPair_1_fq = os.path.join(outfolder, trimAdapPair_1_fq)
trimAdapPair_2_fq = args.sample_name + '_npa_2.fastq.gz'
pathAdapTrimPair_2_fq = os.path.join(outfolder, trimAdapPair_2_fq)

if(os.path.isfile(pathAdapTrimPair_1_fq) & os.path.isfile(pathAdapTrimPair_2_fq  )):
	pm.timestamp('### FASTQ files already exist!')
else:
	cmd = 'java -Xmx2g -jar  /cm/shared/apps/picard-tools/2.6.0/picard.jar SamToFastq INPUT=' + args.data_path + ' FASTQ=' + pathAdapTrimPair_1_fq + '  SECOND_END_FASTQ=' + pathAdapTrimPair_2_fq
	pm.run(cmd, pathAdapTrimPair_1_fq)



# ################
# ### Repair pairing ###
# ################

# # in the previous process some reads might have gotten discarded and the pairs are not corresponding. I will use this bbtools repai function to repair the pairs

# pm.timestamp('### Fix the pairs correspondance: ')
# repairPair_1_fq = args.sample_name + '_repair_1.fq'
# pathRepairPair_1_fq = os.path.join(outfolder, repairPair_1_fq)
# repairPair_2_fq = args.sample_name + '_repair_2.fq'
# pathRepairPair_2_fq = os.path.join(outfolder, repairPair_2_fq)
# repairSingle = args.sample_name + '_single.fq'
# pathRepairSingle = os.path.join(outfolder, repairSingle)

# cmd = "/home/apopa/Desktop/work/software/bbmap/repair.sh in1=" + pathFirstRead_fq + " in2=" + pathSecondRead_fq + " out1=" + pathRepairPair_1_fq + " out2=" + pathRepairPair_2_fq + " outsingle=" + pathRepairSingle
# pm.run(cmd, pathRepairSingle)


################
### Error correct the overlaping pair regions ###
################
pm.timestamp('### Error correct overlap: ')

trimCorrOver_1_fq = args.sample_name + '_ecco_1.fq'
pathCorrOver_1_fq = os.path.join(outfolder, trimCorrOver_1_fq)
trimCorrOver_2_fq = args.sample_name + '_ecco_2.fq'
pathCorrOver_2_fq = os.path.join(outfolder, trimCorrOver_2_fq)

cmd = BBMERGE+" in1=" + pathAdapTrimPair_1_fq+ " in2="+pathAdapTrimPair_2_fq+ " out="+pathCorrOver_1_fq +" out2="+pathCorrOver_2_fq+" mininsert=200 ecco mix ordered=t qtrim2=r "
pm.run(cmd, pathCorrOver_2_fq)

################
### BWA-MEM ###
################
pm.timestamp('### BWA-MEM: ')

bwaOutputBAM = args.sample_name + '_bwa.bam'
pathBWAOut = os.path.join(outfolder, bwaOutputBAM)
bwaSortedBAM = args.sample_name + '_sorted.bam'
pathSortedBAM = os.path.join(outfolder, bwaSortedBAM)
bwaErrorLog = args.sample_name + '_bwa.err.log'
pathBWAErrorLog = os.path.join(outfolder, bwaErrorLog)
tmpBamBWA = args.sample_name + 'tmpBWA'
pathBWATmp = os.path.join(outfolder, tmpBamBWA)

### make the read group for the bwa call

# ID = args.flowcell + '_' + args.sample_name + '_npa_npp'
# LB = args.flowcell
# SM = args.sample_name.split('_')[0] + args.sample_name.split('_')[1]
# PL = "ILLUMINA"
# PU = "HL37TBBXX.5"

cmd = "bwa mem -k 17 -r 1.25 -M -t 17 " + refGENOME + " " + pathCorrOver_1_fq + " " + pathCorrOver_2_fq + " > " + pathBWAOut

pm.timestamp(cmd)
pm.run(cmd, pathBWAOut)

### samtools to sort the bwa file

cmd = " samtools view -Shb " + pathBWAOut +  "| samtools sort -T " + pathBWATmp + " > " + pathSortedBAM

pm.timestamp(cmd)
pm.run(cmd, pathSortedBAM)

################
### Index bam file ###
################
pm.timestamp('### Index Bam file: ')

indexSortedBAM = args.sample_name + '_sorted.bam.bai'
pathIndexSortedBAM = os.path.join(outfolder, indexSortedBAM)
cmd = "samtools index " + pathSortedBAM
pm.timestamp(cmd)
pm.run(cmd, pathIndexSortedBAM)

################
### Stats mapping ###
################
pm.timestamp('### Stats mapping: ')

flagstatBAM = args.sample_name + '.flagstat'
pathFlagStat = os.path.join(outfolder, flagstatBAM)
cmd = "samtools flagstat " + pathSortedBAM + " >  " + pathFlagStat
pm.timestamp(cmd)
#pm.run(cmd, pathFlagStat)

idxstatsBAM = args.sample_name + '.idxstats'
pathIdxStat = os.path.join(outfolder, idxstatsBAM)
cmd = "samtools idxstats " + pathSortedBAM + " >  " + pathIdxStat
pm.timestamp(cmd)
pm.run(cmd, pathIdxStat)

fullStatsBAM = args.sample_name + '.stats'
pathFullStat = os.path.join(outfolder, fullStatsBAM)
cmd = "samtools stats " + pathSortedBAM + " >  " + pathFullStat
pm.timestamp(cmd)
pm.run(cmd, pathFullStat)


################
### Extract virus sequences ###
################
pm.timestamp('### Extract virus sequences: ')
virusBAM = args.sample_name + '_sorted_viral.bam'
pathVirusBAM = os.path.join(outfolder, virusBAM)
cmd = "samtools view -bh -f 2 -F 256 " + pathSortedBAM + " NC_045512.2 > " + pathVirusBAM
pm.timestamp(cmd)
pm.run(cmd, pathVirusBAM)

indexVirusBAM = args.sample_name + '_sorted_viral.bam.bai'
pathIndexVirusBAM = os.path.join(outfolder, indexVirusBAM)
cmd = "samtools index " + pathVirusBAM
pm.timestamp(cmd)
pm.run(cmd, pathIndexVirusBAM)

################
### Reheader ###
################
pm.timestamp('### Reheader virus bam: ')

virusReheaderBAM = args.sample_name + '_sorted_viral_rh.bam'
pathVirusReheaderBAM = os.path.join(outfolder, virusReheaderBAM)
cmd = "samtools view -h " + pathVirusBAM + " | grep -Ev \'@SQ.*SN:chr[0-9]|Un|chrX|chrY|chrM|chrEBV|NT_|NW_|NC_000|NC_012\' | samtools view -b > " + pathVirusReheaderBAM
pm.timestamp(cmd)
pm.run(cmd, pathVirusReheaderBAM)
### I don't know why this command is not working and giving a truncated bam file
#samtools view -h SRR11314339_sorted.bam | grep -Ev '@SQ.*SN:chr[0-9]|Un|chrX|chrY|chrM|chrEBV' > test.bam

indexVirusReheaderBAM = args.sample_name + '_sorted_viral_rh.bam.bai'
pathIndexVirusReheaderBAM = os.path.join(outfolder, indexVirusReheaderBAM)
cmd = "samtools index " + pathVirusReheaderBAM  
pm.timestamp(cmd)
pm.run(cmd, pathIndexVirusReheaderBAM)

################
### iVar for primer masking ###
################
pm.timestamp('### Trim primers: ')

virusTrimmBAM = args.sample_name + '_trim.bam'
pathVirusTrimmBAM = os.path.join(outfolder, virusTrimmBAM)
pathTmpMv = os.path.join(outfolder, virusTrimmBAM)


if os.path.isfile(pathVirusTrimmBAM) :
    pm.timestamp('### Trim file already exists exist!')
else:
	cmd = iVar+" trim -b "+PRIMERSBED+" -s 4 -q 0 -m 20 -i "+pathVirusReheaderBAM+" -p " + args.sample_name + "_trim"
	pm.timestamp(cmd)
	pm.run(cmd, pathVirusTrimmBAM)
	cmd = "mv " + virusTrimmBAM + "  " + pathVirusTrimmBAM
	pm.run(cmd, pathBWAErrorLog)

tmpBamTrim = args.sample_name + 'tmpTrim'
pathTrimTmp = os.path.join(outfolder, tmpBamTrim)
virusTrimmSortBAM = args.sample_name + '_sorted_viral_rh_trim.bam'
pathVirusTrimmSortBAM = os.path.join(outfolder, virusTrimmSortBAM)
cmd = " samtools view -Shb " + pathVirusTrimmBAM +  " | samtools sort -T " + pathTrimTmp + " > " + pathVirusTrimmSortBAM
pm.timestamp(cmd)
pm.run(cmd, pathVirusTrimmSortBAM)

indexVirusTrimmSortBAM = args.sample_name + '_sorted_viral_rh_trim.bam.bai'
pathIndexVirusTrimmSortBAM = os.path.join(outfolder, indexVirusTrimmSortBAM)
cmd = "samtools index " + pathVirusTrimmSortBAM 
pm.timestamp(cmd)
pm.run(cmd, pathIndexVirusTrimmSortBAM)


################
### Stats mapping Virus ###
################
pm.timestamp('### Stats mapping virus: ')
flagstatVirusBAM = args.sample_name + '_virus.flagstat'
pathFlagStatVirus = os.path.join(outfolder, flagstatVirusBAM)
cmd = "samtools flagstat " + pathVirusTrimmSortBAM  + " >  " + pathFlagStatVirus
pm.timestamp(cmd)
pm.run(cmd, pathFlagStatVirus)

idxstatsVirusBAM = args.sample_name + '_virus.idxstats'
pathIdxStatVirus = os.path.join(outfolder, idxstatsVirusBAM)
cmd = "samtools idxstats " + pathVirusTrimmSortBAM + " >  " + pathIdxStatVirus
pm.timestamp(cmd)
pm.run(cmd, pathIdxStatVirus)

fullStatsVirusBAM = args.sample_name + '_virus.stats'
pathFullStatVirus = os.path.join(outfolder, fullStatsVirusBAM)
cmd = "samtools stats " + pathVirusTrimmSortBAM + " >  " + pathFullStatVirus
pm.timestamp(cmd)
pm.run(cmd, pathFullStatVirus)

################
### Get Coverages ###
################

pm.timestamp('### Get Coverages: ')
coverageFileBw = args.sample_name + '.bw'
pathCoverageBw = os.path.join(outfolder, coverageFileBw)
coverageFileBed = args.sample_name + '.bed'
pathCoverageBed = os.path.join(outfolder, coverageFileBed)


#this bamCoverage from deeptools creates a nice bwa coverage file
cmd = "bamCoverage --bam " + pathVirusTrimmSortBAM + " -o "  + pathCoverageBw + " --binSize 10"
pm.timestamp(cmd)
pm.run(cmd, pathCoverageBw)
cmd = "bamCoverage --bam " + pathVirusTrimmSortBAM + " -o "  + pathCoverageBed + " --outFileFormat bedgraph --binSize 10"
pm.timestamp(cmd)
pm.run(cmd, pathCoverageBed)

################
### Clip overlap ###
################
# when the fragment size and the reads lengths are very close

pm.timestamp('### Clip overlaps: ')
clipViralBam = args.sample_name + '_clip_viral.bam'
pathClipVirus = os.path.join(outfolder, clipViralBam)
statsClip = args.sample_name + '_clip_viral_stats.txt'
pathStatsClip = os.path.join(outfolder, statsClip)
cmd = bamUTIL + " clipOverlap --in " + pathVirusTrimmSortBAM + " --out " + pathClipVirus + " --stats " + pathStatsClip
pm.timestamp(cmd)
pm.run(cmd, pathClipVirus)

##index
clipIndexFile = args.sample_name + '_clip_viral.bam.bai'
pathIdxClip = os.path.join(outfolder, clipIndexFile)
cmd = "samtools idxstats " + pathClipVirus
pm.timestamp(cmd)
pm.run(cmd, pathIdxClip)


################
### Create FASTA consensus from BAM ###
################
pm.timestamp('### FastaConsensus: ')

maskBed = args.sample_name + '_mask.bed'
pathMaskBed = os.path.join(outfolder, maskBed)
maskFasta = args.sample_name + '_ref_mask.fasta'
pathMaskFasta = os.path.join(outfolder, maskFasta )
fqConsensus = args.sample_name + '_cons.fq'
pathFQConsensus = os.path.join(outfolder, fqConsensus)

# cmd = 'awk \'{if($4 < 100){print $0}}\' ' + pathCoverageBed + " > " + pathMaskBed
# pm.timestamp(cmd)
# pm.run(cmd, pathMaskBed)

# cmd = "bedtools maskfasta -fi " + refSARSCoV2 + " -bed " + pathMaskBed + " -fo " + pathMaskFasta
# pm.timestamp(cmd)
# pm.run(cmd, pathMaskFasta)

#cmd= "samtools mpileup -d 0 -f " + refSARSCoV2 + " " +  pathVirusTrimmSortBAM + " | bcftools call -c | vcfutils.pl vcf2fq > " + pathFQConsensus
#cmd= "bcftools mpileup -d 100000 -q 20 -Q 10 -f " + refSARSCoV2 + " " +  pathVirusTrimmSortBAM + " | bcftools call -c | vcfutils.pl vcf2fq > " + pathFQConsensus
#cmd= "bcftools mpileup -d 200000 -q 20 -Q 10 -f " + pathMaskFasta + " " +  pathVirusTrimmSortBAM + " | bcftools call -c | vcfutils.pl vcf2fq > " + pathFQConsensus
cmd= "samtools mpileup -d 0 -uf " + refSARSCoV2 + " " +  pathClipVirus + " | bcftools call -c | vcfutils.pl vcf2fq > " + pathFQConsensus
pm.timestamp(cmd)
pm.run(cmd, pathFQConsensus)

#bcftools mpileup -d 50000 -f refSARSCoV2 -L 10000 -q 20 -Q 10 pathVirusTrimmSortBAM | bcftools call -c | bcftools consensus -f refSARSCoV2 > 
#bcftools mpileup -d 50000 -f /scratch/lab_bergthaler/genomes/SARS-CoV-2/RefSeq_sequence_Wuhan-Hu1.fa -L 10000 -q 20 -Q 10 ../results_pipeline/CoV_023_S63673_S120/CoV_023_S63673_S120_sorted_viral_rh_trim.bam | bcftools call -c | bcftools consensus -f /scratch/lab_bergthaler/genomes/SARS-CoV-2/RefSeq_sequence_Wuhan-Hu1.fa > test.fasta
#bcftools mpileup -d 50000 -f /scratch/lab_bergthaler/genomes/SARS-CoV-2/RefSeq_sequence_Wuhan-Hu1.fa -L 10000 -q 20 -Q 10 ../results_pipeline/CoV_023_S63673_S120/CoV_023_S63673_S120_sorted_viral_rh_trim.bam | bcftools call -c --ploidy 1 -O z -o test.vcf.gz
#bcftools consensus -e "QUAL < 50" -I -f /scratch/lab_bergthaler/genomes/SARS-CoV-2/RefSeq_sequence_Wuhan-Hu1.fa -o test.cons.fasta test.vcf.gz

#samtools mpileup -uf -d 0 /scratch/lab_bergthaler/genomes/SARS-CoV-2/RefSeq_sequence_Wuhan-Hu1.fa 

#https://www.biostars.org/p/367626/
fastaConsensus = args.sample_name + '_cons.fasta'
pathFastaConsensus = os.path.join(outfolder, fastaConsensus)

#[apopa@n001 statsFasta]$ tabix -f -p vcf test.vcf.gz 
#[apopa@n001 statsFasta]$ bcftools consensus -e "QUAL < 50" -I -f /scratch/lab_bergthaler/genomes/SARS-CoV-2/RefSeq_sequence_Wuhan-Hu1.fa -o test.cons.fasta test.vcf.gz

cmd= SEQTK + " seq -aQ64 -q20 -n N " + pathFQConsensus + " > " +  pathFastaConsensus
pm.timestamp(cmd)
pm.run(cmd, pathFastaConsensus)



################
### LoFreq - realignment ###
################

pm.timestamp('### LoFreq viterbi: ')
loFreqViterbiFile = args.sample_name + '_real_viterbi.bam'
pathViterbiFile = os.path.join(outfolder, loFreqViterbiFile)
tmpViterbi = "tmp_" + args.sample_name + '.bam'
pathTmpViterbi = os.path.join(outfolder, tmpViterbi)
cmd = LOFREQ + " viterbi -f " +  refSARSCoV2 + " " + pathClipVirus + " | samtools sort -T " + pathTmpViterbi + " - > " + pathViterbiFile
pm.timestamp(cmd)
pm.run(cmd, pathViterbiFile)

##index
viterbiIndexFile = args.sample_name + '_real_viterbi.bam.bai'
pathIdxViterbi = os.path.join(outfolder, viterbiIndexFile)
cmd = "samtools idxstats " + pathViterbiFile + " >  " + pathIdxViterbi
pm.timestamp(cmd)
pm.run(cmd, pathIdxViterbi)


################
### LoFreq - mark indel qualities ###
################

pm.timestamp('### LoFreq indelqual: ')
loFreqIndelQualFile = args.sample_name + '_real_viterbi_indelQual.bam'
pathIndelQualFile = os.path.join(outfolder, loFreqIndelQualFile)
cmd = LOFREQ + " indelqual --dindel -f " +  refSARSCoV2 + " -o " + pathIndelQualFile + " " + pathViterbiFile
pm.timestamp(cmd)
pm.run(cmd, pathIndelQualFile)

##index
indelQualIndexFile = args.sample_name + '_real_viterbi_indelQual.bam.bai'
pathIdxIndelQual = os.path.join(outfolder, indelQualIndexFile)
cmd = "samtools idxstats " + pathIndelQualFile + " >  " + pathIdxIndelQual
pm.timestamp(cmd)
pm.run(cmd, pathIdxViterbi)



################
### LoFreq - call variants ###
################
pm.timestamp('### LoFreq: ')
loFreqFile = args.sample_name + '_lofreq.vcf'
pathLoFreqFile = os.path.join(outfolder, loFreqFile)
loFreqFilter = args.sample_name + '_lofreq_filter.vcf'
pathLoFreqFilter = os.path.join(outfolder, loFreqFilter)
loFreqComp = args.sample_name + '_lofreq_filter.vcf.gz'
pathLoFreqComp = os.path.join(outfolder, loFreqComp)
loFreqIdx = args.sample_name + '_lofreq_filter.vcf.gz.tbi'
pathLoFreqIdx = os.path.join(outfolder, loFreqIdx)

#- q Skip any base with baseQ smaller
#-Q Skip alternate bases with baseQ smaller
#-m Skip reads with mapping qualirty smaller
#-C Test only positions having at least this coverage
#-a P-Value cutoff / significance level
#cmd = LOFREQ + " call -f " +refSARSCoV2 +  " -q 20 -Q 20 -m 20 -C 75 -a 0.05 --call-indels -o " + pathLoFreqFile + " " + pathIndelQualFile 
#cmd = LOFREQ + " call -f " +refSARSCoV2 +  " -q 0 -Q 15 -B -m 20 -C 75 -a 0.05 --call-indels -o " + pathLoFreqFile + " " + pathIndelQualFile 
cmd = LOFREQ + " call -f " +refSARSCoV2 +  "  -C 75  --no-default-filter  --call-indels -o " + pathLoFreqFile + " " + pathIndelQualFile 
pm.timestamp(cmd)
pm.run(cmd, pathLoFreqFile)


#filter on B (maximum phred-value allowed)
#cmd = LOFREQ + " filter -i " + pathLoFreqFile + " -B 30 --print-all -o " + pathLoFreqFilter
cmd = LOFREQ + " filter -i " + pathLoFreqFile + " -v 75 -a 0.01 --no-defaults -Q 90 --print-all | bcftools filter -m + -s \"HRUN_gt_3\" -e \"HRUN > 4\" > " + pathLoFreqFilter
#cmd = LOFREQ + " filter -i " + pathLoFreqFile + " -v 75 -a 0.01 --no-defaults -Q 90 | bcftools filter -m + -s \"HRUN_gt_3\" -e \"HRUN > 4\" > " + pathLoFreqFilter
pm.timestamp(cmd)
pm.run(cmd, pathLoFreqFilter)

# cmd = "bgzip -f " + pathLoFreqComp
# pm.run(cmd, pathLoFreqComp)

# cmd = "tabix -f -p vcf " + pathLoFreqIdx
# pm.run(cmd, pathLoFreqIdx)

################
### Variants to BED formating ###
################
pm.timestamp('### Prepare output files: ')

#splitSampName = args.sample_name.split("_", 2)
#sampShort = splitSampName[0]+splitSampName[1]
loFreqRename = args.sample_name + "_samp.lofreq.bed.vcf"
pathLoFreqRename = os.path.join(outfolder, loFreqRename)
loFreqRenameZip = args.sample_name + "_samp.lofreq.bed.vcf.gz"
pathLoFreqRenameZip = os.path.join(outfolder, loFreqRenameZip)
loFreqRenIx = args.sample_name + "_samp.lofreq.bed.vcf.gz.tbi"
pathLoFreqRenIx = os.path.join(outfolder, loFreqRenIx)

cmd = "awk -v OFS=\"\t\" \'/^\\#\\#[^I]/ {print} /^\\#\\#INFO/ {sub(\"AF,Number=1\",\"AF,Number=A\",$0); print $0; sub(\"INFO\",\"FORMAT\",$0); print $0; print \"##FORMAT=<ID=PQ,Number=1,Type=Integer,Description=\\\"Phred-scaled variant call P value\\\">\" } /\\#CH/ {print $0,\"FORMAT\",\"\'"+ args.sample_name + "\'\"} !/^\\#/ {form=$8;  gsub(/=[^A-Z]+/,\":\",form); gsub(/;/,\":\",form); sub(/:$/,\"\",form); sub(/INDEL:/,\"\",form);  samp=$8; gsub(/[A-Z4]+=/,\"\",samp); gsub(/;/,\":\",samp); sub(/INDEL:/,\"\",samp); print $0,form\":PQ\",samp\":\"$6}\' " + pathLoFreqFilter + " > " + pathLoFreqRename
#+" | bgzip -c >" +pathLoFreqRename 
pm.timestamp(cmd)
pm.run(cmd, pathLoFreqRename)
cmd = "bgzip -f " + pathLoFreqRename
pm.timestamp(cmd)
pm.run(cmd, pathLoFreqRenameZip)
cmd = "tabix -f -p vcf " + pathLoFreqRenameZip
pm.timestamp(cmd)
pm.run(cmd, pathLoFreqRenIx)

################
### DOWNsample bam file ###
################
pm.timestamp('### Downsample: ')
downSampleBam = args.sample_name + '_real_viterbi_indelQual_max_cov_50000.bam'
pathDownSample = os.path.join(outfolder, downSampleBam)
pathDownSampleIx = os.path.join(outfolder, args.sample_name + '_real_viterbi_indelQual_max_cov_50000.bam.bai')

cmd = 'python ./software/downsample_to_cov.py -b ' + pathIndelQualFile +' -c 50000 --bedtools ./software/bedtools2/bin/bedtools'
pm.timestamp(cmd)
pm.run(cmd, pathDownSample)
cmd = 'samtools index '+ pathDownSample
pm.run(cmd, pathDownSampleIx)

pm.stop_pipeline()
