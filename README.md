# SARSCoV2_Code
## Mapping, Genome Reconstruction and Low Frequency Variant Calling
authors: Alexandra Popa and Lukas Endler
This pipeline is implemented with looper (https://looper.databio.org/en/latest/) and pyPiper (https://pypi.org/project/pyPiper/). 
### SarsVirSeq_VariantCaller.py
This pipeline starts with raw reads, filters and maps them on the hg38/SARSCoV2 genome, extracts the SARSCoV2 sequences, reconstructs the viral genome, and calls for low frequency variants.

<ol>
<li>BBMerge (do not remove primers)</li>  
<li>FASTQC</li>  
<li>Repair pairing if needed</li> 
<li>BWA-MEM pipeline</li> 
<li>MultiQC stats</li> 
<li>Extract viral sequences</li> 
<li>Get Coverages</li> 
<li>iVar to correct for the primers</li> 
<li>Realignment Viterbi</li> 
<li>Mark indel qualities</li> 
<li>Call the variants with LoFreq</li> 
</ol>

### SarsVirSeq_VariantCaller.yaml
It contains the configuration file of the pipeline specifying the submission criteria to the server

### SARS_LowFreqScript.sh
Script that combines the individual VCF files to call and annotate the low frequency variants.

## Phylogenetic Analysis
Scripts for the reconstruction of global phylogeny (Global), phylogeny of Austrian strains (OnlyAustrian) and phylogeny of early European strains (EarlyEuropean) can be found in respective folders

### Phyloscript.sh
<ol>
<li>Multiple sequence alignment of selected strains</li> 
<li>Building of raw phylogenetic tree</li> 
<li>Refinement of phylogenetic tree</li> 
<li>Identify synonymous and non-synonymous mutations compared to reference sequence</li> 
<li>Create JSON file for visualisation</li> 
</ol>

### subsampling.py
<ol>
<li>Creates a random subset of given sequence selection</li> 
</ol>

### Reformat_fasta_files.py
<ol>
<li>Adjusts formats for Phyloscript.sh</li> 
</ol>
