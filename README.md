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
