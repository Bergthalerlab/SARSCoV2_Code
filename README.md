# SARSCoV2_Code
## Mapping, Genome Reconstruction and Low Frequency Variant Calling
This pipeline is implemented with looper (https://looper.databio.org/en/latest/) and pyPiper (https://pypi.org/project/pyPiper/). 
### SarsVirSeq_VariantCaller.py

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
