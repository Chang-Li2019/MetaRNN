# MetaRNN for pathogenicity prediction of nsSNVs and nfINDELs

## Quick start guide
1. Download file package (file size: 13G):
    ```
    wget -c https://usf.box.com/shared/static/c4656df0rz10vw6z7ra11mkw69ulr9d3 -O MetaRNN.tar.gz 
    ```
    OR open the save link in any web browser: https://usf.box.com/shared/static/c4656df0rz10vw6z7ra11mkw69ulr9d3
    
    OR get it from [figshare](https://figshare.com/articles/software/MetaRNN_Differentiating_Rare_Pathogenic_and_Rare_Benign_Missense_SNVs_and_InDels_Using_Deep_Learning/19742503).
  
2. Extract downloaded folder:
    ```
    tar -xf MetaRNN.tar.gz
    cd MetaRNN
    ```
  
3. Create conda environment with required packages and activate the virtual environment:
    ```
    conda env create -f environment.yml
    conda activate MetaRNN
    ```
5. Run example VCF file in hg38 genome assembly:
    ```
    python ./MetaRNN.py hg38 test.vcf
    ```

## Notes

Supplementary Tables and model analyses associated with the publication are available under folders Data and Notebooks.

