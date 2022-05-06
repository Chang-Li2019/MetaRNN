# MetaRNN for pathogenicity prediction of nsSNVs and nfINDELs

## Quick start guide
1. Download file package (file size: 13G):
    ```
    wget -c https://usf.box.com/shared/static/46fsvir38fd7uxcnsdwayde2mga34hl1 -O MetaRNN.tar.gz 
    ```
    OR open the save link in any web browser: https://usf.box.com/shared/static/46fsvir38fd7uxcnsdwayde2mga34hl1
  
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
