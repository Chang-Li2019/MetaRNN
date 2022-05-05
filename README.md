# MetaRNN for pathogenicity prediction of nsSNVs and nfINDELs

## Installation
1. Download file package:
  ```
  wget -c -O MetaRNN.tar.gz 
  ```
  OR open the link in any web browser
  
2. Extract downloaded folder:
  ```
  tar -xf MetaRNN.tar.gz
  cd MetaRNN
  ```
  
3. create conda environment with required packages and activate the virtual environment:
  ```
  conda env create -f environment.yml
  conda activate MetaRNN
  ```
5. Run example VCF file in hg38 genome assembly:
  ```
  python ./MetaRNN hg38 test.vcf
  ```
