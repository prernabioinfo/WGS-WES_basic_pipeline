### This is a Variant Calling Pipeline 
#### Usage
```
conda env create -f variant_analysis.yml
conda activate variant_analysis
python variant_calling_pipeline.py --help
python variant_calling_pipeline.py -1 /path/to/directory/sample1_1.fastq.gz -2 /path/to/directory/sample1_2.fastq.gz -g /path/to/the/genome/genome.fna -o /path/to/output/directory/directory_name --threads 16
or
python variant_calling_pipeline.py -sample-sheet sample_list.tsv -g /path/to/the/genome/genome.fna -o /path/to/output/directory/directory_name --threads 16

```
#### Input file format
File should contain 2 columns, if you are giving file as input. The first column should contain path to read1 file and the read1 file name and the second column should contain rea2 path along with the file name. 

```
samples.tsv (tab/space/comma separated, col1=R1, col2=R2):
/data/nili_ravi_1.fastq.gz    /data/nili_ravi_2.fastq.gz
/data/murrah_1.fastq.gz     /data/murrah_2.fastq.gz

```
