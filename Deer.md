# Detecting Multiple Paternity in White-tailed Deer

##### A 2bRAD data analysis pipeline walkthrough

Adam Koller - Luther College '24 - 2022

---


```python
# Put packages and versions in notebook
```

### Importing libraries


```python
import rpy2
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import sys
import scipy.stats
```


```python
%load_ext rpy2.ipython
%R .libPaths("~/Programs/R/x86_64-pc-linux-gnu-library/4.1")
```





<span>StrVector with 4 elements.</span>
<table>
<tbody>
  <tr>

    <td>
    '/home/LC/...
    </td>

    <td>
    '/usr/loca...
    </td>

    <td>
    '/usr/lib/...
    </td>

    <td>
    '/usr/lib/...
    </td>

  </tr>
</tbody>
</table>





```python
%%capture importLibraries
%%R
library(ggplot2)
library(ggrepel)
library(sets)
library(dplyr)
```


```python
# Changing working directory (set equal to project folder)
os.environ['WORKDIR'] = './DeerProject'
```

### Adding software to path


```python
path = '/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/opt/anaconda3/bin' # <- insert current paths into this string
path += ':/home/LC/kollad01/Programs/2bRAD_denovo' # Enter path to directory with scripts
path += ':/home/LC/kollad01/.local/bin'
path += ':home/LC/kollad01/Programs/bowtie2-2.4.5-mingw-x86_64'
path += ':home/LC/kollad01/Programs/samtools-1.15.1'
path += ':/home/LC/kollad01/Programs/angsd'
path += ':/home/LC/kollad01/Programs/angsd/NgsRelate'
path += ':/home/LC/kollad01/Programs/htslib'
path += ':/home/LC/kollad01/Programs/R/x86_64-pc-linux-gnu-library/4.1'
os.environ['PATH'] = path
```


```python
print(os.environ['PATH'])
```

    /bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/opt/anaconda3/bin:/home/LC/kollad01/Programs/2bRAD_denovo:/home/LC/kollad01/.local/bin:home/LC/kollad01/Programs/bowtie2-2.4.5-mingw-x86_64:home/LC/kollad01/Programs/samtools-1.15.1:/home/LC/kollad01/Programs/angsd:/home/LC/kollad01/Programs/angsd/NgsRelate:/home/LC/kollad01/Programs/htslib:/home/LC/kollad01/Programs/R/x86_64-pc-linux-gnu-library/4.1


### Understanding Data Structure

200 samples (195 unique, 5 replicates) split among 17 .fastq files.


```bash
%%bash
head Data/ILL-RAD01_S1_L001_R1_001.fastq
```

    @J00102:28:HTWWLBBXX:1:1101:2747:1050 1:N:0:NTCACG
    NGAACCGGCACAACCCAGCGAGGCACCTGCCTTTGGCCAGCCTCACAGAT
    +
    #AAFFJJJJJJJJJJJJJJJJJJAJJJJJJJJFJJJJJJ<JJJJJJJJAJ
    @J00102:28:HTWWLBBXX:1:1101:3072:1050 1:N:0:NTCACG
    NGAACCAGGACGCCGATGGCAACATCGTCGAACTGCGTTGCTGTGAAGAT
    +
    #AA<FFFFJFJJJJJJJJJJAFJJJJJJJJJJJAJJJJFJF-FJJJFJFF
    @J00102:28:HTWWLBBXX:1:1101:3579:1050 1:N:0:NTCACG
    NGAACCGGCTGGCCAAGGGCAGGGGCCTCGCTGGGTTGTGCCTCAGAGAT


50bp chunks: N(18) CGA N(6) TGC N(12) 4bp barcode

### Counting raw reads


```python
barcodes = ['GTGT','AGAC','ACCA','AGTG','CATC','GTGA','TCAG','GCTT','CTAC','TGTC','TCAC','GACT'] # <- Our barcodes
```


```python
# Function to count reads for each barcode w/wo bcg1 enzyme site. Input directory path with .fastq files, list of barcodes, 
# and bcg1 TRUE/FALSE
def count(directoryName, barcodes, bcg1Bool):  
    barCounts = {'Unmatched':0} # initialize dictionary
    for fileName in os.listdir(directoryName): # iterate through list of file paths
        if fileName.endswith(".fastq"): # ignore files not ending in .fastq
            fileID = fileName.split("_")[1] # ID contains well info used to distinguish samples
            file = open(directoryName+fileName, 'r') # open file
            count = 0 # count variable to ready every 4 lines (lines with sequence)
            for line in file: # iterate through lines of file
                if count%4 == 1: # "every 4th line"
                    enzymeSite = line[18:30] # This is the index of the bcg1 enzyme site
                    if bcg1Bool: # If we restrict counting to sequences containing bcg1 site
                        if (enzymeSite[0:3] == 'CGA' and enzymeSite[9:12] == 'TGC') \ # bcg1 sequence
                        or (enzymeSite[0:3] == 'GCA' and enzymeSite[9:12] == 'TCG'): # or its reverse complement
                            barID = fileID + "_" + line[42:46] #42:46 is the position of our 4bp barcode
                            try:
                                barCounts[barID] += 1 # increase count in dictionary with that barcode
                            except:
                                if line[42:46] in barcodes: # make sure it is a valid barcode                      
                                    barCounts[barID] = 1 # Otherwise add barcode to dictionary with count 1
                                else:
                                    barCounts['Unmatched'] += 1 # Invalid barcodes are thrown in unmatched category
                                    
                    # Same thing as above but for case where we are not counting those with bcg1 site
                    else:
                        barID = fileID + "_" + line[42:46] 
                        try:
                            barCounts[barID] += 1
                        except:
                            if line[42:46] in barcodes:                       
                                barCounts[barID] = 1
                            else:
                                barCounts['Unmatched'] += 1
                count += 1
    def writeCountstoFile(outputName):  # function to write count data to csv file
        reads = []
        for ID in barCounts.keys():
            reads.append({"ID":ID,outputName+"_count":barCounts[ID]})
        reads = pd.DataFrame.from_dict(reads)
        reads.to_csv('./Output/'+outputName+'_counts.csv', index=False)
        
    if bcg1Bool: # Naming csv files properly 
        writeCountstoFile('bcg1')
    else:
        writeCountstoFile('raw')
```


```python
%%time
count('Data/',barcodes, False)
```

    ILL-RAD09_S9_L001_R1_001.fastq
    ILL-RAD06_S6_L001_R1_001.fastq
    ILL-RAD17_S17_L001_R1_001.fastq
    ILL-RAD11_S11_L001_R1_001.fastq
    ILL-RAD10_S10_L001_R1_001.fastq
    ILL-RAD02_S2_L001_R1_001.fastq
    ILL-RAD13_S13_L001_R1_001.fastq
    ILL-RAD03_S3_L001_R1_001.fastq
    ILL-RAD15_S15_L001_R1_001.fastq
    ILL-RAD08_S8_L001_R1_001.fastq
    ILL-RAD04_S4_L001_R1_001.fastq
    ILL-RAD16_S16_L001_R1_001.fastq
    ILL-RAD05_S5_L001_R1_001.fastq
    ILL-RAD07_S7_L001_R1_001.fastq
    ILL-RAD01_S1_L001_R1_001.fastq
    ILL-RAD14_S14_L001_R1_001.fastq
    ILL-RAD12_S12_L001_R1_001.fastq
    CPU times: user 5min 57s, sys: 40.5 s, total: 6min 38s
    Wall time: 6min 4s



```python
raw_reads = pd.read_csv('Output/raw_counts.csv')
print(sum(raw_reads['raw_count']), 'total reads (including unmatched)')
```

    343721394 total reads (including unmatched)



```python
def describeCounts(fileName): # Function to give basic statistics for count data
    raw_count = pd.read_csv(fileName) # Read in file
    omitted = ['S17_CTAC','S17_TGTC','S17_TCAC','S17_GACT','Unmatched'] # Count function created barcodes that are not samples
    raw_count.drop(raw_count[raw_count['ID'].isin(omitted)].index, inplace =True) # Exclude previous barcodes
    sns.histplot(data=raw_count, x = raw_count.columns[1], fill=False) # Create histogram with seaborn
    plt.axvline(raw_count.iloc[:,1].median(), color='green', linestyle='dashed') # median vertical line
    plt.axvline(raw_count.iloc[:,1].mean(), color='red', linestyle='dashed') # mean vertical line
    print(raw_count.describe().round()) # show stats
```


```python
describeCounts("Output/raw_counts.csv")
```

           raw_count
    count      200.0
    mean   1659439.0
    std    1431810.0
    min      12717.0
    25%     739973.0
    50%    1311156.0
    75%    1976969.0
    max    9650018.0



    
![png](output_23_1.png)
    



```python
%%time
count('./Data/',barcodes, True)
```

    ILL-RAD09_S9_L001_R1_001.fastq
    ILL-RAD06_S6_L001_R1_001.fastq
    ILL-RAD17_S17_L001_R1_001.fastq
    ILL-RAD11_S11_L001_R1_001.fastq
    ILL-RAD10_S10_L001_R1_001.fastq
    ILL-RAD02_S2_L001_R1_001.fastq
    ILL-RAD13_S13_L001_R1_001.fastq
    ILL-RAD03_S3_L001_R1_001.fastq
    ILL-RAD15_S15_L001_R1_001.fastq
    ILL-RAD08_S8_L001_R1_001.fastq
    ILL-RAD04_S4_L001_R1_001.fastq
    ILL-RAD16_S16_L001_R1_001.fastq
    ILL-RAD05_S5_L001_R1_001.fastq
    ILL-RAD07_S7_L001_R1_001.fastq
    ILL-RAD01_S1_L001_R1_001.fastq
    ILL-RAD14_S14_L001_R1_001.fastq
    ILL-RAD12_S12_L001_R1_001.fastq
    CPU times: user 7min 23s, sys: 46.6 s, total: 8min 10s
    Wall time: 7min 30s



```python
describeCounts("Output/bcg1_counts.csv")
```

           bcg1_count
    count       200.0
    mean    1621326.0
    std     1401299.0
    min       11924.0
    25%      716104.0
    50%     1284612.0
    75%     1914270.0
    max     9452617.0



    
![png](output_25_1.png)
    


### Splitting by barcode / Removing PCR duplicates

**note.** We changed the variable *minBCcount* from 100000 to 100 in the *trim2bRAD_2barcodes_dedup_N.pl* script found in Mikhael Matz's 2b_RAD Denovo scripts cloned from github. We did this in order to preserve our samples which contained less than 100000 reads.


```python
%%time
%%capture dedupLog
%%bash
cd Data
2bRAD_trim_launch_dedup_N.pl .fastq > trims
bash trims
```

    CPU times: user 40.2 ms, sys: 25.7 ms, total: 65.9 ms
    Wall time: 30min 22s



```python
# Saving output from previous cell to file
f = open("./Output/dedupLog.txt", "w")
for line in dedupLog.stdout.split('\n'):
    f.write(line+'\n')
f.close()
```

Since we changed *minBCcount*, we created demultiplexed files for barcodes not associated with one of our samples (junk files). The next three cells are removing those files.


```python
def removeJunkFiles(directory, Barcodes):
    files = os.listdir(directory)
    for fileName in files:
        try:
            ext = fileName.split("_")[5]
        except:
            pass
        else:
            bcode = ext[0:4]
            if bcode not in Barcodes:
                os.remove(directory + fileName)
```


```python
removeJunkFiles('Data/', barcodes)
```

Removing files that the previous step missed...


```bash
%%bash
cd Data
rm ILL-RAD17_S17_L001_R1_001_CTAC.tr0
rm ILL-RAD17_S17_L001_R1_001_TGTC.tr0
rm ILL-RAD17_S17_L001_R1_001_TCAC.tr0
rm ILL-RAD17_S17_L001_R1_001_GACT.tr0
```

Check that the number of .tr0 files matches number of samples.


```bash
%%bash
cd Data
ls *.tr0 | wc -l
```

    200


### Counting reads after removal of PCR duplicates


```python
%%time
%%bash
>./Output/post_dedup_counts.csv
echo "ID,post_dedup_count" >> ./Output/post_dedup_counts.csv;
cd Data
for FILE in *.tr0;
do
    COUNT=`wc -l $FILE`;
    IFS=' '
    read -a strarr <<< "$COUNT"
    COUNT=${strarr[0]}
    let COUNT=$COUNT/4
    IFS='_' read -ra strarr2 <<< $FILE
    ID1=${strarr2[1]}
    ID2=${strarr2[-1]}
    ID2="${ID2:0:4}" 
    OUT=$ID1"_"$ID2","$COUNT
    echo $OUT >> ../Output/post_dedup_counts.csv;
done
    
```

    CPU times: user 4.79 ms, sys: 1.26 ms, total: 6.05 ms
    Wall time: 4.13 s



```python
describeCounts('Output/post_dedup_counts.csv')
```

           post_dedup_count
    count             200.0
    mean           296931.0
    std            166819.0
    min              5285.0
    25%            190888.0
    50%            262810.0
    75%            382353.0
    max           1358318.0



    
![png](output_39_1.png)
    


### Quality filtering

Trimming low-quality base pairs at the end of reads with Cutadapt.


```python
%%time
%%bash
cd Data
> trims
for file in *.tr0; do
echo "cutadapt -q 15,15 -m 25 -o ${file/.tr0/}.trim $file > ${file}.trimlog.txt" >> qFilt; done
bash qFilt
```

    CPU times: user 5.07 ms, sys: 5.53 ms, total: 10.6 ms
    Wall time: 3min 16s


Check that we have expected number of .trim files.


```bash
%%bash
cd Data
ls *.trim | wc -l
```

    200


### Counting reads after quality filtering


```python
%%time
%%bash
>./Output/post_filt_counts.csv
echo "ID,post_filt_count" >> ./Output/post_filt_counts.csv;
cd Data
for FILE in *.trim;
do
    COUNT=`wc -l $FILE`;
    IFS=' '
    read -a strarr <<< "$COUNT"
    COUNT=${strarr[0]}
    let COUNT=$COUNT/4
    IFS='_' read -ra strarr2 <<< $FILE
    ID1=${strarr2[1]}
    ID2=${strarr2[-1]}
    ID2="${ID2:0:4}" 
    OUT=$ID1"_"$ID2","$COUNT
    echo $OUT >> ../Output/post_filt_counts.csv;
done
```

    CPU times: user 3.31 ms, sys: 2.45 ms, total: 5.75 ms
    Wall time: 4.13 s



```python
describeCounts('Output/post_filt_counts.csv')
```

           post_filt_count
    count            200.0
    mean          296931.0
    std           166819.0
    min             5285.0
    25%           190888.0
    50%           262810.0
    75%           382353.0
    max          1358318.0



    
![png](output_47_1.png)
    



```bash
%%bash
head Data/ILL-RAD17_S17_L001_R1_001_GCTT.trim
```

    @J00102:28:HTWWLBBXX:1:1101:10216:1086 1:N:0:NTAGAG bcd=GCTT
    AGCGCGTAAACCGCATCGGTATCGTATAGCGTGGCG
    +
    <FJJAJJJJJJJJJJJJJJAJJJJJJJAFFJJJF<J
    @J00102:28:HTWWLBBXX:1:1101:19532:1086 1:N:0:NTAGAG bcd=GCTT
    CAGATCCAGACTGCAAACCTATCGTACTCGTTGGCA
    +
    FFAFJJJJJJJJFJJJJJFFJJJJJJJJJJJJJ<JF
    @J00102:28:HTWWLBBXX:1:1101:26991:1103 1:N:0:NTAGAG bcd=GCTT
    GACAGCTGGCTGGCAGGCGCCTCGGCCTCCGGTTCG


### Building reference genome

*Odocoileus virginianus* reference genome from the University of Illinois: https://www.ncbi.nlm.nih.gov/data-hub/genome/GCA_023699985.1/

Copy all of the .fna sequence files (one for each chromosome) to a new directory 'genome'


```bash
%%bash
mkdir genome 
cd ncbi_dataset/data/GCA_023699985.1
cp *.fna ../../../genome
```

    mkdir: cannot create directory ‘genome’: File exists


    2bRAD_denovo
    Data
    Deer.ipynb
    genome
    ncbi_dataset
    Output


Remove sex chromosomes and unplaced scaffolds from the genome directory


```bash
%%bash
rm genome/chrX.fna
rm genome/chrY.fna
rm genome/unplaced.scaf.fna
```

Concatenate chromosomes together into a single genome file. Remove chromosome files.


```bash
%%bash
shopt -s extglob
cd genome
cat *.fna > genome.fasta
rm !("genome.fasta")
```

Use bowtie2, samtools, and picard to create genome index files. See more on these packages: http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml, http://www.htslib.org/doc/#manual-pages, https://broadinstitute.github.io/picard/.


```python
%%time
%%capture genomeBuild
%%bash
cd genome
bowtie2-build genome.fasta genome.fasta
```

    CPU times: user 124 ms, sys: 53.4 ms, total: 178 ms
    Wall time: 1h 33min 49s



```python
%%time
%%bash
cd genome
samtools faidx genome.fasta
```

    CPU times: user 6.26 ms, sys: 430 µs, total: 6.69 ms
    Wall time: 9.25 s



```python
%%time
%%capture picardOut
%%bash
# NEED JAVA INSTALLED
#java -jar picard.jar CreateSequenceDictionary R=genome/genome.fasta O=genome/genome.dict
```

    CPU times: user 8.03 ms, sys: 0 ns, total: 8.03 ms
    Wall time: 7.34 ms


### Aligning reads to the reference genome

Aligning reads to the reference genome. %%capture keeps record of the standard output which contains useful data on allignment rates.


```python
%%time
%%capture allignmentRates
%%bash
cd Data
2bRAD_bowtie2_launch.pl '\.trim$' ../genome/genome.fasta > maps
bash maps
```

    CPU times: user 157 ms, sys: 76.6 ms, total: 233 ms
    Wall time: 1h 21min 3s


Saving standard output from read alligning to a file.


```python
f = open("./Output/allignmentRates.txt", "w")
for line in allignmentRates.stderr.split('\n'):
    f.write(line+'\n')
f.close()
```

This function takes the data from the file just created, and stores it into a csv file.


```python
# Allignment rate function
def getAllignmentRates(allignFile, output, fileOrder):
    allignFile = open(allignFile, "r")
    orderFile = open(fileOrder ,'r')
    
    order = []
    allignment = []
    readCount = []
    for line in orderFile:
        line = line.split(" ")[-3]
        line = line.split("_")
        ida = line[1]
        idb = line[5].split(".")[0]
        id = ida + "_" + idb
        order.append(id)
        
    count = 0
    for line in allignFile:
        if count%6 ==0:
            reads = line.strip().split(" ")[0]
            readCount.append(reads)
        
        if count % 6 == 5:
            line = line.strip()
            percent = line.split(" ")[0][0:-1]
            allignment.append(percent) 
        count += 1    
    df_allign = []
    
    for i in range(len(order)):
        df_allign.append({'ID':order[i], 'reads':readCount[i],'allignment_rate':allignment[i]})
        
    df_allign = pd.DataFrame.from_dict(df_allign)
    df_allign.to_csv('./Output/allignmentRates'+output+'.csv', index=False)

```


```python
getAllignmentRates('./Output/allignmentRates.txt','','./Data/maps')
```


```python
#Importing files with pandas
countRaw = pd.read_csv("./Output/raw_counts.csv")
countBcg1 = pd.read_csv("./Output/bcg1_counts.csv")
countPostDedup = pd.read_csv("./Output/post_dedup_counts.csv")
countPostFiltering = pd.read_csv("./Output/post_filt_counts.csv")
allignmentRate = pd.read_csv("./Output/allignmentRates.csv")
sampleWell = pd.read_csv("./SampleData/well_plate.csv")
sampleWell = sampleWell.rename(columns={"Sample ID":"ID"})

# Merging dataframes
countsMerged = countRaw.merge(countBcg1, how='right')
countsMerged = countsMerged.merge(countPostDedup, how='right')
countsMerged = countsMerged.merge(countPostFiltering, how='right')
countsMerged = countsMerged.merge(allignmentRate, how = 'right')
countsMerged = countsMerged.merge(sampleWell, how = 'right')

# Calculating zscores, changing col names, and merging to main dataframe
numeric_cols = countsMerged.select_dtypes(include=[np.number]).columns
zscores = countsMerged[numeric_cols].apply(scipy.stats.zscore)
zscores = zscores.rename(columns=lambda x: x+"_Z")
countsMerged = pd.concat([countsMerged, zscores], axis = 1)

countsMerged.to_csv('./Output/readCountSummary.csv', index=False)
display(countsMerged)

print('Weighted average allignment rate: ',round(sum(countsMerged['allignment_rate'] * countsMerged['reads'])/sum(countsMerged['reads']),2),"%",sep='')
print('Non-weighted average allignment rate: ', round(sum(countsMerged['allignment_rate'])/len(countsMerged['allignment_rate']),2),"%", sep='')
```


<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>ID</th>
      <th>raw_count</th>
      <th>bcg1_count</th>
      <th>post_dedup_count</th>
      <th>post_filt_count</th>
      <th>reads</th>
      <th>allignment_rate</th>
      <th>Sample Number</th>
      <th>raw_count_Z</th>
      <th>bcg1_count_Z</th>
      <th>post_dedup_count_Z</th>
      <th>post_filt_count_Z</th>
      <th>reads_Z</th>
      <th>allignment_rate_Z</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>S7_TGTC</td>
      <td>4324452</td>
      <td>4235342</td>
      <td>517398</td>
      <td>517398</td>
      <td>517398</td>
      <td>28.83</td>
      <td>1</td>
      <td>1.865960</td>
      <td>1.870104</td>
      <td>1.324911</td>
      <td>1.324911</td>
      <td>1.324911</td>
      <td>-0.848647</td>
    </tr>
    <tr>
      <th>1</th>
      <td>S10_AGTG</td>
      <td>2511035</td>
      <td>2444053</td>
      <td>496764</td>
      <td>496764</td>
      <td>496764</td>
      <td>4.17</td>
      <td>2</td>
      <td>0.596262</td>
      <td>0.588590</td>
      <td>1.200909</td>
      <td>1.200909</td>
      <td>1.200909</td>
      <td>-1.712575</td>
    </tr>
    <tr>
      <th>2</th>
      <td>S13_CTAC</td>
      <td>2195954</td>
      <td>2166437</td>
      <td>411489</td>
      <td>411489</td>
      <td>411489</td>
      <td>0.45</td>
      <td>3</td>
      <td>0.375652</td>
      <td>0.389980</td>
      <td>0.688443</td>
      <td>0.688443</td>
      <td>0.688443</td>
      <td>-1.842900</td>
    </tr>
    <tr>
      <th>3</th>
      <td>S16_GTGT</td>
      <td>2137286</td>
      <td>2100111</td>
      <td>404740</td>
      <td>404740</td>
      <td>404740</td>
      <td>11.74</td>
      <td>4</td>
      <td>0.334574</td>
      <td>0.342530</td>
      <td>0.647884</td>
      <td>0.647884</td>
      <td>0.647884</td>
      <td>-1.447371</td>
    </tr>
    <tr>
      <th>4</th>
      <td>S14_TGTC</td>
      <td>247758</td>
      <td>242169</td>
      <td>63079</td>
      <td>63079</td>
      <td>63079</td>
      <td>76.53</td>
      <td>5</td>
      <td>-0.988415</td>
      <td>-0.986668</td>
      <td>-1.405353</td>
      <td>-1.405353</td>
      <td>-1.405353</td>
      <td>0.822456</td>
    </tr>
    <tr>
      <th>...</th>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
    </tr>
    <tr>
      <th>195</th>
      <td>S6_GACT</td>
      <td>2782969</td>
      <td>2706626</td>
      <td>552994</td>
      <td>552994</td>
      <td>552994</td>
      <td>57.12</td>
      <td>68B</td>
      <td>0.786661</td>
      <td>0.776439</td>
      <td>1.538827</td>
      <td>1.538828</td>
      <td>1.538828</td>
      <td>0.142454</td>
    </tr>
    <tr>
      <th>196</th>
      <td>S4_TGTC</td>
      <td>2108754</td>
      <td>2070261</td>
      <td>330555</td>
      <td>330555</td>
      <td>330555</td>
      <td>81.67</td>
      <td>8A</td>
      <td>0.314597</td>
      <td>0.321174</td>
      <td>0.202064</td>
      <td>0.202064</td>
      <td>0.202064</td>
      <td>1.002528</td>
    </tr>
    <tr>
      <th>197</th>
      <td>S4_TCAC</td>
      <td>1601901</td>
      <td>1563440</td>
      <td>355262</td>
      <td>355262</td>
      <td>355262</td>
      <td>82.77</td>
      <td>8B</td>
      <td>-0.040286</td>
      <td>-0.041412</td>
      <td>0.350542</td>
      <td>0.350542</td>
      <td>0.350542</td>
      <td>1.041065</td>
    </tr>
    <tr>
      <th>198</th>
      <td>S9_GTGT</td>
      <td>283729</td>
      <td>276837</td>
      <td>91338</td>
      <td>91338</td>
      <td>91338</td>
      <td>80.27</td>
      <td>9A</td>
      <td>-0.963230</td>
      <td>-0.961866</td>
      <td>-1.235529</td>
      <td>-1.235529</td>
      <td>-1.235529</td>
      <td>0.953481</td>
    </tr>
    <tr>
      <th>199</th>
      <td>S17_GTGT</td>
      <td>12717</td>
      <td>11924</td>
      <td>5285</td>
      <td>5285</td>
      <td>5285</td>
      <td>49.78</td>
      <td>9B</td>
      <td>-1.152984</td>
      <td>-1.151389</td>
      <td>-1.752671</td>
      <td>-1.752671</td>
      <td>-1.752671</td>
      <td>-0.114693</td>
    </tr>
  </tbody>
</table>
<p>200 rows × 14 columns</p>
</div>


    Weighted average allignment rate: 45.81%
    Non-weighted average allignment rate: 53.05%



```r
%%R
readSummary <- read.csv("./Output/readCountSummary.csv")
par (mfrow = c(2,2))
hist(readSummary$raw_count, main = 'raw_count')
hist(readSummary$bcg1_count, main = 'bcg1_count')
hist(readSummary$post_dedup_count, main = 'post_dedup_count')
hist(readSummary$allignment_rate, main = 'allignment_rate')
```


    
![png](output_70_0.png)
    


Convert sam files to bam files as input for ANGSD using Samtools


```python
%%time
%%bash
cd Data
>s2b
for file in *.sam; do
echo "samtools sort -O bam -o ${file/.sam/}.bam $file && samtools index ${file/.sam/}.bam">>s2b;
done
bash s2b
```

    CPU times: user 6.31 ms, sys: 3.61 ms, total: 9.92 ms
    Wall time: 1min 47s


### Genotype calling with ANGSD


```python
%%time
%%capture angsd_geno_likely
%%bash
# Fuzzy genotyping
cd Data
ls *.bam > bams
angsd -b bams -gl 2 -domajorminor 1 -snp_pval 1e-6 -domaf 1 -minmaf 0.05 -doGlf 3 -uniqueOnly 1 -remove_bads 1 -minMapQ 20 -minQ 25 -minInd 150 -out geno_likely
mv geno_likely* ../Output
```

    CPU times: user 28 ms, sys: 7.77 ms, total: 35.7 ms
    Wall time: 2min 14s



```python
%%time
%%capture angsd_geno_call
%%bash
# Calling genotypes
cd Data
angsd -b bams -gl 2 -domajorminor 1 -snp_pval 1e-6 -domaf 1 -minmaf 0.05 -doGlf 3 -uniqueOnly 1 -remove_bads 1 -minMapQ 20 -minQ 25 -minInd 150 -doCounts 1 -doDepth 1 -dumpCounts 2 -doPost 1 -doGeno 4 -postCutoff 0.6 -out geno_call
mv geno_call* ../Output
```

    CPU times: user 32.6 ms, sys: 293 µs, total: 32.9 ms
    Wall time: 2min 26s



```python
def getANGSDresults(fileName, excludeZero):
    df = pd.read_table(fileName, header=None)
    df.drop(101,axis=1,inplace=True)
    if excludeZero:
        df.drop(0, axis=1, inplace=True)
    df['count_depth'] = df[list(df.columns)].sum(axis=1)
    df['sumProduct'] = 0
    for i in list(df.columns):
        if type(i) == int:
            df['sumProduct'] += df[i].multiply(i)
    df['average_depth'] = df['sumProduct'] / df['count_depth']
    
    def median(row):
        row = row.iloc[0:100]
        row = list(row)
        n = sum(row)
        position = 0
        i = 0
        while position < n//2:
            position += row[i]
            i += 1
        return i 
   
    df['median_depth'] = df.apply(lambda row: median(row), axis=1)

    
    print(df['average_depth'].describe().round(2))
    weightedAverage = sum(df['average_depth'] * df['count_depth']) / sum(df['count_depth'])
    print('Weighted average: ', weightedAverage)
    sns.histplot(data=df, x = 'average_depth', fill=False)
    plt.axvline(df['average_depth'].median(), color='green', linestyle='dashed')
    plt.axvline(df['average_depth'].mean(), color='red', linestyle='dashed')
    display(df)
    df[['count_depth','average_depth','median_depth']].to_csv('./Output/coverageSummary.csv', index=True)
```


```python
getANGSDresults("./Output/geno_call.depthSample", True)
```

    count    200.00
    mean       8.29
    std        4.62
    min        1.84
    25%        4.61
    50%        7.88
    75%       10.73
    max       30.62
    Name: average_depth, dtype: float64
    Weighted average:  9.10367442449871



<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>1</th>
      <th>2</th>
      <th>3</th>
      <th>4</th>
      <th>5</th>
      <th>6</th>
      <th>7</th>
      <th>8</th>
      <th>9</th>
      <th>10</th>
      <th>...</th>
      <th>95</th>
      <th>96</th>
      <th>97</th>
      <th>98</th>
      <th>99</th>
      <th>100</th>
      <th>count_depth</th>
      <th>sumProduct</th>
      <th>average_depth</th>
      <th>median_depth</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>42</td>
      <td>64</td>
      <td>96</td>
      <td>141</td>
      <td>139</td>
      <td>162</td>
      <td>166</td>
      <td>143</td>
      <td>158</td>
      <td>164</td>
      <td>...</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>40</td>
      <td>3376</td>
      <td>51975</td>
      <td>15.395438</td>
      <td>13</td>
    </tr>
    <tr>
      <th>1</th>
      <td>93</td>
      <td>145</td>
      <td>200</td>
      <td>228</td>
      <td>248</td>
      <td>238</td>
      <td>205</td>
      <td>200</td>
      <td>216</td>
      <td>186</td>
      <td>...</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>28</td>
      <td>3339</td>
      <td>38315</td>
      <td>11.474993</td>
      <td>9</td>
    </tr>
    <tr>
      <th>2</th>
      <td>425</td>
      <td>425</td>
      <td>442</td>
      <td>353</td>
      <td>317</td>
      <td>207</td>
      <td>173</td>
      <td>140</td>
      <td>114</td>
      <td>98</td>
      <td>...</td>
      <td>1</td>
      <td>1</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>13</td>
      <td>3147</td>
      <td>19840</td>
      <td>6.304417</td>
      <td>4</td>
    </tr>
    <tr>
      <th>3</th>
      <td>81</td>
      <td>131</td>
      <td>127</td>
      <td>169</td>
      <td>160</td>
      <td>172</td>
      <td>175</td>
      <td>191</td>
      <td>159</td>
      <td>148</td>
      <td>...</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>26</td>
      <td>3377</td>
      <td>49906</td>
      <td>14.778206</td>
      <td>12</td>
    </tr>
    <tr>
      <th>4</th>
      <td>27</td>
      <td>28</td>
      <td>58</td>
      <td>110</td>
      <td>125</td>
      <td>140</td>
      <td>174</td>
      <td>146</td>
      <td>155</td>
      <td>152</td>
      <td>...</td>
      <td>1</td>
      <td>2</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>26</td>
      <td>3381</td>
      <td>57567</td>
      <td>17.026619</td>
      <td>14</td>
    </tr>
    <tr>
      <th>...</th>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
    </tr>
    <tr>
      <th>195</th>
      <td>506</td>
      <td>469</td>
      <td>419</td>
      <td>296</td>
      <td>236</td>
      <td>175</td>
      <td>166</td>
      <td>110</td>
      <td>98</td>
      <td>69</td>
      <td>...</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>15</td>
      <td>2962</td>
      <td>18729</td>
      <td>6.323093</td>
      <td>4</td>
    </tr>
    <tr>
      <th>196</th>
      <td>78</td>
      <td>3</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>5</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>...</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>1</td>
      <td>87</td>
      <td>214</td>
      <td>2.459770</td>
      <td>1</td>
    </tr>
    <tr>
      <th>197</th>
      <td>768</td>
      <td>595</td>
      <td>380</td>
      <td>288</td>
      <td>167</td>
      <td>124</td>
      <td>59</td>
      <td>74</td>
      <td>39</td>
      <td>38</td>
      <td>...</td>
      <td>0</td>
      <td>1</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>11</td>
      <td>2714</td>
      <td>13267</td>
      <td>4.888357</td>
      <td>2</td>
    </tr>
    <tr>
      <th>198</th>
      <td>126</td>
      <td>23</td>
      <td>9</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>1</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>...</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>1</td>
      <td>165</td>
      <td>694</td>
      <td>4.206061</td>
      <td>1</td>
    </tr>
    <tr>
      <th>199</th>
      <td>171</td>
      <td>294</td>
      <td>363</td>
      <td>318</td>
      <td>347</td>
      <td>307</td>
      <td>243</td>
      <td>214</td>
      <td>179</td>
      <td>134</td>
      <td>...</td>
      <td>1</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>11</td>
      <td>3309</td>
      <td>27012</td>
      <td>8.163191</td>
      <td>6</td>
    </tr>
  </tbody>
</table>
<p>200 rows × 104 columns</p>
</div>



    
![png](output_77_2.png)
    


Mean depth = 8.29 (weighted = 9.1). 

### NGSRelate


```bash
%%bash
zcat Output/geno_likely.mafs.gz | cut -f5 | sed 1d>Output/allele_freq
```


```python
%%time
%%capture ngsRelate
%%bash
cd Data
ngsRelate -f ../Output/allele_freq -g ../Output/geno_likely.glf.gz -n 200 -z bams -O ../Output/ngsOutput.res
```

    CPU times: user 9.39 ms, sys: 1.97 ms, total: 11.4 ms
    Wall time: 39.2 s


The following function formats the NgsRelate output. It creates identifier columns (sample ID, sample number, family, mother?, child?, and replicate?) to help with easier subsetting and analysis down the road.


```r
%%R
formatNGS <- function(outputFileName, wellPlateFileName){
    ngs_out <- read.table(outputFileName, header=TRUE)
    
    well <- read.csv(wellPlateFileName, header=TRUE)
    colnames(ngs_out)[3] <- "file_a"
    colnames(ngs_out)[4] <- "file_b"
    
    
    
    # Keep only columns we care about
    ngs_out <- subset(ngs_out, select = c(1:5, 15:20, 25, 31:33))
   
    # Create sampleID columns from file name for easier identification 
    sample.ID.A <- unlist(lapply(strsplit(ngs_out$file_a, "_" ), function(x) {paste(x[[2]], substr(x[[6]],1,4), sep ="_")}))
    ngs_out$sample.ID.A <- sample.ID.A
    
    sample.ID.B <- unlist(lapply(strsplit(ngs_out$file_b, "_" ), function(x) {paste(x[[2]], substr(x[[6]],1,4), sep ="_")}))
    ngs_out$sample.ID.B <- sample.ID.B
    
    # Joining sample ID with sample number
    ngs_out <- merge(ngs_out, well, by.x ="sample.ID.A", by.y="Sample.ID")
    colnames(ngs_out)[length(ngs_out)] <- "sample.num.A"
    
    ngs_out <- merge(ngs_out, well, by.x ="sample.ID.B", by.y="Sample.ID")
    colnames(ngs_out)[length(ngs_out)] <- "sample.num.B"
    
    
    
    # Reinstate ordering
    ngs_out <- ngs_out %>% select(sample.ID.A, everything())
    ngs_out <- ngs_out[order(ngs_out$a, ngs_out$b),]
    
    # Adding family columns
    findFamily <- function(sampleID){
        family = ''
        for (char in strsplit(sampleID,"")[[1]]){
            if (!is.na(as.numeric(char))){
                family = paste(family, char, sep = '')}
            else{break}}
        return (family)}
        
    ngs_out$family.A <- unlist(lapply(ngs_out$sample.num.A, function(x) {findFamily(x)}))    
    ngs_out$family.B <- unlist(lapply(ngs_out$sample.num.B, function(x) {findFamily(x)}))
    
    # Adding mother columns
    findMother <- function(x){
        for (char in strsplit(x, "")[[1]]){
            if (is.na(as.numeric(char)) & !grepl(char, ' replic.')) {
                return (FALSE)}}
        return (TRUE) }
    
    ngs_out$mother.A <- unlist(lapply(ngs_out$sample.num.A, function(x) {findMother(x)}))
    ngs_out$mother.B <- unlist(lapply(ngs_out$sample.num.B, function(x) {findMother(x)}))
    
    # Adding child columns
    ngs_out$child.A <- unlist(lapply(ngs_out$mother.A, function(x) {!x}))
    ngs_out$child.B <- unlist(lapply(ngs_out$mother.B, function(x) {!x}))
    
    # Adding replicate columns
    ngs_out$replic.A <- unlist(lapply(ngs_out$sample.num.A, function(x){grepl(' replic.', x)}))
    ngs_out$replic.B <- unlist(lapply(ngs_out$sample.num.B, function(x){grepl(' replic.', x)}))
    
    
    # Adding coverage depth median/average statistics
    cov <- read.csv('./Output/coverageSummary.csv',header=TRUE)
    ngs_out <- merge(ngs_out, cov, by.x='a',by.y=c(1))
    ngs_out <- merge(ngs_out, cov, by.x='b',by.y=c(1))
    
    return (ngs_out)

  
}
```


```r
%%R
ngs <- formatNGS('Output/ngsOutput.res','SampleData/well_plate.csv')
str(ngs)
```

    'data.frame':	19900 obs. of  33 variables:
     $ b                     : int  1 2 2 3 3 3 4 4 4 4 ...
     $ a                     : int  0 0 1 0 2 1 1 0 2 3 ...
     $ sample.ID.A           : chr  "S1_ACCA" "S1_ACCA" "S1_AGAC" "S1_ACCA" ...
     $ sample.ID.B           : chr  "S1_AGAC" "S1_AGTG" "S1_AGTG" "S1_CATC" ...
     $ file_a                : chr  "ILL-RAD01_S1_L001_R1_001_ACCA.trim.bt2.bam" "ILL-RAD01_S1_L001_R1_001_ACCA.trim.bt2.bam" "ILL-RAD01_S1_L001_R1_001_AGAC.trim.bt2.bam" "ILL-RAD01_S1_L001_R1_001_ACCA.trim.bt2.bam" ...
     $ file_b                : chr  "ILL-RAD01_S1_L001_R1_001_AGAC.trim.bt2.bam" "ILL-RAD01_S1_L001_R1_001_AGTG.trim.bt2.bam" "ILL-RAD01_S1_L001_R1_001_AGTG.trim.bt2.bam" "ILL-RAD01_S1_L001_R1_001_CATC.trim.bt2.bam" ...
     $ nSites                : int  3643 3408 3373 3680 3413 3638 3645 3684 3413 3686 ...
     $ rab                   : num  0 0.00232 0.008608 0.000004 0.298146 ...
     $ Fa                    : num  0.0201 0.0169 0.0141 0.0224 0 ...
     $ Fb                    : num  0.021755 0.000005 0.000001 0 0 ...
     $ theta                 : num  0 0.001162 0.005739 0.000002 0.149073 ...
     $ inbred_relatedness_1_2: num  0 0.000005 0.005739 0 0 ...
     $ inbred_relatedness_2_1: num  0e+00 5e-06 1e-06 0e+00 0e+00 ...
     $ F_diff_a_b            : num  -0.000814 0.008451 0.001301 0.011188 0 ...
     $ R0                    : num  0.5598 0.4068 0.3611 0.3213 0.0811 ...
     $ R1                    : num  0.235 0.269 0.254 0.265 0.574 ...
     $ KING                  : num  -0.021 0.0351 0.0499 0.0655 0.2288 ...
     $ sample.num.A          : chr  "43A replic." "43A replic." "65B" "43A replic." ...
     $ sample.num.B          : chr  "65B" "45A" "45A" "45B" ...
     $ family.A              : chr  "43" "43" "65" "43" ...
     $ family.B              : chr  "65" "45" "45" "45" ...
     $ mother.A              : logi  FALSE FALSE FALSE FALSE FALSE FALSE ...
     $ mother.B              : logi  FALSE FALSE FALSE FALSE FALSE FALSE ...
     $ child.A               : logi  TRUE TRUE TRUE TRUE TRUE TRUE ...
     $ child.B               : logi  TRUE TRUE TRUE TRUE TRUE TRUE ...
     $ replic.A              : logi  TRUE TRUE FALSE TRUE FALSE FALSE ...
     $ replic.B              : logi  FALSE FALSE FALSE FALSE FALSE FALSE ...
     $ count_depth.x         : int  3376 3376 3339 3376 3147 3339 3339 3376 3147 3377 ...
     $ average_depth.x       : num  15.4 15.4 11.5 15.4 6.3 ...
     $ median_depth.x        : int  13 13 9 13 4 9 9 13 4 12 ...
     $ count_depth.y         : int  3339 3147 3147 3377 3377 3377 3381 3381 3381 3381 ...
     $ average_depth.y       : num  11.5 6.3 6.3 14.8 14.8 ...
     $ median_depth.y        : int  9 4 4 12 12 12 14 14 14 14 ...



```r
%%R
ggplot(ngs, aes(x=nSites)) +
    geom_histogram(bins = 60, fill = 'cornflowerblue') +
    theme_classic() +
    ggtitle("SNP count distribution") 
```


    
![png](output_85_0.png)
    


Skewed left with a spike at around 200 SNPs. 


```r
%%R
unrelated <- ngs %>% filter(family.A != family.B)
outliers <- unrelated %>% filter(rab > 0.16)
ggplot(unrelated, aes(y=rab)) + 
    geom_boxplot() +
    geom_text_repel(data=outliers, aes(x=0,label=nSites)) +
    theme_bw() +
    ggtitle("Relatedness values for unrelated pairs")
```


    
![png](output_87_0.png)
    


Median is around 0, as we would expect, but many outliers with high relatedness (all with relatively low quantity of site comparisons).


```r
%%R
ggplot(ngs %>% filter(family.A == family.B & mother.A + mother.B == 1 & replic.A + replic.B == 0), 
       aes(x=rab)) +
    geom_histogram(bins=50) +
    theme_bw() +
    ggtitle("Relatedness values for mother offspring pairs")
```


    
![png](output_89_0.png)
    



```r
%%R
outliers <- ngs %>% filter(family.A == family.B & mother.A + mother.B == 1 & replic.A + replic.B == 0) %>% filter(rab<0.38 | rab > 0.6)
ggplot(ngs %>% filter(family.A == family.B & mother.A + mother.B == 1 & replic.A + replic.B == 0), 
       aes(y=rab)) + 
    geom_boxplot() +
    geom_text_repel(data=outliers, aes(x=0,label=nSites)) +
    theme_bw() +
    ggtitle("Relatedness values for mother offspring pairs")
```


    
![png](output_90_0.png)
    


Median ~0.5 as we expect, but large number of pairs with low relatedness, all with relatively low quantity of sites (<500)


```r
%%R
sibPairs <- ngs %>% filter(family.A==family.B & child.A + child.B ==2 & replic.A+replic.B ==0)
ggplot(sibPairs, aes(x=rab)) + 
    geom_histogram(bins=40) + 
    theme_classic() + 
    geom_vline(xintercept = 0.16, linetype = 'dashed', color='red') +
    geom_vline(xintercept = 0.34, linetype = 'dashed', color='red') +
    geom_vline(xintercept = 0.41, linetype = 'dashed', color='blue') +
    geom_vline(xintercept = 0.59, linetype = 'dashed', color='blue') +
    ggtitle("Relatedness distribution for sibling pairs") + 
    labs(caption = "Red region marks expected relatedness values for half-siblings, blue for full-siblings")
    
```


    
![png](output_92_0.png)
    



```r
%%R
outliers <- sibPairs %>% filter(rab<0.16 | rab > 0.59)
ggplot(sibPairs, aes(y=rab)) + 
    geom_boxplot() +
    geom_text_repel(data=outliers, aes(x=0,label=nSites)) +
    theme_classic() +
    ggtitle("Relatedness values for sibling pairs")
```


    
![png](output_93_0.png)
    


### Outlier investigation: 17A / 17B

By filtering to nSites equals 143, we find that the outlier with low relatedness are supposed siblings 17A and 17B.


```r
%%R
print(sibPairs %>% filter (nSites == 143))
```

       b  a sample.ID.A sample.ID.B                                     file_a
    1 98 43     S4_GTGA     S9_AGTG ILL-RAD04_S4_L001_R1_001_GTGA.trim.bt2.bam
                                          file_b nSites rab Fa      Fb theta
    1 ILL-RAD09_S9_L001_R1_001_AGTG.trim.bt2.bam    143   0  0 6.3e-05     0
      inbred_relatedness_1_2 inbred_relatedness_2_1 F_diff_a_b       R0       R1
    1                      0                      0   -3.1e-05 0.292035 0.345635
          KING sample.num.A sample.num.B family.A family.B mother.A mother.B
    1 0.090396          17B          17A       17       17    FALSE    FALSE
      child.A child.B replic.A replic.B count_depth.x average_depth.x
    1    TRUE    TRUE    FALSE    FALSE          3359        11.73415
      median_depth.x count_depth.y average_depth.y median_depth.y
    1              9           142        2.387324              1



```r
%%R
pairs.17A <- ngs %>% filter(sample.num.A == '17A' | sample.num.B == '17A')
ggplot(pairs.17A, aes(y=rab)) + 
    geom_boxplot() + 
    geom_text_repel(data=pairs.17A %>% filter(rab > 0.25), aes(x=0,label=nSites)) + theme_bw()
```


    
![png](output_97_0.png)
    


Many inviduals with high relatedness to 17A, all have relatively low SNP count (16-132). 17A likely has poor sequence quality.


```r
%%R
sampleQuality <- read.csv("./Output/readCountSummary.csv")
sample17A.quality <- sampleQuality %>% filter(Sample.Number == '17A')
sample17A.quality
```

           ID raw_count bcg1_count post_dedup_count post_filt_count  reads
    1 S9_AGTG   1916960    1879736           200420          200420 200420
      allignment_rate Sample.Number raw_count_Z bcg1_count_Z post_dedup_count_Z
    1            0.88           17A   0.1803085    0.1848701         -0.5799924
      post_filt_count_Z    reads_Z allignment_rate_Z
    1        -0.5799923 -0.5799923         -1.827836


17A has below average reads/allignment rates (Z-scores: Z = -0.5799923, Z = -1.827836)

Perhaps 17A was mislabeled and belongs to another family. Checking to see if there are any other potential families that sample 17A reasonably could fall into where the relatedness between 17A and potential siblings is greater than 0.2, and the relatedness between 17A and potential mother is greater than 0.4.


```r
%%R
potentialSibs.17A <- pairs.17A %>% filter(child.A + child.B ==2 & rab > 0.2)
potentialMothers.17A <- pairs.17A %>% filter(child.A + child.B ==1 & rab > 0.4)

possibleFams.17A <- intersect(union(potentialSibs.17A$family.A,potentialSibs.17A$family.B),
                                union(potentialMothers.17A$family.A, potentialMothers.17A$family.B))

print(possibleFams.17A)

```

    [1] "17" "28"



```r
%%R
# Only potential family identified is 28.
potentialFams.17A <- pairs.17A %>% filter ((family.A %in% c("17","28") & 
                                            (family.B %in% c("17","28"))))
potFamily <- function(x,y){
    if (x=='17'){return (y)}else{return (x)}}

otherID <- function(x,y){
    if (x=='17A'){return(y)} else{return(x)}}
potentialFams.17A$pot.fam <- mapply(potFamily, potentialFams.17A$family.A, potentialFams.17A$family.B )
potentialFams.17A$otherID <- mapply(otherID, potentialFams.17A$sample.num.A, potentialFams.17A$sample.num.B)

ggplot(data = potentialFams.17A, aes(y = rab, x = pot.fam)) + 
    theme_bw() +
    geom_point() +
    geom_text_repel(aes(x=pot.fam, label=nSites), nudge_x=0.2) + 
    geom_text_repel(aes(x=pot.fam, label=otherID), nudge_y=-0.025) + 
    ylab("Relatedness") + 
    xlab("Potential Family")
```


    
![png](output_103_0.png)
    


Family 28 is the only family that meets these criteria. Lets assume 17A and 28A could be full siblings and 17A and 28B are half siblings (their rab is about 0.16). This implies that 28A and and 28B are also half siblings. Is this likely?


```r
%%R
sibs.28 <- sibPairs %>% filter(sample.num.A == '28B' & sample.num.B == '28A')
print(sibs.28$rab)
```

    [1] 0.456836


28A and 28B are likely full siblings (rab = 0.46). Our conclusion is that 17A was not misabled but rather suffers from poor read quality.


```r
%%R
# A quick look at 17B
pairs.17B <- ngs %>% filter(sample.num.A == '17B' | sample.num.B == '17B')
ggplot(pairs.17B, aes(y=rab)) + 
    geom_boxplot() + 
    geom_text(data=pairs.17B %>% filter(rab > 0.05), aes(x=0,label=nSites), nudge_x=-0.05) + 
    geom_text(data=pairs.17B %>% filter(rab > 0.05), aes(x=0,label=sample.num.B), nudge_x=0.05) +
    theme_bw()
```


    
![png](output_107_0.png)
    


17B has expected relatedness with mother, and we are very confident in this value (# SNPs = 2856). It is very likely that 17B was not mislabeled.


```r
%%R
pairs.28B <- ngs %>% filter(sample.num.A == '28B' | sample.num.B == '28B')
ggplot(pairs.28B, aes(y=rab)) + 
    geom_boxplot() + 
    geom_text(data=pairs.28B %>% filter(rab > 0.05), aes(x=0,label=nSites), nudge_x=-0.05) + 
    geom_text(data=pairs.28B %>% filter(rab > 0.05), aes(x=0,label=sample.num.B), nudge_x=0.05) +
    theme_bw()
```


    
![png](output_109_0.png)
    



```r
%%R
pairs.28A <- ngs %>% filter(sample.num.A == '28A' | sample.num.B == '28A')
ggplot(pairs.28A, aes(y=rab)) + 
    geom_boxplot() + 
    geom_text(data=pairs.28A %>% filter(rab > 0.05), aes(x=0,label=nSites), nudge_x=-0.05) + 
    geom_text(data=pairs.28A %>% filter(rab > 0.05), aes(x=0,label=sample.num.A), nudge_x=0.05) +
    theme_bw()
```


    
![png](output_110_0.png)
    


28B is highly related to its mother 28, so it too is likely correctly labeled.

We conclude that 17A suffers from poor sequencing/genotyping quality and was not mislabeled.

### Bootstrap Confidence Intervals


```r
%%R
pairs <- data.frame(sibPairs$a, sibPairs$b)
write.table(pairs, file = 'CI/pairs.csv',sep=',',col.names=FALSE,row.names=FALSE)
```


```python
def makeCIlaunch(frequencyFile, ANGSDglfFile, bamsFile):
    launch_pairs = open("./CI/ngs_CI_launch.sh","w")
    pairCombos = open("./CI/pairs.csv", "r")

    for line in pairCombos:
        line = line.strip()
        line=line.split(",")
        a = line[0]
        b = line[1]

        launch_pairs.write('ngsRelate -f '+frequencyFile+' -g '+ ANGSDglfFile+' -n 200 -z '+bamsFile+' -O CI/'+a+'_'+b+'_CI_output.res -a '+a+ ' -b '+b+ ' -B 10000')
        launch_pairs.write('\n')


    launch_pairs.close()
```


```python
makeCIlaunch("Output/allele_freq",'Output/geno_likely.glf.gz','Data/bams')
```


```python
%%time
%%capture ngsBootstrap
%%bash
bash CI/ngs_CI_launch.sh
```

    CPU times: user 127 ms, sys: 66.8 ms, total: 193 ms
    Wall time: 34min 40s



```bash
%%bash
cd CI
> CI_R_launch
for file in *.res; do
echo "Rscript get_CI.R $file ../Output/CI_output.csv" >> CI_R_launch;
done
> ../Output/CI_output.csv
echo "a,b,rab_lower,rab_upper,R0_lower,R0_upper,R1_lower,R1_upper" >> ../Output/CI_output.csv
bash CI_R_launch
```


```r
%%R
bootstrapCI <- read.csv('Output/CI_output.csv') 
sibPairs <- merge(sibPairs, bootstrapCI)
ggplot(sibPairs, aes(x=reorder(family.A, rab), y=rab)) +
    geom_point() + 
    theme_bw() + theme(axis.text.x = element_blank(), panel.grid.major = element_blank(), 
                       panel.grid.minor = element_blank() ) +
    geom_pointrange(aes(ymin = rab_lower, ymax = rab_upper),width=0.1, alpha = 0.5) +
    geom_text_repel(data=sibPairs%>%filter(family.A=='28'),aes(x=reorder(family.A, rab), label=family.A))

```


    
![png](output_119_0.png)
    



```r
%%R
str(ngs)
```

    'data.frame':	19900 obs. of  33 variables:
     $ b                     : int  1 2 2 3 3 3 4 4 4 4 ...
     $ a                     : int  0 0 1 0 2 1 1 0 2 3 ...
     $ sample.ID.A           : chr  "S1_ACCA" "S1_ACCA" "S1_AGAC" "S1_ACCA" ...
     $ sample.ID.B           : chr  "S1_AGAC" "S1_AGTG" "S1_AGTG" "S1_CATC" ...
     $ file_a                : chr  "ILL-RAD01_S1_L001_R1_001_ACCA.trim.bt2.bam" "ILL-RAD01_S1_L001_R1_001_ACCA.trim.bt2.bam" "ILL-RAD01_S1_L001_R1_001_AGAC.trim.bt2.bam" "ILL-RAD01_S1_L001_R1_001_ACCA.trim.bt2.bam" ...
     $ file_b                : chr  "ILL-RAD01_S1_L001_R1_001_AGAC.trim.bt2.bam" "ILL-RAD01_S1_L001_R1_001_AGTG.trim.bt2.bam" "ILL-RAD01_S1_L001_R1_001_AGTG.trim.bt2.bam" "ILL-RAD01_S1_L001_R1_001_CATC.trim.bt2.bam" ...
     $ nSites                : int  3643 3408 3373 3680 3413 3638 3645 3684 3413 3686 ...
     $ rab                   : num  0 0.00232 0.008608 0.000004 0.298146 ...
     $ Fa                    : num  0.0201 0.0169 0.0141 0.0224 0 ...
     $ Fb                    : num  0.021755 0.000005 0.000001 0 0 ...
     $ theta                 : num  0 0.001162 0.005739 0.000002 0.149073 ...
     $ inbred_relatedness_1_2: num  0 0.000005 0.005739 0 0 ...
     $ inbred_relatedness_2_1: num  0e+00 5e-06 1e-06 0e+00 0e+00 ...
     $ F_diff_a_b            : num  -0.000814 0.008451 0.001301 0.011188 0 ...
     $ R0                    : num  0.5598 0.4068 0.3611 0.3213 0.0811 ...
     $ R1                    : num  0.235 0.269 0.254 0.265 0.574 ...
     $ KING                  : num  -0.021 0.0351 0.0499 0.0655 0.2288 ...
     $ sample.num.A          : chr  "43A replic." "43A replic." "65B" "43A replic." ...
     $ sample.num.B          : chr  "65B" "45A" "45A" "45B" ...
     $ family.A              : chr  "43" "43" "65" "43" ...
     $ family.B              : chr  "65" "45" "45" "45" ...
     $ mother.A              : logi  FALSE FALSE FALSE FALSE FALSE FALSE ...
     $ mother.B              : logi  FALSE FALSE FALSE FALSE FALSE FALSE ...
     $ child.A               : logi  TRUE TRUE TRUE TRUE TRUE TRUE ...
     $ child.B               : logi  TRUE TRUE TRUE TRUE TRUE TRUE ...
     $ replic.A              : logi  TRUE TRUE FALSE TRUE FALSE FALSE ...
     $ replic.B              : logi  FALSE FALSE FALSE FALSE FALSE FALSE ...
     $ count_depth.x         : int  3376 3376 3339 3376 3147 3339 3339 3376 3147 3377 ...
     $ average_depth.x       : num  15.4 15.4 11.5 15.4 6.3 ...
     $ median_depth.x        : int  13 13 9 13 4 9 9 13 4 12 ...
     $ count_depth.y         : int  3339 3147 3147 3377 3377 3377 3381 3381 3381 3381 ...
     $ average_depth.y       : num  11.5 6.3 6.3 14.8 14.8 ...
     $ median_depth.y        : int  9 4 4 12 12 12 14 14 14 14 ...


### Colony

##### Preparing Input


```r
%%R
# Creating table containing information on each sample (File, ngsID, sample.num, family, ... )
sample_A <- subset(ngs, select=c("sample.ID.A",'a','sample.num.A','family.A','mother.A','child.A','replic.A'))#, 'Fa')) # Select columns from ngs we are interested in (a)
sample_A <- sample_A[!duplicated(sample_A),] # Remove duplicate records
colnames(sample_A) <- c('sample.ID','ID','sample.num','family','mother','child','replic')#, 'internal_relatedness') # Change col names

sample_B <- subset(ngs, select=c("sample.ID.B",'b','sample.num.B','family.B','mother.B','child.B','replic.B'))#, 'Fb')) # Repeat the same process for the 'b' records
sample_B <- sample_B[!duplicated(sample_B),]
colnames(sample_B) <- c('sample.ID','ID','sample.num','family','mother','child','replic')#, 'internal_relatedness')

sample_info <- rbind(sample_A, sample_B) # Add their rows together (stack on top of each other)
sample_info <- sample_info[!duplicated(sample_info),] # Remove duplicates
write.table(sample_info, 'SampleData/sample_info.txt', sep="\t", row.names=FALSE, col.names=TRUE, quote = FALSE ) # Write to file
```


```r
%%R
getColonyInput <- function(ngs, genoFile, sample_info){
    ### MARKER INFO
    geno <- read.delim(genoFile, header=FALSE) # Read geno file
    geno <- subset(geno, select=c(3:202)) # Get rid of columns we don't need
    geno <- data.frame(t(geno)) # Transpose dataframe (switch columns/rows)
    rownames(geno) <- 1:200 # Change row names
 
    baseToInt <- function(x){ # Function converts nucleotide to integer
        if (x=='C'){return(1)} else if(x=='G'){return(2)}
            else if(x=='T'){return(3)} else if(x=='A'){return(4)}
                else{return(0)}}
    
    geno_new <- data.frame('index'=c(0:199)) # Create new dataframe
    for (col in colnames(geno)){ # Iterate through columns of geno dataframe
        
        col.1 = paste('mk',substr(col, 2,nchar(col)),'-1',sep='') # Make new col name ex. mk1-1
        geno_new[col.1] <- unlist(lapply(geno[,col], function(x) {baseToInt(substr(x,1,1))})) # Add column with 1st allele
        
        col.2 = paste('mk',substr(col, 2,nchar(col)),'-2',sep='') # Make other col name ex. mk1-2
        geno_new[col.2] <- unlist(lapply(geno[,col], function(x) {baseToInt(substr(x,2,2))})) # Add column with 2nd allele
    }
          
    geno <- merge(sample_info, geno_new, by.x='ID', by.y="index") # Merge with sample information
    write.table(geno, 'ColonyInput/all_geno.txt', sep="\t", row.names=FALSE, col.names=FALSE, quote = FALSE ) # write to table
           
    # Offspring subset
    geno.child <- geno %>% filter(child==TRUE, replic==FALSE) %>% subset(select=c(3, 8:6813)) # Subset
    write.table(geno.child, 'ColonyInput/offspring_geno.txt', sep="\t", row.names=FALSE, col.names=FALSE, quote = FALSE ) # write to table
    
    # Maternal subset
    geno.maternal <- geno %>% filter(mother==TRUE, replic==FALSE) %>% subset(select=c(3, 8:6813)) # Subset
    write.table(geno.maternal, 'ColonyInput/maternal_geno.txt', sep="\t", row.names=FALSE, col.names=FALSE, quote = FALSE ) # write to table
    
    # Replicates subset    
    geno.replicates <- geno %>% filter(replic==TRUE) %>% subset(select=c(3, 8:6813)) # Subset
    geno.replicates$sample.num <- unlist(lapply(geno.replicates$sample.num, function(x)  # Replace space with underscore  
        {paste(unlist(strsplit(x, " "))[1], unlist(strsplit(x, " "))[2], sep="_")}))
    sampleNumReplic <- unlist(lapply(geno.replicates$sample.num, function(x) {unlist(strsplit(x, "_"))[1]} )) # Ex. getting "43A" from "43A replic."
    geno.replicates <- rbind(geno.replicates, geno %>% filter(sample.num %in% sampleNumReplic)  %>% subset(select=c(3, 8:6813))) # Adding replicate rows.
    write.table(geno.replicates, 'ColonyInput/replicates_geno.txt', sep="\t", row.names=FALSE, col.names=FALSE, quote = FALSE ) # write to table
 
        
    ### Creating Family Structure Table (rows look like: 1 -> 1a -> 1b) i.e. mother ID followed by offsrping.
    # *NOTE: we have one instance (siblings 63A and 63B) where there is no data on the mother.
    nrows = length(geno.maternal$sample.num)
    fam_str <- data.frame(rep(NA,nrows), rep(NA,nrows), rep(NA,nrows)) # Create empty dataframe
    rownames(fam_str) <- geno.maternal$sample.num # set index equal to maternal sample numbers
    colnames(fam_str) <- c(1,2,3) # Change column names (allows for easy indexing)
    offspring <- sample_info%>%filter(mother==FALSE, replic==FALSE)%>%subset(select=c(3,4)) #Filter out mothers/replicates
        
    for (i in 1:length(offspring$family)){ # Iterate through each row in offspring table
        family = offspring$family[i] # ex. family = '1'
        sample.num = offspring$sample.num[i] #ex. sample.num = '1a'
        if (family %in% rownames(fam_str) & is.na(fam_str[family,1])){ # "If family is a valid mother ID and 1st cell in row with family is empty"
            fam_str[family,1] = sample.num # Put in sample number
        } else if (family %in% rownames(fam_str) & is.na(fam_str[family,2])){ # and so on for the next cells...
            fam_str[family,2] = sample.num
        } else if (family %in% rownames(fam_str) & is.na(fam_str[family,3])){
            fam_str[family,3] = sample.num}
    } 
     
    write.table(fam_str, "ColonyInput/family_structure.txt", sep="\t", row.names=TRUE, col.names=FALSE, na="", quote = FALSE ) # write to table       
}    
```


```r
%%R
getColonyInput(ngs, './Output/geno_call.geno.gz', sample_info)
```

##### Running COLONY with replicates to approximate genotyping error rates


```python
%%time
%%capture
%%bash
cd ../Programs/COLONY
./colony2s.ifort.out IFN:../../DeerProject/ColonyInput/replicates.Dat
```

    CPU times: user 95.3 ms, sys: 7.75 ms, total: 103 ms
    Wall time: 43.6 s



```bash
%%bash
mv ../Programs/COLONY/deerReplicates.* ColonyOutput/deerReplicates
```


```python
replicErrorRates = pd.read_csv('ColonyOutput/deerReplicates/deerReplicates.ErrorRate')
replicErrorRates.describe().round(4)
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>StartDropRate</th>
      <th>DropRateEst</th>
      <th>DropRateCI95LB</th>
      <th>DropRateCI95UB</th>
      <th>StartOtherErrorRate</th>
      <th>OtherErrorRateEst</th>
      <th>OtherErrorRateCI95LB</th>
      <th>OtherErrorRateCI95UB</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>count</th>
      <td>2582.00</td>
      <td>2582.0000</td>
      <td>2582.0000</td>
      <td>2582.0000</td>
      <td>2582.0000</td>
      <td>2582.0000</td>
      <td>2582.0000</td>
      <td>2582.0000</td>
    </tr>
    <tr>
      <th>mean</th>
      <td>0.05</td>
      <td>0.0519</td>
      <td>0.0012</td>
      <td>0.4291</td>
      <td>0.0001</td>
      <td>0.0073</td>
      <td>0.0000</td>
      <td>0.4827</td>
    </tr>
    <tr>
      <th>std</th>
      <td>0.00</td>
      <td>0.1926</td>
      <td>0.0133</td>
      <td>0.2237</td>
      <td>0.0000</td>
      <td>0.0250</td>
      <td>0.0000</td>
      <td>0.2154</td>
    </tr>
    <tr>
      <th>min</th>
      <td>0.05</td>
      <td>0.0000</td>
      <td>0.0000</td>
      <td>0.0926</td>
      <td>0.0001</td>
      <td>0.0000</td>
      <td>0.0000</td>
      <td>0.0098</td>
    </tr>
    <tr>
      <th>25%</th>
      <td>0.05</td>
      <td>0.0000</td>
      <td>0.0000</td>
      <td>0.2430</td>
      <td>0.0001</td>
      <td>0.0000</td>
      <td>0.0000</td>
      <td>0.3197</td>
    </tr>
    <tr>
      <th>50%</th>
      <td>0.05</td>
      <td>0.0000</td>
      <td>0.0000</td>
      <td>0.4220</td>
      <td>0.0001</td>
      <td>0.0000</td>
      <td>0.0000</td>
      <td>0.4274</td>
    </tr>
    <tr>
      <th>75%</th>
      <td>0.05</td>
      <td>0.0000</td>
      <td>0.0000</td>
      <td>0.5062</td>
      <td>0.0001</td>
      <td>0.0000</td>
      <td>0.0000</td>
      <td>0.5877</td>
    </tr>
    <tr>
      <th>max</th>
      <td>0.05</td>
      <td>1.0000</td>
      <td>0.1675</td>
      <td>1.0000</td>
      <td>0.0001</td>
      <td>0.2612</td>
      <td>0.0002</td>
      <td>1.0000</td>
    </tr>
  </tbody>
</table>
</div>




```python
with open('ColonyInput/sibsMarkerInfo.txt','w') as markerFile:
    markerFile.write('mk@\n')
    markerFile.write('0@\n')
    markerFile.write(replicErrorRates['DropRateEst'].mean().round(4).astype('str') + '@\n')
    markerFile.write(replicErrorRates['OtherErrorRateEst'].mean().round(4).astype('str') + '@\n')
```

COLONY run for sibship inferences using estimated average alleic dropout rate and other error rate from replicates.


```python
%%capture ColonySibs1
%%time
%%bash
cd ../Programs/COLONY
./colony2s.ifort.out IFN:../../DeerProject/ColonyInput/sibs.Dat
```


```bash
%%bash
mv ../Programs/COLONY/deerSibs.* ColonyOutput/deerSibs
```

COLONY run for sibship inferences using loci specific error rates from inferred pedigree.


```r
%%R
error_rates <- read.csv("./ColonyOutput/deerSibs/deerSibs.ErrorRate")
error_rates$dominance = 0
error_rates$OtherErrorRateEst <- sapply(error_rates$OtherErrorRateEst, function(x){if (x==0){return(0.0001)}else{return(x)}})
error_rates <- t(error_rates[,c(1,10,3,7)])
write.table(error_rates, "ColonyInput/sibs2MarkerInfo.txt", sep="\t", row.names=FALSE, col.names=FALSE, quote = FALSE )
```


```python
%%time
%%capture ColonySibs2
%%bash
cd ../Programs/COLONY
./colony2s.ifort.out IFN:../../DeerProject/ColonyInput/sibs2.Dat
```

    CPU times: user 27.3 s, sys: 4.02 s, total: 31.3 s
    Wall time: 2h 38min 16s



```bash
%%bash
mv ../Programs/COLONY/deerSibs2.* ColonyOutput/deerSibs2
```

### Results


```r
%%R
halfSibDyad <- read.csv('ColonyOutput/deerSibs2/deerSibs2.HalfSibDyad') # Read in COLONY output for half sibs
colnames(halfSibDyad)[3] <- 'Prob.Half' # Change column name
fullSibDyad <- read.csv('ColonyOutput/deerSibs2/deerSibs2.FullSibDyad') # Read in COLONY output for full sibs
colnames(fullSibDyad)[3] <- 'Prob.Full' # Change column name
# Merge fullSib and halfSib dataframes 
Sib.Pairs <- merge(halfSibDyad, fullSibDyad, by.x= c('OffspringID1','OffspringID2'), by.y = c('OffspringID1','OffspringID2'), all.x =TRUE, all.y = TRUE )
Sib.Pairs[is.na(Sib.Pairs)] <- 0 # Replace NAs with 0s
sibAssignFunc <- function(x,y){if (x>0.95){return("Half")} else if(y>0.95){return("Full")} else{return('Unknown')}} # 0.95 is our probability cutoff
Sib.Pairs$assignment <- mapply(sibAssignFunc, Sib.Pairs$Prob.Half, Sib.Pairs$Prob.Full ) # Apply previous function to our half/full probabilities

findFamily <- function(sampleID){ # Function that inputs sample ID and outputs family. Ex. '10A' -> '10'
        family = ''
        for (char in strsplit(sampleID,"")[[1]]){
            if (!is.na(as.numeric(char))){
                family = paste(family, char, sep = '')}
            else{break}}
        return (family)}


removeNonFam <- function(x,y){return(findFamily(x)==findFamily(y))} # Returns true if familes are the same. Ex. (10A, 10B) -> TRUE
Sib.Pairs$Same.Family <- mapply(removeNonFam, Sib.Pairs$OffspringID1, Sib.Pairs$OffspringID2) # Applies previous function to our dataframe
Sib.Pairs <- Sib.Pairs[Sib.Pairs$Same.Family,] # Remove entries that are not of the same family
head(Sib.Pairs) # Display results
```

      OffspringID1 OffspringID2 Prob.Half Prob.Full assignment Same.Family
    1          10A          10B     0.000     1.000       Full        TRUE
    2          11A          11B     0.000     1.000       Full        TRUE
    4          12B          12A     0.693     0.307    Unknown        TRUE
    6          13A          13B     1.000     0.000       Half        TRUE
    7          14B          14A     1.000     0.000       Half        TRUE
    8          15A          15B     1.000     0.000       Half        TRUE



```r
%%R
sibPairs <- left_join(sibPairs, Sib.Pairs, by = c('sample.num.A' = 'OffspringID1', 'sample.num.B' = 'OffspringID2'))
```


```python
%R str(sibPairs)
```

    'data.frame':	82 obs. of  47 variables:
     $ b                     : int  106 107 108 11 113 115 116 118 121 121 ...
     $ a                     : int  85 17 61 9 105 111 101 114 27 32 ...
     $ sample.ID.A           : chr  "S8_AGAC" "S2_GACT" "S6_AGAC" "S1_TCAC" ...
     $ sample.ID.B           : chr  "S9_TCAG" "S9_TGTC" "S10_ACCA" "S1_TGTC" ...
     $ file_a                : chr  "ILL-RAD08_S8_L001_R1_001_AGAC.trim.bt2.bam" "ILL-RAD02_S2_L001_R1_001_GACT.trim.bt2.bam" "ILL-RAD06_S6_L001_R1_001_AGAC.trim.bt2.bam" "ILL-RAD01_S1_L001_R1_001_TCAC.trim.bt2.bam" ...
     $ file_b                : chr  "ILL-RAD09_S9_L001_R1_001_TCAG.trim.bt2.bam" "ILL-RAD09_S9_L001_R1_001_TGTC.trim.bt2.bam" "ILL-RAD10_S10_L001_R1_001_ACCA.trim.bt2.bam" "ILL-RAD01_S1_L001_R1_001_TGTC.trim.bt2.bam" ...
     $ nSites                : int  1863 2393 3462 3576 1766 32 30 3172 3439 3531 ...
     $ rab                   : num  0.307 0.455 0.432 0.483 0.468 ...
     $ Fa                    : num  0.00001 0.0303 0.00936 0.03562 0.01789 ...
     $ Fb                    : num  0.01711 0.01468 0.01294 0.00511 0.00561 ...
     $ theta                 : num  0.156 0.233 0.219 0.246 0.235 ...
     $ inbred_relatedness_1_2: num  0 0.015148 0.004675 0.018289 0.000003 ...
     $ inbred_relatedness_2_1: num  0.00855 0.007337 0.006278 0.000981 0.002807 ...
     $ F_diff_a_b            : num  -0.000002 -0.000004 -0.000186 -0.002054 0.008944 ...
     $ R0                    : num  0.1162 0.1135 0.1163 0.0895 0.081 ...
     $ R1                    : num  0.435 0.628 0.625 0.623 0.625 ...
     $ KING                  : num  0.183 0.222 0.22 0.233 0.238 ...
     $ sample.num.A          : chr  "36B" "55B" "59B" "47B" ...
     $ sample.num.B          : chr  "36A" "55A" "59A" "47A" ...
     $ family.A              : chr  "36" "55" "59" "47" ...
     $ family.B              : chr  "36" "55" "59" "47" ...
     $ mother.A              : logi  FALSE FALSE FALSE FALSE FALSE FALSE ...
     $ mother.B              : logi  FALSE FALSE FALSE FALSE FALSE FALSE ...
     $ child.A               : logi  TRUE TRUE TRUE TRUE TRUE TRUE ...
     $ child.B               : logi  TRUE TRUE TRUE TRUE TRUE TRUE ...
     $ replic.A              : logi  FALSE FALSE FALSE FALSE FALSE FALSE ...
     $ replic.B              : logi  FALSE FALSE FALSE FALSE FALSE FALSE ...
     $ count_depth.x         : int  3355 3111 3361 3372 1719 187 197 3272 3280 3350 ...
     $ average_depth.x       : num  11.39 7.28 13.7 16.18 2.75 ...
     $ median_depth.x        : int  9 4 11 13 1 1 1 5 7 10 ...
     $ count_depth.y         : int  1748 2408 3208 3269 3216 158 200 2999 3275 3275 ...
     $ average_depth.y       : num  3.33 3 6.73 9.24 6.71 ...
     $ median_depth.y        : int  2 2 5 7 5 1 1 3 6 6 ...
     $ rab_lower             : num  0.239 0.393 0.393 0.447 0.4 ...
     $ rab_upper             : num  0.373 0.512 0.471 0.519 0.539 ...
     $ R0_lower              : num  0.0338 0.0538 0.0795 0.0593 0.0215 ...
     $ R0_upper              : num  0.214 0.182 0.156 0.123 0.16 ...
     $ R1_lower              : num  0.334 0.5 0.551 0.551 0.468 ...
     $ R1_upper              : num  0.556 0.777 0.706 0.7 0.825 ...
     $ min_med_depth         : int  2 2 5 7 1 1 1 3 6 6 ...
     $ min_avg_depth         : num  3.33 3 6.73 9.24 2.75 ...
     $ avg_med_depth         : num  5.5 3 8 10 3 1 1 4 6.5 8 ...
     $ avg_avg_depth         : num  7.36 5.14 10.21 12.71 4.73 ...
     $ Prob.Half             : num  1 0 0 0 0 0.219 0.246 0 0 0 ...
     $ Prob.Full             : num  0 1 1 1 1 0.781 0.754 1 1 1 ...
     $ assignment            : chr  "Half" "Full" "Full" "Full" ...
     $ Same.Family           : logi  TRUE TRUE TRUE TRUE TRUE TRUE ...



```r
%%R
p1 <- ggplot(sibPairs, aes(x=reorder(family.A, rab), y=rab, color=assignment)) +
    geom_point() + 
    theme_bw() + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"), 
          axis.text.x = element_blank(),
          plot.caption.position = "plot",
          plot.caption = element_text(hjust = 0)) +
    geom_pointrange(aes(ymin = rab_lower, ymax = rab_upper),width=0.1, alpha = 0.5) +
    labs(color='assignment') +  
    xlab("Litter (ordered by relatedness)") +
    ylab("Relatedness") +
    scale_color_manual('Sibship',values=c("#EE8866", "#44BB99", "#BBCC33")) +
    theme(legend.key.width = unit(2, 'cm'), legend.title = element_text(size=20), legend.text = element_text(size=20),
         axis.title.x = element_text(size = 16), axis.text.y = element_text(size = 16), axis.title.y = element_text(size=16),
          legend.position = c(.95, .95),
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.margin = margin(6, 6, 6, 6))


p2 <- p1 + scale_color_manual('Sibship',values=c("grey70", "grey10", "grey40"))
print(p1)
print(p2)
```

    R[write to console]: Scale for 'colour' is already present. Adding another scale for 'colour',
    which will replace the existing scale.
    



    
![png](output_142_1.png)
    



    
![png](output_142_2.png)
    



```r
%%R
family <- unlist(lapply(Sib.Pairs$OffspringID1, function(x){findFamily(x)}))
familyUnique <- unique(family)
familyAssignments <- data.frame(rep(0,length(familyUnique)))
colnames(familyAssignments) <- c("Evd.Of.Mult.Pat")
rownames(familyAssignments) <- familyUnique
for (i in 1:length(Sib.Pairs$OffspringID1)){
    family = findFamily(Sib.Pairs$OffspringID1[i])
    if (Sib.Pairs$assignment[i] == 'Half'){
        familyAssignments[family, 1] = 1
    }
    else if (Sib.Pairs$assignment[i] == 'Unknown'){
        familyAssignments[family, 1] = 'Unknown'
    }
}
```

% with multiple paternity


```r
%%R
known <- familyAssignments %>% filter(Evd.Of.Mult.Pat != 'Unknown')
mean(as.numeric(known$Evd.Of.Mult.Pat))
```

    [1] 0.2241379



```r
%%R
familyAssignments
```

       Evd.Of.Mult.Pat
    10               0
    11               0
    12         Unknown
    13               1
    14               1
    15               1
    16               0
    17               1
    1          Unknown
    20               0
    21               0
    22               0
    24               0
    26               0
    27               0
    28               0
    29               1
    2                1
    30               0
    31               0
    32               0
    33               0
    34               0
    35               0
    36               1
    37               1
    38               0
    39               0
    3          Unknown
    40               0
    41               0
    42               0
    43               0
    44               1
    45               1
    46               0
    47               0
    48               0
    49               0
    4                0
    50               0
    51               1
    52               0
    53               0
    54               0
    55               0
    56               0
    57               0
    58               0
    59               0
    5                0
    60               0
    61               0
    62               0
    63               1
    64               0
    65               0
    66               0
    67               1
    68               0
    8                0
    9          Unknown



```r
%%R
table(familyAssignments$Evd.Of.Mult.Pat)
```

    
          0       1 Unknown 
         45      13       4 



```r
%%R
pieDist <- data.frame(group=c('Unknown','Single','Multiple'), value=c(4,45,13))
ggplot(pieDist, aes(x='', y=value, fill=group)) + 
    geom_bar(stat="identity",width=1,color='white') + 
    coord_polar('y',start=0) +
    theme_void() +
    #scale_fill_manual("Paternity",values=c("#44BB99","#EE8866",  "#BBCC33")) +
    #geom_text(aes(label=value), position = position_stack(vjust=0.5), size =10) +
    theme(legend.key.size = unit(2, 'cm')) + 
    theme(legend.key.width = unit(1, 'cm'), legend.title = element_text(size=20), legend.text = element_text(size=20)) +
    scale_fill_manual('Paternity',values=c("grey10", "grey70", "grey40"))

    
```


    
![png](output_148_0.png)
    


### SNP distribution for sibling pairs


```r
%%R
ggplot(sibPairs, aes(x=nSites, fill = assignment, color = assignment)) +
    geom_histogram(bins = 20, position='identity', alpha=0.5) +
    theme_bw() +
    ggtitle("SNP count distribution")
```


    
![png](output_150_0.png)
    



```r
%%R
ggplot(sibPairs, aes(x=assignment, y=nSites, fill=assignment)) +
  geom_boxplot() +
  theme_classic()

```


    
![png](output_151_0.png)
    


### Coverage depth for sibling pairs


```r
%%R 
sibPairs$min_med_depth <- mapply(function(x,y){return(min(x,y))},sibPairs$median_depth.x, sibPairs$median_depth.y)
sibPairs$min_avg_depth <- mapply(function(x,y){return(min(x,y))},sibPairs$average_depth.x, sibPairs$average_depth.y)
sibPairs$avg_med_depth <- mapply(function(x,y){return(mean(c(x,y)))},sibPairs$median_depth.x, sibPairs$median_depth.y)
sibPairs$avg_avg_depth <- mapply(function(x,y){return(mean(c(x,y)))},sibPairs$average_depth.x, sibPairs$average_depth.y)
```


```r
%%R
ggplot(sibPairs, aes(x=assignment, y=min_med_depth, fill=assignment)) + 
    geom_boxplot() +
    theme_classic()
```


    
![png](output_154_0.png)
    



```r
%%R
ggplot(sibPairs, aes(x=assignment, y=min_avg_depth, fill=assignment)) + 
    geom_boxplot() +
    theme_classic()
```


    
![png](output_155_0.png)
    



```r
%%R
ggplot(sibPairs, aes(x=assignment, y=avg_med_depth, fill=assignment)) + 
    geom_boxplot() +
    theme_classic()
```


    
![png](output_156_0.png)
    



```r
%%R
ggplot(sibPairs, aes(x=assignment, y=avg_avg_depth, fill=assignment)) + 
    geom_boxplot() +
    theme_classic()
```


    
![png](output_157_0.png)
    



```r
%%R
sibPairsKnown <- sibPairs %>% filter(assignment != 'Unknown')
readQualityLogit <- glm(factor(assignment) ~ avg_med_depth + nSites + nSites * avg_med_depth, data= sibPairsKnown, family='binomial')
summary(readQualityLogit)
```

    
    Call:
    glm(formula = factor(assignment) ~ avg_med_depth + nSites + nSites * 
        avg_med_depth, family = "binomial", data = sibPairsKnown)
    
    Deviance Residuals: 
        Min       1Q   Median       3Q      Max  
    -1.4677  -0.6703  -0.5222  -0.2983   2.1541  
    
    Coefficients:
                           Estimate Std. Error z value Pr(>|z|)
    (Intercept)          -3.1793366  3.1268759  -1.017    0.309
    avg_med_depth         0.7870571  0.6002302   1.311    0.190
    nSites                0.0008185  0.0010185   0.804    0.422
    avg_med_depth:nSites -0.0002847  0.0001873  -1.520    0.129
    
    (Dispersion parameter for binomial family taken to be 1)
    
        Null deviance: 76.370  on 77  degrees of freedom
    Residual deviance: 68.601  on 74  degrees of freedom
    AIC: 76.601
    
    Number of Fisher Scoring iterations: 5
    



```r
%%R
predpr = predict(readQualityLogit, type='response')
predlogit = predict(readQualityLogit)
plot(jitter(ifelse(sibPairsKnown$assignment=='Half',1,0),0.1) ~ predlogit, 
xlab='predicted logit', ylab="Probability of Half siblings")
pred.ord = order(predlogit)
lines(predlogit[pred.ord], predpr[pred.ord])
lines (lowess(predlogit[pred.ord],predpr [pred.ord]), col='red',lty=2)
```


    
![png](output_159_0.png)
    


### Relatedness for sibling pairs


```r
%%R
ggplot(sibPairs, aes(x=assignment, y=rab, fill=assignment)) + 
    geom_boxplot() +
    theme_classic()
```


    
![png](output_161_0.png)
    



```r
%%R
t.test(rab ~ assignment, data=sibPairsKnown)
```

    
    	Welch Two Sample t-test
    
    data:  rab by assignment
    t = 10.327, df = 15.843, p-value = 1.93e-08
    alternative hypothesis: true difference in means between group Full and group Half is not equal to 0
    95 percent confidence interval:
     0.1800790 0.2732057
    sample estimates:
    mean in group Full mean in group Half 
             0.4769129          0.2502705 
    


### Comparing Mother, Offspring relatedness between mult. pat and non mult. pat groups


```r
%%R
motherOffspring <- ngs %>% filter(family.A == family.B & mother.A + mother.B == 1 & replic.A + replic.B == 0)
motherOffspring <- merge(motherOffspring, familyAssignments, by.x='family.A', by.y='row.names')
motherOffspring <- motherOffspring %>% filter(nSites > 500)
str(motherOffspring)
```

    'data.frame':	112 obs. of  34 variables:
     $ family.A              : chr  "10" "10" "11" "11" ...
     $ b                     : int  97 122 160 160 142 142 72 126 91 94 ...
     $ a                     : int  25 97 123 127 74 130 14 72 87 87 ...
     $ sample.ID.A           : chr  "S3_AGAC" "S9_AGAC" "S11_CATC" "S11_GTGA" ...
     $ sample.ID.B           : chr  "S9_AGAC" "S11_AGTG" "S14_CTAC" "S14_CTAC" ...
     $ file_a                : chr  "ILL-RAD03_S3_L001_R1_001_AGAC.trim.bt2.bam" "ILL-RAD09_S9_L001_R1_001_AGAC.trim.bt2.bam" "ILL-RAD11_S11_L001_R1_001_CATC.trim.bt2.bam" "ILL-RAD11_S11_L001_R1_001_GTGA.trim.bt2.bam" ...
     $ file_b                : chr  "ILL-RAD09_S9_L001_R1_001_AGAC.trim.bt2.bam" "ILL-RAD11_S11_L001_R1_001_AGTG.trim.bt2.bam" "ILL-RAD14_S14_L001_R1_001_CTAC.trim.bt2.bam" "ILL-RAD14_S14_L001_R1_001_CTAC.trim.bt2.bam" ...
     $ nSites                : int  3361 2991 3457 3462 3280 3091 3570 3552 3389 3359 ...
     $ rab                   : num  0.474 0.494 0.458 0.458 0.511 ...
     $ Fa                    : num  0.037363 0.004812 0.044335 0.063107 0.000028 ...
     $ Fb                    : num  0.00287 0.05004 0.01457 0.00757 0 ...
     $ theta                 : num  0.238 0.255 0.237 0.239 0.256 ...
     $ inbred_relatedness_1_2: num  0.002821 0.004807 0.017647 0.035338 0.000014 ...
     $ inbred_relatedness_2_1: num  0.00285 0.02742 0.01421 0.00757 0 ...
     $ F_diff_a_b            : num  0.0173 0 0.0114 0 0 ...
     $ R0                    : num  0.0582 0.0177 0.06265 0.05325 0.00343 ...
     $ R1                    : num  0.515 0.522 0.48 0.488 0.594 ...
     $ KING                  : num  0.227 0.247 0.217 0.224 0.27 ...
     $ sample.num.A          : chr  "10A" "10" "11A" "11B" ...
     $ sample.num.B          : chr  "10" "10B" "11" "11" ...
     $ family.B              : chr  "10" "10" "11" "11" ...
     $ mother.A              : logi  FALSE TRUE FALSE FALSE FALSE FALSE ...
     $ mother.B              : logi  TRUE FALSE TRUE TRUE TRUE TRUE ...
     $ child.A               : logi  TRUE FALSE TRUE TRUE TRUE TRUE ...
     $ child.B               : logi  FALSE TRUE FALSE FALSE FALSE FALSE ...
     $ replic.A              : logi  FALSE FALSE FALSE FALSE FALSE FALSE ...
     $ replic.B              : logi  FALSE FALSE FALSE FALSE FALSE FALSE ...
     $ count_depth.x         : int  3374 2653 3009 3317 3325 3339 3338 3328 3184 3184 ...
     $ average_depth.x       : num  18.21 3.29 5.3 8.61 10.74 ...
     $ median_depth.x        : int  16 2 3 7 8 7 8 8 4 4 ...
     $ count_depth.y         : int  2653 3275 3287 3287 1481 1481 3328 3287 3247 3187 ...
     $ average_depth.y       : num  3.29 8.1 10.09 10.09 2.49 ...
     $ median_depth.y        : int  2 6 8 8 1 1 8 9 4 4 ...
     $ Evd.Of.Mult.Pat       : chr  "0" "0" "0" "0" ...



```r
%%R
ggplot(data=motherOffspring, aes(x=rab, fill=Evd.Of.Mult.Pat, color=Evd.Of.Mult.Pat)) + 
    geom_histogram(position='identity', alpha=0.5,bins=30) +
    geom_vline(xintercept=mean((motherOffspring %>% filter(Evd.Of.Mult.Pat == 1))$rab, na.rm=TRUE), color='green') +
    geom_vline(xintercept=mean((motherOffspring %>% filter(Evd.Of.Mult.Pat == 0))$rab, na.rm=TRUE), color='#F8766D') +
    theme_classic()
```


    
![png](output_165_0.png)
    



```r
%%R
ggplot(motherOffspring, aes(x=Evd.Of.Mult.Pat, y=rab, fill=Evd.Of.Mult.Pat)) + 
    geom_boxplot() +
    theme_classic()
```


    
![png](output_166_0.png)
    



```r
%%R
qplot(sample = rab, data= motherOffspring, color=Evd.Of.Mult.Pat) +
    theme_classic() + 
    stat_qq_line()
```


    
![png](output_167_0.png)
    


QQ plots show good evidence of normality in the non-mult. pat. group, and some deviance from normality in the mult. pat. group.


```r
%%R
rab0 <- as.numeric(unlist(motherOffspring %>% filter(Evd.Of.Mult.Pat == 0) %>% select(rab)))
rab1 <- as.numeric(unlist(motherOffspring %>% filter(Evd.Of.Mult.Pat == 1) %>% select(rab)))
shapiro.test(rab0) # Very normal
```

    
    	Shapiro-Wilk normality test
    
    data:  rab0
    W = 0.99561, p-value = 0.9921
    



```r
%%R
shapiro.test(rab1) # Statistically significant deviation from normality
```

    
    	Shapiro-Wilk normality test
    
    data:  rab1
    W = 0.90138, p-value = 0.04378
    



```r
%%R
var.test(rab ~ Evd.Of.Mult.Pat, data = (motherOffspring %>% filter(Evd.Of.Mult.Pat != 'Unknown')), alternative = 'two.sided')
```

    
    	F test to compare two variances
    
    data:  rab by Evd.Of.Mult.Pat
    F = 0.5851, num df = 90, denom df = 19, p-value = 0.09783
    alternative hypothesis: true ratio of variances is not equal to 1
    95 percent confidence interval:
     0.2628871 1.1023075
    sample estimates:
    ratio of variances 
             0.5851013 
    


Variances are not significantly different (*P = 0.09783)*


```r
%%R
t.test(rab0, rab1)
```

    
    	Welch Two Sample t-test
    
    data:  rab0 and rab1
    t = -2.2736, df = 24.117, p-value = 0.03218
    alternative hypothesis: true difference in means is not equal to 0
    95 percent confidence interval:
     -0.035620648 -0.001726647
    sample estimates:
    mean of x mean of y 
    0.4729008 0.4915744 
    


Mean relatedness between mothers and their offspring are significantly higher in mult. pat. litters than in non. mult. pat. (*P = 0.03218*)

Confirm with non-parametric test (Mann-Whitney U test), since there are outliers in the multiple paternity group. 


```r
%%R
rabManWU <-wilcox.test(rab ~ Evd.Of.Mult.Pat, data = (motherOffspring %>% filter(Evd.Of.Mult.Pat != 'Unknown')), na.rm=TRUE, paired=FALSE, exact=FALSE, conf.int=TRUE)
print(rabManWU)
```

    
    	Wilcoxon rank sum test with continuity correction
    
    data:  rab by Evd.Of.Mult.Pat
    W = 543, p-value = 0.004923
    alternative hypothesis: true location shift is not equal to 0
    95 percent confidence interval:
     -0.031959317 -0.005912661
    sample estimates:
    difference in location 
               -0.01849531 
    


Difference is more significant in non-parametric test (*P = 0.004923*)

### Comparing inbreeding coefficient


```r
%%R
str(motherOffspring)
```

    'data.frame':	112 obs. of  34 variables:
     $ family.A              : chr  "10" "10" "11" "11" ...
     $ b                     : int  97 122 160 160 142 142 72 126 91 94 ...
     $ a                     : int  25 97 123 127 74 130 14 72 87 87 ...
     $ sample.ID.A           : chr  "S3_AGAC" "S9_AGAC" "S11_CATC" "S11_GTGA" ...
     $ sample.ID.B           : chr  "S9_AGAC" "S11_AGTG" "S14_CTAC" "S14_CTAC" ...
     $ file_a                : chr  "ILL-RAD03_S3_L001_R1_001_AGAC.trim.bt2.bam" "ILL-RAD09_S9_L001_R1_001_AGAC.trim.bt2.bam" "ILL-RAD11_S11_L001_R1_001_CATC.trim.bt2.bam" "ILL-RAD11_S11_L001_R1_001_GTGA.trim.bt2.bam" ...
     $ file_b                : chr  "ILL-RAD09_S9_L001_R1_001_AGAC.trim.bt2.bam" "ILL-RAD11_S11_L001_R1_001_AGTG.trim.bt2.bam" "ILL-RAD14_S14_L001_R1_001_CTAC.trim.bt2.bam" "ILL-RAD14_S14_L001_R1_001_CTAC.trim.bt2.bam" ...
     $ nSites                : int  3361 2991 3457 3462 3280 3091 3570 3552 3389 3359 ...
     $ rab                   : num  0.474 0.494 0.458 0.458 0.511 ...
     $ Fa                    : num  0.037363 0.004812 0.044335 0.063107 0.000028 ...
     $ Fb                    : num  0.00287 0.05004 0.01457 0.00757 0 ...
     $ theta                 : num  0.238 0.255 0.237 0.239 0.256 ...
     $ inbred_relatedness_1_2: num  0.002821 0.004807 0.017647 0.035338 0.000014 ...
     $ inbred_relatedness_2_1: num  0.00285 0.02742 0.01421 0.00757 0 ...
     $ F_diff_a_b            : num  0.0173 0 0.0114 0 0 ...
     $ R0                    : num  0.0582 0.0177 0.06265 0.05325 0.00343 ...
     $ R1                    : num  0.515 0.522 0.48 0.488 0.594 ...
     $ KING                  : num  0.227 0.247 0.217 0.224 0.27 ...
     $ sample.num.A          : chr  "10A" "10" "11A" "11B" ...
     $ sample.num.B          : chr  "10" "10B" "11" "11" ...
     $ family.B              : chr  "10" "10" "11" "11" ...
     $ mother.A              : logi  FALSE TRUE FALSE FALSE FALSE FALSE ...
     $ mother.B              : logi  TRUE FALSE TRUE TRUE TRUE TRUE ...
     $ child.A               : logi  TRUE FALSE TRUE TRUE TRUE TRUE ...
     $ child.B               : logi  FALSE TRUE FALSE FALSE FALSE FALSE ...
     $ replic.A              : logi  FALSE FALSE FALSE FALSE FALSE FALSE ...
     $ replic.B              : logi  FALSE FALSE FALSE FALSE FALSE FALSE ...
     $ count_depth.x         : int  3374 2653 3009 3317 3325 3339 3338 3328 3184 3184 ...
     $ average_depth.x       : num  18.21 3.29 5.3 8.61 10.74 ...
     $ median_depth.x        : int  16 2 3 7 8 7 8 8 4 4 ...
     $ count_depth.y         : int  2653 3275 3287 3287 1481 1481 3328 3287 3247 3187 ...
     $ average_depth.y       : num  3.29 8.1 10.09 10.09 2.49 ...
     $ median_depth.y        : int  2 6 8 8 1 1 8 9 4 4 ...
     $ Evd.Of.Mult.Pat       : chr  "0" "0" "0" "0" ...



```r
%%R
# Creating table containing information on each sample (File, ngsID, sample.num, family, ... )
sample_A <- subset(ngs, select=c("sample.ID.A",'a','sample.num.A','family.A','mother.A','child.A','replic.A', 'Fa')) # Select columns from ngs we are interested in (a)
sample_A <- sample_A[!duplicated(sample_A),] # Remove duplicate records
colnames(sample_A) <- c('sample.ID','ID','sample.num','family','mother','child','replic', 'internal_relatedness') # Change col names

sample_B <- subset(ngs, select=c("sample.ID.B",'b','sample.num.B','family.B','mother.B','child.B','replic.B', 'Fb')) # Repeat the same process for the 'b' records
sample_B <- sample_B[!duplicated(sample_B),]
colnames(sample_B) <- c('sample.ID','ID','sample.num','family','mother','child','replic', 'internal_relatedness')

fGrouped <- rbind(sample_A, sample_B) # Add their rows together (stack on top of each other)
#fGrouped %>% group_by('ID') %>% summarise(across(internal_relatedness, mean, na.rm = TRUE))
fGrouped <- aggregate(fGrouped[, 8], list(fGrouped$ID), mean)
colnames(fGrouped) <- c('ID','mean.F')
head(fGrouped)

```

      ID      mean.F
    1  0 0.024127309
    2  1 0.021983203
    3  2 0.001309133
    4  3 0.000419800
    5  4 0.039986274
    6  5 0.035172513



```r
%%R
ggplot(motherOffspring, aes(x=Evd.Of.Mult.Pat, y=F_diff_a_b, fill=Evd.Of.Mult.Pat)) + 
    geom_boxplot() +
    theme_classic()
```


    
![png](output_181_0.png)
    



```r
%%R
chkdata <-
function(genotypes) {

   genotypes <- as.matrix(genotypes)

   if (ncol(genotypes) %% 2 == 1) {
      return("Odd number of columns in the input.")
   }

   individuals <- nrow(genotypes)
   loci <- ncol(genotypes) / 2
   n_alleles <- array(loci)

   for (l in 1:loci) {
      g <- 2 * l - 1
      h <- 2 * l
      frequencies <- list(table(genotypes[, g:h]))
      n_alleles[l] <- length(frequencies[[1]])
   }

   if (sum(n_alleles < 2) > 0) {
      loci_list <- ""
      for (i in 1:loci) {
         if (n_alleles[i] < 2) {
            if (loci_list == "") {
               loci_list <- i
            }
            else {
               loci_list <- paste(loci_list, i, sep = ", ")
            }
         }
      }
      return(paste("Only one allele in loci:", loci_list))
   }

   write(n_alleles, file="number_of_alleles.txt", ncolumns=1)

   k <- 0
   ind_loci_list <- "\n"

   for (i in 1:individuals) {

      for (l in 1:loci) {

         g <- 2 * l - 1
         h <- 2 * l

         if (xor(!is.na(genotypes[i, g]), !is.na(genotypes[i, h]))) {
            ind_loci_list <- paste(ind_loci_list, i, l, "\n")
            k <- k + 1
         }
      }
   }

   if (k > 0) {
      return(paste("One or more individuals are missing one allele in one or more loci:\nIndividual  Locus", ind_loci_list))
   }
    else{
        return('good to go')
    }

}


```


```r
%%R
hl <-
function(genotypes) {

   genotypes <- as.matrix(genotypes)

   individuals <- nrow(genotypes)
   loci <- ncol(genotypes) / 2
   hl <- array(NA, dim=c(individuals, 1))
   E <- array(loci)
   frequencies <- array(loci)

   for (l in 1:loci) {
      E[l] <- 1
      g <- 2 * l - 1
      h <- 2 * l
      frequencies[l] <- list(table(genotypes[, g:h]))
      E[l] <- 1 - sum((frequencies[[l]] / sum(frequencies[[l]]))^2)
   }

   for (i in 1:individuals) {

      sum.Eh <- 0
      sum.Ej <- 0

      for (l in 1:loci) {

         g <- 2 * l - 1
         h <- 2 * l

         if ((!is.na(genotypes[i, g])) && (!is.na(genotypes[i, h]))) {
            if (genotypes[i, g] == genotypes[i, h]) {
               sum.Eh <- sum.Eh + E[l]
            }
            else {
               sum.Ej <- sum.Ej + E[l]
            }
         }
      }

      hl[i] <- sum.Eh / (sum.Eh + sum.Ej)

   }

   return(hl)

}


```

### Hl coefficient


```bash
%%bash 
gunzip Output/geno_call.mafs.gz
```


```r
%%R
alleleFreq <- read.delim('Output/geno_call.mafs')
head(alleleFreq)
```

          chromo position major minor  knownEM pK.EM nInd
    1 CM042924.1   384438     G     C 0.063544     0  183
    2 CM042924.1   400489     C     T 0.050972     0  176
    3 CM042924.1  1567263     A     G 0.180832     0  178
    4 CM042924.1  1608587     C     T 0.143123     0  170
    5 CM042924.1  1608591     T     C 0.137723     0  172
    6 CM042924.1  1737830     C     T 0.167873     0  153



```r
%%R
genoCall <- read.delim('./Output/geno_call.geno.gz')
genoCall <- subset(genoCall, select=c(3:202)) # Get rid of columns we don't need
genoCall <- data.frame(t(genoCall)) # Transpose dataframe (switch columns/rows)
rownames(genoCall) <- 1:200 # Change row names
```

Subsetting genoCall to find allele frequencies. Taking one offspring representative from each litter, so that our population allele frequencies are unbiased.


```r
%%R
genoCallMothers <-  merge(sample_info, genoCall, by.x='ID', by.y='row.names')
genoCallMothers <- genoCallUR %>% filter(mother == TRUE)
genoCallMothers <- subset(genoCallMothers, select=c(8:3409))
```


```r
%%R
HLResult <- data.frame(hl(genoCall))
```


```python
%R familyAssignments
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>Evd.Of.Mult.Pat</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>10</th>
      <td>0</td>
    </tr>
    <tr>
      <th>11</th>
      <td>0</td>
    </tr>
    <tr>
      <th>12</th>
      <td>Unknown</td>
    </tr>
    <tr>
      <th>13</th>
      <td>1</td>
    </tr>
    <tr>
      <th>14</th>
      <td>1</td>
    </tr>
    <tr>
      <th>...</th>
      <td>...</td>
    </tr>
    <tr>
      <th>66</th>
      <td>0</td>
    </tr>
    <tr>
      <th>67</th>
      <td>1</td>
    </tr>
    <tr>
      <th>68</th>
      <td>0</td>
    </tr>
    <tr>
      <th>8</th>
      <td>0</td>
    </tr>
    <tr>
      <th>9</th>
      <td>Unknown</td>
    </tr>
  </tbody>
</table>
<p>62 rows × 1 columns</p>
</div>




```r
%%R
HLMerged <- merge(sample_info, HLResult, by.x='ID', by.y='row.names')
HLMerged <- merge(HLMerged, familyAssignments, by.x='family', by.y='row.names')
str(HLMerged) # Figure out what happened to one sample
```

    'data.frame':	199 obs. of  9 variables:
     $ family         : chr  "1" "1" "1" "10" ...
     $ ID             : int  83 116 101 25 122 97 127 160 123 143 ...
     $ sample.ID      : chr  "S7_TGTC" "S10_GTGT" "S9_GACT" "S3_AGAC" ...
     $ sample.num     : chr  "1" "1B" "1A" "10A" ...
     $ mother         : logi  TRUE FALSE FALSE FALSE FALSE TRUE ...
     $ child          : logi  FALSE TRUE TRUE TRUE TRUE FALSE ...
     $ replic         : logi  FALSE FALSE FALSE FALSE FALSE FALSE ...
     $ hl.genoCall.   : num  0.132 0.928 0.188 0.14 0.132 ...
     $ Evd.Of.Mult.Pat: chr  "Unknown" "Unknown" "Unknown" "0" ...



```r
%%R
ggplot(HLMerged, aes(x=Evd.Of.Mult.Pat, y=hl.genoCall., fill=Evd.Of.Mult.Pat)) + 
    geom_boxplot() +
    #geom_text_repel(data=outliers, aes(x=0,label=nSites)) +
    theme_classic() #Sequence depth and number of SNPs impact.
```


    
![png](output_193_0.png)
    



```r
%%R
HLManWU <-wilcox.test(hl.genoCall. ~ Evd.Of.Mult.Pat, data = (HLMerged %>% filter(Evd.Of.Mult.Pat != 'Unknown')), na.rm=TRUE, paired=FALSE, exact=FALSE, conf.int=TRUE)
print(HLManWU)
```

    
    	Wilcoxon rank sum test with continuity correction
    
    data:  hl.genoCall. by Evd.Of.Mult.Pat
    W = 2759, p-value = 0.4458
    alternative hypothesis: true location shift is not equal to 0
    95 percent confidence interval:
     -0.006924869  0.002872771
    sample estimates:
    difference in location 
              -0.001885305 
    


Checking internal relatedness (IR)


```r
%%R
source('../Programs/R/x86_64-pc-linux-gnu-library/4.1/Rhh/R/ir.R')
```


```r
%%R
IRResult <- data.frame(ir(genoCall))
IRMerged <- merge(HLMerged, IRResult, by.x='ID', by.y='row.names')
```


```r
%%R
ggplot(IRMerged, aes(x=Evd.Of.Mult.Pat, y=ir.genoCall., fill=Evd.Of.Mult.Pat)) + 
    geom_boxplot() +
    theme_classic()
```


    
![png](output_198_0.png)
    



```r
%%R
plot(IRMerged$ir.genoCall. ~ IRMerged$hl.genoCall.)
```


    
![png](output_199_0.png)
    


Use allele frequencies from hard genotype calls, geno_call.mafs.gz, higher in mult. pat. litters ?


```python
# Poster dimensions: 42 x 36. Posters due Dec. 8. Finalized week before
# See if IR/HL are correlated with mother offspring relatedness
```

### Conception date between single paternity and multiple paternity groups


```r
%%R
sampleData <- read.delim("SampleData/sampleData.tab")
sampleData$Luther.Mother.ID. <- sapply(sampleData$Luther.Mother.ID., function(x) {return(as.numeric(substr(x,4,5)))} )
sampleData <- merge(sampleData, familyAssignments, by.x="Luther.Mother.ID.", by.y="row.names", all.x=TRUE)
sampleDataNoNA <- sampleData%>%filter(Evd.Of.Mult.Pat != 'Unknown')
sampleDataNoNA$Evd.Of.Mult.Pat <- sapply(sampleDataNoNA$Evd.Of.Mult.Pat, function(x) {if (x == "1"){return("Multiple")} else{return("Single")}})
```


```python
%R str(sampleData)
```

    'data.frame':	68 obs. of  19 variables:
     $ Luther.Mother.ID.            : num  1 2 3 4 5 6 7 8 9 10 ...
     $ Master.Deer.ID..             : int  2 3 4 5 8 9 10 13 16 29 ...
     $ Sex..M.F.                    : logi  FALSE FALSE FALSE FALSE FALSE FALSE ...
     $ Age.Estimate..previous.fall. : num  3.5 1.5 2.5 4.5 4.5 0.5 8.5 2.5 3.5 4.5 ...
     $ Tooth.taken.                 : chr  "" "x" "x" "x" ...
     $ K.Age                        : num  3.5 1.5 2.5 6.5 4.5 0.5 9.5 2.5 3.5 6.5 ...
     $ Date..mm.dd.yy.              : chr  "03/04/13" "04/08/13" "04/08/13" "04/09/13" ...
     $ County                       : chr  "Lucas" "Jones" "Dubuque" "Clinton" ...
     $ CO.Num                       : int  59 53 31 23 93 88 80 91 68 59 ...
     $ T                            : chr  "71N" "84N" "87N" "80N" ...
     $ R                            : chr  "23W" "2W" "1W" "3E" ...
     $ S                            : int  2 8 26 1 10 NA 1 20 25 32 ...
     $ Conception.Date..Julian.Days.: int  312 NA 328 319 328 NA 301 327 NA 318 ...
     $ Number.of.Fetuses            : int  2 2 2 2 3 1 2 2 2 2 ...
     $ Sex..1                       : chr  "M" "F" "M" "M" ...
     $ Sex..2                       : chr  "M" "M" "F" "M" ...
     $ Sex..3                       : chr  "" "" "" "" ...
     $ X                            : logi  NA NA NA NA NA NA ...
     $ Evd.Of.Mult.Pat              : chr  "Unknown" "1" "Unknown" "0" ...



```r
%%R
conceptionPlot <- ggplot(data=sampleDataNoNA, aes(x=Conception.Date..Julian.Days., fill=Evd.Of.Mult.Pat, color=Evd.Of.Mult.Pat)) +
    geom_histogram(position='identity',alpha=0.5,bins=12) +
    theme_classic() + 
    scale_fill_manual('Paternity',values=c("#44BB99", "#EE8866")) +
    scale_colour_manual('Paternity',values=c("#44BB99","#EE8866")) +
    theme(legend.key.width = unit(1, 'cm'), legend.title = element_text(size=20), legend.text = element_text(size=20),
         axis.title.x = element_text(size = 16), axis.text.y = element_text(size = 16), axis.title.y = element_text(size=16),
         legend.position = c(.95, .95),
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.margin = margin(6, 6, 6, 6)) + 
    xlab('Conception date') + 
    ylab('Frequency')
print(conceptionPlot)
conceptionPlot + scale_fill_manual(values=c("grey70", "grey10")) + scale_color_manual(values=c("grey70", "grey10"))
```

    R[write to console]: Scale for 'fill' is already present. Adding another scale for 'fill', which
    will replace the existing scale.
    
    R[write to console]: Scale for 'colour' is already present. Adding another scale for 'colour',
    which will replace the existing scale.
    



    
![png](output_205_1.png)
    



    
![png](output_205_2.png)
    



```r
%%R
qplot(sample=Conception.Date..Julian.Days., data=sampleDataNoNA, color=Evd.Of.Mult.Pat) + theme_classic() + stat_qq_line()
```


    
![png](output_206_0.png)
    



```r
%%R 
sampleDataNoNA.noOutliers <- sampleDataNoNA %>% filter(Conception.Date..Julian.Days. < 360)
```


```r
%%R
qplot(data=sampleDataNoNA.noOutliers, sample=Conception.Date..Julian.Days., color=Evd.Of.Mult.Pat) + theme_classic() + stat_qq_line()
```


    
![png](output_208_0.png)
    



```r
%%R
shapiro.test(as.numeric(unlist(sampleDataNoNA.noOutliers %>% filter(Evd.Of.Mult.Pat == 0) %>% select(Conception.Date..Julian.Days.))))
```

    
    	Shapiro-Wilk normality test
    
    data:  as.numeric(unlist(sampleDataNoNA.noOutliers %>% filter(Evd.Of.Mult.Pat == 0) %>% select(Conception.Date..Julian.Days.)))
    W = 0.97404, p-value = 0.4166
    



```r
%%R
shapiro.test(as.numeric(unlist(sampleDataNoNA.noOutliers %>% filter(Evd.Of.Mult.Pat == 1) %>% select(Conception.Date..Julian.Days.))))
```

    
    	Shapiro-Wilk normality test
    
    data:  as.numeric(unlist(sampleDataNoNA.noOutliers %>% filter(Evd.Of.Mult.Pat == 1) %>% select(Conception.Date..Julian.Days.)))
    W = 0.94617, p-value = 0.6481
    



```r
%%R
var.test(Conception.Date..Julian.Days. ~ Evd.Of.Mult.Pat, sampleDataNoNA.noOutliers, alternative = 'two.sided')
```

    
    	F test to compare two variances
    
    data:  Conception.Date..Julian.Days. by Evd.Of.Mult.Pat
    F = 5.8196, num df = 43, denom df = 8, p-value = 0.01285
    alternative hypothesis: true ratio of variances is not equal to 1
    95 percent confidence interval:
      1.520158 14.571795
    sample estimates:
    ratio of variances 
              5.819607 
    

