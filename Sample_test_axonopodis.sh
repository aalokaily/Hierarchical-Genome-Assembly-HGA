#!/bin/bash

curr=`pwd`


script_dir=$(dirname "${0}")
cd $script_dir
script_dir=`pwd`

cd $curr
mkdir Axonopodis_test
cd Axonopodis_test
curr=`pwd`


echo -----------Download the data from GAGE-B-----------  | tee -a $curr/run.log

wget http://ccb.jhu.edu/gage_b/datasets/X_axonopodis_HiSeq.tar.gz  | tee -a $curr/run.log

echo -----------unzip the dataset----------- | tee -a $curr/run.log

tar -xvf X_axonopodis_HiSeq.tar.gz | tee -a $curr/run.log

echo -----------Download interleave_fastq.py to interlace paired fastq files----------- | tee -a $curr/run.log

wget https://gist.githubusercontent.com/ngcrawford/2232505/raw/338758f7fcca8ad24340730b96ba645c46fa1b0e/interleave_fastq.py | tee -a $curr/run.log

echo -----------Interlace paired fastq files-----------

python interleave_fastq.py $curr/trimmed/insert_400_1__cov250x.fastq $curr/trimmed/insert_400_2__cov250x.fastq data.fastq | tee -a $curr/run.log

rm -rf ./trimmed | tee -a $curr/run.log
rm -rf ./raw | tee -a $curr/run.log
rm X_axonopodis_HiSeq.tar.gz | tee -a $curr/run.log
rm interleave_fastq.py | tee -a $curr/run.log


echo -----------install Velvet----------- | tee -a $curr/run.log
wget https://www.ebi.ac.uk/~zerbino/velvet/velvet_1.2.10.tgz | tee -a $curr/run.log
tar zxvf velvet_1.2.10.tgz | tee -a $curr/run.log
cd velvet_1.2.10 
make 'MAXKMERLENGTH=111' 'LONGSEQUENCES=1' | tee -a $curr/run.log

cd $curr 
rm velvet_1.2.10.tgz | tee -a $curr/run.log


echo -----------install SPAdes----------- | tee -a $curr/run.log

wget http://spades.bioinf.spbau.ru/release3.0.0/SPAdes-3.0.0.tar.gz | tee -a $curr/run.log
tar -xzf SPAdes-3.0.0.tar.gz | tee -a $curr/run.log
cd SPAdes-3.0.0 
./spades_compile.sh | tee -a $curr/run.log

cd $curr 
rm SPAdes-3.0.0.tar.gz | tee -a $curr/run.log

echo -----------Run HGA----------- | tee -a $curr/run.log
python $script_dir/HGA.py -velvet $curr/velvet_1.2.10 -spades $curr/SPAdes-3.0.0/bin -PA velvet  -P12  $curr/data.fastq  -R12 $curr/data.fastq -ins 400 -std 40 -P 16 -Pkmer 31 -Rkmer 51 -t 1 -out  ./HGA | tee -a $curr/run.log

cd $curr 
rm -rf $curr/SPAdes-3.0.0 | tee -a $curr/run.log
rm -rf $curr/velvet_1.2.10  | tee -a $curr/run.log
rm data.fastq | tee -a $curr/run.log


