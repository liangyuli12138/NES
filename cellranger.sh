


#!/bin/bash
#PBS -l nodes=1:ppn=10
#PBS -d /public/workspace/lily/PS/Response/
#PBS -N CWH_1
#PBS -o /public/workspace/lily/log_qusb/CWH_1.out
#PBS -e /public/workspace/lily/log_qusb/CWH_1.err

module load cellranger-3.0.2
echo "start ..."
cellranger count --id=CWH_1 \
    --sample=CWH_1-1 \
    --fastqs=/public/workspace/DATA/sorted_by_date/2022-05-02/caoyachan0426036/X101SC22041880-Z01-J001-B1-1_10X_release_20220422/Rawdata/CWH_1 \
    --transcriptome=/public/workspace/lily/REF/INDEX-hg19/refdata-cellranger-hg19-1.2.0 \
    --chemistry=auto \
    --localcores=10 --localmem=60 --nosecondary



echo "done ..."




#!/bin/bash
#PBS -l nodes=1:ppn=10
#PBS -d /public/workspace/lily/PS/Response/
#PBS -N CWH_1
#PBS -o /public/workspace/lily/log_qusb/CWH_1.out
#PBS -e /public/workspace/lily/log_qusb/CWH_1.err

cd /public/workspace/lily/PS/Response/

echo "start ..."
/public/workspace/lily/software/cellranger-3.0.2/cellranger count --id=CWH_1 \
    --sample=CWH_1-1 \
    --fastqs=/public/workspace/lily/gliomas/caoyachan0426036/X101SC22041880-Z01-J001-B1-1_10X_release_20220422/Rawdata/CWH_1 \
    --transcriptome=/public/workspace/lily/REF/INDEX-hg19/refdata-cellranger-hg19-1.2.0 \
    --chemistry=auto \
    --localcores=10 --localmem=60 --nosecondary


echo "done ..."



# for i in `ls -d */ | sed -e 's/\/$//' | grep -v "res" `
# do 
# 	rm -r ${i}
# done




#!/bin/bash
#PBS -l nodes=1:ppn=10
#PBS -d /public/workspace/lily/PS/Response/
#PBS -N CWH_2
#PBS -o /public/workspace/lily/log_qusb/CWH_2.out
#PBS -e /public/workspace/lily/log_qusb/CWH_2.err


module load cellranger-3.0.2
echo "start ..."
cellranger count --id=CWH_2 \
    --sample=CWH_2-1 \
    --fastqs=/public/workspace/DATA/sorted_by_date/2022-05-02/caoyachan0426036/X101SC22041880-Z01-J001-B1-1_10X_release_20220422/Rawdata/CWH_2 \
    --transcriptome=/public/workspace/lily/REF/INDEX-hg19/refdata-cellranger-hg19-1.2.0 \
    --chemistry=auto \
    --localcores=10 --localmem=60 --nosecondary



echo "done ..."


echo "start ..."
/public/workspace/lily/software/cellranger-3.0.2/cellranger count --id=CWH_2 \
    --sample=CWH_2-1 \
    --fastqs=/public/workspace/lily/gliomas/caoyachan0426036/X101SC22041880-Z01-J001-B1-1_10X_release_20220422/Rawdata/CWH_2 \
    --transcriptome=/public/workspace/lily/REF/INDEX-hg19/refdata-cellranger-hg19-1.2.0 \
    --chemistry=auto \
    --localcores=10 --localmem=60 --nosecondary


echo "done ..."






#!/bin/bash
#PBS -l nodes=1:ppn=10
#PBS -d /public/workspace/lily/PS/Response/
#PBS -N CWH_3
#PBS -o /public/workspace/lily/log_qusb/CWH_3.out
#PBS -e /public/workspace/lily/log_qusb/CWH_3.err


module load cellranger-3.0.2

echo "start ..."
cellranger count --id=CWH_3 \
    --sample=CWH_3-1 \
    --fastqs=/public/workspace/DATA/sorted_by_date/2022-05-02/caoyachan0426036/X101SC22041880-Z01-J001-B1-1_10X_release_20220422/Rawdata/CWH_3 \
    --transcriptome=/public/workspace/lily/REF/INDEX-hg19/refdata-cellranger-hg19-1.2.0 \
    --chemistry=auto \
    --localcores=10 --localmem=60 --nosecondary



echo "done ..."


echo "start ..."
/public/workspace/lily/software/cellranger-3.0.2/cellranger count --id=CWH_3 \
    --sample=CWH_3-1 \
    --fastqs=/public/workspace/lily/gliomas/caoyachan0426036/X101SC22041880-Z01-J001-B1-1_10X_release_20220422/Rawdata/CWH_3 \
    --transcriptome=/public/workspace/lily/REF/INDEX-hg19/refdata-cellranger-hg19-1.2.0 \
    --chemistry=auto \
    --localcores=10 --localmem=60 --nosecondary


echo "done ..."







