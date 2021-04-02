#!/usr/bin/env bash

cd /home/arc78/scratch60/covid-placenta/data/00_fastq

mv AL1/AL-1_HHT_S1_L004_I1_001.fastq.gz COVID_4__decidua__HHT_S1_L004_I1_001.fastq.gz
mv AL1/AL-1_HHT_S1_L004_R1_001.fastq.gz COVID_4__decidua__HHT_S1_L004_R1_001.fastq.gz
mv AL1/AL-1_HHT_S1_L004_R2_001.fastq.gz COVID_4__decidua__HHT_S1_L004_R2_001.fastq.gz
mv AL2/AL-2_HHT_S2_L004_I1_001.fastq.gz COVID_4__villi__HHT_S2_L004_I1_001.fastq.gz
mv AL2/AL-2_HHT_S2_L004_R1_001.fastq.gz COVID_4__villi__HHT_S2_L004_R1_001.fastq.gz
mv AL2/AL-2_HHT_S2_L004_R2_001.fastq.gz COVID_4__villi__HHT_S2_L004_R2_001.fastq.gz
mv AL3/AL-3_HHT_S3_L004_I1_001.fastq.gz COVID_5__decidua__HHT_S3_L004_I1_001.fastq.gz
mv AL3/AL-3_HHT_S3_L004_R1_001.fastq.gz COVID_5__decidua__HHT_S3_L004_R1_001.fastq.gz
mv AL3/AL-3_HHT_S3_L004_R2_001.fastq.gz COVID_5__decidua__HHT_S3_L004_R2_001.fastq.gz
mv AL4/AL-4_HHT_S4_L004_I1_001.fastq.gz COVID_5__villi__HHT_S4_L004_I1_001.fastq.gz
mv AL4/AL-4_HHT_S4_L004_R1_001.fastq.gz COVID_5__villi__HHT_S4_L004_R1_001.fastq.gz
mv AL4/AL-4_HHT_S4_L004_R2_001.fastq.gz COVID_5__villi__HHT_S4_L004_R2_001.fastq.gz
mv AL5/AL5_HCR_S11_L002_I1_001.fastq.gz CNTRL_1__decidua__HCR_S11_L002_I1_001.fastq.gz
mv AL5/AL5_HCR_S11_L002_R1_001.fastq.gz CNTRL_1__decidua__HCR_S11_L002_R1_001.fastq.gz
mv AL5/AL5_HCR_S11_L002_R2_001.fastq.gz CNTRL_1__decidua__HCR_S11_L002_R2_001.fastq.gz
mv AL6/AL6_HCR_S12_L002_I1_001.fastq.gz CNTRL_1__villi__HCR_S12_L002_I1_001.fastq.gz
mv AL6/AL6_HCR_S12_L002_R1_001.fastq.gz CNTRL_1__villi__HCR_S12_L002_R1_001.fastq.gz
mv AL6/AL6_HCR_S12_L002_R2_001.fastq.gz CNTRL_1__villi__HCR_S12_L002_R2_001.fastq.gz
mv AL8/AL8_HCR_S5_L003_I1_001.fastq.gz CNTRL_2__villi__HCR_S5_L003_I1_001.fastq.gz
mv AL8/AL8_HCR_S5_L003_R1_001.fastq.gz CNTRL_2__villi__HCR_S5_L003_R1_001.fastq.gz
mv AL8/AL8_HCR_S5_L003_R2_001.fastq.gz CNTRL_2__villi__HCR_S5_L003_R2_001.fastq.gz
mv AL9/AL9_HCR_S6_L003_I1_001.fastq.gz CNTRL_3__decidua__HCR_S6_L003_I1_001.fastq.gz
mv AL9/AL9_HCR_S6_L003_R1_001.fastq.gz CNTRL_3__decidua__HCR_S6_L003_R1_001.fastq.gz
mv AL9/AL9_HCR_S6_L003_R2_001.fastq.gz CNTRL_3__decidua__HCR_S6_L003_R2_001.fastq.gz
mv AL10/AL10_HCR_S7_L003_I1_001.fastq.gz CNTRL_3__villi__HCR_S7_L003_I1_001.fastq.gz
mv AL10/AL10_HCR_S7_L003_R1_001.fastq.gz CNTRL_3__villi__HCR_S7_L003_R1_001.fastq.gz
mv AL10/AL10_HCR_S7_L003_R2_001.fastq.gz CNTRL_3__villi__HCR_S7_L003_R2_001.fastq.gz

rmdir AL1
rmdir AL2
rmdir AL3
rmdir AL4
rmdir AL5
rmdir AL6
rmdir AL7
rmdir AL8
rmdir AL9
rmdir AL10

