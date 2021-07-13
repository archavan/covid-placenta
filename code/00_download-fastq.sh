#!/usr/bin/env bash

# download 10x data: fastq files ==============================================

# samples 1-4
for sample in AL-1 AL-2 AL-3 AL-4
do
	cd /home/arc78/scratch60/covid-placenta/data/00_fastq
	wget -r -nH --cut-dirs=3 -np --accept ".gz" "http://fcb.ycga.yale.edu:3010/SArOwtq7STkIpVS9ct8hDFjIusWPK/042020/"$sample"_HHT_fqs/"
done

# samples 5, 6, 8. we did not use AL7
for sample in AL5 AL6 AL8
do
	cd /home/arc78/scratch60/covid-placenta/data/00_fastq
	wget -r -nH --cut-dirs=3 -np --accept ".gz" "http://fcb.ycga.yale.edu:3010/fyiZHAf0T61ABnyAN7GnonRlVkOa9/060120/"$sample"_HCR_fqs/"
done

# samples 8-10
for sample in AL8 AL9 AL10 # we did not use AL7
do
	cd /home/arc78/scratch60/covid-placenta/data/00_fastq
	wget -r -nH --cut-dirs=3 -np --accept ".gz" "http://fcb.ycga.yale.edu:3010/kSQvzi5anisrelWEbQXa6cCfEa7oN/060520/"$sample"_HCR_fqs/"
done


# end =========================================================================
	
	
