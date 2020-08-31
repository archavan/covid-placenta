# download 10x data: cellranger output

# samples 1-4
for sample in AL-1 AL-2 AL-3 AL-4
do
	cd ~/research/covid-placenta/data/01_cellranger
	mkdir $sample
	cd $sample
	wget -r -nH --cut-dirs=3 -np --accept ".gz" "http://fcb.ycga.yale.edu:3010/SArOwtq7STkIpVS9ct8hDFjIusWPK/042020/"$sample"_HHT_cellranger/filtered_feature_bc_matrix/"
	wget "http://fcb.ycga.yale.edu:3010/SArOwtq7STkIpVS9ct8hDFjIusWPK/042020/"$sample"_HHT_cellranger/metrics_summary.csv"
	wget "http://fcb.ycga.yale.edu:3010/SArOwtq7STkIpVS9ct8hDFjIusWPK/042020/"$sample"_HHT_cellranger/web_summary.html"
done

# samples 5-8
for sample in AL5 AL6 AL7 AL8
do
	cd ~/research/covid-placenta/data/01_cellranger
	mkdir $sample
	cd $sample
	wget -r -nH --cut-dirs=3 -np --accept ".gz" "http://fcb.ycga.yale.edu:3010/fyiZHAf0T61ABnyAN7GnonRlVkOa9/060120/"$sample"_HCR_cellranger/filtered_feature_bc_matrix/"
	wget "http://fcb.ycga.yale.edu:3010/fyiZHAf0T61ABnyAN7GnonRlVkOa9/060120/"$sample"_HCR_cellranger/metrics_summary.csv"
	wget "http://fcb.ycga.yale.edu:3010/fyiZHAf0T61ABnyAN7GnonRlVkOa9/060120/"$sample"_HCR_cellranger/web_summary.html"
done

## rename directories to make their names consistent.
cd ~/research/covid-placenta/data/01_cellranger

mv AL-1 AL1
mv AL-2 AL2
mv AL-3 AL3
mv AL-4 AL4

