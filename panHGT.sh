#echo "Welcome to the PanGenome Pipeline!"

#bash download.sh
#python3 fnaRename.py
#bash prokka.sh
#python3 gffFilter.py
#bash panaroo.sh

#bash cloudExtract.sh
#bash hgtTests.sh
bash calcScoreHGT.sh
python3 categorization.py
mkdir -p ./HGT_evidence/charts/bar ./HGT_evidence/charts/pie

python3 visualizations.py
python3 hgtSummary.py

echo "Pipeline complete!"
echo "Please check the output files in the current directory."