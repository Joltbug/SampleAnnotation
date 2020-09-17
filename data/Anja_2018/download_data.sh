wget https://ndownloader.figshare.com/files/12853469 -O uATAC.PBMC.bed
wget https://ndownloader.figshare.com/files/12853382 -O uATAC.CD8.bed
wget https://ndownloader.figshare.com/files/12853376 -O uATAC.T.bed
wget https://ndownloader.figshare.com/files/12853373 -O uATAC.CD4.bed
wget https://ndownloader.figshare.com/files/12853220 -O uATAC.B.bed
wget https://ndownloader.figshare.com/files/12853442 -O uATAC.M.bed


for f in ./uATAC*.bed
do
awk {'printf ("%s\t%s\t%s\n", $1, $2, $3)'} $f > tmp.bed
sort -k1,1 -k2,2n tmp.bed > $f
rm -f tmp.bed
done
gzip *.bed