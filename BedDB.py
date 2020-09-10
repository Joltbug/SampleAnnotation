## data
## - folder1
## - folder2
## - folder3
## - meta_data: (list of folders and some metadata of each folder)

## in each folder
## -script to download data
## -bed1, bed2, bed3
## -background bed
## -info.yaml

## info.yaml
## key: background, value:background file
## key: cell_type, value: {celltype: file to a celltype} 
## key: source GSEXXXX

## BedDB.py
## def list_dbs():
### read meta_data or os.listdir('./data')

## def read_db("foldername"):
### read beds,  info.yaml in foldername
### return {ct: ct peaks}, background peaks

## class BedDBLoader():