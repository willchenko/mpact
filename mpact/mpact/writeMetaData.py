# write metadata to file
import csv
import os

#writeMetaData(metadata,'metadata_table.csv',folder_path)

def writeMetaData(metadata,filename,folder_path):
    os.chdir(folder_path)
    with open(filename,'a',newline='') as csvfile: 
    # creating a csv writer object 
        csvwriter = csv.writer(csvfile)
        csvwriter.writerow(['id','target','precursor','number of rxns','overall_string']) 
        for path in metadata['pathways']:
            csvwriter.writerow([path,metadata[path]['target'],metadata[path]['precursor'],metadata[path]['number_of_rxns'],metadata[path]['ov_string']]) 
    return 

