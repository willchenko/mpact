

import os
import csv
import re
import numpy as np
import matplotlib.pyplot as plt
os.chdir(r"C:\Users\Owner\Documents\PhD_stuff\research")
from mpact import parseRxnString


os.chdir(r"C:\Users\Owner\Documents\PhD_stuff\research\ModCell2-master\ModCell2-master\problems\ecoli-gem-total")
csv_file = 'metadata_table.csv'
[prods, nodes, nrs, ov_strings] = read_metadata_table(csv_file)

cofactor = 'atp_c'
req = get_cofactor_requirements(ov_strings,cofactor)
labels = list(np.unique(req))
count = []
for label in labels:
    count.append(len([i for i, x in enumerate(req) if x == label]))


plt.bar([str(i) for i in labels],count)
plt.xlabel("Moles of ATP per module") 
plt.ylabel("No. Production Modules") 
#plt.title("Students enrolled in different courses") 
plt.show()

def read_metadata_table(csv_file):
    prods = []
    nodes = []
    nrs = []
    ov_strings = []
    with open(csv_file,newline='') as csvfile:
            reader = csv.reader(csvfile)
            for row in reader:
                if row[0]!= 'id':
                    prods.append(row[0])
                    nodes.append(row[2])
                    nrs.append(row[3])
                    ov_strings.append(row[4])
    return prods, nodes, nrs, ov_strings

def get_cofactor_requirements(ov_strings,cofactor):
    req = []
    for string in ov_strings:
        [mets,stoich,rev] = parseRxnString(string)
        if cofactor in mets:
            req.append(stoich[mets.index(cofactor)])
        else:
            req.append(0)
    return req
        
