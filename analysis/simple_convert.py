#!/usr/bin/env python3
import sys
sys.path.append('./venv/lib/python3.12/site-packages')

import openpyxl
import csv

print("Loading Excel file with openpyxl...")
wb = openpyxl.load_workbook('GSE131005_RNASeq.xlsx', read_only=True, data_only=True)
ws = wb.active

print("Converting to CSV...")
with open('GSE131005_RNASeq.csv', 'w', newline='') as f:
    writer = csv.writer(f)
    row_count = 0
    for row in ws.iter_rows(values_only=True):
        writer.writerow(row)
        row_count += 1
        if row_count % 1000 == 0:
            print(f"Processed {row_count} rows...")
    
print(f"Done! Total rows: {row_count}")
wb.close()