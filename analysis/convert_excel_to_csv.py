#!/usr/bin/env python3

import sys
import subprocess

packages = ['openpyxl', 'pandas']
for package in packages:
    try:
        __import__(package)
    except ImportError:
        print(f"Installing {package}...")
        subprocess.check_call([sys.executable, '-m', 'pip', 'install', package])

import pandas as pd

print("Loading Excel file...")
df = pd.read_excel('GSE131005_RNASeq.xlsx')

print(f"Data shape: {df.shape}")
print(f"Columns: {list(df.columns)}")

print("\nFirst few rows:")
print(df.head())

print("\nConverting to CSV...")
df.to_csv('GSE131005_RNASeq.csv', index=False)
print("Done! File saved as GSE131005_RNASeq.csv")