#!/usr/bin/env python
'''
    Python For Extracting metadata
    1. Extract acronyms
    2. Extract Manufacturers and their models
    3. Extract MRI sequences based on TE and TR
    
'''
################################################################################
# Libraries
import numpy as np
import pandas as pd
import os
import json

################################################################################
''' 1. Exctract Acronyms '''
################################################################################
# Import TSV files

Dir_to_TSVs = 'C:/Users/agolbaf/Desktop/University_of_Plymouth_Biomedical_Research_Laboratories/Plymouth_Data/Files_NIfTI/TSV_Files'
TSV_Subs = (next(os.walk(Dir_to_TSVs))[1])
TSV_Subs.sort()

List_Total = []
for ID in TSV_Subs:
    Path = os.path.join(Dir_to_TSVs,ID)
    TSV_Files = os.listdir(Path)
    TSV_Files.sort()

    for file in TSV_Files:
        Index = file.find('Edited')
        if Index == 0:
            Path_to_File = os.path.join(Path,file)
            Data = pd.read_csv(Path_to_File , sep='\t')
            
            # Drop First Column
            Data = Data[Data.columns[1:]]
            Seq_Column = Data.loc[:,'sequence_name']
            
            for i in range(len(Seq_Column)):
                cell = Seq_Column[i]
                
                if cell not in List_Total and pd.isnull(cell) == False:
                    List_Total.append(cell)    

df = pd.DataFrame(List_Total, columns=["sequence_name"])
df.to_csv('Acronyms.csv', index=False)

################################################################################
'''2. Extract Manufacturers and their models'''
################################################################################
# Import JSON files
# Manufacturer and Manufacturermodel names

Dir_to_Nifti = 'C:/Users/agolbaf/Desktop/University_of_Plymouth_Biomedical_Research_Laboratories/Plymouth_Data/Files_NIfTI/Nifti'
Subs = (next(os.walk(Dir_to_TSVs))[1])
Subs.sort()

List_Total_M_M = []
for ID_Sub in Subs:
    Path_to_Sess = os.path.join(Dir_to_Nifti,ID_Sub)
    Sess = os.listdir(Path_to_Sess)
    Sess.sort()
    for ID_Ses in Sess:
        Path_to_Jsons = os.path.join(Path_to_Sess, ID_Ses + "/anat/")
        Json_Files = os.listdir(Path_to_Jsons)
        Json_Files.sort()
        for file in Json_Files:
            if file.endswith(".json"):
                # Psth to Json files
                Path_to_File = os.path.join(Path_to_Jsons,file)
                
                # Opening JSON file
                # Check file size to avaid empty files
                Size = os.path.getsize(Path_to_File)
                
                if Size > 10:
                    f = open(Path_to_File, 'r', encoding='utf-8')
                    # returns JSON object as a dictionary
                    Data = json.load(f)
                    Manufacturer = ""
                    Manufacturermodel = ""
                    if "Manufacturer" in Data:
                        Manufacturer = Data["Manufacturer"]
                    if "ManufacturersModelName" in Data:
                        Manufacturermodel = Data["ManufacturersModelName"]
                    
                    Raw = [Manufacturer,Manufacturermodel]
                    
                    if Raw not in List_Total_M_M:
                        List_Total_M_M.append(Raw) 
                        
df = pd.DataFrame(List_Total_M_M, columns=["Manufacturer","Manufacturermodel"])
df.to_csv('Manufacturer.csv', index=False)              

################################################################################
'''3. Extract MRI sequences based on TE and TR'''
################################################################################
# Time dimension: ms
# There is no any kinds of overlaps among ranges here so in the code we do not need to apply barriers
# T1 Charactristics
def T1_Check(TR, TE):
    if TR < 800 and TE < 30:
        return("T1")
    else:
        return("T1_Unknown")
    
# T2 Charactristics   
def T2_Check(TR, TE):
    if 3000 >TR > 2000 and TE > 80:
        return("T2")
    else:
        return("T2_Unknown")

# Flair Charactristics
def Flair_Check(TR, TE):
    if TR > 3000 and TE > 80:
        return("Flair")
    else:
        return("Flair_Unknown")
    

# Read TR and TE from dataset
Dir_to_TSVs = 'C:/Users/agolbaf/Desktop/University_of_Plymouth_Biomedical_Research_Laboratories/Plymouth_Data/Files_NIfTI/TSV_Files'
TSV_Subs = (next(os.walk(Dir_to_TSVs))[1])
TSV_Subs.sort()

Sequence_List = []

Sub = 0

for ID in TSV_Subs:
    
    Path = os.path.join(Dir_to_TSVs,ID)
    TSV_Files = os.listdir(Path)
    TSV_Files.sort()
    Sub = Sub + 1
    Ses = 0
    for file in TSV_Files:
        
        Index = file.find('Edited')
        if Index == 0:
            Ses = Ses + 1 
            Path_to_File = os.path.join(Path,file)
            Data = pd.read_csv(Path_to_File , sep='\t')
            
            # Drop First Column
            Data = Data[Data.columns[1:]]
            
            Series_Column = Data.loc[:,'series_id']
            TR_Column = Data.loc[:,'TR']
            TE_Column = Data.loc[:,'TE']
            
            for i in range(len(Series_Column)):
                Series_id = Series_Column[i]
                TR = float(TR_Column[i]) * 1000 # Turn to ms
                TE = float(TE_Column[i]) * 1000 # Turn to ms
                
                # Sequnece Check
                T1 = T1_Check(TR,TE)
                T2 = T2_Check(TR, TE)
                Flair = Flair_Check(TR, TE)
                
                Raw = [Sub, Ses, Series_id, T1, T2, Flair]
                Sequence_List.append(Raw)

# Statistical analysis
T1_Statistices = 0
T2_Statistices = 0
Flair_Statistices = 0

for i in range (len(Sequence_List)):
    T1 = Sequence_List[i][3]
    T2 = Sequence_List[i][4]
    Flair = Sequence_List[i][5]
    
    if T1 == "T1":
        T1_Statistices+=1
        
    if T2 == "T2":
        T2_Statistices+=1
        
    if Flair == "Flair":
        Flair_Statistices+=1
    
        
Number_of_T1 = T1_Statistices
Number_of_T2 = T2_Statistices
Number_of_Flair = Flair_Statistices
    
Raw_01 = ["", "", "", "Number_of_T1", "Number_of_T2", "Number_of_Flair"]    
Raw_02 = ["", "", "", Number_of_T1, Number_of_T2, Number_of_Flair]  

Sequence_List.append(Raw_01)
Sequence_List.append(Raw_02)
              
df = pd.DataFrame(Sequence_List, columns=["Subject","Session", "Series_id", "T1_Check", "T2_Check", "Flair_Check"])
df.to_csv('Sequence_List.csv', index=False)   
    
    
    
            
    
    
    
    


