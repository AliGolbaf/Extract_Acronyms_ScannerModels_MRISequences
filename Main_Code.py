#!/usr/bin/env python
'''
    Python For Extracting metadata
    1. Extract acronyms
    2. Extract Manufacturers and their models
    3. Extract MRI sequences based on TE and TR and acronyms extracted from: 
        https://github.com/rordenlab/dcm2niix/blob/master/console/nii_dicom_batch.cpp (lines: 1176-1197)
    
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
'''3. 
    3.1 Extract MRI sequences based on TE and TR
    3.2 Extract MRI sequences based on acronyms extracted from: 
        https://github.com/rordenlab/dcm2niix/blob/master/console/nii_dicom_batch.cpp (lines: 1176-1197)
    '''
################################################################################
# Time dimension: ms
# There is no any kinds of overlaps among ranges here so in the code we do not need to apply barriers
# T1 Charactristics
def T1_Check(TR, TE, Acronym):
    
    var = []
    # tr and te
    if TR < 800 and TE < 30:
        var.append("T1")
    else:
        var.append("T1_Unknown")
    
    # Acronym
    if Acronym in ["tfl3d", "tfl_me3d5_16ns"]:
        var.append("T1")
        
    else:
        var.append("T1_Unknown")
    
    return(var)
        
# T2 Charactristics   
def T2_Check(TR, TE, Acronym):
    
    var = []
    # tr and te
    if 3000 >TR > 2000 and TE > 80:
        var.append("T2")
    else:
        var.append("T2_Unknown")
        
    # Acronym
    if Acronym in ["spc3d", "tse2d", "tse3d"]:
        var.append("T2")
        
    else:
        var.append("T2_Unknown")
    
    return(var)

# Flair Charactristics
def Flair_Check(TR, TE, Acronym):
    
    var = []
    # tr and te
    if TR > 3000 and TE > 80:
        
        var.append("Flair")
    else:
        var.append("Flair_Unknown")
    
    # Acronym
    if Acronym in ["tir2d", "spcir"]:
        var.append("Flair")
        
    else:
        var.append("Flair_Unknown")
    return(var)
    
def DWI_Check(Acronym):
    
    var = []
    # Acronym
    if Acronym in ["ep_b", "epse2d" ]:
        var.append("DWI")
    else:
        var.append("DWI_Unknown")
    
    return var
    
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
            Acronym_Column = Data.loc[:,'sequence_name']
            for i in range(len(Series_Column)):
                Series_id = Series_Column[i]
                TR = float(TR_Column[i]) * 1000 # Turn to ms
                TE = float(TE_Column[i]) # Currently is in ms
                
                Acronym = (Acronym_Column[i])
                if pd.isnull(Acronym):
                    acronym = 'Nan'
                else:  
                    acronym = Acronym.replace('*', '')
                                    
                # Sequnece Check
                T1 = T1_Check(TR,TE, acronym)
                T2 = T2_Check(TR, TE, acronym)
                Flair = Flair_Check(TR, TE, acronym)
                DWI = DWI_Check(acronym)
                
                
                Raw = [Sub, Ses, Series_id, T1[0],T1[1], T2[0],T2[1], Flair[0], Flair[1], DWI[0]]
                Sequence_List.append(Raw)

# # Statistical analysis
# T1_Statistices    = 0
# T2_Statistices    = 0
# Flair_Statistices = 0
# DWI_Statistics    = 0

# for i in range (len(Sequence_List)):
#     T1 = Sequence_List[i][3]
#     T2 = Sequence_List[i][4]
#     Flair = Sequence_List[i][5]
    
#     if T1 == "T1":
#         T1_Statistices+=1
        
#     if T2 == "T2":
#         T2_Statistices+=1
        
#     if Flair == "Flair":
#         Flair_Statistices+=1
    
        
# Number_of_T1 = T1_Statistices
# Number_of_T2 = T2_Statistices
# Number_of_Flair = Flair_Statistices
    
# Raw_01 = ["", "", "", "Number_of_T1", "Number_of_T2", "Number_of_Flair"]    
# Raw_02 = ["", "", "", Number_of_T1, Number_of_T2, Number_of_Flair]  

# Sequence_List.append(Raw_01)
# Sequence_List.append(Raw_02)
              
df = pd.DataFrame(Sequence_List, columns=["Subject","Session", "Series_id","T1_Check_TR_TE",
                                          "T1_Check_Acronym", "T2_Check_TR_TE", "T2_Check_Acronym", 
                                          "Flair_Check_TR_TE", "Flair_Check_Acronym", "DWI_Check_Acronym"])
df.to_csv('Sequence_List.csv', index=False)   
    
    
    
            
    
    
    
    


