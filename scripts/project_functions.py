import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pandas_profiling import ProfileReport

# df = pd.read_csv("C:/Users/nitchakan/Desktop/DATA301/project-group45-project/data/raw/GSE113690_Autism_16S_rRNA_OTU_assignment_and_abundance.csv")

def load_and_process(path):
    
    def TD_ASD(row):
        row['TD_ASDdiff'] = int(row['TDmean']) - int(row['ASDmean'])
        return row

# df = (
#     pd.read_csv("../../data/raw/GSE113690_Autism_16S_rRNA_OTU_assignment_and_abundance.csv")
#     )

    df = (
        pd.read_csv("/Users/nitchakan/Desktop/DATA301/project-group45-project/data/raw/GSE113690_Autism_16S_rRNA_OTU_assignment_and_abundance.csv")
        )

    df1 = (
        df.assign(AllAvg=df.iloc[:,2:254].mean(axis='columns'))
        .drop(columns = ['OTU'])
        .query('AllAvg > 1')
        .assign(TDmean=df.loc[:,['A1', 'A10', 'A100', 'A101', 'A102', 'A104', 'A105', 'A106', 'A108', 'A109', 'A11', 'A110', 'A111', 
                                 'A112', 'A113', 'A114', 'A115', 'A116', 'A117', 'A118', 'A119', 'A12', 'A120', 'A121', 'A122', 'A123', 
                                 'A124', 'A125', 'A126', 'A128', 'A129', 'A130', 'A131', 'A132', 'A133', 'A135', 'A136', 'A137', 'A139', 
                                 'A14', 'A140', 'A141', 'A142', 'A143', 'A144', 'A145', 'A146', 'A147', 'A148', 'A149', 'A150', 'A151', 
                                 'A152', 'A153', 'A154', 'A155','A156', 'A157', 'A158', 'A159', 'A16', 'A160', 'A161', 'A162', 'A163', 'A164', 
                                 'A165', 'A166', 'A167', 'A168', 'A169', 'A17', 'A18', 'A19', 'A2', 'A20', 'A21', 'A22', 'A23', 'A25', 'A26', 
                                 'A27', 'A28', 'A29', 'A3', 'A31', 'A32', 'A33', 'A34', 'A35', 'A40', 'A43', 'A45', 'A46', 'A47', 'A49', 'A5', 
                                 'A50', 'A51', 'A52', 'A53', 'A54', 'A55', 'A57', 'A58', 'A59', 'A6', 'A60', 'A61', 'A62', 'A63', 'A65', 'A66', 
                                 'A67', 'A68', 'A69', 'A7', 'A70', 'A71', 'A72', 'A73', 'A75', 'A76', 'A77', 'A78', 'A8', 'A80', 'A82', 'A83', 
                                 'A84', 'A85', 'A86', 'A87', 'A88', 'A89', 'A9', 'A90', 'A91', 'A92', 'A93', 'A95', 'A97', 'A99'
                                ]].mean(axis='columns'))
        .assign(ASDmean=df.loc[:, [ 'B1', 'B10', 'B100', 'B101', 'B103', 'B104', 'B105', 'B106', 'B107', 'B108', 'B109', 'B11', 'B111', 'B112', 
                                   'B113', 'B114', 'B115', 'B116', 'B117', 'B119', 'B12', 'B120', 'B122', 'B123', 'B124', 'B125', 'B126', 'B127', 
                                   'B128', 'B129', 'B13', 'B131', 'B132', 'B133', 'B134', 'B135', 'B136', 'B137', 'B138', 'B139', 'B14', 'B141', 
                                   'B142', 'B143', 'B144', 'B145', 'B147', 'B148', 'B149', 'B15', 'B150', 'B151', 'B152', 'B153', 'B154', 'B155', 
                                   'B156', 'B157', 'B158', 'B159', 'B16', 'B160', 'B161', 'B162', 'B163', 'B164', 'B165', 'B166', 'B17', 'B18', 
                                   'B19', 'B2', 'B20', 'B21', 'B22', 'B23', 'B24', 'B25', 'B26', 'B27', 'B28', 'B29', 'B3', 'B30', 'B31', 'B34', 
                                   'B35', 'B36', 'B37', 'B39', 'B4', 'B40', 'B42', 'B44', 'B45', 'B46', 'B48', 'B49', 'B5', 'B50', 'B51', 'B52', 
                                   'B54', 'B55', 'B56', 'B57', 'B58', 'B59', 'B6', 'B60', 'B61']].mean(axis = 'columns'))
        .apply(lambda x: TD_ASD(x), axis='columns')
        .query('TD_ASDdiff != 0')
        .join(df['taxonomy'].str.split(';',8, expand=True).rename(columns={0:'Domain', 1:'Kingdom', 2:'Phylum', 
                                                                           3:'Class', 4:'Order', 5:'Family', 6:'Genus', 
                                                                           7:'Species'}))

        )
    return df1
