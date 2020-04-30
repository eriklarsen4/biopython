# -*- coding: utf-8 -*-
"""
Created on Sat Apr 11 16:33:27 2020

@author: Erik
"""

## Clear the environment!
## Import the package with the pandas function
## Import the package with the numpy function
## Import the Counter, defaultdict modules from the collections package

import pandas as pd
from pandas import DataFrame
import numpy as np
from collections import Counter
from collections import defaultdict


## Import the Panther output ##
Panther_Analysis = pd.read_csv(r'M:\Erik\Data\Omics\pant2016_enriched.csv', index_col = 0, comment = ";")
## Create a variable out of the column that includes the DEGs identified by Panther to be associated with the term in the first column
Genes = Panther_Analysis.loc[:,'Panther_2016.Genes']
## Remove the semicolons from the list
Genes = Genes.str.split(pat = ";")
## Convert the new "series" to a dataframe and add the column title, "DEGs"
Genes = Genes.to_frame(name = "DEGs")

#Genes = Genes.iloc[0]

## Add the "Genes" dataframe to the Panther output dataframe as a column, titled "DEGs"
Panther_Analysis["DEGs"] = pd.DataFrame(Genes, index = Panther_Analysis.index)
## Remove unused/unnecessary columns
Panther_Analysis = Panther_Analysis.drop(['Panther_2016.Old.P.value', 'Panther_2016.Old.Adjusted.P.value', 'Panther_2016.Odds.Ratio', 'Panther_2016.Genes'], axis = 1)
## Rename the columns
Panther_Analysis = Panther_Analysis.rename(columns = {'Panther_2016.Term': 'Panther Term', 'Panther_2016.Overlap' : 'DEGs Identified in Panther Pathway', 'Panther_2016.P.value' : 'P-value', 'Panther_2016.Adjusted.P.value' : 'Adj. P-value', 'Panther_2016.Combined.Score': 'Panther Combined Score'})
## Subset the dataframe for export, including the first 11 Panther Gene Ontology terms and the associated data
Panther_Analysis2 = Panther_Analysis[:11]



## Export the new dataframe to the server
#Panther_Analysis2.to_csv('M://Omics//e13 Panther Pathway Analysis.csv')
## Export the bigger dataframe to the server
#Panther_Analysis.to_csv('M://Omics//Complete e13 Panther Pathway Analysis.csv')



## Re-order the dataframe by highest Combined Score
Panther_Analysis = Panther_Analysis.sort_values(by = "Panther Combined Score", ascending = False)
## Make a new subset for clustergram
Panther_Analysis3 = Panther_Analysis[:41]

## Export this dataframe for Bar plot generation in ggplot
Panther_Analysis3.to_csv('M://Omics//Panther Pathway Analysis e13 GT.csv')

## Remove columns besides Ontology terms and Genes
Panther_Analysis3 = Panther_Analysis3[["Panther Term", "DEGs"]]
## Re-set the newest Panther Analysis dataframe's indeces and remove the old index column
Panther_Analysis3 = Panther_Analysis3.reset_index()
Panther_Analysis3 = Panther_Analysis3.drop(["index"], axis = 1)
## Slice the appropriate rows (4, 18, 26, 33, 39, 41) as indexed by excel! -1!!!
## Create a list of the row indeces of those pathways
a = [3,17,24,32,39,40]
## Slice the data accordingly
Panther_Analysis3 = Panther_Analysis3.iloc[a]
Panther_Analysis3 = Panther_Analysis3.rename(columns = {0 : "DEGs"})
## Re-set the newest Panther Analysis dataframe's indeces and remove the old index column
Panther_Analysis3 = Panther_Analysis3.reset_index()
Panther_Analysis3 = Panther_Analysis3.drop(["index"], axis = 1)

## Unpack the DEGs into a list of DEGs within the selected pathways identified by Panther
Genes_Wnt_Subset = []
for i in Panther_Analysis3["DEGs"]:
    if i not in Genes_Wnt_Subset:
        Genes_Wnt_Subset.append(i)

## Linearize all the DEGs
Genes_Wnt_Subset = [item for elem in Genes_Wnt_Subset for item in elem]
## Turn the list into a data frame and rename the column of genes
Genes_Wnt_Subset = pd.DataFrame(Genes_Wnt_Subset)
Genes_Wnt_Subset = Genes_Wnt_Subset.rename(columns = {0 : "DEGs"})

## Identify genes that come up in multiple pathways
dup_genes = [k for k,v in Counter(Genes_Wnt_Subset["DEGs"]).items() if v>1]
dup_gene_indeces = defaultdict(list)

## Find indeces of duplicated genes
for i,item in enumerate(Genes_Wnt_Subset["DEGs"]):
    dup_gene_indeces[item].append(i)
dup_gene_indeces = {k:v for k,v in dup_gene_indeces.items() if len(v)>1}

## Remove duplicate genes
seen = set()
duplicates = []
for item in Genes_Wnt_Subset["DEGs"]:
    if item not in seen:
        seen.add(item)
        duplicates.append(item)

## Rename the sliced data to include only one DEG for all the relevant pathways
Genes_Wnt_Subset = duplicates
## Turn it into a dataframe
Genes_Wnt_Subset = pd.DataFrame(Genes_Wnt_Subset)
## Re-orient the dataframe
Genes_Wnt_Subset = Genes_Wnt_Subset.rename(columns = {0 : "DEGs"})

## Create the Clustergram dataframe first as an array, filling it with the appropriate dimensions of zeros, to be re-filled
Clustergram = np.zeros(shape = (len(Genes_Wnt_Subset), len(a)))

## Check to see if boolean evaluations can be performed
Genes_Wnt_Subset["DEGs"][0] in Panther_Analysis3["DEGs"][0]
#Panther_Analysis3["Panther Term"][0]
#Genes_Wnt_Subset["DEGs"][0]

## Re-populate the Clustergram array with values where DEGs (indeces) along the y axis are in Panther Pathways (indeces) across the x axis
for i in range(Genes_Wnt_Subset.shape[0]):
    for j in range(Panther_Analysis3.shape[0]):
        if Genes_Wnt_Subset["DEGs"][i] in Panther_Analysis3["DEGs"][j]:
            Clustergram[i, j] = 1
        else:
            Clustergram[i, j] = 0

## Convert the array to a dataframe, re-populating the indeces with the DEGs' and Pathways' names            
Clustergram = pd.DataFrame(data = Clustergram, index = Genes_Wnt_Subset["DEGs"], columns = Panther_Analysis3["Panther Term"])

## Remove the list
del a

## Check to see if dataframe has the correct number of DEGs overlapping with Panther-identified genes in the pathway/process
## (compare to excel file)
Clustergram[Panther_Analysis3["Panther Term"][1]].value_counts()

## Export the dataframe for ggplotting
Clustergram.to_csv(r'M:\Omics\Custom Python Enrichr Clustergram e13 GT.csv')


----------------------------------------------------------------------------------
