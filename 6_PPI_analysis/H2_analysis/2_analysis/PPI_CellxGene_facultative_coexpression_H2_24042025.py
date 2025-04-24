#!/usr/bin/env python
# coding: utf-8

# <p style="font-size:18px; font-weight:bold;"> 2025 Olivia Debnath</p>
# <p style="font-size:14px;">Dana-Farber Cancer Institute & Harvard Medical School</p>

# **Hypothesis 2 (H₂): Facilitated co-expression & functional dependence**
#     
# If a lowly expressed gene (POI: Protein of Interest) maintains stable but widespread co-expression (≥10%) with multiple interactors across diverse cell types, its function is cellularly indispensable despite its low abundance.
# 
# 🔹 Why this hypothesis?
# Some genes are not highly expressed, but they are still critical for cellular homeostasis and operate in a facilitated interaction manner, where multiple weakly co-expressed partners compensate for expression variability.
# 
# For example:
# 
# MLH1 (MutL Homolog 1) is a key DNA mismatch repair (MMR) gene but does not show strong expression peaks in any single cell type.
# Instead, MLH1 and its partners (e.g., PMS2, MSH2, MSH6) maintain steady but weak co-expression across multiple proliferating cell states—suggesting functional compensation rather than strong enrichment.
# 
# Technical pipeline (adjusted based on previous exploratory analyses):
# 
# If a lowly expressed POI (≤ Median% expression across clusters in a tissue) is consistently co-expressed with one or more interactors (≥30%), it remains functionally relevant, as its interactors might compensate for its deficiency. We can reduce the cut-off to 20% for now to pick up a few more cases for MLH1. 
# 
# Certain genes—like MLH1 (Mismatch Repair)— are weakly expressed (often <10% across most clusters or even <5%) but still cell-essential due to their role in highly conserved pathways. Instead of relying on strong co-expression (Hypothesis 1), these genes depend on their highly expressed interactors to maintain functionality.
# 
# 🔍 H₂ filtering strategy: 
# 
# 1. POI expression threshold:
# 
# - Identify cases where POI expression ≥ median % expression of the gene across clusters within the same tissue (adaptive thresholding).
# - This allows the detection of functional low-expression genes while ignoring highly abundant ones.
#  
# 2. Interactor rescue requirements:
# 
# - At least one interactor must be expressed at ≥30% in the same cell type/state.
# - If multiple interactors are highly expressed, it strengthens the functional dependence hypothesis.
# 
# 
# 3. Rescue strength annotation: (new column: "Rescue Category")
# 
#     - Weak Rescue (10-30%) → Interaction present, but not strong.
#     - Moderate Rescue (30-50%) → Likely functionally relevant.
#     - High Rescue (50-80%) → Strong compensatory effect.
#     - Robust Rescue (80-100%) → Near-total functional compensation.

# In[1]:


import os
print(os.path.exists("./results_11032025/Jess_PPI_21032025/PPI_preprocessed_15042025/Level2/filtered_S0/"))


# In[2]:


#Define input directory: 
input_dir = "./results_11032025/Jess_PPI_21032025/PPI_preprocessed_15042025/Level2/filtered_S0/"  

#Dynamically find all relevant input files
input_files = sorted([f for f in os.listdir(input_dir) if f.endswith("_filtered_S0_17042025.csv")]) 
print(input_files) 

#Print total count of files
print(f"\n✅ Total Input Files Found: {len(input_files)}")


# In[3]:


#Specify output directory: 
output_dir = "./results_11032025/Jess_PPI_21032025/PPI_preprocessed_15042025/PPI_contextualization/filtered_Step2_H2/"
os.makedirs(output_dir, exist_ok=True)  #Ensure output directory exists


# In[4]:


import pandas as pd
import numpy as np


# 1). POI expression cutoff: - POI must be ≥ median expression % (across clusters in that tissue). - We use the tissue-specific median rather than a fixed 30% cutoff, because some genes (e.g., MLH1) have consistently low but functional expression.
# 
# 2). Interactor expression cutoff: - At least one interactor must be ≥20% expression (instead of the 30% threshold used in H₁).

# In[5]:


#Applied code from PPI_CellxGene_H2_corr_27022025.ipynb on all PPI cases: 

#**Define Expression Cutoffs for H₂**
WEAK_COMPENSATION = 20   #At least 20% interactor expression
MODERATE_COMPENSATION = 30  #At least 30% interactor expression
STRONG_COMPENSATION = 50  #At least 50% interactor expression


# In[6]:


import pandas as pd
import numpy as np
from scipy.stats import zscore
import matplotlib.pyplot as plt


# We need to remove cases where both the POI & the interactor are ≥30% expression to prevent redundancy with H3 (strong co-expression & functional disruption). However, if the POI is ≥30% but the interactor is <30%, we still keep those cases.

# In[7]:


def process_H2_file(df, protein_of_interest):
    """ 
    Filters for H₂: Facilitated Co-Operation
    - POI must be expressed at least at the **median expression %** (within its tissue).
    - At least one interactor must be highly expressed (≥20%) in the same cell type.
    - **POI expression is explicitly retained & displayed.**
    """

    print(f"\n🔍STEP-1: Computing tissue median expression for {protein_of_interest}...")

    #Step-1: Compute POI-Specific Median Expression Cutoff
    #Computes median expression % of POI across all cell types in each tissue
    #If POI is absent in a tissue, it won't be processed
    tissue_median_expr = df[df["Gene Symbol"] == protein_of_interest].groupby("Tissue")["%Cells Expressing Gene"].median()

    if tissue_median_expr.empty:
        print(f"⚠️WARNING: {protein_of_interest} does not meet median expression in any tissue. Skipping file.")
        return None

    print(f"📊Tissue-specific median expression for {protein_of_interest}:\n{tissue_median_expr}\n")

    #Merge computed tissue-specific median expression into df for reference
    df = df.merge(tissue_median_expr.rename("Tissue_Median_Exp"), on="Tissue", how="left")

    #Explicitly store POI expression per cell type
    #This ensures that even after filtering, POI expression remains in output
    poi_expression_dict = df[df["Gene Symbol"] == protein_of_interest].set_index(["Tissue", "Cell Type"])["%Cells Expressing Gene"].to_dict()
    df["POI_Expression_%"] = df.apply(lambda row: poi_expression_dict.get((row["Tissue"], row["Cell Type"]), np.nan), axis=1)

    #Ensure POI_Above_Median is explicitly written (True/False instead of NaN)
    # - If POI expression in a cell type is ≥ tissue median, it's marked as True
    df["POI_Above_Median"] = df["POI_Expression_%"] >= df["Tissue_Median_Exp"]
    df["POI_Above_Median"] = df["POI_Above_Median"].fillna(False).astype(bool) #Convert NaN to False where applicable

    #STEP-2: Extract only POI rows that pass the median expression cutoff
    #Ensures POI_Expression_% is stored per tissue-cell type.
    #Ensures POI_Above_Median is explicitly marked True or False
    #Since Jupyter prints full df before filtering, some False values will still show in output logs
    main_protein_present = df[df["POI_Above_Median"]]

    print(f"✅STEP-2: Found {len(main_protein_present)} clusters where {protein_of_interest} is expressed above median.")
    print(main_protein_present.head(), "\n")

    if main_protein_present.empty:
        print(f"⚠️WARNING: {protein_of_interest} does not surpass median expression threshold in any tissue. Skipping file.")
        return None

    #STEP-3: Identify interactors that are **highly expressed** (≥20%)
    #Filters only interactors (excluding POI itself) where %Cells Expressing Gene is ≥ 20%
    interactors_present = df[(df["Gene Symbol"] != protein_of_interest) & (df["%Cells Expressing Gene"] >= WEAK_COMPENSATION)]

    print(f"✅STEP-3: Found {len(interactors_present)} interactors with ≥20% expression.")
    print(interactors_present.head(), "\n")

    if interactors_present.empty:
        print(f"⚠️WARNING: No interactors reach ≥20% expression. Skipping file.")
        return None

    #STEP-4: Identify valid interactions (POI + at least one interactor in the same cluster)
    #Merges clusters where POI is above median with clusters where an interactor is ≥ 20%
    valid_interactions = pd.merge(
        main_protein_present[["Tissue", "Cell Type"]],
        interactors_present[["Tissue", "Cell Type"]],
        on=["Tissue", "Cell Type"], how="inner"
    ).drop_duplicates()

    print(f"✅STEP-4: Found {len(valid_interactions)} valid POI-interactor pairs.")
    print(valid_interactions.head(), "\n")

    if valid_interactions.empty:
        print(f"⚠️WARNING: No valid POI-interactor pairs found for {protein_of_interest}. Skipping file.")
        return None

    #STEP-5: Merge valid interactions back into the dataset
    #Ensures only rows with a POI-interactor pair are retained
    df = df.merge(valid_interactions, on=["Tissue", "Cell Type"], how="inner")

    print(f"✅STEP-5: Filtered dataset now contains {len(df)} rows.")
    print(df.head(), "\n")

    #STEP-6: Define Compensation Category Based on Interactor Expression
    #Categorizes interactions based on interactor expression level
    conditions = [
        df["%Cells Expressing Gene"] >= STRONG_COMPENSATION,  #Robust Compensation (≥50%)
        df["%Cells Expressing Gene"] >= 30,  #Moderate Compensation (30-49%)  <-- ✅ FIXED RANGE
        df["%Cells Expressing Gene"] >= 20  #Weak Compensation (20-29%)  <-- ✅ FIXED RANGE
    ]
    categories = ["Robust Compensation (≥50%)", "Moderate Compensation (30-49%)", "Weak Compensation (20-29%)"]

    df["Compensation Category"] = np.select(conditions, categories, default="No Compensation (<20%)")

    #STEP-7: Remove Lowly Expressed Interactors (<20%)
    #Interactions that fail the ≥20% expression cutoff are removed
    df = df[df["Compensation Category"] != "No Compensation (<20%)"]

    print(f"✅STEP-7: Final dataset contains {len(df)} rows after removing non-compensating interactors.")

    return df


# In[8]:


output_dir = "./results_11032025/Jess_PPI_21032025/PPI_preprocessed_15042025/PPI_contextualization/filtered_Step2_H2/"


for file in input_files:
    input_path = os.path.join(input_dir, file)
    
    print(f"\n📂 Processing file: {file}")

    #Extract POI from filename
    protein_of_interest = file.split("_")[0]
    print(f"🧬 Identified POI: {protein_of_interest}")

    #Read CSV (Fix encoding issues)
    df = pd.read_csv(input_path, encoding="utf-8")  
    print(f"✅ Loaded {file} | Shape: {df.shape}")

    if df.empty:
        print(f"⚠️ Skipping {file} as it's empty.")
        continue

    #Process file using H2 criteria
    processed_df = process_H2_file(df, protein_of_interest)

    if processed_df is None or processed_df.empty:
        print(f"⚠️ WARNING: No valid data after processing {file}. Skipping...")
        continue

    #Save output in Excel format (Fix encoding issues)**
    output_file = os.path.join(output_dir, file.replace("_S0_17042025", "_S2_H2_24042025").replace(".csv", ".xlsx")) 
    processed_df.to_excel(output_file, index=False, engine= "openpyxl")  

    print(f"✅ Successfully saved H2 output: {output_file}")


# In[ ]:





# Key confirmations:
# 
# 1. All POI rows are first filtered for POI_Above_Median.
# 
# df["POI_Above_Median"] = df["POI_Expression_%"] >= df["Tissue_Median_Exp"]
# Only True rows move to the next step (Step 2).
# 
# 
# 2. Only interactors with ≥20% expression are considered.
# 
# interactors_present = df[(df["Gene Symbol"] != protein_of_interest) & (df["%Cells Expressing Gene"] >= WEAK_COMPENSATION)]
# If no interactors qualify, the file is skipped (Step 3).
# 
# 
# 3. Final filtering merges valid POI & interactors, ensuring that both conditions are met.
# 
# df = df.merge(valid_interactions, on=["Tissue", "Cell Type"], how="inner")
# Only valid POI-interactor pairs are retained (Step 5).
# 
# 
# 4. The output will only contain PPIs where both conditions hold (Steps 6-7).
# 
# POI expression is shown in the final file for reference (POI_Expression_%).
# Rows failing either condition are removed before saving.

# In[10]:


#Just print the moderate & robust compensations in the notebook:

#Find all processed H2 files
h2_files = [f for f in os.listdir(output_dir) if f.endswith("_S2_H2_24042025.xlsx")]
len(h2_files)


# - We need to ensure that while filtering for Moderate & Robust Compensation, we retain weakly expressed interactors within the same cluster for context.
# 
# - Also, instead of completely removing POI as an interactor, let's explicitly flag self-interactions so that they can be reviewed manually.
# 
# - Updates for handling self-interactions:
#     - If POI is present as an interactor, mark it as "Self-Interaction" in a new column.
#     - Keep all cases in the output but differentiate them so you can review them later.
#     - Sorting remains unchanged, but self-interactions will be labeled clearly.
#     - Filters out H3 redundant cases at the end (Step-9) without affecting H2
#     - Renames %Cells Expressing Gene → %Cells Expressing Interactor
#     - Renames POI_Expression_% → %Cells Expressing POI
#     - Ensures sorting remains correct (R → M → W)

# In[11]:


#let's retain all relevant columns, including Cell Count, so you have full flexibility
#Self_Interaction == "Yes" rows should be removed just before saving the Excel file.
#Refer to line-16 of PPI_CellxGene_H2_corr_27022025.ipynb 

#Process each file
for file in h2_files:
    file_path = os.path.join(output_dir, file)
    
    #Load the processed file
    df = pd.read_excel(file_path, engine="openpyxl")
    
    print(f"\n📂 Processing file: {file}")
    print(f"✅ Loaded {file} | Shape: {df.shape}")

    #Extract POI from filename
    protein_of_interest = file.split("_")[0]

    #Step 1: Identify clusters where at least one interactor is Moderate or Robust
    valid_clusters = df[df["Compensation Category"].isin(["Robust Compensation (≥50%)", "Moderate Compensation (30-49%)"])][["Tissue", "Cell Type"]].drop_duplicates()

    if valid_clusters.empty:
        print(f"⚠️ No valid interactor-driven Robust or Moderate Compensation cases found in {file}. Skipping...\n")
        continue

    #Step 2: Retain all interactors (including weak ones) in these clusters
    filtered_df = df.merge(valid_clusters, on=["Tissue", "Cell Type"], how="inner")

    print(f"✅ Found {len(filtered_df)} total interactors in clusters with at least one Moderate/Robust Compensation.")

    #Step 3: Construct the PPI column (POI-Interactor Pair)
    filtered_df["PPI"] = protein_of_interest + "-" + filtered_df["Gene Symbol"]  #Format: POI-Interactor (e.g., STXBP1-STX5)

    #Step 4: Assign W, M, R for Compensation Category (Fixing Mapping Issue)
    category_mapping = {
        "Robust Compensation (≥50%)": "R",
        "Moderate Compensation (30-49%)": "M",
        "Weak Compensation (20-29%)": "W"
    }

    #Standardize category column before mapping
    filtered_df["Compensation Category"] = filtered_df["Compensation Category"].astype(str).str.strip()

    #Map category names
    filtered_df["Comp_Category"] = filtered_df["Compensation Category"].map(category_mapping)

    #Handle missing values (should not happen, but just in case)
    filtered_df["Comp_Category"] = filtered_df["Comp_Category"].fillna("Unknown")

    #Step 5: Flag Self-Interactions
    filtered_df["Self_Interaction"] = filtered_df["Gene Symbol"].apply(lambda x: "Yes" if x == protein_of_interest else "No")

    #Step 6: Rename columns for consistency
    filtered_df.rename(columns={
        "%Cells Expressing Gene": "%Cells Expressing Interactor",
        "POI_Expression_%": "%Cells Expressing POI"
    }, inplace=True)

    #Step 7: Select only required columns (PPI-related & Cell Count)
    output_cols = ["Tissue", "Cell Type", "Cell Count", "PPI", "%Cells Expressing POI", "%Cells Expressing Interactor", "Comp_Category", "Self_Interaction"]
    final_df = filtered_df[output_cols]  #Retain only relevant columns

    #Step 8: Sort by Tissue → Cell Type → Compensation (R → M → W)
    comp_rank = {"R": 1, "M": 2, "W": 3}
    final_df["Comp_Rank"] = final_df["Comp_Category"].map(comp_rank)

    sorted_df = final_df.sort_values(by=["Tissue", "Cell Type", "Comp_Rank"],
                                     ascending=[True, True, True]).drop(columns=["Comp_Rank"])

    #Step 9: **Exclude H3 redundant cases (where both POI & Interactor are ≥30%)** => very important to avoid redundancy 
    sorted_df = sorted_df[~((sorted_df["%Cells Expressing POI"] >= 30) & 
                            (sorted_df["%Cells Expressing Interactor"] >= 30))]

    print(f"✅ Removed H3 redundant cases. Final dataset contains {len(sorted_df)} rows.")

    #Step 10: Print first 10 rows for quick verification
    print(sorted_df.head(10)) 

    #Step 11: Remove Self-Interactions before saving the files
    sorted_df = sorted_df[sorted_df["Self_Interaction"] != "Yes"]

    #Step 12: Save the sorted version
    sorted_file_path = os.path.join(output_dir, file.replace("_S2_H2_24042025.xlsx", "_Sorted_S2_H2_24042025.xlsx"))
    sorted_df.to_excel(sorted_file_path, index=False, engine="openpyxl")

    print(f"✅ Saved sorted file (Only PPI + Cell Count, Self-Interactions Removed): {sorted_file_path}\n")
