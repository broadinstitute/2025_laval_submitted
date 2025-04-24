A). In the H₂ analysis (1_PPI_CellxGene_facultative_coexpression_H2_analysis.ipynb), we investigated whether a lowly expressed protein of interest (POI) could still be functionally relevant if it exhibited stable co-expression with its interactors across diverse cell types. For instance, we've proteins like MLH1 that are not cell type specific but is broadly expressed across tissues/clusters for their fundamental roles (like DNA repair). First, we identified lowly expressed POIs based on an adaptive expression threshold (≤ median % expression across clusters within a tissue). Then, we assessed whether these POIs were "rescued" by their interactors. Rescue was defined by the presence of at least one interactor expressed above a certain threshold (≥30%) in the same cell type.

Based on the expression levels of these interactors, we assigned a "Compensation category" (Weak, Moderate, Robust) to represent the strength of functional compensation. We renamed the "%Cells Expressing Gene" to "%Cells Expressing Interactor" and renamed "POI_Expression_%" to "%Cells Expressing POI." Finally, the pipeline filters out redundant H3 cases at the end without affecting H2 & ensures correct sorting (R -> M -> W).

Broad review:

File Processing: Iterates through files, loads data, and extracts the POI name.
Moderate/Robust Interactors: Identifies clusters (Tissue/Cell Type combinations) where at least one interactor shows moderate or robust compensation.
Filtering: Retains all interactors (weak, moderate, robust) within these identified clusters.
PPI Construction: Creates the "POI-Interactor" pair string.
Compensation category mapping: Converts compensation categories into simplified "R," "M," "W" codes.
Self-Interaction flagging: Marks PPIs where the interactor is the POI itself.
Column renaming/selection: As described in the notebook
Sorting: Sorts by Tissue, Cell Type & Compensation category (R > M > W).
H3 Exclusion: Removes cases where both the POI & interactor have >= 30% expression (to avoid redundant data with H1 results)
Self-Interaction exclusion: Removes self-interaction rows

B). Then, we also have a script to prioritize PPIs  (2_PPI_priortization_analysis_H2.ipynb)

PPI prioritization script is designed to assess the impact of the H2 filtering process by quantifying how many interactors are retained or lost per Protein Of Interest (POI). First, the script extracts the complete list of unique interactors for each POI from the pre-filtered (Step-0) data. Then, it extracts the list of interactors that passed the H2 criteria from the post-filtered data. By comparing these two lists, the script calculates the number and percentage of interactors lost during the H2 filtering. The script then flags POIs where more than 50% of their interactors were lost. Finally, the script generates a summary bar plot visualizing the number of retained interactions per POI, with a reference line indicating the median retention rate. This allows for an overview of which POIs retained most of their interacting partners after the H2 analysis.
