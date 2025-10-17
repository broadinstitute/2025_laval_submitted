import seaborn as sns

clinvar_category = ['1_Pathogenic', '2_Benign', '3_Conflicting', '4_VUS', '5_Others', '6_No_ClinVar']
clinvar_palette_clinvar = sns.color_palette("Set2")
clinvar_palette_clinvar[0], clinvar_palette_clinvar[1], clinvar_palette_clinvar[2], \
    clinvar_palette_clinvar[3], clinvar_palette_clinvar[4], clinvar_palette_clinvar[5], clinvar_palette_clinvar[7] = \
clinvar_palette_clinvar[1], clinvar_palette_clinvar[4], clinvar_palette_clinvar[3], \
    clinvar_palette_clinvar[2], clinvar_palette_clinvar[0], clinvar_palette_clinvar[7], clinvar_palette_clinvar[5]

palette_dict = {
    "clinvar_clnsig_clean_oneperc": dict(
        zip(clinvar_category, clinvar_palette_clinvar[:7])
    ),
    "clinvar_clnsig_clean":{
        '1_Pathogenic': "#CA7682",
        '2_Benign': "#1D7AAB",
        '3_Conflicting': "#505050", #"grey", or lighter red "#E6B1B8"
        '4_VUS': "#A0A0A0", ## Pillar default
        '5_Others': "#E0E0E0"
        # '6_No_ClinVar'
    }
}

## From Malvika
# colors_custom = {
#     "Pathogenic": "#CA7682",
#     "Pathogenic*": "#CA7682",
#     "Likely Pathogenic": "#E6B1B8",
#     "Likely Pathogenic*": "#E6B1B8",
#     "Likely pathogenic": "#E6B1B8",
#     "Benign": "#1D7AAB",
#     "Benign*": "#1D7AAB",
#     "Uncertain significance": "#A0A0A0",
#     "Combined LR is not <0.5 or >2": "#A0A0A0",
#     "VUS": "#A0A0A0",
#     "No Classification": "#A0A0A0",
#     "Model not applied": "#A0A0A0",
#     "Likely Benign": "#63A1C4",
#     "Likely Benign*": "#63A1C4",
#     "Likely benign": "#63A1C4",
#     '-12': "#176082",  # deepest blue
#     '-11': "#1D7AAB",  # deep steel blue
#     '-10': "#4B91A6",  # medium steel blue
#     '-9': "#63A1C4",   # sky blue
#     '-8': "#7AB5D1",   # soft blue-gray
#     '-7': "#99C8DC",   # light blue-gray
#     '-6': "#B8DCE8",   # very light blue
#     '-5': "#D0E8F0",   # icy pale blue
#     '-4': "#E4F1F6",   # faint pale blue
#     '-3': "#EDF6FA",   # hint of blue
#     '-2': "#F4F9FC",   # nearly white blue
#     '-1': "#F9FCFE",   # even paler blue
#     '0': "#E0E0E0",    # true neutral
#     '1': "#F7E4E7",    # very pale pink
#     '2': "#F0D0D5",    # pale dusty pink
#     '3': "#E6B1B8",    # blush pink
#     '4': "#D68F99",    # dusty rose
#     '5': "#CA7682",    # muted rose
#     '6': "#B85C6B",    # deeper rose
#     '7': "#A84957",    # muted wine
#     '8': "#943744",    # rich berry red
#     '9': "#7F2936",    # burgundy
#     '10': "#671B28",   # dark wine red
#     '11': "#520F1C",   # deeper red
#     '12': "#3A060D",   # near black red
#     'Pathogenic to Benign': "#CA7682",
#     'Benign to Pathogenic': "#7AB5D1",
#     'Pathogenic to Pathogenic': "#B85F6E",
#     'Benign to Benign': "#5C9EBD",
#     'Benign to Uncertain': "#6FA6C5",
#     'Pathogenic to Uncertain': "#C27B87",
#     'Benign to Conflicting': "#80B3CC",
#     'Pathogenic to Conflicting': "#CC8A94",
#     'Benign to No Evidence': "#92C1D6",
#     'Pathogenic to No Evidence': "#D29AA3",
#     'No evidence': "#E0E0E0",        # base gray
#     'No Classification': "#A0A0A0",  # silver gray
#     'Conflicting evidence': "grey"   # neutral gray
# }