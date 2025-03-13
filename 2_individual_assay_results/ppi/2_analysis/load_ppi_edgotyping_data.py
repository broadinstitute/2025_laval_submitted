#TO DO:
# include hsPRS-RRSv2
#

import pandas as pd
import numpy as np
from TableToDict import GetDFv2, GetDF, DisplayColumns

from load_pdest_stability import load_pdest_dataset
from load_n2h_data import get_mn2h_allele_df, get_yN2H_allele_df, get_mn2h_normalized_df, get_yn2h_nlr_df

# def load_raw_results():
#     q = 'select disease_orf_id, swim_seq_validated from vus.wt_disease_nodes_summary;'
#     col = ['disease_orf_id', 'swim_seq_validated']

#     wt_disease_nodes_summary = GetDFv2(q,col)
#     wt_validated_d = dict(zip(wt_disease_nodes_summary.disease_orf_id, wt_disease_nodes_summary.swim_seq_validated))

#     q = 'select ad_orf_id, assay_id, swim_seq_validated from vus.vus_ad_interactors;'
#     col = ['ad_orf_id', 'assay_id', 'swim_seq_validated']

#     vus_ad_interactors = GetDFv2(q,col)
#     vus_ad_interactors_d = dict(zip(vus_ad_interactors.ad_orf_id, vus_ad_interactors.swim_seq_validated))

#     q = 'select orf_id, mut_id, nt_change, swim_seq_validated from vus.vus_db_alleles;'
#     col = ['orf_id', 'mut_id', 'nt_change', 'swim_seq_validated']

#     vus_db_alleles = GetDFv2(q,col)
#     vus_db_alleles_d = {(o,m):s for o,m,s in zip(vus_db_alleles.orf_id, vus_db_alleles.mut_id, vus_db_alleles.swim_seq_validated)}

#     q = 'select ad_orf_id, db_orf_id, db_mut_id, retest_pla, retest_pos, standard_batch, assay_id, db_source_plate, db_source_well, ad_source_plate, ad_source_well, LW, LWH1, LWH10, LWH25, LWA, LWAH1, db_validated, ad_validated from vus.retest where standard_batch in ("VUSAPWT1B1","VUSAPWT1B2","VUSAPWT2B1","VUSAPWT6B1") and db_orf_id is not null and ad_orf_id is not null;'

#     col = ['ad_orf_id', 'db_orf_id', 'db_mut_id', 'retest_pla', 'retest_pos', 'standard_batch', 'assay_id', 'db_source_plate', 'db_source_well', 'ad_source_plate', 'ad_source_well', 'LW', 'LWH1', 'LWH10', 'LWH25', 'LWA', 'LWAH1', 'db_validated', 'ad_validated']

#     retest_df = GetDFv2(q,col)
#     retest_df['db_source_validated'] = [vus_db_alleles_d[db,m] if (db,m) in vus_db_alleles_d else wt_validated_d[db] for db,m in zip(retest_df.db_orf_id, retest_df.db_mut_id)]
#     retest_df['ad_source_validated'] = retest_df.ad_orf_id.map(vus_ad_interactors_d)

#     retest_df['seq_validated'] = [(dbs == 1 or dbd == 1) and (ads == 1 or add == 1) for dbs, ads, dbd, add in zip(retest_df.db_source_validated, retest_df.ad_source_validated, retest_df.db_validated, retest_df.ad_validated)]

#     return retest_df

def load_raw_results():
    import os
    if 'georgescoppin' not in os.getcwd():
        q = 'select disease_orf_id, swim_seq_validated from vus.wt_disease_nodes_summary;'
        col = ['disease_orf_id', 'swim_seq_validated']

        wt_disease_nodes_summary = GetDFv2(q,col)
        wt_validated_d = dict(zip(wt_disease_nodes_summary.disease_orf_id, wt_disease_nodes_summary.swim_seq_validated))

        q = 'select ad_orf_id, assay_id, swim_seq_validated from vus.vus_ad_interactors;'
        col = ['ad_orf_id', 'assay_id', 'swim_seq_validated']

        vus_ad_interactors = GetDFv2(q,col)
        vus_ad_interactors_d = dict(zip(vus_ad_interactors.ad_orf_id, vus_ad_interactors.swim_seq_validated))

        q = 'select orf_id, mut_id, nt_change, swim_seq_validated from vus.vus_db_alleles;'
        col = ['orf_id', 'mut_id', 'nt_change', 'swim_seq_validated']

        vus_db_alleles = GetDFv2(q,col)
        vus_db_alleles_d = {(o,m):s for o,m,s in zip(vus_db_alleles.orf_id, vus_db_alleles.mut_id, vus_db_alleles.swim_seq_validated)}

        q = 'select ad_orf_id, db_orf_id, db_mut_id, retest_pla, retest_pos, standard_batch, assay_id, db_source_plate, db_source_well, ad_source_plate, ad_source_well, LW, LWH1, LWH10, LWH25, LWA, LWAH1, db_validated, ad_validated from vus.retest where standard_batch in ("VUSAPWT1B1","VUSAPWT1B2","VUSAPWT2B1","VUSAPWT6B1") and db_orf_id is not null and ad_orf_id is not null;'

        col = ['ad_orf_id', 'db_orf_id', 'db_mut_id', 'retest_pla', 'retest_pos', 'standard_batch', 'assay_id', 'db_source_plate', 'db_source_well', 'ad_source_plate', 'ad_source_well', 'LW', 'LWH1', 'LWH10', 'LWH25', 'LWA', 'LWAH1', 'db_validated', 'ad_validated']

        retest_df = GetDFv2(q,col)
        retest_df['db_source_validated'] = [vus_db_alleles_d[db,m] if (db,m) in vus_db_alleles_d else wt_validated_d[db] for db,m in zip(retest_df.db_orf_id, retest_df.db_mut_id)]
        retest_df['ad_source_validated'] = retest_df.ad_orf_id.map(vus_ad_interactors_d)

        retest_df['seq_validated'] = [(dbs == 1 or dbd == 1) and (ads == 1 or add == 1) for dbs, ads, dbd, add in zip(retest_df.db_source_validated, retest_df.ad_source_validated, retest_df.db_validated, retest_df.ad_validated)]
    else:
        retest_df = load_raw_results_from_dropbox_backup()

    return retest_df

def load_raw_results_from_dropbox_backup():
    dropbox_path = '/Users/georgescoppin/Dropbox/edgotyping_dropbox/large_data_files/ccsb_database_useful_tables'
    wt_disease_nodes_summary = pd.read_csv(f'{dropbox_path}/vus/wt_disease_nodes_summary.tsv', sep='\t')
    wt_validated_d = dict(zip(wt_disease_nodes_summary.disease_orf_id, wt_disease_nodes_summary.swim_seq_validated))

    vus_ad_interactors = pd.read_csv(f'{dropbox_path}/vus/vus_ad_interactors.tsv', sep='\t')
    vus_ad_interactors_d = dict(zip(vus_ad_interactors.ad_orf_id, vus_ad_interactors.swim_seq_validated))

    vus_db_alleles = pd.read_csv(f'{dropbox_path}/vus/vus_db_alleles.tsv', sep='\t')
    vus_db_alleles_d = {(o,m):s for o,m,s in zip(vus_db_alleles.orf_id, vus_db_alleles.mut_id, vus_db_alleles.swim_seq_validated)}

    retest_df = pd.read_csv(f'{dropbox_path}/vus/retest.tsv', sep='\t', low_memory=False)
    col = ['ad_orf_id', 'db_orf_id', 'db_mut_id', 'retest_pla', 'retest_pos', 'standard_batch', 'assay_id', 'db_source_plate', 'db_source_well', 'ad_source_plate', 'ad_source_well', 'LW', 'LWH1', 'LWH10', 'LWH25', 'LWA', 'LWAH1', 'db_validated', 'ad_validated']

    retest_df = retest_df[(retest_df.standard_batch.isin(["VUSAPWT1B1","VUSAPWT1B2","VUSAPWT2B1","VUSAPWT6B1"]))&
                          (~retest_df.db_orf_id.isnull())&(~retest_df.ad_orf_id.isnull())
                          ][col].copy()
    
    retest_df['db_source_validated'] = [vus_db_alleles_d[db,m] if (db,m) in vus_db_alleles_d else wt_validated_d[db] for db,m in zip(retest_df.db_orf_id, retest_df.db_mut_id)]
    retest_df['ad_source_validated'] = retest_df.ad_orf_id.map(vus_ad_interactors_d)

    retest_df['seq_validated'] = [(dbs == 1 or dbd == 1) and (ads == 1 or add == 1) for dbs, ads, dbd, add in zip(retest_df.db_source_validated, retest_df.ad_source_validated, retest_df.db_validated, retest_df.ad_validated)]

    return retest_df


def get_norm_score_df(retest_df, output='alleles'):
    '''
    if output = 'all': output contains PRS/RRS/WT/Allele scores
    if output = 'alleles': output restricted to valid alleles (wt positive and LW etc)
    '''
    
    df = retest_df.query("~db_orf_id.isnull() and ad_orf_id == 0 and (db_validated == 1 or db_source_validated == 1)")
    EAD_d = {(db_orf_id,db_mut_id):{'LW':{},'LWA':{}, 'LWH1':{}, 'LWH10':{}, 'LWH25':{}, 'LWAH1':{}} for db_orf_id,db_mut_id in  zip(df.db_orf_id, df.db_mut_id)}

    for db_orf_id, db_mut_id, LW, LWA, LWH1, LWH10, LWH25, LWAH1 in zip(df.db_orf_id, df.db_mut_id, df.LW, df.LWA,df.LWH1,df.LWH10,df.LWH25,df.LWAH1):
        for media, score in zip(['LW','LWA', 'LWH1', 'LWH10', 'LWH25', 'LWAH1'],[LW,LWA, LWH1, LWH10, LWH25, LWAH1]):
            EAD_d[(db_orf_id,db_mut_id)][media] = score

    col = ['ad_orf_id', 'db_orf_id', 'db_mut_id','standard_batch', 'assay_id', 'LW', 'LWH1', 'LWH10', 'LWH25','LWA', 'LWAH1']
    norm_score_df = retest_df.query("ad_orf_id > 0 and seq_validated and LW > 0")[col].copy()
    
    def normalize_score(lw_l,score_l):
        lw_test,lw_ead,lw = lw_l
        test,ead_ctrl = score_l
        
        res = np.nan 
        if lw and test == test and ead_ctrl == ead_ctrl:
            
            if test > 1 and ead_ctrl <= 1:
                res = 1
            elif test > 1 and ead_ctrl > 1:
                res = 'AA'
            elif test <= 1:
                res = 0
                
        return [test,ead_ctrl,res]
    
    for media in ['LW','LWH1', 'LWH10', 'LWH25', 'LWA', 'LWAH1']:
        norm_score_df[media] = [[score, EAD_d[(db_orf_id, db_mut_id)][media]] if (db_orf_id, db_mut_id) in EAD_d else [score,np.nan] for db_orf_id, db_mut_id, score in zip(norm_score_df.db_orf_id, norm_score_df.db_mut_id, norm_score_df[media])]
        if media == 'LW':
            norm_score_df[media] = norm_score_df[media].apply(lambda x: x+[x[0] > 1 and x[1] > 1])
        else:
            norm_score_df[media] = norm_score_df.apply(lambda row: normalize_score(row['LW'], row[media]),axis=1)
            
    norm_score_df['sum_score'] = [sum([e[-1] for e in [lwh1,lwh10,lwh25,lwa,lwah1] if e[-1] in [0,1]]) if lw[-1] else np.nan for lw,lwh1,lwh10,lwh25,lwa,lwah1 in zip(norm_score_df.LW, norm_score_df.LWH1, norm_score_df.LWH10, norm_score_df.LWH25, norm_score_df.LWA, norm_score_df.LWAH1)]
    norm_score_df['n_condition'] = [sum([1 for e in [lwh1,lwh10,lwh25,lwa,lwah1] if e[-1] in [0,1]]) if lw[-1] else np.nan for lw,lwh1,lwh10,lwh25,lwa,lwah1 in zip(norm_score_df.LW, norm_score_df.LWH1, norm_score_df.LWH10, norm_score_df.LWH25, norm_score_df.LWA, norm_score_df.LWAH1)]
    norm_score_df['score_norm'] =  norm_score_df.sum_score / norm_score_df.n_condition

    df = norm_score_df.query("db_mut_id == 0")
    wt_scored_d =  {(db,ad,assay):s for db,ad,assay,s in zip(df.db_orf_id, df.ad_orf_id, df.assay_id,df.score_norm)}
    
    norm_score_df['wt_score_norm'] = [wt_scored_d[db,ad,assay] if (db,ad,assay) in wt_scored_d else np.nan for db,ad,assay in zip(norm_score_df.db_orf_id, norm_score_df.ad_orf_id, norm_score_df.assay_id)]
    norm_score_df.rename(columns = {'score_norm':'allele_score_norm'}, inplace=True)
    norm_score_df['valid'] = (norm_score_df.wt_score_norm>0)&(~norm_score_df.allele_score_norm.isnull())&(norm_score_df.n_condition > 1)

    return norm_score_df


def calculate_similarity(mut_profile1, mut_profile2):
    '''
    - return 1 if same as WT
    - less than 1 if 
    '''
    
    keys = {k for k in mut_profile1.keys()&mut_profile2.keys() if mut_profile1[k] == mut_profile1[k] and mut_profile2[k] == mut_profile2[k]}
    res = np.mean([1 - (abs(mut_profile1[k] - mut_profile2[k])/max([mut_profile1[k], mut_profile2[k]])) if mut_profile1[k] + mut_profile2[k] > 0 else 1 for k in keys])
    
    return res


def calculate_similarity_euclid_max_score_norm(mut_profile1, mut_profile2, wt_profile, how='euclid'):
    '''calculate similarity between two profiles
    
    return a score between 0 and 1
    distance can be calculated using euclidean 
    score = distance / max distance 
    how: euclid or taxicab (default is euclid)
    
    '''
    
    method = 2 if how == 'euclid' else 1
    
    keys = {k for k in mut_profile1.keys()&mut_profile2.keys() if mut_profile1[k] == mut_profile1[k] and mut_profile2[k] == mut_profile2[k]}
    
    dist = sum(abs(mut_profile1[k] - mut_profile2[k])**method for k in keys)**(1/method) 
    max_dist = sum(max(mut_profile1[k], mut_profile2[k], wt_profile[k])**method for k in keys)**(1/method) 

    if len(keys) > 0:
        try:
            score = 1-(dist/max_dist)
        except:
            print(dist,max_dist,mut_profile1,mut_profile2)
            raise
    else:
        score = np.nan

    return score

def get_edg_profile_df(normalized_df):
    df = normalized_df.query("db_mut_id != 999999")
    
    edg_d = dict()
    for i,(ad_orf_id,db_orf_id,db_mut_id,standard_batch,assay_id,LW,LWH1,LWH10,LWH25,LWA,LWAH1,sum_score,n_condition,allele_score_norm,wt_score_norm,valid) in df.iterrows():
        if (db_orf_id, db_mut_id) not in edg_d:
            edg_d[(db_orf_id, db_mut_id)] =  {'ad_profile':{ad_orf_id:allele_score_norm}}
        else:
            edg_d[(db_orf_id, db_mut_id)]['ad_profile'][ad_orf_id] = allele_score_norm
            
    all_d = {orf_id:dict() for orf_id in df.db_orf_id}
    for i,(ad_orf_id,db_orf_id,db_mut_id,standard_batch,assay_id,LW,LWH1,LWH10,LWH25,LWA,LWAH1,sum_score,n_condition,allele_score_norm,wt_score_norm,valid) in df.iterrows():
        if db_mut_id not in all_d[db_orf_id]:
            all_d[db_orf_id][db_mut_id] =  {ad_orf_id:allele_score_norm}
        else:
            all_d[db_orf_id][db_mut_id][ad_orf_id] = allele_score_norm


    edg_df = pd.DataFrame(edg_d).T
    edg_df.reset_index(inplace=True)
    edg_df.rename(columns={'level_0':'db_orf_id', 'level_1':'db_mut_id'},inplace=True)

    wt_d = dict(zip(edg_df[edg_df.db_mut_id == 0].db_orf_id, edg_df[edg_df.db_mut_id == 0].ad_profile))
    
    edg_df['wt_score'] = edg_df.db_orf_id.map(wt_d)
    edg_df['n_allele'] = edg_df.db_orf_id.apply(lambda x: dict(edg_df.db_orf_id.value_counts())[x] - 1)
    edg_df['degree_tested'] = edg_df.ad_profile.apply(len)
    edg_df.sort_values(by=['db_orf_id','db_mut_id'],inplace=True)
    
    
    edg_df['wt_similarity'] = [calculate_similarity_euclid_max_score_norm(allele_profile, wt_profile, wt_profile) for allele_profile, wt_profile in zip(edg_df.ad_profile, edg_df.wt_score)]
    edg_df['similarity'] = [{db_mut_id:calculate_similarity_euclid_max_score_norm(allele_profile,all_d[db_orf_id][db_mut_id], wt_profile) for db_mut_id in  all_d[db_orf_id]} for allele_profile, db_orf_id, wt_profile in zip(edg_df.ad_profile, edg_df.db_orf_id, edg_df.wt_score)]
    #edg_df['similarity'] = [{db_mut_id:(allele_profile,all_d[db_orf_id][db_mut_id], wt_profile) for db_mut_id in  all_d[db_orf_id]} for allele_profile, db_orf_id, wt_profile in zip(edg_df.ad_profile, edg_df.db_orf_id, edg_df.wt_score)]

    return edg_df
    


##################################################################
# main function

def load_qY2H(dataset='normalized_df'):
    '''
    return qY2H data for variants and PRS RRSv2

    datasets:
    - retest_df: raw data
    - normalized_df: normalized data (based on empty AD controls), includes hsPRSRRSv2 (using db_mut_id = 999999)
    - edgotyping_profile_df: profile similarity, WT and allele similarity, does not include hsPRSRRSv2

    All data has been sequence confirmed (ORF level post pairwise test, variant level post mutagenesis)
    '''
    
    if dataset == 'retest_df':
        res = load_raw_results()

    elif dataset == 'normalized_df':
        retest_df = load_raw_results()
        res = get_norm_score_df(retest_df)

    elif dataset == 'edgotyping_profile_df':
        retest_df = load_raw_results()
        normalized_df = get_norm_score_df(retest_df)
        res = get_edg_profile_df(normalized_df.query("valid == True"))

    else:
        res = None
        print(f'error: wrong dataset name, {dataset} does not exist')

    return res

if __name__ == "__main__":      
    x = load_qY2H(dataset='normalized_df', include_hsPRSRRSv2=False)
    print(x.query("db_mut_id == 999999").head())


