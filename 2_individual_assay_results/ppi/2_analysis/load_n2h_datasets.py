#script adapted for dropbox

from TableToDict import GetDFv2
import glob
import xlrd
import pandas as pd
import numpy as np

# def get_path():
#     '''
#     get path to edgotyping_shared folder path
#     '''
#     from os import getcwd
#     current = getcwd()
#     def get_idx():
#         for i,e in enumerate(current.split('/')):
#             if e == 'Edgotyping_shared':
#                 idx_drive = i
#         return idx_drive
#     idx_drive = get_idx()
#     path = '/'.join(e for i,e in enumerate(current.split('/')) if i <= idx_drive)
#     return(path)


#YEAST data
def load_yN2H_apwt(dropbox_path):
    res = {}
    path = get_path()
    #for f in glob.glob('../experiment_data_processing/yN2H_experiment/20210331_yN2H_APWT/20210329_VUS_yN2Hapwt/*.xls'):
    for f in glob.glob(f'{dropbox_path}/experiment_data_processing/yN2H_experiment/20210331_yN2H_APWT/20210329_VUS_yN2Hapwt/*.xls'):
        

        retest_pla = int(f.split('_')[-4].split('Plate')[-1])

        if  retest_pla < 38:
            replicate = 1
        else:
            retest_pla -= 37
            replicate = 2

        workbook = xlrd.open_workbook(f, encoding_override='iso-8859-1') 
        df = pd.read_excel(workbook, skiprows=10, skipfooter=35)
        df.rename(columns={c:i+1 for i,c in enumerate(df.columns)}, inplace=True)

        value_d = dict()
        for l, (i, row) in zip(['A','B','C','D','E','F','G','H'], df.iterrows()):
            for i, value in enumerate(row):
                key = (retest_pla,replicate,f'{l}{str(i+1).zfill(2)}')
                res[key] = value
    return res

def get_yn2h_apwt_layout(dropbox_path):
    import pandas as pd
    import numpy as np
    import os


    apwt_df = pd.read_csv(f'{dropbox_path}/large_data_files/ccsb_database_useful_tables/vus/retest.tsv', sep='\t', low_memory=False)
    apwt_df = apwt_df[apwt_df.standard_batch == 'yN2Hapwt'].copy()

    rc4_mut_ids = pd.read_csv(f'{dropbox_path}/large_data_files/ccsb_database_useful_tables/vus/N2H_validation_nodes.tsv', sep='\t')
    rc4_mut_ids = rc4_mut_ids[rc4_mut_ids.collection == 'rc4'].copy()
    rc4_mut_ids = set(zip(rc4_mut_ids.orf_id, rc4_mut_ids.mut_id))


    col = ['ad_orf_id', 'ad_symbol', 'db_orf_id','db_symbol', 'db_mut_id', 'variant_id','retest_pla','retest_pos',
           'N2H_raw_r01', 'N2H_raw_r02']

    apwt_df = apwt_df[col].copy()

    empty_well_df = pd.read_excel(f'{dropbox_path}/experiment_data_processing/yN2H_experiment/20210331_yN2H_APWT/20210328_yN2H_empty_wells.xlsx',sheet_name='Sheet2')
    empty_well_df['replicate'] = empty_well_df.Plates.apply(lambda x: int(x.split('_')[2][-1]))
    empty_well_df['plate'] = empty_well_df.Plates.apply(lambda x: int(x.split('_')[-1].replace('Plate','')))
    empty_well_df['empty_wells'] = empty_well_df['Empty wells'].apply(lambda x: {f"{e.replace(' ','')[0]}{e.replace(' ','')[1:].zfill(2)}" for e in x.split(',')})
    empty_well_set =  {(plate, rep, well) for plate, rep, wells in zip(empty_well_df.plate, empty_well_df.replicate, empty_well_df.empty_wells) for well in wells }

    apwt_df['empty_well_rep1'] = [(db_orf_id, db_mut_id) in rc4_mut_ids or (plate,1,retest_pos) in empty_well_set for db_orf_id, db_mut_id, plate, retest_pos in zip(apwt_df.db_orf_id, apwt_df.db_mut_id, apwt_df.retest_pla, apwt_df.retest_pos)]
    apwt_df['empty_well_rep2'] = [(db_orf_id, db_mut_id) in rc4_mut_ids or (plate,2,retest_pos) in empty_well_set for db_orf_id, db_mut_id, plate, retest_pos in zip(apwt_df.db_orf_id, apwt_df.db_mut_id, apwt_df.retest_pla, apwt_df.retest_pos)]

    return apwt_df

def get_yn2h_nlr_df(dropbox_path):
    import pandas as pd
    import numpy as np
    apwt_df = get_yn2h_apwt_layout(dropbox_path)

    empty_N1_rep1_d = {(db,db_mut_id) if db_mut_id == db_mut_id else (db,999999):score for db, db_mut_id, score in zip(apwt_df[(apwt_df.variant_id == 'empty_N1')&(apwt_df.empty_well_rep1 == False)].db_orf_id, 
                                                                                                                       apwt_df[(apwt_df.variant_id == 'empty_N1')&(apwt_df.empty_well_rep1 == False)].db_mut_id, 
                                                                                                                       apwt_df[(apwt_df.variant_id == 'empty_N1')&(apwt_df.empty_well_rep1 == False)].N2H_raw_r01)}
    empty_N2_rep1_d = dict(zip(apwt_df[(apwt_df.variant_id == 'empty_N2')&(apwt_df.empty_well_rep1 == False)].ad_orf_id, apwt_df[(apwt_df.variant_id == 'empty_N2')&(apwt_df.empty_well_rep1 == False)].N2H_raw_r01))

    empty_N1_rep2_d = {(db,db_mut_id) if db_mut_id == db_mut_id else (db,999999):score for db, db_mut_id, score in zip(apwt_df[(apwt_df.variant_id == 'empty_N1')&(apwt_df.empty_well_rep2 == False)].db_orf_id, 
                                                                                                                       apwt_df[(apwt_df.variant_id == 'empty_N1')&(apwt_df.empty_well_rep2 == False)].db_mut_id, 
                                                                                                                       apwt_df[(apwt_df.variant_id == 'empty_N1')&(apwt_df.empty_well_rep2 == False)].N2H_raw_r02)}
    empty_N2_rep2_d = dict(zip(apwt_df[(apwt_df.variant_id == 'empty_N2')&(apwt_df.empty_well_rep2 == False)].ad_orf_id, apwt_df[(apwt_df.variant_id == 'empty_N2')&(apwt_df.empty_well_rep2 == False)].N2H_raw_r02))
     

    yN2H_normalized = apwt_df[apwt_df.variant_id.isin(['wt', 'allele', 'RRS','PRS'])].copy()

    yN2H_normalized['db_orf_mut'] = [(orf_id,mut_id) if mut_id == mut_id else (orf_id,999999) for orf_id,mut_id in zip(yN2H_normalized.db_orf_id, yN2H_normalized.db_mut_id)]
    yN2H_normalized['empty_N1_rep1'] = yN2H_normalized.db_orf_mut.map(empty_N1_rep1_d)
    yN2H_normalized['empty_N2_rep1'] = yN2H_normalized.ad_orf_id.map(empty_N2_rep1_d)
    yN2H_normalized['empty_N1_rep2'] = yN2H_normalized.db_orf_mut.map(empty_N1_rep2_d)
    yN2H_normalized['empty_N2_rep2'] = yN2H_normalized.ad_orf_id.map(empty_N2_rep2_d)
    yN2H_normalized.drop(columns='db_orf_mut',inplace=True)

    yN2H_normalized['max_ctrl_r1'] = [max([en1,en2]) for en1, en2 in zip(yN2H_normalized.empty_N1_rep1, yN2H_normalized.empty_N2_rep1)]
    yN2H_normalized['max_ctrl_r2'] = [max([en1,en2]) for en1, en2 in zip(yN2H_normalized.empty_N1_rep2, yN2H_normalized.empty_N2_rep2)]

    yN2H_normalized['valid_r1'] = [ew == False and  en1 == en1 and en2 == en2 for ew, en1, en2 in zip(yN2H_normalized.empty_well_rep1, yN2H_normalized.empty_N1_rep1, yN2H_normalized.empty_N2_rep1)]
    yN2H_normalized['valid_r2'] = [ew == False and  en1 == en1 and en2 == en2 for ew, en1, en2 in zip(yN2H_normalized.empty_well_rep2, yN2H_normalized.empty_N1_rep2, yN2H_normalized.empty_N2_rep2)]

    yN2H_normalized['NLR1'] = [np.nan if not valid_r1 else N2H_raw_r01/max_ctrl_r1 for N2H_raw_r01,valid_r1,max_ctrl_r1 in zip(yN2H_normalized.N2H_raw_r01, yN2H_normalized.valid_r1, yN2H_normalized.max_ctrl_r1)]
    yN2H_normalized['NLR2'] = [np.nan if not valid_r2 else N2H_raw_r02/max_ctrl_r2 for N2H_raw_r02,valid_r2,max_ctrl_r2 in zip(yN2H_normalized.N2H_raw_r02, yN2H_normalized.valid_r2, yN2H_normalized.max_ctrl_r2)]

    mean_nlr_rrs_r1 = np.mean(yN2H_normalized[(yN2H_normalized.valid_r1)&(yN2H_normalized.variant_id == 'RRS')].NLR1)
    mean_nlr_rrs_r2 = np.mean(yN2H_normalized[(yN2H_normalized.valid_r1)&(yN2H_normalized.variant_id == 'RRS')].NLR2)
    std_nlr_rrs_r1 = np.std(yN2H_normalized[(yN2H_normalized.valid_r1)&(yN2H_normalized.variant_id == 'RRS')].NLR1)
    std_nlr_rrs_r2 = np.std(yN2H_normalized[(yN2H_normalized.valid_r1)&(yN2H_normalized.variant_id == 'RRS')].NLR2)

    d = {0:1,1:0}
    yN2H_normalized['repeat_ratio'] = [[nlr1,nlr2][i%2]/[nlr1,nlr2][d[(i%2)]] for i,(nlr1, nlr2) in enumerate(zip(yN2H_normalized.NLR1, yN2H_normalized.NLR2))]
    yN2H_normalized['repeat_ratio_log2'] = yN2H_normalized.repeat_ratio.apply(np.log2)

    #yN2H_normalized['avg_NLR'] = [np.mean([e for e in [nlr1, nlr2] if e == e]) if nlr1 == nlr1 or nlr2 == nlr2 else np.nan for nlr1,nlr2 in zip(yN2H_normalized.NLR1, yN2H_normalized.NLR2)]
    yN2H_normalized['avg_NLR'] = [np.mean([nlr1, nlr2]) for nlr1,nlr2 in zip(yN2H_normalized.NLR1, yN2H_normalized.NLR2)]

    return yN2H_normalized


def get_yN2H_allele_df(dropbox_path):
    import pandas as pd
    import numpy as np
    yN2H_normalized = get_yn2h_nlr_df(dropbox_path)
    allele_df = yN2H_normalized.query("variant_id == 'allele' and ~avg_NLR.isnull()").copy()
    df = yN2H_normalized.query("variant_id == 'wt'")

    wt_NLR_d = {(r.db_orf_id, r.ad_orf_id):r.avg_NLR for i,r in df.iterrows()}
    allele_df['avg_WT_NLR'] = [wt_NLR_d[N2,N1] for N2,N1 in zip(allele_df.db_orf_id, allele_df.ad_orf_id)]

    allele_df = allele_df[~allele_df.avg_WT_NLR.isnull()].copy() #remove only one allele-WT
    allele_df['AWT_ratio'] = allele_df.avg_NLR / allele_df.avg_WT_NLR

    df = yN2H_normalized.query("variant_id == 'RRS'")
    l = [np.mean([nlr1,nlr2]) for nlr1,nlr2 in zip(df.NLR1, df.NLR2) if nlr1 == nlr1 and nlr2 == nlr2]
    std_rrs, mean_rrs = np.std(l), np.mean(l)
    allele_df['z_WT_NLR'] = allele_df.avg_WT_NLR.apply(lambda x: (x-mean_rrs)/std_rrs)
    allele_df['WT_detected'] = allele_df.z_WT_NLR.apply(lambda x: x >= 2.23)

    std_ratio = np.std(yN2H_normalized.repeat_ratio_log2)
    mean_ratio = np.mean(yN2H_normalized.repeat_ratio_log2)
    allele_df['z_perturbation'] = allele_df.AWT_ratio.apply(lambda x: (np.log2(x)-mean_ratio)/std_ratio)
    allele_df['z_cat'] = ['not_detected' if not wt_detected else 'increase_ppi' if zp>=2 else 'decrease_ppi' if zp<=-2 else 'no_perturbation' for zp,wt_detected in zip(allele_df.z_perturbation, allele_df.WT_detected)] 

    return allele_df


#MAMMALIAN N2H
def get_mn2h_normalized_df(dropbox_path):
    xp1_mapping_d = {1:(1,1),   #plate number (1-69) (reader): retest_pla (1-23), replicate number (1-2-3)
                     2:(1,2),
                     3:(1,3),
                     4:(2,1),
                     5:(2,2),
                     6:(2,3),
                     7:(3,1),
                     8:(3,2),
                     9:(3,3),
                     10:(4,1),
                     11:(4,2),
                     12:(4,3),
                     13:(5,1),
                     14:(5,2),
                     15:(5,3),
                     16:(6,1),
                     17:(6,2),
                     18:(6,3),
                     19:(7,1),
                     20:(7,2),
                     21:(7,3),
                     22:(8,1),
                     23:(8,2),
                     24:(8,3),
                     25:(9,1),
                     26:(9,2),
                     27:(9,3),
                     28:(10,1),
                     29:(10,2),
                     30:(10,3),
                     31:(11,1),
                     32:(11,2),
                     33:(11,3),
                     34:(12,1),
                     35:(12,2),
                     36:(12,3),
                     37:(13,1),
                     38:(13,2),
                     39:(13,3),
                     40:(14,1),
                     41:(14,2),
                     42:(14,3),
                     43:(15,1),
                     44:(15,2),
                     45:(15,3),
                     46:(16,1),
                     47:(16,2),
                     48:(16,3),
                     49:(17,1),
                     50:(17,2),
                     51:(17,3),
                     52:(18,1),
                     53:(18,2),
                     54:(18,3),
                     55:(19,1),
                     56:(19,2),
                     57:(19,3),
                     58:(20,1),
                     59:(20,2),
                     60:(20,3),
                     61:(21,1),
                     62:(21,2),
                     63:(21,3),
                     64:(22,1),
                     65:(22,2),
                     66:(22,3),
                     67:(23,1),
                     68:(23,2),
                     69:(23,3)
                    }

    xp2_mapping_d = {1:(1,1),
                     2:(1,2),
                     3:(1,3),
                     4:(2,1),
                     5:(2,2),
                     6:(2,3),
                     7:(3,1),
                     8:(3,2),
                     9:(3,3),
                     10:(4,1),
                     11:(4,2),
                     12:(4,3),
                     13:(7,1),
                     14:(7,2),
                     15:(7,3),
                     16:(6,1),
                     17:(6,2),
                     18:(6,3),
                     19:(5,1),
                     20:(5,2),
                     21:(5,3),
                     22:(9,1),
                     23:(9,2),
                     24:(9,3),
                     25:(10,1),
                     26:(10,2),
                     27:(10,3),
                     28:(11,1),
                     29:(11,2),
                     30:(11,3),
                     31:(12,1),
                     32:(12,2),
                     33:(12,3),
                     34:(13,1),
                     35:(13,2),
                     36:(13,3),
                     37:(14,1),
                     38:(14,2),
                     39:(14,3),
                     40:(15,1),
                     41:(15,2),
                     42:(15,3),
                     43:(16,1),
                     44:(16,2),
                     45:(16,3),
                     46:(17,1),
                     47:(17,2),
                     48:(17,3),
                     49:(18,1),
                     50:(18,2),
                     51:(18,3),
                     52:(19,1),
                     53:(19,2),
                     54:(19,3),
                     55:(20,1),
                     56:(20,2),
                     57:(20,3),
                     58:(21,1),
                     59:(21,2),
                     60:(21,3),
                     61:(22,1), #plate is missing
                     62:(22,2),
                     63:(22,3),
                     64:(23,1),
                     65:(23,2),
                     66:(23,3),
                     67:(8,1),
                     68:(8,2),
                     69:(8,3)
                    }

    #xp1_mapping_d pla:pla_id,repeat
    #xp2_mapping_d pla:pla_id,repeat

    res = {'xp1':{1:{},2:{},3:{}},
           'xp2':{1:{},2:{},3:{}}
          }

    for rep_id in ['1','2']:
        for i,elem in enumerate(glob.glob(f'{dropbox_path}/experiment_data_processing/mN2H_experiment/20210823_vus_apwt/2022_APWT_data/N2H_raw_replicat0{rep_id}/*')):
            pla = int(elem.split('-')[-1][1:-4])
            retest_pla, repeat = xp1_mapping_d[pla] if rep_id == '1' else xp2_mapping_d[pla]

            wb = xlrd.open_workbook(elem, encoding_override='ascii')  #ascii  #'iso-8859-1'
            df = pd.read_excel(wb, skiprows=24, skipfooter=0)
            if len(df) > 8:
                df = pd.read_excel(wb, skiprows=35, skipfooter=10)

            df.rename(columns={c:i+1 for i,c in enumerate(df.columns)}, inplace=True)

            value_d = dict()
            for l, (i, row) in zip(['A','B','C','D','E','F','G','H'], df.iterrows()):
                for i, value in enumerate(row):
                    key = (retest_pla,f'{l}{str(i+1).zfill(2)}')

                    res[f'xp{rep_id}'][repeat][key] = value

    col = ['ad_orf_id', 'ad_symbol', 'db_orf_id',
           'db_symbol', 'db_mut_id','variant_id', 'retest_pla', 'retest_pos', 
           'standard_batch']


    retest_df = pd.read_csv(f'{dropbox_path}/large_data_files/ccsb_database_useful_tables/vus/retest.tsv', sep='\t', low_memory=False)
    retest_df = retest_df[retest_df.standard_batch == 'mN2Hapwtv2'].copy()
    retest_df.rename(columns={c:c.replace('ad','N1').replace('db','N2') for c in retest_df.columns},inplace=True)

    for xp in ['xp1','xp2']:
        for rep in [1,2,3]:
            retest_df[f'mN2H_{xp}r{rep}'] = [res[xp][rep][pla,pos] if (pla,pos) in res[xp][rep] else np.nan for pla,pos in zip(retest_df.retest_pla, retest_df.retest_pos)]

    #discard measure if signal is 20% of those of 2 other replicates (it's never 5 fold more than the 2 other measures)    
    retest_df['mN2H_xp1r1'] = [np.nan if r1/r2<0.2 and r1/r3<0.2 else r1 for r1,r2,r3 in zip(retest_df.mN2H_xp1r1, retest_df.mN2H_xp1r2, retest_df.mN2H_xp1r3)]
    retest_df['mN2H_xp1r2'] = [np.nan if r2/r1<0.2 and r2/r3<0.2 else r2 for r1,r2,r3 in zip(retest_df.mN2H_xp1r1, retest_df.mN2H_xp1r2, retest_df.mN2H_xp1r3)]
    retest_df['mN2H_xp1r3'] = [np.nan if r3/r2<0.2 and r3/r1<0.2 else r3 for r1,r2,r3 in zip(retest_df.mN2H_xp1r1, retest_df.mN2H_xp1r2, retest_df.mN2H_xp1r3)]
    retest_df['mN2H_xp2r1'] = [np.nan if r1/r2<0.2 and r1/r3<0.2 else r1 for r1,r2,r3 in zip(retest_df.mN2H_xp2r1, retest_df.mN2H_xp2r2, retest_df.mN2H_xp2r3)]
    retest_df['mN2H_xp2r2'] = [np.nan if r2/r1<0.2 and r2/r3<0.2 else r2 for r1,r2,r3 in zip(retest_df.mN2H_xp2r1, retest_df.mN2H_xp2r2, retest_df.mN2H_xp2r3)]
    retest_df['mN2H_xp2r3'] = [np.nan if r3/r2<0.2 and r3/r1<0.2 else r3 for r1,r2,r3 in zip(retest_df.mN2H_xp2r1, retest_df.mN2H_xp2r2, retest_df.mN2H_xp2r3)]

    xp_errors = {1:{(5,'E03')},
                 2:{(5,'E03'),
                    (13,'A05'),
                    (13,'A06'),
                    (13,'D06'),
                    (14,'C01'),
                    (14,'B04'),
                    (14,'E06'),
                    (14,'H06'),
                    (14,'B08'),
                    (14,'A10'),
                   }
                }

    retest_df['error1'] = [(p,w) in [(pla,well) for pla,well in xp_errors[1]] for p,w in zip(retest_df.retest_pla, retest_df.retest_pos)]
    retest_df['error2'] = [(p,w) in [(pla,well) for pla,well in xp_errors[2]] for p,w in zip(retest_df.retest_pla, retest_df.retest_pos)]

    retest_df['xp1'] = [np.mean([e for e in [r1,r2,r3] if e == e]) for r1,r2,r3 in zip(retest_df.mN2H_xp1r1, retest_df.mN2H_xp1r2, retest_df.mN2H_xp1r3)]
    retest_df['xp2'] = [np.mean([e for e in [r1,r2,r3] if e == e]) for r1,r2,r3 in zip(retest_df.mN2H_xp2r1, retest_df.mN2H_xp2r2, retest_df.mN2H_xp2r3)]

    retest_df['xp1'] = [xp1 if not error else np.nan for xp1,error in zip(retest_df.xp1, retest_df.error1)]
    retest_df['xp2'] = [xp2 if not error else np.nan for xp2,error in zip(retest_df.xp2, retest_df.error2)]

    empty_N1_xp1_d = dict(zip(zip(retest_df[retest_df.variant_id == 'empty_N1'].N2_orf_id,retest_df[retest_df.variant_id == 'empty_N1'].N2_mut_id),retest_df[retest_df.variant_id == 'empty_N1'].xp1))
    empty_N2_xp1_d = dict(zip(retest_df[retest_df.variant_id == 'empty_N2'].N1_orf_id, retest_df[retest_df.variant_id == 'empty_N2'].xp1))

    empty_N1_xp2_d = dict(zip(zip(retest_df[retest_df.variant_id == 'empty_N1'].N2_orf_id,retest_df[retest_df.variant_id == 'empty_N1'].N2_mut_id),retest_df[retest_df.variant_id == 'empty_N1'].xp2))
    empty_N2_xp2_d = dict(zip(retest_df[retest_df.variant_id == 'empty_N2'].N1_orf_id, retest_df[retest_df.variant_id == 'empty_N2'].xp2))

    normalized_df = retest_df[retest_df.variant_id.isin(['wt','RRS','PRS', 'allele'])].copy()

    normalized_df.drop(columns=['error1','error2','mN2H_xp1r1','mN2H_xp1r2','mN2H_xp1r3','mN2H_xp2r1','mN2H_xp2r2','mN2H_xp2r3'], inplace=True)

    normalized_df['empty_N1_xp1'] = normalized_df.apply(lambda x: empty_N1_xp1_d[x.N2_orf_id,x.N2_mut_id], axis=1)
    normalized_df['empty_N2_xp1'] = normalized_df.N1_orf_id.map(empty_N2_xp1_d)

    normalized_df['empty_N1_xp2'] = normalized_df.apply(lambda x: empty_N1_xp2_d[x.N2_orf_id,x.N2_mut_id], axis=1)
    normalized_df['empty_N2_xp2'] = normalized_df.N1_orf_id.map(empty_N2_xp2_d)

    normalized_df['max_ctrl_r1'] = [max([en1,en2]) for en1, en2 in zip(normalized_df.empty_N1_xp1, normalized_df.empty_N2_xp1)]
    normalized_df['max_ctrl_r2'] = [max([en1,en2]) for en1, en2 in zip(normalized_df.empty_N1_xp2, normalized_df.empty_N2_xp2)]

    normalized_df['NLR1'] = [N1N2/max([N1F2,N2F1]) for N1N2, N1F2, N2F1 in zip(normalized_df.xp1, normalized_df.empty_N1_xp1, normalized_df.empty_N2_xp1)]
    normalized_df['NLR2'] = [N1N2/max([N1F2,N2F1]) for N1N2, N1F2, N2F1 in zip(normalized_df.xp2, normalized_df.empty_N1_xp2, normalized_df.empty_N2_xp2)]

    d = {0:1,1:0}
    normalized_df['repeat_ratio'] = [[nlr1,nlr2][i%2]/[nlr1,nlr2][d[(i%2)]] for i,(nlr1, nlr2) in enumerate(zip(normalized_df.NLR1, normalized_df.NLR2))]
    normalized_df['repeat_ratio_log2'] = normalized_df.repeat_ratio.apply(np.log2)

    normalized_df['avg_NLR'] = [np.mean([nlr1, nlr2]) for nlr1,nlr2 in zip(normalized_df.NLR1, normalized_df.NLR2)]


    return normalized_df


def get_mn2h_allele_df(dropbox_path):
    normalized_df = get_mn2h_normalized_df(dropbox_path)
    #mean NLR
    #allele-wt ratio
    #z-score perturbation

    #two metrics: detection (above RRS, noise) + perturbation quantification allele/wt zscore
    #NA if not 


    allele_df = normalized_df.query("variant_id == 'allele' and ~avg_NLR.isnull()").copy()
    #allele_df['nt_change'] = [mut_d[orf_id,mut_id][0] for orf_id, mut_id in zip(allele_df.N1_orf_id, allele_df.N1_mut_id)]
    #allele_df['aa_change'] = [mut_d[orf_id,mut_id][1] for orf_id, mut_id in zip(allele_df.N1_orf_id, allele_df.N1_mut_id)]
    #allele_df['clinical_significance'] = [mut_d[orf_id,mut_id][2] for orf_id, mut_id in zip(allele_df.N1_orf_id, allele_df.N1_mut_id)]

    df = normalized_df.query("variant_id == 'wt'")
    wt_NLR_d = {(r.N2_orf_id, r.N1_orf_id):r.avg_NLR for i,r in df.iterrows()}
    
    allele_df['avg_WT_NLR'] = [wt_NLR_d[N2,N1] for N2,N1 in zip(allele_df.N2_orf_id, allele_df.N1_orf_id)]
    allele_df = allele_df[~allele_df.avg_WT_NLR.isnull()].copy() #remove only one allele-WT
    allele_df['AWT_ratio'] = allele_df.avg_NLR / allele_df.avg_WT_NLR

    df = normalized_df.query("variant_id == 'RRS' and avg_NLR < 800")
    l = [np.mean([nlr1,nlr2]) for nlr1,nlr2 in zip(df.NLR1, df.NLR2) if nlr1 == nlr1 and nlr2 == nlr2]
    std_rrs, mean_rrs = np.std(l), np.mean(l)
    allele_df['z_WT_NLR'] = allele_df.avg_WT_NLR.apply(lambda x: (x-mean_rrs)/std_rrs)
    allele_df['WT_detected'] = allele_df.z_WT_NLR.apply(lambda x: x >= 2.23)

    std_ratio = np.std(normalized_df.repeat_ratio_log2)
    mean_ratio = np.mean(normalized_df.repeat_ratio_log2)
    allele_df['z_perturbation'] = allele_df.AWT_ratio.apply(lambda x: (np.log2(x)-mean_ratio)/std_ratio)
    #allele_df['z_cat'] = allele_df.z_perturbation.apply(lambda x: 'increase_ppi' if x>=3 else 'decrease_ppi' if x<=-3 else 'no_perturbation')
    allele_df['z_cat'] = ['not_detected' if not wt_detected else 'increase_ppi' if zp>=2 else 'decrease_ppi' if zp<=-2 else 'no_perturbation' for zp,wt_detected in zip(allele_df.z_perturbation, allele_df.WT_detected)] 

    return allele_df

if __name__ == "__main__":   
    res = load_yN2H_apwt(dropbox_path)
    yN2H_normalized = get_yn2h_nlr_df(dropbox_path)
    allele_df = get_yN2H_allele_df(dropbox_path)

    #normalized_df = get_mn2h_normalized_df(dropbox_path)
    #allele_df = get_mn2h_allele_df()

    print(yN2H_normalized[(yN2H_normalized.db_mut_id==0)&(yN2H_normalized.db_symbol != 'UBQLN2')].columns)







