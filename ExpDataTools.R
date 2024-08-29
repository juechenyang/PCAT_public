# Title     : TODO
# Objective : TODO
# Created by: JasonY
# Created on: 2019-07-25

chooseDB <- function(name){
    if(name=='target'){
        DB <- list(name='target',exp_db_name='target_expr',
            clinical_db_name='target_clin', clinical_id_name='patient_id',biospec_id_name='sample_id',
            disease_tag='Disease', cor_db_name='target_cor', mut_db_name='target_mut_bbk2', cpn_db_name='target_cpn', mut_info_db='target_mut_info_bbk',
            hover_string = c("sample_id","Variant_Type", "HGVSc", "HGVSp", "HGVSp_Short"), fusion_db='target_fusion_bbk', donor_name='donor_gene', acceptor_name='acceptor_gene',
            fusion_search_db = 'target_fusion_search_db'
        )
    }else if(name=='pptc'){
        DB <- list(name='pptc',exp_db_name='pptc_expr', clinical_db_name='pdx_clin', drug_db='pptc_drug_bk', drug_anno_db='pptc_drug_anno',
            biospec_id_name='sample_id', disease_tag='Histology', mut_db_name='pptc_mut_bbk', mut_info_db='pptc_mut_info_bbk',cpn_db_name='pptc_cpn_bk',
            hover_string = c("sample_id","Variant_Type","HGVSp_Short", "t_alt_count", "t_ref_count"), protein_change_identifier='HGVSp_Short',fusion_db='pptc_fusion',
            donor_name='donor_gene', acceptor_name='acceptor_gene', fusion_search_db = 'pptc_fusion_search_db'
        )
    }else if(name=='pbta'){
        DB = list(
            name='pbta',exp_db_name='pbta_expr', clinical_db_name='pbta_clin', biospec_id_name='sample_id', disease_tag='short_histology'
        )
    }
    return(DB)
}
get_sub_disease_clinical_features <- function(sub_disease){
    if(sub_disease=='TARGET-OS'){
        return(c('MetaSite','PrimTumorSite', 'Gender','Race','Ethnicity','AgeDxDays','DxStatus','SpecificTumorSite','SpecificTumorRegion','PrimSiteProgression','InitRelapseSite','HistResponse','NecrosisPerc','RelapseType'))
    }else if(sub_disease=='TARGET-RT'){
        return(c('AgeDxDays','StageCls','INI1Mut', 'Gender','Race','Ethnicity'))
    }else if(sub_disease=='TARGET-WT'){
        return(c('StageCls','Histology', 'Gender','Race','AgeDxDays','Ethnicity'))
    }else if(sub_disease=='TARGET-NBL'){
        return(c('StageCls','MYCNAmp','Histology','Grade','MKI', 'Gender','Race','Ethnicity',"AgeDxDays",'DxCat','RiskGrp','Ploidy'))
    }else if(sub_disease=='TARGET-ALL-P2'){
        return(c('Gender', 'Race', 'Ethnicity', 'AgeDxDays', 'WBCDx', 'CNS.Status.at.Diagnosis', 'Bone.Marrow.Site.of.Relapse', 'CNS.Site.of.Relapse', 'Testes.Site.of.Relapse',
        'Other.Site.of.Relapse', 'ETV6_RUNX1_Stat', 'Tris_4_10_Stat', 'MLL_Stat', 'TCF3_PBX1_Stat', 'BCR_ABL1_Stat', 'Down.Syndrome', 'Cell.of.Origin'))
    }else if(sub_disease=='TARGET-ALL-P3'){
        return(c('AgeDxDays','InitTherapy', 'ALALPres','WHO_ALALCls','WHO_Dx', 'Gender','WBCDx'))
    }else if(sub_disease=='TARGET-AML'){
        return(c('RiskGrp', 'GeneFusCls','Gender','Race','Ethnicity','AgeDxDays','WBCDx','BM_LeukBlast_Perc','Peripheral_Blast_Perc','CNSDisease',
        'Chloroma','FABCat','PrimCytoCode',
        'FLT3.ITD','FLT3_PM','NPMMut','CEBPAMut','WT1Mut','cKitEx8Mut','cKitEx17Mut','MRDCour1',
        'MRDCour2','CRStatCour1','CRStatCour2','SCT_in_1stCR','BoneMarrow_SiteRelapse.IndFail',
        'CNS_SiteRelapse.IndFail','Chloroma_SiteRelapse.IndFail','Cytogenetic_SiteRelapse.IndFail','Other_SiteRelapse.IndFail'))
    }else if(sub_disease=='pan-cancer'){
        return(c('Disease', 'Gender','Race','Ethnicity','AgeDxDays'))
    }
}
process_exp_data = function(DB, gene_name, mode){
    query_db <- DB$exp_db_name
    if(mode=='single'){
        exp_data <- QueryOnRow(query_db, gene_name)
    }else{
        exp_data <- QueryOnRow(query_db, gene_name)
    }
    wide_exp_data = tidyr::spread(exp_data, key='gene_name', value='expr')
    return(wide_exp_data)
}

