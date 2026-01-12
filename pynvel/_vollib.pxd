"""
External references to the NVEL/VOLLIB library
"""

import numpy as np
cimport numpy as np

cdef extern from "vollib.h":
    # vernum.f
    void vernum_(int *v)

    void scrib_(float *dia, float *len, char* cor, float *vol, int cl)

    # Return the default equation for a species
    # voleqdef.f
    void voleqdef_(char* var, int* region, char* forest, char* district
            , int* species, char* product, char* vol_eq, int* err_flag
            , int vl, int fl, int dl, int pl, int el)

    # Return the FIA default equation for a species
    # voleqdef.f
    void fiavoleqdef_(char* var, int *region, char* forest, char* district
            , int *species, char* vol_eq, int *err_flag
            , int vl, int fl, int dl, int el)

    # Call the NSVB equation lookup
            # Implemented in nsvb.f
    void nvb_defaulteq_(int *region, char* forest, char* district
            , int *species, char* vol_eq, int *err_flag
            , int fl, int dl, int el)

    # vollibcs.f
    void vollibc2_(
            int *regn, char* forsti, char* voleqi, float *mtopp, float *mtops
            , float *stump, float *dbhob, float *drcob, char* httypei
            , float *httot, int *htlog, float *ht1prd, float *ht2prd
            , float *upsht1, float *upsht2, float *upsd1, float *upsd2
            , int *htref, float *avgz1, float *avgz2, int *fclass
            , float *dbtbh, float *btr
            , int *i3, int *i7, int *i15, int *i20, int *i21
            , float *vol, float *logvoli, float *logdiai
            , float *loglen, float *bolht, int *tlogs
            , float *nologp, float *nologs
            , int *cutflg, int *bfpflg, int *cupflg, int *cdpflg, int *spflg
            , char* conspeci, char* prodi, int *httfll, char* livei
            , int *ba, int *si, char* ctypei, int *errflag, int *indeb
            , int *pmtflg, merchrules_ *merrules
            , int forsti_len, int voleqi_len, int httypei_len, int conspeci_len
            , int prodi_len, int livei_len, int ctypei_len
            )

    # Calculate tree volume according to user defined merchandizing rules.
    # volinit2.f
    void volinit2_(
            int *regn, char* forsti, char* voleqi, float *mtopp, float *mtops
            , float *stump, float *dbhob, float *drcob, char* httypei
            , float *httot, int *htlog, float *ht1prd, float *ht2prd
            , float *upsht1, float *upsht2, float *upsd1, float *upsd2
            , int *htref, float *avgz1, float *avgz2, int *fclass
            , float *dbtbh, float *btr
            , int *i3, int *i7, int *i15, int *i20, int *i21
            , float *vol, float *logvoli, float *logdiai
            , float *loglen, float *bolht, int *tlogs
            , float *nologp
            , float *nologs
            , int *cutflg
            , int *bfpflg
            , int *cupflg
            , int *cdpflg
            , int *spflg
            , char* conspeci
            , char* prodi
            , int *httfll
            , char* livei
            , int *ba
            , int *si
            , char* ctypei
            , int *errflag
            , merchrules_ *merrules
            , int *idist
            , int forsti_len, int voleqi_len, int httypei_len, int conspeci_len
            , int prodi_len, int livei_len, int ctypei_len
            )

    ## TODO: Implement VOLINITNVB, it is the new top level NVEL subroutine
    ##       It accomodates merch rules and hands off to VOLINIT for non-FIA equations
    # VOLINITNVB(REGN,FORST,VOLEQ,MTOPP,MTOPS,STUMP,DBHOB,
    #  +    DRCOB,HTTYPE,HTTOT,HTLOG,HT1PRD,HT2PRD,UPSHT1,UPSHT2,UPSD1,
    #  +    UPSD2,HTREF,AVGZ1,AVGZ2,FCLASS,DBTBH,BTR,CR,CULL,DECAYCD,
    #  +    VOL,LOGVOL,LOGDIA,LOGLEN,BOLHT,TLOGS,NOLOGP,NOLOGS,CUTFLG,
    #  +    BFPFLG,CUPFLG,CDPFLG,SPFLG,CONSPEC,PROD,HTTFLL,LIVE,
    #  +    BA,SI,CTYPE,ERRFLAG,IDIST,BRKHT,BRKHTD,FIASPCD,DRYBIO,
    #  +    GRNBIO,MRULEFLG,MERRULES)
    void volinitnvb_(
            int *region, char* forest, char* vol_eq, float *min_top_prim, float *min_top_sec
            , float *stump_ht, float *dbh_ob, float *drc_ob, char* ht_type
            , float *total_ht, int *ht_log, float *ht_prim, float *ht_sec
            , float *upper_ht1, float *upper_ht2, float *upper_diam1, float *upper_diam2
            , int *ht_ref, float *avg_z1, float *avg_z2, int *form_class
            , float *bark_2thick, float *bark_ratio

            , float *crown, float *cull, int *decaycd

            , float *volume, float *log_volume, float *log_diam
            , float *log_len, float *bol_ht, int *num_logs
            , float *num_logs_prim
            , float *num_logs_sec
            , int *cubic_total_flag
            , int *bdft_prim_flag
            , int *cubic_prim_flag
            , int *cord_prim_flag
            , int *sec_vol_flag
            , char* con_spp
            , char* prod_code
            , int *ht_1st_limb
            , char* live
            , int *basal_area
            , int *site_index
            , char* calc_type
            , int *error_flag
            , int *district
            , float *brkht
            , float *brkhtd
            , int *fiaspcd
            , float *drybio
            , float *grnbio
            , int *mruleflag
            , merchrules_ *merrules

            , int forest_len, int voleq_len, int httype_len, int conspec_len
            , int prod_len, int live_len, int ctype_len
            )

    # Calculate the height to desired top dib
    # ht2topd.f
    void ht2topd_(
            int *regn, char* forsti, char* voleqi
            , float *dbhob, float *httot, float *ht1prd, float *ht2prd
            , float *upsht1, float *upsht2, float *upsd1, float *upsd2
            , float *avgz1, float *avgz2, int *htref, float *dbtbh, float *btr
            , int *fclass, float *stemdib, float *stemht, int *errflag)

    # calcdia.f
    void calcdia2_(
            int *regn, char* forsti, char* voleqi, float *stump_ht
            , float *dbhob, float *drcob, float *httot
            , float *upsht1, float *upsht2, float *upsd1, float *upsd2
            , int *htref, float *avgz1, float *avgz2, int *fclass
            , float *dbtbh, float *btr
            , float *stemht, float *stemdib, float *stemdob
            , int *errflag)

    # Compute biomass using the Jenkins component ratio methods
    # jenkins.f
    void jenkins(int *spec, float *dbhob, float *bioms)

# A struct to represent the merchrules defined type
cdef extern from "vollib.h":
    ctypedef struct merchrules_:
        int evod
        int opt
        float maxlen
        float minlen
        float minlent
        float merchl
        float mtopp
        float mtops
        float stump
        float trim
        float btr
        float dbtbh
        float minbfd
        char* cor #NOTE: **Must be single char, not pointer. Unsure how to pass a len arg in a struct
