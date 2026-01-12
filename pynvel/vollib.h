typedef struct {
  int evod;
  int opt;
  float maxlen;
  float minlen;
  float minlent;
  float merchl;
  float mtopp;
  float mtops;
  float stump;
  float trim;
  float btr;
  float dbtbh;
  float minbfd;
  char * cor;
} merchrules_;

extern void vernum_(int *v);

extern void voleqdef_(
  char* var, int* region, char* forest, char* district
  , int* species, char* product, char* vol_eq, int* err_flag
  , int vl, int fl, int dl, int pl, int el
);

extern void fiavoleqdef_(
  char* var, int *region, char* forest, char* district
  , int *species, char* vol_eq, int *err_flag
  , int vl, int fl, int dl, int el
);

extern void nvb_defaulteq_(
  int *region, char* forest, char* district
  , int *species, char* vol_eq, int *err_flag
  , int fl, int dl, int el
);

extern void ht2topd_(
  int *regn, char* forsti, char* voleqi
  , float *dbhob, float *httot, float *ht1prd, float *ht2prd
  , float *upsht1, float *upsht2, float *upsd1, float *upsd2
  , float *avgz1, float *avgz2, int *htref, float *dbtbh, float *btr
  , int *fclass, float *stemdib, float *stemht, int *errflag
);

extern void calcdia2_(
  int *regn, char* forsti, char* voleqi, float *stump_ht
  , float *dbhob, float *drcob, float *httot
  , float *upsht1, float *upsht2, float *upsd1, float *upsd2
  , int *htref, float *avgz1, float *avgz2, int *fclass
  , float *dbtbh, float *btr
  , float *stemht, float *stemdib, float *stemdob
  , int *errflag
);

extern void scrib_();

extern void volinitnvb_(
  int *region
  , char* forest
  , char* vol_eq
  , float *min_top_prim
  , float *min_top_sec
  , float *stump_ht
  , float *dbh_ob
  , float *drc_ob
  , char* ht_type
  , float *total_ht
  , int *ht_log
  , float *ht_prim
  , float *ht_sec
  , float *upper_ht1
  , float *upper_ht2
  , float *upper_diam1
  , float *upper_diam2
  , int *ht_ref
  , float *avg_z1
  , float *avg_z2
  , int *form_class
  , float *bark_2thick
  , float *bark_ratio

  , float *crown
  , float *cull
  , int *decaycd

  , float *volume
  , float *log_volume
  , float *log_diam
  , float *log_len
  , float *bol_ht
  , int *num_logs
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
  , int prod_len, int live_len, int ctype_len);