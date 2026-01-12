data alzheimer25;
    set "~/Project 3 LDA/alzheimer25.sas7bdat"; 
run;

proc means data=alzheimer25 noprint;
    var AGE BMI;
    output out=stats mean=mean_age mean_bmi std=sd_age sd_bmi;
run;

proc mi data=alzheimer25 seed=1234 out=alzheimer25_mi simple nimpute=20 round=1;
var age bmi sex adl taupet0 taupet1 taupet2 taupet3 taupet4 taupet5 taupet6
abpet0 abpet1 abpet2 abpet3 abpet4 abpet5 abpet6
cdrsb0 cdrsb1 cdrsb2 cdrsb3 cdrsb4 cdrsb5 cdrsb6;
/*by patid;*/
run;

data alzheimer_long;
    set alzheimer25_mi; 
    
    array cdrsb_arr[0:6] cdrsb0-cdrsb6;
    array bprs_arr[0:6] bprs0-bprs6;
    array abpet_arr[0:6] abpet0-abpet6;
    array taupet_arr[0:6] taupet0-taupet6;
    
    do TIME = 0 to 6;
        CDRSB  = cdrsb_arr[TIME];
        BPRS   = bprs_arr[TIME];
        ABPET  = abpet_arr[TIME];
        TAUPET = taupet_arr[TIME];

        if missing(CDRSB) then R = 0; 
        else R = 1;

        if R = 1 then do;
            if CDRSB < 10 then CDRSB_CAT = 0;
            else CDRSB_CAT = 1;
        end;
        else CDRSB_CAT = .; 

        TIMECLSS = put(TIME, 1.);
        
        
        output; 
    end;
    
    drop cdrsb0 cdrsb1 cdrsb2 cdrsb3 cdrsb4 cdrsb5 cdrsb6
    abpet0 abpet1 abpet2 abpet3 abpet4 abpet5 abpet6
    taupet0 taupet1 taupet2 taupet3 taupet4 taupet5 taupet6;
run;

proc means data=ALZHEIMER_LONG noprint;
    var ABPET TAUPET AGE BMI; 
    output out=bio_stats mean=m_ab m_tau m_age m_bmi std=s_ab s_tau s_age s_bmi;
run;

data alzheimer_long_centered;
    if _n_=1 then set bio_stats;
    set alzheimer_long;
    
    ABPET_STD = (ABPET - m_ab) / s_ab;
    TAUPET_STD = (TAUPET - m_tau) / s_tau;
    AGE_STD = (AGE - m_age) / s_age;
    BMI_STD = (BMI - m_bmi) / s_bmi;
run;

ods trace on;

proc genmod data = alzheimer_long_centered descending;
    class PATID SEX(ref = '1') TIMECLSS / param = ref;
    by _imputation_;
    model CDRSB_CAT = TIME TIME * TIME SEX AGE_STD BMI_STD ADL TAUPET_STD ABPET_STD 
    TIME * SEX TIME * AGE_STD TIME * BMI_STD TIME * ADL TIME * TAUPET_STD TIME * ABPET_STD 
    TIME * TIME * SEX TIME * TIME * AGE_STD TIME * TIME * BMI_STD TIME * TIME * ADL 
    TIME * TIME * TAUPET_STD TIME * TIME * ABPET_STD / 
    dist=binomial link=logit;
    repeated subject=PATID / withinsubject=TIMECLSS type=UN modelse COVB;
    ods output GEEEmpPEst = gmparms 
               ParmInfo   = gmpinfo 
               GEERCov    = gmcovb;
run;


data gmpinfo_clean;
    set gmpinfo;
    if parameter = 'Prm5' then delete; 
run;

/* 3. Pooling senza Intercept */
proc mianalyze parms=gmparms covb=gmcovb parminfo=gmpinfo wcov bcov tcov;
    modeleffects Intercept TIME TIME * TIME SEX AGE_STD BMI_STD ADL TAUPET_STD ABPET_STD 
    TIME * SEX TIME * AGE_STD TIME * BMI_STD TIME * ADL TIME * TAUPET_STD TIME * ABPET_STD 
    TIME * TIME * SEX TIME * TIME * AGE_STD TIME * TIME * BMI_STD TIME * TIME * ADL 
    TIME * TIME * TAUPET_STD TIME * TIME * ABPET_STD ;
    ods output ParameterEstimates=mianalyze_results;
run;
