"""
sfh_depth_builder_empirical.py
Builds sparc_sfh_depth_catalog.csv from sparc_z0mgs_crossmatch.csv.
Layer 1: z0MGS empirical SFR/Mstar where z0mgs_SFR_Msun_yr exists.
Layer 2: delayed-tau fallback where empirical SFR is missing.
"""
import numpy as np
import pandas as pd
INPUT = "sparc_z0mgs_crossmatch.csv"
OUTPUT = "sparc_sfh_depth_catalog.csv"
EPS=1e-30; F_REM=0.15; T_OBS_GYR=12.0; N_GRID=800; TAU_MIN=1.0; TAU_MAX=9.0; TAU0_KERNEL=2.0; P_KERNEL=0.5

def delayed_tau_depth(Mstar, Mgas):
    Mstar=max(float(Mstar), EPS); Mgas=max(float(Mgas),0.0); fgas=Mgas/max(Mstar+Mgas,EPS)
    tau=TAU_MIN+(TAU_MAX-TAU_MIN)*fgas; t=np.linspace(0.0,T_OBS_GYR,N_GRID); dt=t[1]-t[0]
    shape=t*np.exp(-t/max(tau,EPS)); formed=np.trapezoid(shape, dx=dt)
    if formed<=EPS: return Mstar*F_REM
    sfr=shape*(Mstar/formed); lookback=T_OBS_GYR-t; kernel=(1.0+lookback/TAU0_KERNEL)**(-P_KERNEL)
    return float(np.trapezoid(sfr*kernel, dx=dt)*F_REM)

def empirical_depth(Mstar_emp, sfr_emp):
    Mstar_emp=max(float(Mstar_emp), EPS); sfr_emp=max(float(sfr_emp), EPS); ssfr=sfr_emp/Mstar_emp
    SSFR_REF=1e-10; persistence=1.0+(SSFR_REF/(ssfr+SSFR_REF))
    return float(Mstar_emp*F_REM*persistence)

df=pd.read_csv(INPUT)
for c in ["z0mgs_SFR_Msun_yr","z0mgs_Mstar_Msun"]:
    if c not in df.columns: df[c]=np.nan
D=[]; src=[]
for _, row in df.iterrows():
    emp=(pd.notna(row.get('z0mgs_SFR_Msun_yr')) and row.get('z0mgs_SFR_Msun_yr')>0 and pd.notna(row.get('z0mgs_Mstar_Msun')) and row.get('z0mgs_Mstar_Msun')>0)
    if emp:
        D.append(empirical_depth(row['z0mgs_Mstar_Msun'], row['z0mgs_SFR_Msun_yr'])); src.append('z0MGS_Empirical')
    else:
        D.append(delayed_tau_depth(row['Mstar'], row['Mgas'])); src.append('Delayed_Tau_Fallback')
df['D_hist_sfh']=D; df['D_hist_source']=src; df['SFH_source']=np.where(df['D_hist_source']=='z0MGS_Empirical','z0MGS_empirical','Delayed_Tau_Fallback')
mask=df['D_hist_source']=='z0MGS_Empirical'
df.loc[mask,'final_Mstar_for_SFH_Msun']=df.loc[mask,'z0mgs_Mstar_Msun']
df.loc[mask,'final_SFR_for_SFH_Msun_yr']=df.loc[mask,'z0mgs_SFR_Msun_yr']
df.loc[mask,'final_sSFR_for_SFH_yr']=df.loc[mask,'z0mgs_SFR_Msun_yr']/np.maximum(df.loc[mask,'z0mgs_Mstar_Msun'], EPS)
df.to_csv(OUTPUT,index=False)
print('Layer 1 z0MGS:', int((df['D_hist_source']=='z0MGS_Empirical').sum()))
print('Layer 2 fallback:', int((df['D_hist_source']=='Delayed_Tau_Fallback').sum()))
