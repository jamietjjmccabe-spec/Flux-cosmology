"""Coordinate-match SPARC rows to z0MGS after sparc_ra_deg/sparc_dec_deg have been resolved."""
import numpy as np, pandas as pd
from astropy.coordinates import SkyCoord
import astropy.units as u
SPARC_PATH='sparc_z0mgs_crossmatch.csv'; Z0_PATH='z0mgs_table4.csv'; OUT_PATH='sparc_z0mgs_crossmatch.csv'; MATCH_RADIUS_ARCMIN=3.0
sparc=pd.read_csv(SPARC_PATH); z0=pd.read_csv(Z0_PATH)
z=z0.dropna(subset=['ra','dec','logMstar','logSFR']).copy()
zc=SkyCoord(ra=z['ra'].values*u.deg, dec=z['dec'].values*u.deg)
for col in ['z0mgs_PGC','z0mgs_logMstar','z0mgs_logSFR','z0mgs_Mstar_Msun','z0mgs_SFR_Msun_yr','z0mgs_sSFR_yr','z0mgs_Dist_Mpc','z0mgs_match_sep_arcmin']:
    sparc[col]=np.nan
sparc['z0mgs_match_status']='unmatched'; sparc['z0mgs_match_method']=''
for idx,row in sparc.iterrows():
    if pd.isna(row.get('sparc_ra_deg')) or pd.isna(row.get('sparc_dec_deg')): continue
    c=SkyCoord(row['sparc_ra_deg']*u.deg, row['sparc_dec_deg']*u.deg); sep=c.separation(zc); j=int(np.argmin(sep.arcmin))
    if float(sep.arcmin[j])<=MATCH_RADIUS_ARCMIN:
        zr=z.iloc[j]; logm=float(zr['logMstar']); logsfr=float(zr['logSFR']); m=10**logm; sfr=10**logsfr
        sparc.loc[idx,'z0mgs_PGC']=zr['PGC']; sparc.loc[idx,'z0mgs_logMstar']=logm; sparc.loc[idx,'z0mgs_logSFR']=logsfr
        sparc.loc[idx,'z0mgs_Mstar_Msun']=m; sparc.loc[idx,'z0mgs_SFR_Msun_yr']=sfr; sparc.loc[idx,'z0mgs_sSFR_yr']=sfr/max(m,1e-30)
        sparc.loc[idx,'z0mgs_Dist_Mpc']=zr['Dist'] if 'Dist' in zr.index else np.nan; sparc.loc[idx,'z0mgs_match_sep_arcmin']=float(sep.arcmin[j])
        sparc.loc[idx,'z0mgs_match_status']='matched'; sparc.loc[idx,'z0mgs_match_method']='coordinate'
matched=sparc['z0mgs_SFR_Msun_yr'].notna() & (sparc['z0mgs_SFR_Msun_yr']>0)
sparc['SFH_source']=np.where(matched,'z0MGS_empirical','Delayed_Tau_Fallback')
sparc['D_hist_source']=np.where(matched,'z0MGS_Empirical','Delayed_Tau_Fallback')
sparc.to_csv(OUT_PATH,index=False)
print(sparc['D_hist_source'].value_counts(dropna=False))
