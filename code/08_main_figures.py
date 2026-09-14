"""08_main_figures.py  Main figures 2 to 4 (Fig. 1 is drawn in Illustrator; its data panels come from 09_fig1_panels.py).

Reads the ledger and the participant-level frames from 02_survival_cohort.py.
Output: figures/main/Fig2_model_and_validation, Fig3_OIS_clinical, Fig4_what_the_residual_is (.pdf with editable Type 1 text, .png)
"""
import json, pickle, sys, numpy as np, pandas as pd
import matplotlib as mpl; mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import Ridge
from sklearn.model_selection import KFold, cross_val_predict
from sklearn.pipeline import make_pipeline
from lifelines import KaplanMeierFitter
from lifelines.utils import concordance_index
from lifelines.statistics import multivariate_logrank_test
from scipy import stats
import statsmodels.api as sm

mpl.rcParams.update({'pdf.fonttype':42,'ps.fonttype':42,'axes.unicode_minus':False,'figure.dpi':120,
    'axes.spines.top':False,'axes.spines.right':False,'axes.linewidth':0.9,'font.size':8.5,
    'axes.titlesize':9.5,'axes.labelsize':9,'xtick.labelsize':8,'ytick.labelsize':8,'legend.fontsize':7.8})
LANG = sys.argv[1] if len(sys.argv) > 1 else 'en'          # figures are produced in English
OUT = sys.argv[2] if len(sys.argv) > 2 else 'figures/main'
mpl.rcParams['font.sans-serif']=['Noto Sans CJK JP','Noto Sans CJK SC','DejaVu Sans']
mpl.rcParams['font.serif']=['Times','Times New Roman','Nimbus Roman','Liberation Serif','DejaVu Serif']
# 英文图件一律 Times New Roman(2026-09-11 用户定),中文版仍走 CJK sans。
mpl.rcParams['font.family']=('sans-serif' if LANG=='zh' else 'serif')
if LANG=='en':
    mpl.rcParams.update({'pdf.use14corefonts':True,'font.serif':['Times','Times New Roman','Nimbus Roman','DejaVu Serif']})

# ---- Illustrator-editable PDF for the English build: base-14 Type 1 fonts (WinAnsi), so
# ---- every non-WinAnsi glyph (Greek, arrows, minus sign, comparison signs) is spelled out.
import matplotlib.text as _mtext
_SAN = [('ΔC', 'delta C'), ('Δ', 'delta '), ('A\u03b2' + '42', 'A-beta-42'), ('A\u00df' + '42', 'A-beta-42'), ('α', 'alpha'), ('β', 'beta'), ('ρ', 'rho'), ('χ²', 'chi-square'),
        ('χ', 'chi'), ('≥', '>='), ('≤', '<='), ('≈', '~'), ('−', '-'), ('→', '->'), ('↔', '<->'), ('²', '2'), ('⁻', '-'),
        ('⁰', '0'), ('¹', '1'), ('³', '3'), ('⁴', '4'), ('⁵', '5'), ('⁶', '6'), ('⁷', '7'), ('⁸', '8'), ('⁹', '9'), ('…', '...'), ('′', "'")]
def _sanitize(fig):
    import unicodedata as _ud
    def _fallback(ch):
        try: ch.encode('cp1252'); return ch
        except UnicodeEncodeError: pass
        nm = _ud.name(ch, '')
        if 'GREEK SMALL LETTER' in nm: return nm.split()[-1].lower()
        if 'GREEK CAPITAL LETTER' in nm: return nm.split()[-1].capitalize()
        return '?'
    def _san(s0):
        s1 = s0
        for a, b in _SAN: s1 = s1.replace(a, b)
        return ''.join(_fallback(c) for c in s1)
    fig.canvas.draw()   # materialise tick-label strings from their formatters
    for ax in fig.axes:
        for axis in (ax.xaxis, ax.yaxis):
            labs = [t.get_text() for t in axis.get_ticklabels()]
            new = [_san(x) for x in labs]
            if new != labs:
                axis.set_ticks(axis.get_ticklocs()); axis.set_ticklabels(new)
    for t in fig.findobj(_mtext.Text):
        s0 = t.get_text()
        if s0:
            s1 = _san(s0)
            if s1 != s0: t.set_text(s1)

def T(en,zh): return zh if LANG=='zh' else en
# A4 竖版画布(210 x 297 mm),使每张主图就是一块页面大小的画板
A4W, A4H = 8.27, 11.69
TEAL,RED,GREY,DARK,ORANGE,BLUE,PURPLE='#16A085','#C0392B','#7F8C8D','#2C3E50','#E67E22','#2980B9','#8E44AD'
import os
ROOT=os.path.dirname(os.path.dirname(os.path.abspath(__file__))); os.chdir(ROOT); os.makedirs(OUT, exist_ok=True)
J=json.load(open('results/repro/section16_stats.json'))
# participant-level frames written by 02_survival_cohort.py (not distributed)
st=pickle.load(open('results/intermediate/02_cohort.pkl','rb')); df=st['df']; prod=st['prod'].copy(); FEAT=st['FEATURES']
H=[d for d in J['decomposition'] if d['cohort']=='Hyposmia'][0]
rng=np.random.default_rng(42)
from lifelines.statistics import logrank_test as _lrt
from scipy.stats import norm as _norm
def wilson(k,n,z=1.96):
    if n==0: return (0,0)
    ph=k/n; d=1+z*z/n; c=(ph+z*z/(2*n))/d; h=z*np.sqrt(ph*(1-ph)/n+z*z/(4*n*n))/d
    return (max(0,c-h)*100, min(1,c+h)*100)
def r_ci(r,n,z=1.96):
    zf=np.arctanh(r); se=1/np.sqrt(n-3); return np.tanh(zf-z*se), np.tanh(zf+z*se)
def pfmt(pv):
    return 'P < 0.001' if pv<0.001 else f'P = {pv:.3f}' if pv<0.01 else f'P = {pv:.2f}'
def km_band(ax,kmf,col,xmax=2.0):
    ci=kmf.confidence_interval_survival_function_; m=ci.index<=xmax
    ax.fill_between(ci.index[m],ci.iloc[:,0][m],ci.iloc[:,1][m],step='post',color=col,alpha=.10,lw=0)
def pairwise_text(ax,d,tcol,labs,x=.97,y0=.32,dy=.065,fs=6.6):
    pairs=[('T1','T2'),('T2','T3'),('T1','T3')]
    for i,(a_,b_) in enumerate(pairs):
        ga,gb=d[d[tcol]==a_],d[d[tcol]==b_]
        pv=_lrt(ga.time_years,gb.time_years,ga.converted,gb.converted).p_value
        ax.text(x,y0-i*dy,f'{labs[a_]} vs {labs[b_]}: {pfmt(pv)}'.replace(': ',' '),transform=ax.transAxes,ha='right',fontsize=fs,color='#444')
def panel(ax,l,t=None,dx=-0.16,dy=1.06):
    ax.text(dx,dy,l,transform=ax.transAxes,fontsize=12,fontweight='bold',va='top')
    if t: ax.set_title(t,loc='left',fontsize=9.3,fontweight='bold',pad=6)
# 散点一律保持矢量(不用 rasterized),否则 Illustrator 里点是一张位图,放大发虚、也改不了颜色。
# 代价是 PDF 里多几千个矢量圆,体积增加有限。
def save(fig,name,a4=False):
    """a4=True keeps the full canvas so the page stays A4 landscape, which is what
    Illustrator opens as a page-sized artboard. Cropping to the ink with bbox_inches
    would shrink the page back to the bounding box of the drawing."""
    if LANG=='en': _sanitize(fig)
    kw = {} if a4 else dict(bbox_inches='tight')
    fig.savefig(f'{OUT}/{name}.pdf',**kw); fig.savefig(f'{OUT}/{name}.png',dpi=300,**kw)
    plt.close(fig); print('  saved',name,'(A4)' if a4 else '')

pool=df[df.COHORT.isin([1,2])].dropna(subset=FEAT+['upsit']).copy()
pr=df[df.COHORT==4].dropna(subset=['OIS','upsit']).copy()
pr['arm']=np.where(pr.subgroup=='Hyposmia',T('Hyposmia','嗅觉减退队列'),
          np.where(pr.subgroup.astype(str).str.contains('RBD'),T('RBD','RBD 队列'),T('Variant carriers','携带者队列')))
ARM={T('Hyposmia','嗅觉减退队列'):RED,T('RBD','RBD 队列'):BLUE,T('Variant carriers','携带者队列'):PURPLE}

# ============================================================ Fig 2
def fig2():
    fig=plt.figure(figsize=(A4W,A4H)); gs=GridSpec(2,2,figure=fig,hspace=.30,wspace=.32,
                                                  top=.96,bottom=.30,left=.11,right=.97)
    # a 训练 CV
    ax=fig.add_subplot(gs[0,0])
    X,y=pool[FEAT].values,pool.upsit.values
    cv=cross_val_predict(make_pipeline(StandardScaler(),Ridge(alpha=1.0)),X,y,cv=KFold(5,shuffle=True,random_state=41))
    ispd=(pool.COHORT==1).values
    ax.scatter(cv[~ispd],y[~ispd],s=11,alpha=.55,c=TEAL,lw=0,label=T(f'Controls (n={(~ispd).sum()})',f'健康对照 n={(~ispd).sum()}'))
    ax.scatter(cv[ispd],y[ispd],s=9,alpha=.32,c=RED,lw=0,label=T(f'PD (n={ispd.sum()})',f'帕金森病 n={ispd.sum()}'))
    lim=[8,42]; ax.plot(lim,lim,'--',c=GREY,lw=.9,zorder=0)
    sl,ic=np.polyfit(cv,y,1); xs=np.linspace(*lim,50); ax.plot(xs,sl*xs+ic,c=DARK,lw=1.8)
    _r2a,_p2a=stats.pearsonr(cv,y); _lo2a,_hi2a=r_ci(_r2a,len(y))
    ax.text(.04,.96,T(f'5-fold CV  r = {_r2a:.3f} (95% CI {_lo2a:.3f} to {_hi2a:.3f}), P = {_p2a:.1e}',f'五折交叉验证  r = {_r2a:.3f}(95% CI {_lo2a:.3f} 至 {_hi2a:.3f}),P = {_p2a:.1e}'),
            transform=ax.transAxes,va='top',fontsize=8.6,fontweight='bold',
            bbox=dict(boxstyle='round,pad=.3',fc='white',ec=GREY,lw=.7,alpha=.92))
    ax.set_xlim(lim); ax.set_ylim(0,41); ax.legend(frameon=False,loc='lower right',fontsize=7.4)
    ax.set_xlabel(T('Cross-validated OIS','交叉验证 OIS')); ax.set_ylabel(T('Observed UPSIT','实测 UPSIT'))
    panel(ax,'a',T(f'Training cohort (n={len(pool)})',f'训练队列 PD + HC(n={len(pool)})'))
    # b LOSO
    ax=fig.add_subplot(gs[0,1])
    sites=pool.groupby('SITE').size(); sites=sites[sites>=10].index; rs=[]
    for s_ in sites:
        tr,te=pool[pool.SITE!=s_],pool[pool.SITE==s_]
        sc=StandardScaler().fit(tr[FEAT].values); m=Ridge(alpha=1.0).fit(sc.transform(tr[FEAT].values),tr.upsit.values)
        p=m.predict(sc.transform(te[FEAT].values))
        if te.upsit.std()>0 and len(te)>3:
            rs.append(stats.pearsonr(p,te.upsit.values)[0])
    rs=np.sort(np.array(rs))
    # 单色。柱子不按符号也不按构成着色:符号不说明判别力(2026-09-10 用户指出),
    # 而对照占比与 r 只是中等相关(0.521),用深浅编码在图上看着像没规律,反而削弱正文。
    # 构成这条线索由正文与 Supplementary Table 1d 的逐中心数字承担。
    ax.barh(np.arange(len(rs)),rs,color=TEAL,height=.78)
    ax.axvspan(rs.mean()-rs.std(),rs.mean()+rs.std(),color=DARK,alpha=.09,zorder=0)
    ax.axvline(rs.mean(),c=DARK,lw=1.4); ax.axvline(0,c='#34495E',lw=.9)
    _q1,_md,_q3=np.percentile(rs,[25,50,75])
    # 不标「几个为正」:符号本身不说明判别力,要标的是量级与离散
    # 放右下角:柱按 r 升序,最底几行为负值或极小,该区域为空,不会压住柱子
    ax.text(.98,.02,T(f'mean {rs.mean():.3f} ± {rs.std():.3f}\nmedian {_md:.3f} (IQR {_q1:.3f} to {_q3:.3f})',
            f'均值 {rs.mean():.3f} ± {rs.std():.3f}\n中位数 {_md:.3f}(四分位距 {_q1:.3f} 至 {_q3:.3f})'),
            transform=ax.transAxes,ha='right',va='bottom',fontsize=7.4,fontweight='bold',
            bbox=dict(boxstyle='round,pad=.28',fc='white',ec=GREY,lw=.7,alpha=.92))
    ax.set_yticks([]); ax.set_xlabel(T('Held-out r','留出中心上 OIS 与 UPSIT 的相关'))
    ax.set_ylabel(T(f'{len(rs)} PPMI sites','44 个 PPMI 中心'))
    panel(ax,'b',T('Leave-one-site-out','留一中心验证'))
    # c 前驱期校准,按入组队列:十分位均值 + 各队列回归线 + 等值线,斜率与 r 写进图例
    ax=fig.add_subplot(gs[1,:])
    lim=[15,42]; ax.plot(lim,lim,'--',c=GREY,lw=.9,zorder=0,label=T('identity (perfect calibration)','等值线(完美校准)'))
    for a_,g in pr.groupby('arm'):
        r,_pc=stats.pearsonr(g.OIS,g.upsit); s2,i2=np.polyfit(g.OIS,g.upsit,1)
        dec=pd.qcut(g.OIS,10,labels=False,duplicates='drop')
        mx=g.groupby(dec).OIS.mean(); my=g.groupby(dec).upsit.mean()
        ax.scatter(g.OIS,g.upsit,s=5,alpha=.10,c=ARM[a_],lw=0)
        ax.scatter(mx,my,s=34,c=ARM[a_],ec='white',lw=.8,zorder=4)
        xs=np.linspace(g.OIS.quantile(.01),g.OIS.quantile(.99),30)
        ax.plot(xs,s2*xs+i2,c=ARM[a_],lw=2,zorder=3,label=T(f'{a_}  slope {s2:.2f}, r {r:.3f}, P {_pc:.0e}, n {len(g)}',f'{a_}  斜率 {s2:.2f},r {r:.3f},P {_pc:.0e},n {len(g)}'))
    rp=stats.pearsonr(pr.OIS,pr.upsit)[0]
    ax.text(.03,.97,T(f'Pooled r = {rp:.3f}',f'合并队列 r = {rp:.3f}'),transform=ax.transAxes,va='top',
            fontsize=8.8,fontweight='bold',bbox=dict(boxstyle='round,pad=.3',fc='#FDF2E9',ec=ORANGE,lw=.8))
    ax.legend(frameon=False,loc='lower right',fontsize=6.6,handlelength=1.6)
    ax.set_xlim(lim); ax.set_ylim(0,41)
    ax.set_xlabel(T('OIS (predicted UPSIT)','OIS(预测 UPSIT)')); ax.set_ylabel(T('Observed UPSIT','实测 UPSIT'))
    ax.text(.03,.86,T('faint points, individuals. Filled circles, decile means within each arm','淡点为个体,实心圆为各队列内的十分位均值'),transform=ax.transAxes,fontsize=6.8,color='#555')
    panel(ax,'c',T(f'Prodromal calibration (n={len(pr)}), by arm',f'前驱期校准(n={len(pr)}),按入组队列'),dx=-.075)
    # d 分解
    save(fig,'Fig2_model_and_validation',a4=True)
fig2()

# ============================================================ Fig 3
def fig3():
    hh=prod[prod.subgroup=='Hyposmia'].dropna(subset=['OIS','upsit','PUTAMEN_REF_CWM','PUTAMEN_L_REF_CWM',
        'PUTAMEN_R_REF_CWM','time_years','converted','age','sex_num']+st['SBR_COLS']).copy()
    # 0.25 年下限已在上游仅对转化者施加，此处不得对删失者再截断
    hh=hh[(hh.converted==1)|(hh.time_years>0)]
    hh['pct']=hh[['PUTAMEN_L_REF_CWM','PUTAMEN_R_REF_CWM']].min(axis=1)/(1.8547-.00303*hh.age-.2419*hh.sex_num)*100
    hh['age_neg']=-hh['age'].values
    # 仅影像版 OIS：用 33 个 SBR 重训一个岭回归（与 notebook §40 同法）
    from sklearn.linear_model import Ridge as _R; from sklearn.preprocessing import StandardScaler as _S
    _tr=st['train'].dropna(subset=st['SBR_COLS']+['upsit'])
    _sc=_S().fit(_tr[st['SBR_COLS']].values)
    _mo=_R(alpha=1.0).fit(_sc.transform(_tr[st['SBR_COLS']].values), _tr['upsit'].values)
    hh['OIS_sbr']=_mo.predict(_sc.transform(hh[st['SBR_COLS']].values))
    hh['band']=np.where(hh.pct<65,'def',np.where(hh.pct<80,'ind','nor'))
    dfc, ndf = hh[hh.pct<65], hh[hh.pct>=65]
    fig=plt.figure(figsize=(A4W,A4H)); gs=GridSpec(3,6,figure=fig,hspace=.70,wspace=.84,
                                                  top=.965,bottom=.05,left=.10,right=.97,
                                                  height_ratios=[0.72,1.22,1.32])
    kmf=KaplanMeierFitter(); cols={'T1':RED,'T2':ORANGE,'T3':TEAL}
    nm={'T1':T('Lowest OIS','OIS 最低'),'T2':T('Middle','中间'),'T3':T('Highest OIS','OIS 最高')}

    # a 分解:嗅觉减退组内三个分数的转化 C-index 与方差占比(原 Fig 2d,2026-09-07 移入)
    ax=fig.add_subplot(gs[0,:4])
    # C 与 95% CI 直接取 §36(与正文表 2 同一次 bootstrap);方差占比取 decomposition
    _C36={r_['reading']:r_ for r_ in J['cindex_ci'] if r_['stratum']=='hyposmia'}
    rows=[(T('UPSIT','UPSIT 总分'),_C36['UPSIT total']['c'],_C36['UPSIT total']['lo'],_C36['UPSIT total']['hi'],GREY,100.0),
          (T('OIS','OIS 拟合值'),_C36['OIS']['c'],_C36['OIS']['lo'],_C36['OIS']['hi'],TEAL,H['sd_ois']**2/H['sd_upsit']**2*100),
          (T('OMI','OMI 残差'),_C36['OMI']['c'],_C36['OMI']['lo'],_C36['OMI']['hi'],RED,H['sd_omi']**2/H['sd_upsit']**2*100)]
    # 自上而下:UPSIT 总分、OMI、OIS(用户 2026-09-04 定)
    for i,(lab,c,lo,hi,col,var) in enumerate(rows):
        y=2-i
        ax.plot([lo,hi],[y,y],c=col,lw=2.6,solid_capstyle='round')
        ax.scatter([c],[y],s=70,c=col,zorder=3,ec='white',lw=1.1)
        ax.text(hi+.012,y+.16,f'{c:.3f} ({lo:.3f} to {hi:.3f})',va='center',fontsize=7.6,color=col,fontweight='bold')
        if var<99: ax.text(hi+.012,y-.20,T(f'{var:.1f}% of variance',f'占总分方差 {var:.1f}%'),va='center',fontsize=7.2,color=col)
        if lo<0.5: ax.text(hi+.012,y-.42,T('interval includes 0.5','区间含 0.5'),va='center',fontsize=6.8,color=col,style='italic')
    _pre={r_['comparison']:r_ for r_ in J['prespecified_comparisons'] if r_['stratum']=='hyposmia'}
    _pu=_pre['OIS - UPSIT']
    ax.annotate('',xy=(1.30,1),xytext=(1.30,2),arrowprops=dict(arrowstyle='-',color=DARK,lw=1))
    ax.text(1.33,1.5,T(f"OIS vs UPSIT\ndelta C {_pu['dC']:+.3f}\n({_pu['lo']:+.3f} to {_pu['hi']:+.3f})\n{pfmt(_pu['p'])}",
                       f"OIS 对 UPSIT\nΔC {_pu['dC']:+.3f}\n({_pu['lo']:+.3f} 至 {_pu['hi']:+.3f})\n{pfmt(_pu['p'])}"),
            va='center',ha='left',fontsize=6.8,color=DARK)
    ax.axvline(.5,c=GREY,ls='--',lw=.9,zorder=0)
    ax.text(.5,2.52,T('chance','随机水平'),ha='center',fontsize=7.2,color=GREY)
    ax.set_yticks(range(3)); ax.set_yticklabels([r[0] for r in rows][::-1])
    ax.set_xlim(.40,1.70); ax.set_ylim(-.6,2.7); ax.set_xticks([.4,.5,.6,.7,.8,.9,1.0])
    ax.set_xlabel(T('C-index for conversion (95% CI)','转化 C-index(95% CI)'))
    _EV=68
    panel(ax,'a',T(f'Hyposmia group ({H["n"]} participants, {_EV} conversions)',f'嗅觉减退组({H["n"]} 人,{_EV} 例转化)'),dx=-.085)
    _R=H['r_ois_omi']
    fig.text(.5,-.030,T(f'Variance shares exceed 100% because the two components are negatively correlated (r = {_R:+.3f}). The residual is at chance in this group only.',
             f'两成分方差占比之和超过 100%,因二者负相关(r = {_R:+.3f})。残差等同随机仅在本组成立。'),ha='center',fontsize=7.4,style='italic',color='#555')

    # b 时依 AUC(IPCW),§46。与 a 同一组人群,但按时点而非全程给判别力
    ax=fig.add_subplot(gs[0,4:])
    _TD=J['td_auc_ipcw']; _tr=pd.DataFrame(_TD['rows'])
    _ORD=[('UPSIT total',T('UPSIT','UPSIT 总分'),GREY),
          ('putamen SBR',T('Putamen SBR','单区壳核'),DARK),
          ('OIS',T('OIS','OIS'),TEAL)]
    _ts=_TD['times']; _x=np.arange(len(_ts))
    for _k,(_rd,_lab,_col) in enumerate(_ORD):
        _r=_tr[_tr.reading==_rd].set_index('horizon_yr').loc[_ts]
        _off=(_k-1)*0.10
        ax.errorbar(_x+_off,_r.auc.values,
                    yerr=[_r.auc.values-_r.lo.values,_r.hi.values-_r.auc.values],
                    fmt='o-',color=_col,lw=1.6,ms=4.6,capsize=2.4,elinewidth=.9,label=_lab)
    _d=_tr[_tr.reading=='OIS minus UPSIT total'].set_index('horizon_yr').loc[_ts]
    for _i,_t in enumerate(_ts):
        ax.text(_x[_i],1.015,pfmt(float(_d.loc[_t,'p'])),ha='center',fontsize=6.2,color=TEAL,fontweight='bold')
    ax.axhline(.5,c=GREY,ls='--',lw=.9)
    ax.text(.02,.505,T('chance','随机水平'),transform=ax.get_yaxis_transform(),fontsize=6.4,color=GREY)
    ax.set_xticks(_x)
    _n0=_tr[_tr.reading=='OIS'].set_index('horizon_yr').loc[_ts]
    ax.set_xticklabels([f"{_t:.1f}\n{int(_n0.loc[_t,'n_at_risk'])} / {int(_n0.loc[_t,'events_by_t'])}" for _t in _ts],fontsize=7)
    ax.set_xlabel(T('Years    at risk / conversions','年    在险 / 累计转化'),fontsize=7.6)
    ax.set_ylabel(T('Time-dependent AUC (95% CI)','时依 AUC(95% CI)'),fontsize=7.6)
    ax.set_ylim(.42,1.04); ax.set_xlim(-.45,len(_ts)-.55)
    ax.legend(frameon=False,fontsize=6.2,loc='lower right',handlelength=1.3)
    panel(ax,'b',T('Discrimination by horizon','按时点的判别力'),dx=-.30)

    # c 三个影像情形 × 七种读法,数值全部取自 notebook §40(calibration_and_age),与正文两张表同源
    ax=fig.add_subplot(gs[1,:])
    METH=[('single-region putamen',T('Putamen SBR','单区壳核'),DARK),
          ('lowest putamen %expected',T('Lower putamen, % expected','较低侧壳核 %预期'),'#5D6D7E'),
          ('33 SBR + demographics, CV Cox',T('33 SBR + demographics, supervised Cox','33 区 + 人口学,转化标签 Cox'),'#9B8EC4'),
          ('age alone',T('Age','年龄'),'#C8A27F'),
          ('UPSIT total',T('UPSIT total','UPSIT 总分'),GREY),
          ('OIS, imaging only, zero-shot',T('OIS, imaging only','OIS 仅影像版'),'#7FB3C8'),
          ('OIS, full, zero-shot','OIS',TEAL)]
    _CA={(r['stratum'],r['reading']):r for r in J['calibration_and_age'] if 'c' in r and r['c']==r['c']}
    _ST=['hyposmia','hyposmia + DAT deficit','hyposmia + no DAT deficit']
    _SL=[T('Whole hyposmia group','全嗅觉减退组'),T('DAT deficit','DAT 缺损亚组'),T('No DAT deficit','DAT 非缺损亚组')]
    x=np.arange(len(_ST)); w=.115; nm_=len(METH)
    for k2,(key,lab,col) in enumerate(METH):
        vals,los,his=[],[],[]
        for stn in _ST:
            r=_CA[(stn,key)]
            vals.append(r['c']); los.append(r['c']-r['lo']); his.append(r['hi']-r['c'])
        pos=x+(k2-(nm_-1)/2)*w
        ax.bar(pos,vals,w,color=col,label=lab,yerr=[los,his],
               error_kw=dict(ecolor='#566573',lw=.8,capsize=1.6))
        for xi,vv,lo,hi_ in zip(pos,vals,los,his):
            ax.text(xi,vv+hi_+.012,f'{vv:.3f}',ha='center',fontsize=6.2,fontweight='bold',rotation=90,va='bottom')
            if vv-lo < 0.5:
                ax.text(xi,.312,'ns',ha='center',va='bottom',fontsize=6.4,color=RED,fontweight='bold')
    ax.axhline(.5,c=GREY,ls='--',lw=1)
    ax.text(-.5,.503,T('chance','随机水平'),fontsize=7,color=GREY,ha='left',va='bottom')
    ax.set_xticks(x)
    _PP={}
    for r_ in J['prespecified_comparisons']:
        _PP.setdefault(r_['stratum'],[]).append((r_['comparison'],r_['p']))
    def _pline(stn):
        parts=[]
        for cmp_,pv in _PP.get(stn,[]):
            lab=T('OIS vs UPSIT','OIS 对 UPSIT') if 'UPSIT' in cmp_ else T('OIS vs putamen','OIS 对壳核')
            parts.append(f'{lab} {pfmt(pv)}')
        return '\n'.join(parts)      # one comparison per line, the A4 width is too narrow for a single row
    ax.set_xticklabels([f"{sl}\n(n={_CA[(stn,'OIS, full, zero-shot')]['n']}, {_CA[(stn,'OIS, full, zero-shot')]['events']} {T('conv','转化')})\n{_pline(stn)}"
                        for sl,stn in zip(_SL,_ST)],fontsize=7.0)
    ax.set_xlim(-.55,len(_ST)-.45)
    ax.set_ylim(.30,1.20); ax.set_yticks([.3,.4,.5,.6,.7,.8,.9,1.0])
    ax.legend(frameon=False,ncol=4,fontsize=6.4,loc='upper center',bbox_to_anchor=(.5,1.005),columnspacing=1.0,handlelength=1.2)
    ax.set_ylabel(T('C-index [95% CI]','转化 C-index [95% CI]'),fontsize=8.5)
    ax.text(1.0,-.255,T('ns: interval includes 0.5. Bars come from one bootstrap run. P values are the paired comparisons of Table 3.',
                       'ns:区间含 0.5。柱来自单次 bootstrap。P 值为表 3 的配对比较。'),
            transform=ax.transAxes,ha='right',va='top',fontsize=6.4,color='#777',style='italic')
    panel(ax,'c',T('Seven readings in three imaging situations','三种影像情形下七种读法的判别力'),dx=-.09)

    # c 全组 KM,两组:OIS 最低三分位(高危)对其余三分之二(低危),按三分位边界划分,不看结局
    ax=fig.add_subplot(gs[2,:3])
    hh2=hh.copy(); hh2['grp']=np.where(hh2.OIS<=hh2.OIS.median(),'high','low')
    two_cols={'high':RED,'low':TEAL}
    two_nm={'high':T('High risk (OIS at or below median)','高危(OIS 不高于中位数)'),'low':T('Low risk (OIS above median)','低危(OIS 高于中位数)')}
    def _km_two(ax,d):
        for k_ in ['high','low']:
            g=d[d.grp==k_]; kmf.fit(g.time_years,g.converted)
            sf=kmf.survival_function_; m=sf.index<=2
            ax.step(sf.index[m],sf.iloc[:,0][m],where='post',color=two_cols[k_],lw=2.2,
                    label=f'{two_nm[k_]} (n={len(g)}, {int(g.converted.sum())} ev)')
            km_band(ax,kmf,two_cols[k_])
            ci=kmf.confidence_interval_survival_function_; idx=ci.index[ci.index<=2.0][-1]
            r2=(1-float(kmf.predict(2.)))*100; lo2=(1-ci.loc[idx].iloc[1])*100; hi2=(1-ci.loc[idx].iloc[0])*100
            ax.text(2.04,float(kmf.predict(2.)),f'{r2:.1f}%\n({lo2:.1f} to {hi2:.1f})',fontsize=7.2,color=two_cols[k_],va='center',fontweight='bold')
        _a,_b=d[d.grp=='high'],d[d.grp=='low']
        return _lrt(_a.time_years,_b.time_years,_a.converted,_b.converted).p_value
    p2c=_km_two(ax,hh2)
    ax.text(.97,.36,f'log-rank P = {p2c:.1e}',transform=ax.transAxes,ha='right',fontsize=7.6,fontweight='bold')
    ax.text(.97,.29,T(f'split at the OIS median ({hh2.OIS.median():.2f}), outcome-blind',f'按 OIS 中位数({hh2.OIS.median():.2f})划分,不看结局'),transform=ax.transAxes,ha='right',fontsize=6.4,color='#555')
    ax.text(.97,.23,T('bands: 95% CI','带:95% CI'),transform=ax.transAxes,ha='right',fontsize=6.2,color='#777',style='italic')
    ax.set_xlim(0,2.55); ax.set_ylim(.68,1.005); ax.legend(frameon=False,loc='lower left',fontsize=5.5)
    ax.set_xlabel(T('Years','距基线年数'),fontsize=8.5); ax.set_ylabel(T('Free of PD','未确诊比例'),fontsize=8.5)
    panel(ax,'d',T(f'All hyposmic (n={len(hh)})',f'全嗅觉减退组(n={len(hh)})'),dx=-.155)

    # d PARS 三档:人数 vs 事件
    ax=fig.add_subplot(gs[2,3:])
    # e 非缺损亚组内部按 OIS 中位数再分。影像分档的构成只进正文,不占主图版面
    # 下半:非缺损亚组按 OIS 中位数分两组
    nd=ndf.copy(); nd['grp']=np.where(nd.OIS<=nd.OIS.median(),'high','low')
    p2e=_km_two(ax,nd)
    def _med_split(d,col):
        a_,b_=d[d[col]<=d[col].median()],d[d[col]>d[col].median()]
        k1=KaplanMeierFitter().fit(a_.time_years,a_.converted); k2=KaplanMeierFitter().fit(b_.time_years,b_.converted)
        return (1-float(k1.predict(2.)))*100,(1-float(k2.predict(2.)))*100,_lrt(a_.time_years,b_.time_years,a_.converted,b_.converted).p_value
    _u=_med_split(nd,'upsit'); _i=_med_split(nd,'pct')
    ax.text(.97,.36,f'log-rank P = {p2e:.1e}',transform=ax.transAxes,ha='right',fontsize=7.6,fontweight='bold')
    ax.text(.97,.29,T(f'{len(ndf)} non-deficit, split at the OIS median ({nd.OIS.median():.2f}), outcome-blind',
            f'非缺损 {len(ndf)} 人,按 OIS 中位数({nd.OIS.median():.2f})划分,不看结局'),transform=ax.transAxes,ha='right',fontsize=6.4,color='#555')
    ax.text(.97,.23,T('bands: 95% CI','带:95% CI'),transform=ax.transAxes,ha='right',fontsize=6.2,color='#777',style='italic')
    ax.set_xlim(0,2.55); ax.set_ylim(.875,1.003)
    ax.legend(frameon=False,loc='lower left',fontsize=5.8,bbox_to_anchor=(-.01,-.02))
    ax.set_xlabel(T('Years','距基线年数'),fontsize=8.5); ax.set_ylabel(T('Free of PD','未确诊比例'),fontsize=8.5)
    panel(ax,'e',T(f'Within non-deficit (n={len(ndf)})',f'非缺损亚组内部(n={len(ndf)})'),dx=-.155)

    fig.text(.5,-.02,T(f'Curves are shown to 2 years. Panels d and e split each population at its OIS median. In panel e the same split by UPSIT gives {_u[0]:.1f}% against {_u[1]:.1f}% (P = {_u[2]:.1e}) and by lower putamen %expected {_i[0]:.1f}% against {_i[1]:.1f}% (P = {_i[2]:.2f}), so the stratification there comes from the smell information, not from imaging, and that stratum has 22 events. Tertile and permutation-corrected three-band versions are in Supplementary Figs. 4 and 5.',
             f'曲线绘至两年。面板 d 与 e 按各自的 OIS 中位数划分。面板 e 中按 UPSIT 中位数划分为 {_u[0]:.1f}% 对 {_u[1]:.1f}%(P = {_u[2]:.1e}),按较低侧壳核 %预期划分为 {_i[0]:.1f}% 对 {_i[1]:.1f}%(P = {_i[2]:.2f}),该组的分层来自嗅觉信息而非影像,且仅 22 例事件。三分位版与置换校正三档版见 Supplementary Fig 4 与 5。'),
             ha='center',fontsize=7.2,style='italic',color='#555')
    save(fig,'Fig3_OIS_clinical',a4=True)
fig3()

# ============================================================ Fig 4
def fig4():
    fig=plt.figure(figsize=(A4W,A4H)); gs=GridSpec(2,2,figure=fig,hspace=.34,wspace=.36,
                                                  top=.96,bottom=.34,left=.12,right=.97)
    m=J['csf_tau']['main']
    # a tau + NfL,横跨整个上行
    ax=fig.add_subplot(gs[0,:])
    keys=[('Prodromal|pTau181_over_ABeta42_CSF',T('CSF pTau181\n/Aβ42','脑脊液 pTau181\n/Aβ42')),
          ('Prodromal|pTau181_CSF',T('CSF\npTau181','脑脊液\npTau181')),
          ('Prodromal|eMTBR_TAU243_CSF',T('CSF eMTBR\n-tau243','脑脊液 eMTBR\n-tau243')),
          ('Prodromal|NfL_serum',T('serum\nNfL','血清\nNfL')),
          ('Prodromal|NfL_plasma',T('plasma\nNfL','血浆\nNfL'))]
    keys=[(k,l) for k,l in keys if k in m]
    x=np.arange(len(keys)); w=.26; _PROWS={}
    for off,(f_,pf,lab,col) in enumerate([('b_UPSIT','p_UPSIT',T('UPSIT','总分'),GREY),
                                          ('b_OIS','p_OIS',T('OIS','OIS'),TEAL),
                                          ('b_OMI','p_OMI',T('OMI','残差'),RED)]):
        v=[m[k][f_] for k,_ in keys]; ps=[m[k][pf] for k,_ in keys]; pos=x+(off-1)*w
        _se=[abs(b_)/_norm.isf(pp/2) if 0<pp<1 else 0 for b_,pp in zip(v,ps)]
        ax.bar(pos,v,w,color=col,label=lab,yerr=[1.96*e for e in _se],error_kw=dict(ecolor='#566573',lw=.7,capsize=1.5))
        for xi,vv,pp,e in zip(pos,v,ps,_se):
            st_='***' if pp<.001 else ('**' if pp<.01 else ('*' if pp<.05 else ''))
            yy=vv-1.96*e-.004 if vv<0 else vv+1.96*e+.004
            if st_: ax.text(xi,yy,st_,ha='center',va='top' if vv<0 else 'bottom',fontsize=7,color='#333')
        _PROWS[lab]=[f'P={pp:.1e}' if pp<0.001 else f'P={pp:.3f}' if pp<0.01 else f'P={pp:.2f}' for pp in ps]
    ax.axhline(0,c='#34495E',lw=.9); ax.set_xticks(x)
    _labs=[T('UPSIT','总分'),T('OIS','OIS'),T('OMI','残差')]
    ax.set_xticklabels([f"{l}\n(n={m[k]['n']})\n"+"\n".join(f"{lb} {_PROWS[lb][i]}" for lb in _labs) for i,(k,l) in enumerate(keys)],fontsize=5.6)
    ax.set_ylim(-.22,.20); ax.legend(frameon=False,ncol=3,loc='upper left',fontsize=7.4,bbox_to_anchor=(0,1.0))
    ax.text(.99,.02,T('bars: standardised beta with 95% CI','柱:标准化系数与 95% CI'),transform=ax.transAxes,ha='right',fontsize=6,color='#777',style='italic')
    ax.set_ylabel(T('Standardised β','标准化回归系数'))
    ax.axvspan(2.5,4.5,color='#F4F6F6',zorder=0)
    ax.text(3.5,.185,T('negative control','阴性对照'),ha='center',fontsize=7,color='#777',style='italic')
    panel(ax,'a',T('Tau axis with NfL control (prodromal, n = 827)',
                   'tau 轴与 NfL 对照(前驱期,n=827)'),dx=-.055)
    # c 基因型抵消
    _G={(r_['cohort'],r_['measure']):r_ for r_ in J['genotype']}
    for k,(ttl,coh,nn) in enumerate([(T('Non-manifesting carriers','未发病携带者'),'Prodromal','GBA1 161 / LRRK2 166'),
                                     (T('Diagnosed PD','已确诊帕金森病'),'PD','GBA1 69 / LRRK2 126')]):
        ax=fig.add_subplot(gs[1,k])
        vals=[('OIS',_G[(coh,'OIS')]),(T('OMI','残差'),_G[(coh,'OMI')]),(T('UPSIT total','嗅觉总分'),_G[(coh,'UPSIT total')])]
        cols=[TEAL,RED,DARK]
        bb=ax.bar([v[0] for v in vals],[v[1]['beta'] for v in vals],color=cols,width=.5,
                  yerr=[[v[1]['beta']-v[1]['lo'] for v in vals],[v[1]['hi']-v[1]['beta'] for v in vals]],error_kw=dict(ecolor='#566573',lw=.9,capsize=3))
        for r,(lab,g_) in zip(bb,vals):
            v=g_['beta']; edge=g_['hi'] if v>=0 else g_['lo']
            ax.text(r.get_x()+r.get_width()/2, edge+(0.25 if v>=0 else -0.25), f'{v:+.2f}',
                    ha='center',va='bottom' if v>=0 else 'top',fontsize=9,fontweight='bold')
            ax.text(r.get_x()+r.get_width()/2, edge+(1.1 if v>=0 else -1.1), pfmt(g_['p']), ha='center',
                    va='bottom' if v>=0 else 'top', fontsize=7.2, color='#555')
        ax.axhline(0,c='#34495E',lw=1)
        lo,hi=min(v[1]['lo'] for v in vals),max(v[1]['hi'] for v in vals)
        pad=(hi-lo)*.35+.6; ax.set_ylim(lo-pad,hi+pad)
        ax.text(.5,.03,T('OIS + OMI = total; bars are 95% CI of the adjusted difference','OIS + 残差 = 总分;误差线为校正后差值的 95% CI').replace('; ',', ').replace(';',','),transform=ax.transAxes,ha='center',
                fontsize=7,style='italic',color=GREY)
        ax.set_ylabel(T('GBA1 minus LRRK2 (points)','GBA1 减 LRRK2(分)'))
        panel(ax,'b' if k==0 else 'c', f'{ttl}\n({nn})',dx=-.22)
    fig.text(.5,-.015,T('Negative controls: neurofilament light in a, and OIS in the diagnosed cohort in c (P = 0.41).',
             '阴性对照:a 中为神经丝轻链,c 中为已确诊队列内的 OIS(P=0.41)。'),
             ha='center',fontsize=7.6,style='italic',color='#555')
    save(fig,'Fig4_what_the_residual_is',a4=True)
fig4()
