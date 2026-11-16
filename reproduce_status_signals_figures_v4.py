"""Reproduce the figures and statistical reconstructions accompanying
'What Can Status Signals Really Tell Us?' (fourth revision).

Input is transcribed from rounded published plots/tables, NOT individual-level data.
Run: python reproduce_status_signals_figures_v4.py
Dependencies: numpy, scipy, matplotlib. The outputs are written beside this script.
A meta-analysis of these heterogeneous correlations is descriptive, not a substitute
for a systematic synthesis. The stochastic model/interval are illustrative.
"""
from pathlib import Path
from statistics import NormalDist
import csv
import numpy as np
from scipy.stats import norm, t
import matplotlib.pyplot as plt

OUT=Path(__file__).resolve().parent
INK='#183043'; BLUE='#0e7895'; ORANGE='#c76743'; SOFT='#dde5e9'; GRAY='#536579'
plt.rcParams.update({'font.size':12, 'figure.facecolor':'white', 'axes.spines.top':False, 'axes.spines.right':False, 'savefig.facecolor':'white'})

# Forest plot input as read from Alexander's plotted labels (rounded r and CI).
META=[
 ('Blanchard & Clemmensen (1985)', .30,.13,.47),
 ('Latkin et al. (1993)', .20,.16,.23),
 ('Latkin & Vlahov (1998)', .15,-.04,.34),
 ('Meston et al. (1998)', .10,.01,.19),
 ('Gibson et al. (1999)', .10,.01,.19),
 ('Bradley & Wildman (2002)', .20,.10,.30),
 ('Rao et al. (2017)', .38,.31,.45),
 ('Latkin et al. (2017)', .11,.03,.19),
 ('Rasmussen et al. (2018)', .07,.04,.10),
 ('King et al. (2019)', .16,.01,.31),
]
name=[row[0] for row in META]; r=np.array([row[1] for row in META]); low=np.array([row[2] for row in META]); high=np.array([row[3] for row in META]); k=len(r)
z=np.arctanh(r)
# Recover approximate normal SE from Fisher-z-transformed published 95% CI widths.
se=(np.arctanh(high)-np.arctanh(low))/(2*norm.ppf(.975)); var=se**2
wf=1/var; mu_fixed=np.sum(wf*z)/np.sum(wf)
Q=np.sum(wf*(z-mu_fixed)**2); C=np.sum(wf)-np.sum(wf**2)/np.sum(wf)
tau2=max(0,(Q-(k-1))/C) # DerSimonian-Laird moment estimate
wr=1/(var+tau2); weights=wr/np.sum(wr)
mu=np.sum(wr*z)/np.sum(wr); semu=np.sqrt(1/np.sum(wr))
pooled=np.tanh(mu); pooled_CI=np.tanh(mu+np.array([-1,1])*norm.ppf(.975)*semu)
PI=np.tanh(mu+np.array([-1,1])*t.ppf(.975,k-2)*np.sqrt(tau2+semu**2))
I2=max(0,(Q-(k-1))/Q)
print('META pooled r:',round(pooled,4),'CI:',np.round(pooled_CI,4),'Q:',round(Q,3),'I2:',round(I2,4),'tau:',round(tau2**.5,4),'PI:',np.round(PI,4))

with open(OUT/'meta_reconstruction_inputs.csv','w',newline='',encoding='utf-8') as f:
    w=csv.writer(f);w.writerow(['study','published_rounded_r','published_rounded_ci_lower','published_rounded_ci_upper','reconstructed_random_effects_weight_pct'])
    for i,row in enumerate(META):w.writerow([*row,round(weights[i]*100,2)])

# 1 -- WEIGHTED FOREST
fig,ax=plt.subplots(figsize=(12.5,7.5),dpi=180)
y=np.arange(k)[::-1]+2
ax.axvline(0,ls='--',lw=1,color='#a8b4bb')
ax.axvline(pooled,ls=':',lw=1.7,color=ORANGE,alpha=.7)
for i,yy in enumerate(y):
    ax.hlines(yy, low[i],high[i],color=INK,lw=1.8)
    ax.plot([low[i],low[i]],[yy-.07,yy+.07],color=INK,lw=1.5)
    ax.plot([high[i],high[i]],[yy-.07,yy+.07],color=INK,lw=1.5)
    # Marker AREA, not side-length, scales linearly with reconstructed RE weight.
    ax.scatter(r[i],yy,s=10+1000*weights[i],marker='s',color=BLUE,zorder=4)
    ax.text(.55,yy,f'{r[i]:.2f} [{low[i]:+.2f}, {high[i]:+.2f}]   {weights[i]*100:4.1f}%',va='center',fontsize=10,color=GRAY)
ax.fill([pooled_CI[0],pooled,pooled_CI[1],pooled], [1,1.23,1,0.77],facecolor=ORANGE,alpha=.92,edgecolor='#9f3e1f',lw=2.2,zorder=5)
ax.hlines(0,PI[0],PI[1],color=ORANGE,lw=2.3,alpha=.65)
ax.plot(PI, [0,0],'|',color=ORANGE,ms=12,mew=1.7,alpha=.65)
ax.scatter(pooled,0,marker='D',color=ORANGE,s=48,zorder=4)
ax.text(.55,1,f'Pooled {pooled:.2f} [{pooled_CI[0]:+.2f}, {pooled_CI[1]:+.2f}]',va='center',fontsize=10,weight='bold',color=INK)
ax.text(.55,0,f'Prediction [{PI[0]:+.2f}, {PI[1]:+.2f}]',va='center',fontsize=10,weight='bold',color=INK)
ax.set_yticks(list(y)+[1,0],name+['Random-effects pool','95% prediction interval'])
ax.set_ylim(-.7,k+2)
ax.set_xlim(-.15,.53)
ax.set_xlabel('Correlation with social-desirability measure')
ax.set_title('A small average masks variation across studies',fontsize=16,weight='bold',color=INK,pad=24)
ax.grid(axis='x',alpha=.18);ax.set_axisbelow(True)
ax.text(.55,k+1.3,'Correlation [95% CI]     weight',fontsize=10,color=GRAY,weight='bold')
fig.subplots_adjust(left=.31,right=.75,top=.9,bottom=.1)
fig.savefig(OUT/'figure_01_meta_analysis_v4.png',bbox_inches='tight');plt.close(fig)

# 2 -- DIFFERENCES and CIs in Krueger & Forest Study 2, derived from reported mean differences / t / df.
# Perceived relationship quality: difference 1.23, t(217)=7.39.
# Perceived romantic receptiveness: difference -0.60, t(215)=3.35.
DD=[('Perceived high-quality\nrelationship',1.23,7.39,217),('Perceived receptiveness\nto romantic advances',-.60,3.35,215)]
fig,ax=plt.subplots(figsize=(10.5,4.5),dpi=180)
for i,(label,diff,tt,df) in enumerate(DD):
    se_diff=abs(diff/tt);delta=t.ppf(.975,df)*se_diff
    yy=1-i
    ax.hlines(yy,diff-delta,diff+delta,lw=3,color=BLUE)
    ax.plot([diff-delta,diff+delta],[yy,yy],'|',color=BLUE,ms=17,mew=2)
    ax.scatter(diff,yy,s=125,color=INK,zorder=3)
    ax.text(1.83,yy,f'{diff:+.2f}   [{diff-delta:+.2f}, {diff+delta:+.2f}]',va='center',fontsize=12,color=INK)
ax.axvline(0,ls='--',color=ORANGE,lw=1.5)
ax.set_yticks([1,0],[x[0] for x in DD]);ax.set_ylim(-.55,1.6)
ax.set_xlim(-1.15,1.75);ax.grid(axis='x',alpha=.18)
ax.set_xlabel('Difference in perceived rating: displays minus control (1–7 scale)')
ax.set_title('Relationship cues changed observers’ perceptions',fontsize=16,weight='bold',color=INK,pad=22)
ax.text(1.83,1.37,'Difference   [approx. 95% CI]',fontsize=11,color=GRAY,weight='bold')
fig.subplots_adjust(left=.30,right=.74,top=.80,bottom=.2)
fig.savefig(OUT/'figure_03_relationship_display_v4.png',bbox_inches='tight');plt.close(fig)

# 3 -- correlation-scaled results in August belief table. Only these three outcomes
# are correlations/point-biserial r reconstructions; three other rows are mean differences.
B=[('Age',.11,-.02,.23),('Social activities',-.03,-.15,.09),('Extraversion',-.11,-.23,.01)]
fig,ax=plt.subplots(figsize=(10,4.2),dpi=180)
for i,(label,x,lo,hi) in enumerate(B):
    yy=2-i
    ax.hlines(yy,lo,hi,lw=3,color=BLUE)
    ax.plot([lo,hi],[yy,yy],'|',color=BLUE,ms=15,mew=2)
    ax.scatter(x,yy,s=95,color=INK,zorder=3)
    ax.text(.315,yy,f'{x:+.2f}   [{lo:+.2f}, {hi:+.2f}]',va='center',color=GRAY,fontsize=11)
ax.axvline(0,color=ORANGE,ls='--',lw=1.5)
# +/-0.18 illustrates approximately 80% detectability at n=252, NOT an importance threshold.
for xc in (-.18,.18):
    ax.axvline(xc,ls=(0,(3,4)),lw=1.3,color='#8b99a5',alpha=.8)
ax.set_yticks([2,1,0],[d[0] for d in B]);ax.set_xlim(-.3,.30);ax.set_ylim(-.55,2.6)
ax.grid(axis='x',alpha=.2)
ax.set_xlabel('Approximate correlation with red/black-pill belief score')
ax.set_title('Uncertainty in three belief-score associations',fontsize=16,weight='bold',color=INK,pad=22)
ax.text(.315,2.38,'Estimate [95% CI]',fontsize=10.5,color=GRAY,weight='bold')
ax.text(0.5,-.32,'Dashed reference lines: approx. 80% power at |r|≈.18 (n≈252); not a practical-importance threshold',transform=ax.transAxes,ha='center',fontsize=10,color=GRAY)
fig.subplots_adjust(left=.24,right=.7,top=.82,bottom=.29)
fig.savefig(OUT/'figure_04_belief_correlations_v4.png',bbox_inches='tight');plt.close(fig)

# 4 -- complement plotted DIRECTLY: event of no public partner display, hence LR.
# Assumed rates are entirely hypothetical; no empirical base rate implied.
A_success_no=.10;A_unsuccess_no=.90;B_success_no=.80;B_unsuccess_no=.90
LRs=[A_success_no/A_unsuccess_no,B_success_no/B_unsuccess_no]
fig,ax=plt.subplots(figsize=(10.3,5.2),dpi=180)
x=np.arange(2);width=.34
bars1=ax.bar(x-width/2,[10,80],width,label='Genuinely successful (hypothetical)',color=BLUE)
bars2=ax.bar(x+width/2,[90,90],width,label='Unsuccessful (hypothetical)',color='#dfad97')
for b in list(bars1)+list(bars2):
    ax.text(b.get_x()+b.get_width()/2,b.get_height()+2, f'{b.get_height():.0f}%',ha='center',fontsize=12,weight='bold',color=INK)
ax.set_ylim(0,117);ax.set_yticks([0,20,40,60,80,100]);ax.set_ylabel('No public display (%)')
ax.set_xticks(x,[f'Scenario A: strong evidence\nLR(no display) = {LRs[0]:.2f}',f'Scenario B: weak evidence\nLR(no display) = {LRs[1]:.2f}'])
ax.set_title('The same absence can be highly or weakly informative',fontsize=16,weight='bold',color=INK,pad=22)
ax.legend(loc='upper center',ncol=2,frameon=False,bbox_to_anchor=(.5,1.01),fontsize=10)
ax.grid(axis='y',alpha=.19);ax.set_axisbelow(True)
ax.text(.5,-.23,'Likelihood ratio below 1 favors the unsuccessful hypothesis; all rates are illustrative.',transform=ax.transAxes,ha='center',color=ORANGE,fontsize=10.5)
fig.subplots_adjust(bottom=.30,left=.13,right=.97,top=.84)
fig.savefig(OUT/'figure_02_missing_displays_v4.png',bbox_inches='tight');plt.close(fig)

# Pairwise Dunn z-tests reproduced from the JANUARY report's rounded z statistics.
# Conventional two-sided Bonferroni adjustment among three pairwise contrasts.
PAIRS=[('Lifetime partners: Yes vs Other',-2.30,.030),('Lifetime partners: Yes vs No',-.04,1.000),('Lifetime partners: Other vs No',-1.64,.152),('Dating-app success: Other vs No',-3.02,.004),('Dating-app success: Yes vs No',-.22,1.000),('Dating-app success: Yes vs Other',2.67,.011)]
with open(OUT/'dunn_pairwise_recalculation.csv','w',newline='',encoding='utf-8') as f:
    w=csv.writer(f);w.writerow(['comparison','reported_z_rounded','reported_adjusted_p','one_tail_bonferroni_from_z','two_sided_bonferroni_from_z'])
    for label,zval,p in PAIRS:
        one=min(1,3*norm.sf(abs(zval)));two=min(1,6*norm.sf(abs(zval)))
        w.writerow([label,zval,p,round(one,5),round(two,5)])
        print(label,'reported:',p,'one-tail:',round(one,4),'two-sided:',round(two,4))
print('All output saved in',OUT)
