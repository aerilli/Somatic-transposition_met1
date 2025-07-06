#!/usr/bin/env python
# coding: utf-8

# # 

# In[27]:


import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from random import choices # random sampling with replacement
#import matplotlib.backends.backend_pdf


# In[3]:


file="/ebio/abt6_projects/met1_somatic_transpositions/data/3_somatic-events/r_mutations/insertion-bias/bootstrap/allfeat.csv"
df = pd.read_csv(file,delimiter="\t")
df


totreads=4218210
dflen=len(df)
remainingreads=totreads-dflen


data = {'fam': ["pos"]*remainingreads,
        'copy': ["pos"]*remainingreads,
        'V2': ["pos"]*remainingreads,
        'V3': ["pos"]*remainingreads,
        'V4': ["pos"]*remainingreads,
       }


# In[4]:


mdf=pd.concat([pd.DataFrame.from_dict(data), df])


# In[5]:


A= pd.DataFrame.from_dict({'fam': [],
        'copy': [],
        'V2': [],
        'V3': [],
        'V4': [],
    'V5': [],
       })


# In[20]:


#get_ipython().run_cell_magic('time', '', '\ngene=[]\nnonessential=[]\nessential=[]\nintergenic=[]\nTE=[]\nRepeat=[]\n\nfor i in range(1000):\n    mdfs=mdf.sample(n=totreads,replace=True)\n    mdfs=mdfs[mdfs[\'fam\'] != \'pos\']\n\n    \n    gene.append(len(mdfs[mdfs["V3"]==\'gene\']))\n    nonessential.append(len(mdfs[mdfs["V2"]==\'nonessential\']))\n    essential.append(len(mdfs[mdfs["V2"]==\'essential\']))\n    intergenic.append(len(mdfs[mdfs["V3"]==\'intergenic\']))\n    TE.append(len(mdfs[mdfs["V3"]==\'TE\']))\n    Repeat.append(len(mdfs[mdfs["V3"]==\'Repeat\']))\n        \n    #V5={"V5": [i]*int(len(mdfs))}\n    \n    #mdfs=mdfs.assign(**V5)\n\ndata = {\n    "gene" : gene,\n"nonessential" : nonessential,\n"essential" : essential,\n"intergenic" : intergenic,\n"TE" : TE,\n"Repeat" : Repeat,\n}\n\ndf = pd.DataFrame(data)\n')
gene=[]
nonessential=[]
essential=[]
intergenic=[]
TEs=[]
Repeat=[]

for i in range(1000):
    mdfs=mdf.sample(n=totreads,replace=True)
    mdfs=mdfs[mdfs['fam'] != 'pos']
    
    gene.append(len(mdfs[mdfs["V3"]=='gene']))
    nonessential.append(len(mdfs[mdfs["V2"]=='nonessential']))
    essential.append(len(mdfs[mdfs["V2"]=='essential']))
    intergenic.append(len(mdfs[mdfs["V3"]=='intergenic']))
    TEs.append(len(mdfs[mdfs["V3"]=='TE']))
    Repeat.append(len(mdfs[mdfs["V3"]=='Repeat']))
        
    #V5={"V5": [i]*int(len(mdfs))}
    print(i)
    #mdfs=mdfs.assign(**V5)

data = {
    "gene" : gene,
"nonessential" : nonessential,
"essential" : essential,
"intergenic" : intergenic,
"TE" : TEs,
"Repeat" : Repeat,
}

df = pd.DataFrame(data)


# In[42]:


df
df.to_csv('/ebio/abt6_projects/met1_somatic_transpositions/doc/jupyter/bootstrap_insertion-preference.simulated.txt', sep="\t") 


# In[41]:


for i in list(df):
    fig = plt.figure()
    sns.histplot(data=df, x=i, binwidth=5,
                 #palette=sns.color_palette("Paired"),
                 ec='k', stat='probability')
    ax = plt.gca()
    ax.set_ylim([0, 180])
    
    sns.set(style='white')
    plt.show()
    fig.savefig("/ebio/abt6_projects/met1_somatic_transpositions/data/Figures/pltout/prob-distr_feat-bias-insertions.%s.pdf" % i)
    
    


# In[25]:

print("FIN")


