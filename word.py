# -*- coding: utf-8 -*-
"""
Created on Sun Jun  7 22:42:12 2020

@author: zzk
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

import seaborn as sns    
from scipy import stats
from scipy.stats import  norm
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import mean_absolute_error
from sklearn import linear_model, svm, gaussian_process
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import train_test_split

import warnings
warnings.filterwarnings('ignore')



# 蛋白质相互作用DIP数据集：酿酒酵母蛋白质的相互作用数据 该文件一共有24743条相互作用数据，共有4846个不同的蛋白质参与了相互作用。
data_dip = pd.read_csv("./protein_dataset/DIP20101010.txt",sep='\\s+', header=None)
# 基因表达水平数据集 酿酒酵母蛋白质在36个观测时点的基因表达水平值,值越高,说明该蛋白质越活跃
data_act= pd.read_csv("./protein_dataset/gene_express_data.txt",sep='\\s+', header=None,index_col=0)
# 已知的酿酒酵母关键蛋白质数据集 记录了1285个已知关键蛋白质
data_test= pd.read_csv("./protein_dataset/Essential.txt",sep='\\s+', header=None)

# 最终目标 找出关键蛋白质 
# step1 构建特征值  
# 1. 相互作用数 2. 活跃度的平均值  3.（相互作用数 / 活跃数值的平均数）
# 4.活跃的最小值  5.活跃度的最大值 6.相互作用数 * 活跃数值的平均数）
#  7.活跃度的方差 8.（相互作用数 / 活跃度的方差）9. （相互作用数 * 活跃度的方差）

# 1
feature=pd.DataFrame(np.zeros(data_act[1].shape,float),columns=['1'],index=data_act.index)
#feature = pd.concat([feature,tp2],axis=1,join='inner')
count=pd.concat([data_dip[0],data_dip[1]]).value_counts()
for name in count.index.values:
    if(name in feature.index.values):
        feature.loc[name,'1']+=count[name]

# 2, 3,4,5,6,7,8

#feature.loc['Q0085','2']=0.156
index=0
for name in  feature.index.values:
    feature.loc[name,'2'] = data_act.iloc[index].mean()
    feature.loc[name,'3'] = feature.loc[name,'1'] / feature.loc[name,'2']
    feature.loc[name,'4'] = data_act.iloc[index].min()
    feature.loc[name,'5'] = data_act.iloc[index].max()
    feature.loc[name,'6'] = feature.loc[name,'1'] * feature.loc[name,'2']
    feature.loc[name,'7'] = data_act.iloc[index].var()
    feature.loc[name,'8'] = feature.loc[name,'1'] / feature.loc[name,'7']
    feature.loc[name,'9'] = feature.loc[name,'1'] * feature.loc[name,'7']
    index+=1
    
# step2 寻找关键特征值
# 热力图
corrmat= feature.corr()
plt.subplots( figsize=(15,12) )
sns.heatmap(corrmat, vmax=0.9, square=True)

sns.pairplot(feature, height = 2.5)
plt.show();

cols=['2','4','5','6','7']

    
    

