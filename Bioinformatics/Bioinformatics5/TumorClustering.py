from sklearn.linear_model import LogisticRegression
import numpy as np 
import matplotlib.pyplot as plt


# importing data
# -----------------------------------------------
with open(r"C:\Users\tkoti\Desktop\bioinformatics\Bioinformatics5\colon_cancer.txt") as f:
	lines = f.readlines()

	lines = [[float(num) for num in item.split(' ')] for item in lines]

cancer_confirmed = np.array(lines)


with open(r"C:\Users\tkoti\Desktop\bioinformatics\Bioinformatics5\colon_healthy.txt") as f:
	lines = f.readlines()

	lines = [[float(num) for num in item.split(' ')] for item in lines]


healthy_confirmed = np.array(lines)


oracle_data = np.concatenate([cancer_confirmed,healthy_confirmed])


with open(r"C:\Users\tkoti\Desktop\bioinformatics\Bioinformatics5\colon_test.txt") as f:
	lines = f.readlines()

	lines = [[float(num) for num in item.split(' ')] for item in lines]


testdata = np.array(lines)

#Cancer 1, Healthy 0
targets = np.array([1] * 40 + [0] * 21)

#-------------------------------------



#Scaling to obtain homoscedastic parameters for Regression and Machine learning tasks
#-------------------------------------------------
from sklearn import preprocessing

scaler = preprocessing.StandardScaler().fit(oracle_data)

X_scaled = scaler.transform(oracle_data)

testdata_scaled = scaler.transform(testdata)
#-------------------------


'''
First method for determining if individual has Cancer
Logistic Regression
'''
clf = LogisticRegression(random_state = 0).fit(X_scaled, targets)
print(clf.predict(testdata_scaled), "Logistic regression: No Cancer")





'''
Second Method: Shallow MLP
'''

# MLP Layer -> Relu -> MLP Layer -> Relu... -> Softmax -> randomchoice | Backprop
#-------------------------------------------
from sklearn.neural_network import MLPClassifier
from sklearn.model_selection import train_test_split
X_train, X_test, y_train, y_test = train_test_split(X_scaled, targets, stratify=targets,
                                                    random_state=1)

clf2 = MLPClassifier(random_state = 1, max_iter = 200).fit(X_train, y_train)

print(clf2.score(X_test, y_test), "Model Score: 0.8125 (Model Probably Overfit bc too many parameters)")

probs = (clf2.predict_proba(testdata_scaled)).reshape(2,)
print(np.random.choice([0,1], p = probs), "MLP: No Cancer")

#------------------------------------------




'''
Third Method: Support Vector Machine
'''

#----------------------
from sklearn.svm import SVC
from sklearn.pipeline import make_pipeline

clf3 = make_pipeline(SVC(gamma = 'auto'))
clf3.fit(X_scaled,targets)
print(clf3.predict(testdata_scaled), "Support Vector Machine: No Cancer")

#-------------------------




'''
Fourth Method
K Nearest Neighbors and then UMAP dimensionality reduction for visual analysis
'''

#------------------------
import scanpy as sc 
import anndata as ad 


adatas = {}

cancer = ad.AnnData(cancer_confirmed)
cancer.obs_names = [f"C{i}" for i in range(cancer.n_obs)]
healthy = ad.AnnData(healthy_confirmed)
healthy.obs_names = [f"H{i}" for i in range(healthy.n_obs)]

unknown = ad.AnnData(testdata)

adatas['Cancerous'] = ad.AnnData(cancer)
adatas['Healthy'] = ad.AnnData(healthy)
adatas['Unknown'] = unknown
adata = ad.concat(adatas,label = "sample")



sc.pp.normalize_total(adata)
sc.pp.log1p(adata)


sc.pp.neighbors(adata, use_rep = 'X')

sc.tl.umap(adata)
sc.pl.umap(
    adata,
    color="sample",
    size=100,
)


#-----------------------------



'''
Displaying a sample decision Boundary
'''

#------------------

data = adata.obsm['X_umap']
params = adata.uns['umap']

# print(data.shape)

reduced_dim_Cancerous = data[:40,:]
reduced_dim_Healthy = data[40:-1, :]
reduced_dim_Unknown = data[-1,:]

X = reduced_dim_train = data[:-1,:]
y = targets
# print(reduced_dim_Unknown.shape, reduced_dim_Healthy.shape, reduced_dim_Cancerous.shape)






def make_meshgrid(x, y, h=.02):
    x_min, x_max = x.min() - 1, x.max() + 1
    y_min, y_max = y.min() - 1, y.max() + 1
    xx, yy = np.meshgrid(np.arange(x_min, x_max, h), np.arange(y_min, y_max, h))
    return xx, yy

def plot_contours(ax, clf, xx, yy, **params):
    Z = clf.predict(np.c_[xx.ravel(), yy.ravel()])
    Z = Z.reshape(xx.shape)
    out = ax.contourf(xx, yy, Z, **params)
    return out


model = SVC()
clf4 = model.fit(reduced_dim_train, targets)
fig, ax = plt.subplots()
# title for the plots
title = ('Decision surface of SVC ')
# Set-up grid for plotting.

X0, X1 = X[:, 0], X[:, 1]
xx, yy = make_meshgrid(X0, X1)

plot_contours(ax, clf4, xx, yy, cmap=plt.cm.coolwarm, alpha=0.8)
ax.scatter(X0[:40], X1[:40], c='red', s=20, edgecolors='k')
ax.scatter(X0[40:], X1[40:], c='blue', s=20, edgecolors='k')
ax.set_xlabel("UMAP1")
ax.set_ylabel("UMAP2")
ax.set_xticks(())
ax.set_yticks(())
ax.set_title(title)
ax.scatter(reduced_dim_Unknown[0], reduced_dim_Unknown[1], s = 20, edgecolors = 'k', marker = 'v', c = 'black')
ax.legend(['Cancer', 'Healthy', 'Unknown'])
plt.show()

#-------------------------------

#Summary
#Given the four methods, it seems like the unknown patient is healthy and doesn't have colon cancer given the reference datasets.
#The patient is sufficiently far away from the Cancerous Decision boundary