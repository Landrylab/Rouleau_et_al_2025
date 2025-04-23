#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 26 15:24:34 2025

@author: motherfrank
"""
import seaborn as sns
print('Seaborn: ',sns.__version__)
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
print('Matplotlib: ',matplotlib.__version__)
import pandas as pd
print('Pandas: ',pd.__version__)
import sklearn
print('sklearn: ', sklearn.__version__)
import shap
print('SHAP: ',shap.__version__)
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score, balanced_accuracy_score, classification_report, confusion_matrix, PrecisionRecallDisplay, RocCurveDisplay, precision_recall_curve, auc
seed = 999
#%%
# Import data of interest, make list of columns and copy of the raw dataframe
df = pd.read_csv('./Total_dataframe_ML_v12_MTX.csv')
list_cols = list(df.columns)
X = df.copy()

# Drop columns that are no used for the model. These columns are present in the dataframe for subsequent analysis and display
droplist = ['mut_code', 'MTX_resistance', 'Resistant','aa_mut','aa_wt', 'random', 'Unnamed: 0', 'position', 'Function']
X = X.drop(droplist, axis = 1)


# Select column on which you want predictions to be made
y = df['Resistant']

# Splitting the dataset in the train and test sets
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=seed)
X_train_id, X_test_id, y_train_id, y_test_id = train_test_split(df, y, test_size=0.2, random_state=seed)
X_test_id = X_test_id.reset_index()
# z normalization and transformation
scaler = StandardScaler()
scaler.fit(X_train)
scalled_df = pd.DataFrame(scaler.transform(X), columns = X.columns)
X_train_s = pd.DataFrame(scaler.transform(X_train), columns = X_train.columns)
X_test_s = pd.DataFrame(scaler.transform(X_test), columns = X_test.columns)


# Define model using the best found parameters from previous RandomizedSearchCV and GridSearchCV.
# This script contains only the final model. RandomizedSearch results can diverge
# =============================================================================
model = RandomForestClassifier(max_depth = 23, 
                                         min_samples_split = 2, 
                                         n_estimators = 1500,
                                         min_samples_leaf = 20,
                                         max_features = 'sqrt',
                                         bootstrap = False,
                                         class_weight = 'balanced',
                                         random_state = seed).fit(X_train_s, y_train)

# Use these params on training set
train_pred = model.predict(X_train_s) # Predict class of train set
train_pred_prob = model.predict_proba(X_train_s) # Predict class probability of train set

test_pred = model.predict(X_test_s) # Predict class of test set
test_pred_prob = model.predict_proba(X_test_s) # Predict class probability of test set
# =============================================================================

# Display performance metrics
print("Balanced accuracy on train set: %0.1f %%" % (100 * balanced_accuracy_score(y_train, train_pred)))
print("Balanced accuracy on test set: %0.1f %%\n\n" % (100 * balanced_accuracy_score(y_test, test_pred)))
print("Accuracy on train set: %0.1f %%" % (100 * accuracy_score(y_train, train_pred)))
print("Accuracy on test set: %0.1f %%\n\n" % (100 * accuracy_score(y_test, test_pred)))
print(classification_report(y_train, train_pred))
print(classification_report(y_test, test_pred))
cm_train = confusion_matrix(y_train, train_pred)
our_cm = cm_train
TN = our_cm[0, 0]
TP = our_cm[1, 1]
FN = our_cm[1, 0]
FP = our_cm[0, 1]
acc = round((TP + TN) / np.sum(our_cm), 2)  # accuracy
tpr = round(TP / (TP + FN),2)              # true positive rate, sensitivity, recall
tnr = round(TN / (TN + FP),2)              # true negative rate, specificity 
ppv = round(TP / (TP + FP),2)            # positive predictive value, precision
npv = round(TN / (TN + FN),2)             # negative predictive value
print(f"Train set - Accuracy: {acc}, TPR: {tpr}, TNR: {tnr}, PPV: {ppv}, NPV: {npv}")
cm_test = confusion_matrix(y_test, test_pred)
our_cm = cm_test
TN = our_cm[0, 0]
TP = our_cm[1, 1]
FN = our_cm[1, 0]
FP = our_cm[0, 1]
acc = round((TP + TN) / np.sum(our_cm),2)  # accuracy
tpr = round(TP / (TP + FN),2)              # true positive rate, sensitivity, recall
tnr = round(TN / (TN + FP),2)              # true negative rate, specificity 
ppv = round(TP / (TP + FP),2)              # positive predictive value, precision
npv = round(TN / (TN + FN),2)              # negative predictive value
print(f"Test set - Accuracy: {acc}, TPR: {tpr}, TNR: {tnr}, PPV: {ppv}, NPV: {npv}")



#### Generate figures used for model performance display
## Confusion matrices for test and train
categories = ['Sensitive', 'Resistant']
fig, ax = plt.subplots(figsize=(4,4))
sns.set(font_scale=2)
sns.heatmap(cm_train, annot=True, cmap='Blues', fmt='', 
            xticklabels=categories,yticklabels=categories, 
            cbar=False)
plt.ylabel('True label', size = 22)
plt.xlabel('Predicted label', size = 22)
plt.title('Train set', size = 22)
plt.gcf().set_dpi(300)
ax.set_facecolor('white')
plt.savefig('./Figure3A_Test_CM.png', dpi = 300, format = 'png', bbox_inches = 'tight')
plt.show()


categories = ['Sensitive', 'Resistant']
fig, ax = plt.subplots(figsize=(4,4))
sns.set(font_scale=2)
sns.heatmap(cm_test, annot=True, cmap='Blues', fmt='', 
            xticklabels=categories,yticklabels=categories, 
            cbar=False)
plt.ylabel('True label', size = 22)
plt.xlabel('Predicted label', size = 22)
plt.gcf().set_dpi(300)
plt.title('Test set', size = 22)
ax.set_facecolor('white')
plt.savefig('./Figure3A_Train_CM.png', dpi = 300, format = 'png', bbox_inches = 'tight')
plt.show()


## ROC Curves
fig, ax = plt.subplots(figsize=(4,4))
rfc_disp = RocCurveDisplay.from_estimator(model, X_test_s, y_test, ax=ax, alpha=0.8)
plt.plot(ax=ax, alpha=0.8)
plt.title('Test ROC', size = 22)
ax.tick_params(axis = 'both', labelsize = 18)
plt.legend(fontsize = 10, loc = 8, facecolor = 'white', framealpha = 0)
plt.xlabel('False Positive Rate', size = 22)
plt.ylabel('True Positive Rate', size = 22)
ax.grid(False)
ax.set_facecolor('white')
plt.savefig('./Figure3A_Test_ROC.png', dpi = 300, format = 'png', bbox_inches = 'tight')
plt.show()

fig, ax = plt.subplots(figsize=(4,4))
rfc_disp = RocCurveDisplay.from_estimator(model, X_train_s, y_train, ax=ax, alpha=0.8)
plt.plot(ax=ax, alpha=0.8)
plt.title('Train ROC', size = 22)
ax.tick_params(axis = 'both', labelsize = 18)
plt.legend(fontsize = 10, loc = 8, facecolor = 'white', framealpha = 0)
plt.xlabel('False Positive Rate', size = 22)
plt.ylabel('True Positive Rate', size = 22)
ax.grid(False)
ax.set_facecolor('white')
plt.savefig('./Figure3A_Train_ROC.png', dpi = 300, format = 'png', bbox_inches = 'tight')
plt.show()


## PRD Curves
fig, ax = plt.subplots(figsize=(4,4))
rfc_disp = PrecisionRecallDisplay.from_estimator(model, X_test_s, y_test, ax=ax, alpha=0.8)
y_scores = model.predict_proba(X_test_s)[:, 1]
precision_test, recall_test, thresholds_test = precision_recall_curve(y_test, y_scores)
auc_test = auc(recall_test, precision_test)
plt.plot(ax=ax, alpha=0.8)
plt.title('Test PRD ', size = 22)
ax.tick_params(axis = 'both', labelsize = 18)
plt.legend(fontsize = 10, loc = 8, facecolor = 'white', framealpha = 0)
plt.xlabel('Recall', size = 22)
plt.ylabel('Precision', size = 22)
ax.grid(False)
ax.set_facecolor('white')
plt.savefig('./Figure3A_Test_PRD.png', dpi = 300, format = 'png', bbox_inches = 'tight')
plt.show()

fig, ax = plt.subplots(figsize=(4,4))
rfc_disp = PrecisionRecallDisplay.from_estimator(model, X_train_s, y_train, ax=ax, alpha=0.8)
y_scores = model.predict_proba(X_train_s)[:, 1]
precision_train, recall_train, thresholds_train = precision_recall_curve(y_train, y_scores)
auc_train = auc(recall_train, precision_train)
plt.plot(ax=ax, alpha=0.8)
plt.title('Train PRD', size = 22)
ax.tick_params(axis = 'both', labelsize = 18)
plt.legend(fontsize = 10, loc = 8, facecolor = 'white', framealpha = 0)
plt.xlabel('Recall', size = 22)
plt.ylabel('Precision', size = 22)
ax.grid(False)
ax.set_facecolor('white')
plt.savefig('./Figure3A_Train_PRD.png', dpi = 300, format = 'png', bbox_inches = 'tight')
plt.show()

#%%
# Run SHAP's explainer function on the model to see which parameters explain the best model performance
## The base explainer is used, and it will automatically call the tree explainer
# Close all possible previously open figures and reset default parameters
plt.close()
sns.reset_defaults()

# Run explainer based on the test set. This takes a while. 
explainer = shap.Explainer(model.predict, X_test_s)
shap_values = explainer(X_test_s)

#%%
# Display type explainer results. Having it in a different cell allows to modify the figure without the lengthy process of rerunning the explainer 
shap.summary_plot(shap_values, X_test_s, 
                             show = False, 
                             plot_size=[7,10])
plt.xticks(fontsize = 18)
plt.xlabel('SHAP Values', fontsize = 22)
# It is recommended to use the raw column names for the first itterations of the explainer, and to only change them once the 
# final model has been selected, where harmonizing labels can be useful for display.
ylabs = ['Distance to ligand', 'Functionnal', 'Buriedness', 'MutateX score', 'WT sidechain length', 
         'IQR change in distance', 'Allosteric effect', 'GEMME combined', 'GEMME epistasis', 
        'Catalytic', 'Median change in distance', 'Skew change in distance', 'GEMME individual', 'FlexddG total score']
plt.yticks(ticks = [13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0], labels = ylabs, fontsize = 22)
plt.savefig('./SHAP_explainer_final.png', dpi = 300, format = 'png', bbox_inches = 'tight')
plt.close()
#%%
# Re-assign calls made by the model to their original indexes
# When the dataset is split, it shuffles the indexes and resets them. 
# You need to manually go reassign the proper indexes to the calls that have been made. This is why its useful to keep the original dataframe intact
test_FP = pd.DataFrame()
y_test_index = np.array(y_test.index)
y_train_index = np.array(y_train.index)
test_FP['prediction_MTX'] = np.concatenate([test_pred,train_pred])
test_FP['prediction_MTX_0_prob'] = np.concatenate([test_pred_prob[:,0],train_pred_prob[:,0]])
test_FP['prediction_MTX_1_prob'] = np.concatenate([test_pred_prob[:,1],train_pred_prob[:,1]])
test_FP['prediction_index'] = np.concatenate([y_test_index,y_train_index])
test_FP.set_index('prediction_index', drop = True, inplace = True)
test_FP.sort_index(inplace = True)
X_out = df.copy()
X_out['mut_code'] = df['mut_code']
X_out['prediction'] = test_FP['prediction_MTX']
X_out['prediction_MTX_0_prob'] = test_FP['prediction_MTX_0_prob']
X_out['prediction_MTX_1_prob'] = test_FP['prediction_MTX_1_prob']

# Assign TP, FP, TN, FN values to original matrix depending on model performance
X_out['Label'] = ''
for index in X_out.index:
    info = X_out.iloc[index].to_list()
    real = int(info[7])
    call = int(info[23])
    # True negative
    if real == 0 and real == call:
        X_out.at[index,'Label'] = 'True Negative'
    # True positive
    elif real == 1 and real == call:
        X_out.at[index,'Label'] = 'True Positive'
    # False negative
    elif real == 1 and real != call:
        X_out.at[index,'Label'] = "False Negative"
    # False positive
    elif real == 0 and real != call:
        X_out.at[index,'Label'] = "False Positive"          

# Define boolean value for functionality, helps with plotting downstream figures
X_out['Functionnal'] = X_out.apply(lambda row: True if row['Function'] >= -0.1854 else False, axis = 1)


#%%
# Make waterfall plots for misclassified mutants during testing
# Ids are shuffled when sampling, IDs need to be found manually from intermediate dataframes
#G124P = 302
#S37P = 590
#V11E = 101
#V11M = 742
id_list = [101,302,590,742]
for idx in id_list:
    if idx == 302:
        mut = 'G124P'
    elif idx == 590:
        mut = 'S37P'
    elif idx == 101:
        mut = 'V11E'
    elif idx == 742:
        mut = 'V11M'
    
    # Make waterfalls plots for mutants of interest
    plt.rc('font', size=22)
    fig, ax = plt.subplots(figsize = (2,4))
    shap.plots.waterfall(shap_values[idx], max_display=10, show = False)
    ax.grid(False)
    ax.set_facecolor('white')
    plt.savefig('./Waterfall_MTXtest_'+mut+'.png', dpi = 300, format = 'png', bbox_inches = 'tight')
    plt.close()
    

#%%
### Use the generated model to make predictions using data generated using the PjDHFR_TMP complex
# Load data for new drug and pass it through all the same filter as the original dataset
### Very important that the new dataset has the same columns and columns names as the dataset used to train the model
data_TMP = pd.read_csv('./Total_dataframe_ML_v12_TMP.csv')
X_TMP = data_TMP.copy()
X_TMP = X_TMP.drop(droplist, axis = 1)
scaler = StandardScaler()
scaler.fit(X_TMP)
X_TMP = pd.DataFrame(scaler.transform(X_TMP), columns = X_TMP.columns)

# As column names are that same as in the original dataset, the columns are already set
y_TMP_prob = model.predict_proba(X_TMP)
y_TMP = model.predict(X_TMP)

#%%
# As this dataset did not need to be split, we can directly return calls to the dataframe
# If you want to run the TMP exapliner, skip this part, as it will add new columns to the TMP dataset, not recognized by the model
X_TMP['prediction'] = y_TMP
X_TMP['prediction_prob_0'] = y_TMP_prob[:,0]
X_TMP['prediction_prob_1'] = y_TMP_prob[:,1]

# X_TMP now contains all the calls from the model associated with their mutations
# However 

#%%
# Plot Resistant proediction probability with functionality. Interesting plot but did not make it to the final paper. Leaving it here for posterity.
data_TMP_called = pd.read_csv("./TMP_pred_v12_220125.csv")
fig, ax = plt.subplots(figsize = (7,8.45))
sns.scatterplot(x = data_TMP_called["prediction_prob_1"], y = df["Function"], 
                hue = df["min_dist_lig"], palette = "RdYlBu", legend = False, hue_norm = (3.5,20)
                , linewidth = 0, ax = ax)
sns.despine()
plt.xlabel("Resistant prediction probability", fontsize=22)
plt.ylabel("Selection coefficient - 37 °C", fontsize=22)
ax.tick_params(axis = "both", labelsize=18)
plt.axvline(0.5, color = "black", linestyle = ":")
plt.axhline(-0.185, color = "black", linestyle = ":")
norm = plt.Normalize(3.5, 20)
sm = plt.cm.ScalarMappable(cmap="RdYlBu", norm=norm)
cax = fig.add_axes([ax.get_position().x1+0.02,ax.get_position().y0,0.02,ax.get_position().height])
ax.figure.colorbar(sm, cax = cax)
cax.set_ylabel("Distance to folate (Å)", size = 22, labelpad=10)
cax.tick_params(labelsize = 18)
plt.show()
#%%
# Run expaines for TMP, and make waterfall plots for TMP classification of interest
# For overall model, summary plot is the same as with MTX, so not relevant to remake it
explainer = shap.Explainer(model.predict, X_TMP)
shap_values = explainer(X_TMP)

#%%
# Same issue as with MTX. Scalled data is de-indexed and shufled, so mutant IDs need to be recalled manually.
# In this case, data is not de-scaled as with MTX, but is less relevant as the model is already explained. 
#L65P = 2465
#A67V = 1967
#C166Y = 1333
#F36C = 345
for idx in [345,1333,1967,2465]:
    if idx == 2465:
        mut = 'L65P'
    elif idx == 1967:
        mut = 'A67V'
    elif idx == 1333:
        mut = 'C166Y'
    elif idx == 345:
        mut = 'F36C'
    
    plt.rc('font', size=22)
    fig, ax = plt.subplots(figsize = (2,4))
    shap.plots.waterfall(shap_values[idx], max_display=10, show = False)
    ax.grid(False)
    ax.set_facecolor('white')
    plt.savefig('./Waterfall_TMPpred_'+mut+'.png', dpi = 300, format = 'png', bbox_inches = 'tight', facecolor='white')
    plt.close()
    
    


