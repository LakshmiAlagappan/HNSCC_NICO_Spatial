import anndata as ad
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.metrics import classification_report, confusion_matrix
from sklearn.preprocessing import OneHotEncoder, OrdinalEncoder
from sklearn.preprocessing import StandardScaler
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import roc_auc_score, roc_curve
import scanpy as sc
from sklearn.decomposition import PCA
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.metrics import accuracy_score

from sklearn.linear_model import LogisticRegression
import numpy as np
import logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'  # Custom date format
)

def pp_X_quant_train(X):
    logging.info("### Standard Scaling: Quant - Train ###")
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)
    X_scaled = pd.DataFrame(X_scaled, columns=X.columns)
    return X_scaled,scaler

def pp_X_quant_test(X, scaler):
    logging.info("### Standard Scaling: Quant - Test ###")
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)
    X_scaled = pd.DataFrame(X_scaled, columns=X.columns)
    return X_scaled

def pp_X_categ_train(X, type):
    logging.info("### X Encoding: Categ Train ###")
    if (type == "OHE"):
        logging.info("### One Hot Encoder ###")
        encoder = OneHotEncoder(drop='first', handle_unknown='ignore')
        X_encoded = encoder.fit_transform(X).toarray()
        X_encoded = pd.DataFrame(X_encoded,columns=encoder.get_feature_names_out())
    if (type == "OE"):
        logging.info("### Ordinal Encoder ###")
        encoder = OrdinalEncoder(handle_unknown='use_encoded_value', unknown_value=-1)
        X_encoded = encoder.fit_transform(X)
        X_encoded = pd.DataFrame(X_encoded, columns=X.columns)
    return X_encoded, encoder

def pp_X_categ_test(X, encoder, type):
    logging.info("### X Encoding: Categ Test ###")
    if (type == "OHE"):
        logging.info("### One Hot Encoder ###")
        X_encoded = encoder.transform(X).toarray()
        X_encoded = pd.DataFrame(X_encoded,columns=encoder.get_feature_names_out())
    if (type == "OE"):
        logging.info("### Ordinal Encoder ###")
        X_encoded = encoder.transform(X)
        X_encoded = pd.DataFrame(X_encoded,columns=X.columns)
    return X_encoded

def pp_y_train(y):
    logging.info("### Y Encoding: Train ###")
    logging.info("### One Hot Encoder ###")
    encoder = OneHotEncoder(drop='first')
    y_encoded = encoder.fit_transform(y.reshape(-1, 1)).toarray().reshape((len(y),))
    y_encoded = pd.DataFrame(y_encoded, columns=encoder.get_feature_names_out())
    if y_encoded.shape[1] == 2:
        y_encoded = y_encoded[:, 1]  # shape (n_samples,)
    return y_encoded, encoder

def pp_y_test(y, encoder):
    logging.info("### Y Encoding: Test ###")
    y_encoded = encoder.transform(y.reshape(-1, 1)).toarray().reshape((len(y),))
    y_encoded = pd.DataFrame(y_encoded, columns=encoder.get_feature_names_out())
    if y_encoded.shape[1] == 2:
        y_encoded = y_encoded[:, 1]
    return y_encoded


def classify_tcells(tcells_harm, metadata):
    logging.info("### Classifying the double negatives in the T Cells based on Random Forest ###")
    tcells_harm_original = tcells_harm.copy()
    tcells_known = ad.concat([tcells_harm[(tcells_harm.obs['tcell_type'] == "CD4- CD8+")].copy(),
                    tcells_harm[(tcells_harm.obs['tcell_type'] == "CD4+ CD8-")].copy()])
    tcells_known.X[:, tcells_known.var_names.isin(['CD4', 'CD8A', 'CD8B'])] = 0
    tcells_unknown = tcells_harm[(tcells_harm.obs['tcell_type'] == "CD4- CD8-")].copy()

    X = pd.DataFrame(tcells_known.X.toarray(), columns=tcells_known.var_names, index=tcells_known.obs_names)
    for i in metadata:
        if i in tcells_known.obs.columns:
            X[i] = tcells_known.obs[i].values
    y = tcells_known.obs['tcell_type'].values

    quant_cols = X.select_dtypes(include=['number']).columns
    categ_cols = X.select_dtypes(include=['object', 'category']).columns
    print(f"Categ Columns: {categ_cols}")

    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42, stratify=y)
    print(f"\n Data Shapes:")
    print(f"  ├── X_train: {X_train.shape}")
    print(f"  └── X_test : {X_test.shape}")
    print(f"\n Label Distribution (y_train):")
    print(y_train.value_counts())
    print(f"\n Label Distribution (y_test):")
    print(y_test.value_counts())

    logging.info("### Preprocessing ###")
    X_tr_enc, scaler = pp_X_quant_train(X_train[quant_cols])
    X_te_enc = pp_X_quant_test(X_test[quant_cols], scaler)
    X_tr_categ_enc, encoder = pp_X_categ_train(X_train[categ_cols], "OE")
    X_te_categ_enc = pp_X_categ_test(X_test[categ_cols], encoder, "OE")
    X_train_enc = pd.concat([X_tr_enc, X_tr_categ_enc], axis=1)
    X_test_enc = pd.concat([X_te_enc, X_te_categ_enc], axis=1)
    y_train_enc, y_encoder = pp_y_train(y_train)
    y_test_enc = pp_y_test(y_test, y_encoder)

    logging.info("### Random Forest Classifier Built ###")
    rf = RandomForestClassifier(n_estimators=100, 
                                random_state=42)
    rf.fit(X_train_enc, y_train_enc)
    y_pred_rf = rf.predict(X_test_enc)

    from sklearn.metrics import confusion_matrix, classification_report

    print("\n📉 Random Forest Evaluation")
    print("-" * 60)
    print("🔹 Confusion Matrix:\n", confusion_matrix(y_test_enc, y_pred_rf))
    print("\n🔹 Classification Report:\n", classification_report(y_test_enc, y_pred_rf))
    print("=" * 60)

    importances = rf.feature_importances_
    feature_importance_df = pd.DataFrame({
        'Feature': X_train_enc.columns,
        'Importance': importances
    }).sort_values(by='Importance', ascending=False)

    #feature_importance_df.to_csv(res_dir + ex + "_feature_importance.csv", index=False)
    print(feature_importance_df.head(20))

    #plt.figure(figsize=(8, 12))
    #plt.barh(feature_importance_df['Feature'][0:50], feature_importance_df['Importance'][0:50], color='skyblue')
    #plt.xlabel('Importance')
    #plt.ylabel('Features')
    #plt.title('Top 50 Features by Importance')
    #plt.gca().invert_yaxis()  # Invert y-axis for better readability
    #plt.tight_layout()
    #plt.show()

    logging.info("### Retraining the Random Forest Classifier based on top features ###")
    top_n = 200  # Adjust this based on your dataset #100
    top_features = feature_importance_df.head(top_n)['Feature'].values
    X_top = X_train_enc[top_features]
    rf_model_retrained = RandomForestClassifier(n_estimators=1000,#370,
                                                min_samples_split=20,
                                                max_features='log2',
                                                criterion='entropy',
                                                class_weight='balanced',
                                random_state=42)
    rf_model_retrained.fit(X_top, y_train_enc)
    probabilities = rf_model_retrained.predict_proba(X_test_enc[top_features])
    prob = probabilities[:, 1]  # predicted probability of positive class
    auc_score = roc_auc_score(y_test_enc, prob)
    print(f"AUC Score: {auc_score:.2f}")

    fpr, tpr, thresholds = roc_curve(y_test_enc,prob)
    auc_table = pd.DataFrame({
        'False Positive Rate': fpr,
        'True Positive Rate': tpr,
        'Thresholds': thresholds
    })
    #auc_table.to_csv(res_dir + ex + "_auc_table.csv", index=False)

    y_pred_rf = rf_model_retrained.predict(X_test_enc[top_features])


    # best_thresh = 0
    # best_acc = 0

    # for thresh in thresholds:
    #     preds = (prob >= thresh).astype(int)
    #     acc = accuracy_score(y_test_enc, preds)
    #     if acc > best_acc:
    #         best_acc = acc
    #         best_thresh = thresh

    # print("Best threshold:", best_thresh)
    # print("Best accuracy:", best_acc)

    # y_pred_rf = (prob >= best_thresh).astype(int)
    print("\n📉 Retrained Random Forest Evaluation")
    print("-" * 60)
    print("🔹 Confusion Matrix:\n", confusion_matrix(y_test_enc, y_pred_rf))
    print("\n🔹 Classification Report:\n", classification_report(y_test_enc, y_pred_rf))
    print("=" * 60)

    

    logging.info("### Retrieving decision probabilities and AUC values ###")

    #plt.figure(figsize=(8, 6))
    #plt.plot(fpr, tpr, color='blue', label=f'ROC Curve (AUC = {auc_score:.2f})')
    #plt.plot([0, 1], [0, 1], color='red', linestyle='--')  # Diagonal line
    #plt.xlim([0.0, 1.0])
    #plt.ylim([0.0, 1.05])
    #plt.xlabel('False Positive Rate')
    #plt.ylabel('True Positive Rate')
    #plt.title('Receiver Operating Characteristic (ROC) Curve')
    #plt.legend(loc='lower right')
    #plt.grid()
    #plt.show()

    logging.info("### Testing the retrained RF on the DN Unknown ###")
    X = pd.DataFrame(tcells_unknown.X.toarray(), columns=tcells_unknown.var_names, index=tcells_unknown.obs_names)
    for i in metadata:
        if i in tcells_unknown.obs.columns:
            X[i] = tcells_unknown.obs[i].values

    X_enc = pp_X_quant_test(X[quant_cols], scaler)
    X_categ_enc = pp_X_categ_test(X[categ_cols], encoder, "OE")
    X_enc = pd.concat([X_enc, X_categ_enc], axis=1)
    X_enc = X_enc[top_features]
    y_pred_rf_proba = rf_model_retrained.predict_proba(X_enc)
    print(y_pred_rf_proba)
    y_pred_rf = rf_model_retrained.predict(X_enc)
    #y_pred_rf = (y_pred_rf_proba[:,1] >= best_thresh).astype(int)
    y_label = ["CD4- CD8+" if x == 1 else "CD4+ CD8-" for x in y_pred_rf]
    c = 0
    for i in range(len(y_label)):
        if max(y_pred_rf_proba[i])<0.52:
            c+=1
            y_label[i] = "CD4- CD8-"
    tcells_unknown.obs['tcell_type'] = y_label
    logging.info("### Using Threshold 0.52 to retain some as DN ###")

    tcells_harm_annotated = tcells_harm_original.copy()
    if not (tcells_unknown.obs.index == tcells_harm_annotated.obs.index[tcells_harm_annotated.obs['tcell_type'] == 'CD4- CD8-']).all():
        raise ValueError("Cell IDs do not align!")

    tcells_harm_annotated.obs['tcell_type_resolved'] = tcells_harm_annotated.obs['tcell_type'].copy()
    tcells_harm_annotated.obs.loc[tcells_harm_annotated.obs['tcell_type'] == 'CD4- CD8-', 'tcell_type_resolved'] = tcells_unknown.obs['tcell_type']

    #tcells_harm_annotated.write_h5ad(res_dir+"tcells_adata_harm_resolved.h5ad")
    print("\n🔬 T Cell Type Counts (Before Resolving))")
    print("-" * 40)
    print(tcells_harm_annotated.obs['tcell_type'].value_counts())

    print("\n✅ T Cell Type Counts (Resolved)")
    print("-" * 40)
    print(tcells_harm_annotated.obs['tcell_type_resolved'].value_counts())
    print("=" * 40)

    return tcells_harm_annotated

