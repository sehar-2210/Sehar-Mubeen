import os
import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import classification_report, confusion_matrix, accuracy_score
import matplotlib.pyplot as plt
import seaborn as sns

# Step 1: Define folders for ligands
drug_folder = "./drug_compounds"
non_drug_folder = "./non_drug_compounds"

# Step 2: Feature extraction function
def extract_features(mol):
    features = {
        'MolWt': Descriptors.MolWt(mol),
        'MolLogP': Descriptors.MolLogP(mol),
        'NumHDonors': Descriptors.NumHDonors(mol),
        'NumHAcceptors': Descriptors.NumHAcceptors(mol),
        'TPSA': Descriptors.TPSA(mol),
        'NumRotatableBonds': Descriptors.NumRotatableBonds(mol),
        'RingCount': Descriptors.RingCount(mol)
    }
    return features

# Step 3: Process SDF files and extract features
def process_folder(folder, label):
    data = []
    for file in os.listdir(folder):
        if file.endswith(".sdf"):
            suppl = Chem.SDMolSupplier(os.path.join(folder, file))
            for mol in suppl:
                if mol is not None:
                    features = extract_features(mol)
                    features['Label'] = label
                    data.append(features)
    return data

print("Extracting features from drug-like compounds...")
drug_data = process_folder(drug_folder, 1)

print("Extracting features from non-drug-like compounds...")
non_drug_data = process_folder(non_drug_folder, 0)

# Step 4: Combine and save dataset
df = pd.DataFrame(drug_data + non_drug_data)
df.to_csv("ligand_features.csv", index=False)
print("Dataset saved as ligand_features.csv")

# Step 5: Train ML model
X = df.drop(columns=["Label"])
y = df["Label"]

X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

model = RandomForestClassifier(n_estimators=100, random_state=42)
model.fit(X_train, y_train)

# Step 6: Evaluate the model
y_pred = model.predict(X_test)
print("\nClassification Report:")
print(classification_report(y_test, y_pred))
print("Confusion Matrix:")
print(confusion_matrix(y_test, y_pred))
print(f"Accuracy: {accuracy_score(y_test, y_pred)*100:.2f}%")

# Step 7: Visualization
plt.figure(figsize=(6,4))
sns.heatmap(confusion_matrix(y_test, y_pred), annot=True, fmt='d', cmap='Blues')
plt.title("Confusion Matrix")
plt.xlabel("Predicted")
plt.ylabel("Actual")
plt.tight_layout()
plt.savefig("confusion_matrix.png")
plt.show()
