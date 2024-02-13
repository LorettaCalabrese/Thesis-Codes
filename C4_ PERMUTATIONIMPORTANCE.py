

from keras.models import Sequential
from keras.layers import LSTM, Dense, Dropout
from sklearn.metrics import classification_report, accuracy_score
from sklearn.preprocessing import LabelEncoder
from sklearn.model_selection import train_test_split
from sklearn.inspection import permutation_importance
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from keras.optimizers import Adam


# Carica il dataset CSV
file_path = '/Users/L.Calabrese/Desktop/DATABASE_HEELCONTACT_HSP_final.csv'
data = pd.read_csv(file_path)

# Converti i valori della colonna 'Patology' da 0 e 1 a 'NORMATIVO' e 'PARAPARESI'
data['Patology'] = data['Patology'].replace({0: 'NORMATIVO', 1: 'PARAPARESI'})

# Prendi le features e il target
features = ['R_RF', 'R_VL', 'R_VM', 'R_G', 'R_BFCL', 'R_TA', 'R_PL', 'R_SO', 'R_GL', 'R_GMA']
target = 'Patology'

X = data[features].values  # Features
y = data[target].values  # Target

# Codifica il target in valori numerici
label_encoder = LabelEncoder()
encoded_targets = label_encoder.fit_transform(y)

# Dividi il dataset in training set e test set
X_train, X_test, y_train, y_test = train_test_split(X, encoded_targets, test_size=0.2, random_state=42)

# Definisci la funzione per creare il modello LSTM
def create_model():
    model = Sequential()
    model.add(LSTM(units=50, input_shape=(X_train.shape[1], 1)))  # Consideriamo 1 come numero di time steps
    model.add(Dropout(0.2))
    model.add(Dense(1, activation='sigmoid'))
    
    opt = Adam(learning_rate=0.031)
    model.compile(loss='binary_crossentropy', optimizer=opt, metrics=['accuracy'])
    return model

# Crea il modello LSTM
model = create_model()

# Addestra il modello
model.fit(X_train.reshape((X_train.shape[0], X_train.shape[1], 1)), y_train, epochs=20, batch_size=16, verbose=0)

# Valuta il modello sul test set
y_pred_prob = model.predict(X_test.reshape((X_test.shape[0], X_test.shape[1], 1)))
y_pred = (y_pred_prob > 0.5).astype(int)

# Calcola le metriche di performance
target_names = label_encoder.classes_
print("Metriche di performance generali:")
print(classification_report(y_test, y_pred, target_names=target_names))

# Calcola le metriche di performance per le due classi (macro e micro)
for class_idx in range(len(target_names)):
    target_class = label_encoder.inverse_transform([class_idx])[0]
    class_indices = np.where(y_test == class_idx)[0]
    y_test_class = np.zeros_like(y_test)
    y_test_class[class_indices] = 1
    
    y_pred_class = np.zeros_like(y_pred)
    y_pred_class[np.where(y_pred == class_idx)[0]] = 1
    
    print(f"Metriche di performance per la classe '{target_class}':")
    print(classification_report(y_test_class, y_pred_class, target_names=[f'not {target_class}', target_class]))

# Calcola l'accuratezza del modello sul set di test originale
y_pred_original = model.predict(X_test.reshape((X_test.shape[0], X_test.shape[1], 1)))
accuracy_original = accuracy_score(y_test, (y_pred_original > 0.5).astype(int))

# Calcola l'importanza delle caratteristiche
feature_names = features
feature_importance = {}

for feature_idx in range(X_test.shape[1]):
    X_test_shuffled = X_test.copy()
    np.random.shuffle(X_test_shuffled[:, feature_idx])  # Mescola una caratteristica
    
    y_pred_shuffled = model.predict(X_test_shuffled.reshape((X_test_shuffled.shape[0], X_test_shuffled.shape[1], 1)))
    accuracy_shuffled = accuracy_score(y_test, (y_pred_shuffled > 0.5).astype(int))
    
    feature_importance[feature_names[feature_idx]] = accuracy_original - accuracy_shuffled

# Grafico dell'importanza delle caratteristiche
sorted_features = sorted(feature_importance.items(), key=lambda x: x[1], reverse=True)
plt.figure(figsize=(10, 6))
plt.barh(range(len(sorted_features)), [val for _, val in sorted_features], align='center')
plt.yticks(range(len(sorted_features)), [key for key, _ in sorted_features])
plt.xlabel('Importance')
plt.title('Feature Importance (Permutation)')
plt.show()
