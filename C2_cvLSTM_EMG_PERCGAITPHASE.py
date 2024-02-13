
import numpy as np
import pandas as pd
from keras.models import Sequential
from keras.layers import LSTM, Dense, Dropout
from keras.wrappers.scikit_learn import KerasClassifier
from sklearn.model_selection import GridSearchCV
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import LabelEncoder
from sklearn.metrics import accuracy_score
from keras.optimizers import Adam, RMSprop

# Carica il dataset CSV
file_path = '/Users/L.Calabrese/Desktop/PHASES_DATABASE_PARAPARESI.csv'
data = pd.read_csv(file_path)

# Prendi le features e il target
features = ['R_RF_PERC', 'R_VL_PERC', 'R_VM_PERC', 'R_G_PERC', 'R_BFCL_PERC', 'R_TA_PERC', 'R_PL_PERC', 'R_SO_PERC', 'R_GL_PERC', 'R_GMA_PERC']
target = 'patology_phases'

X = data[features].values  # Features
y = data[target].values  # Target

# Codifica il target in valori numerici
label_encoder = LabelEncoder()
encoded_targets = label_encoder.fit_transform(y)

# Dividi il dataset in training set e test set
X_train, X_test, y_train, y_test = train_test_split(X, encoded_targets, test_size=0.2, random_state=42)
print(X_train.shape)

# Reshape di X_train
X_train = X_train.reshape((X_train.shape[0], 1, X_train.shape[1]))

# Definisci la funzione per creare il modello LSTM
def create_model(units=50, dropout_rate=0.2, optimizer='adam', learning_rate=0.001):
    model = Sequential()
    model.add(LSTM(units=units, input_shape=(X_train.shape[1], X_train.shape[2])))  # Consideriamo 1 come numero di time steps
    model.add(Dropout(dropout_rate))

    model.add(Dense(8, activation='softmax'))

    if optimizer == 'adam':
        opt = Adam(learning_rate=learning_rate)
    elif optimizer == 'rmsprop':
        opt = RMSprop(learning_rate=learning_rate)
    else:
        opt = optimizer

    model.compile(loss='categorical_crossentropy', optimizer=opt, metrics=['accuracy'])
    return model

# Crea il modello KerasClassifier
model = KerasClassifier(build_fn=create_model, verbose=0)

# Definisci i parametri da testare nella ricerca in griglia
param_grid = {
    'units': [50, 100],
    'dropout_rate': [0.2, 0.5],
    'optimizer': ['adam', 'rmsprop'],
    'learning_rate': [0.001, 0.01],
    'batch_size': [16, 32],
    'epochs': [5, 10, 20, 30]
}

# Esegui la ricerca in griglia per trovare i migliori iperparametri
grid = GridSearchCV(estimator=model, param_grid=param_grid, cv=3, scoring='accuracy', verbose=2)
grid_result = grid.fit(X_train, y_train)

# Stampa i risultati della ricerca in griglia
print("Miglior accuracy:", grid_result.best_score_)
print("Migliori parametri:", grid_result.best_params_)

# Valuta il modello migliore sul test set
best_model = grid.best_estimator_
y_pred = best_model.predict(X_test)
accuracy = accuracy_score(y_test, y_pred)
print("Accuracy sul test set:", accuracy)

import matplotlib.pyplot as plt

# Esegui l'addestramento del modello e salva la storia dell'addestramento
history = best_model.fit(X_train, y_train, epochs=grid_result.best_params_['epochs'], 
                         batch_size=grid_result.best_params_['batch_size'], verbose=0, validation_data=(X_test, y_test))

# Estrai l'andamento dell'accuracy dal training set e dal test set
train_accuracy = history.history['accuracy']
test_accuracy = history.history['val_accuracy']

# Crea un grafico per visualizzare l'andamento dell'accuracy
plt.plot(train_accuracy, label='Training accuracy')
plt.plot(test_accuracy, label='Test accuracy')
plt.xlabel('Epochs')
plt.ylabel('Accuracy')
plt.title('Training and Test Accuracy')
plt.legend()
plt.show()
