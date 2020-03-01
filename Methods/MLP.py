import numpy as np

from keras.models import Sequential
from keras.layers import Dense, Dropout, Flatten
from keras.regularizers import l2
from keras.callbacks import EarlyStopping, Callback
from keras.layers import Conv2D, MaxPooling2D
from keras import backend as K
from sklearn.model_selection import train_test_split
from sklearn.metrics import r2_score
from scipy.stats import pearsonr


#Eaerly stop condition
early_stopper = EarlyStopping(monitor='val_loss', min_delta=0.1, patience=2, verbose=0, mode='auto')



def compile_model_mlp(geneparam, input_shape):
    # Compile a MLP model based on the given network structure

    # Get our network parameters.
    nb_layers = geneparam['nb_layers']
    nb_neurons = geneparam['nb_neurons']
    activation = geneparam['activation']
    optimizer = geneparam['optimizer']
    dropout = geneparam['dropout']
    weight_decay = geneparam['weight_decay']
    print("Architecture:%d,%s,%s,%d,%.2f%%,%.2f%%" % (nb_neurons, activation, optimizer,
                                                      nb_layers, dropout, weight_decay))

    np.random.seed(1337)
    model = Sequential()

    # Add each layer.
    for i in range(nb_layers):

        # Need input shape for first layer.
        if i == 0:
            if weight_decay > 0:
                model.add(Dense(nb_neurons, activation=activation, input_dim=input_shape,
                                kernel_regularizer=l2(weight_decay)))
            else:
                model.add(Dense(nb_neurons, activation=activation, input_dim=input_shape))

        else:
            if weight_decay > 0:
                model.add(Dense(nb_neurons, activation=activation, kernel_regularizer=l2(weight_decay)))
            else:
                model.add(Dense(nb_neurons, activation=activation))
        if dropout > 0:
            model.add(Dropout(dropout))  # dropout for each layer

    # Output layer.
    model.add(Dense(1))

    model.compile(loss='mse',
                  optimizer=optimizer,
                  metrics=['mae'])

    return model


def read_net_structure(model_file):
# extract the best MLP model structure of a trait
    with open(model_file) as f:
        line = f.readline()
        net_struc = line.split("\t")[0]
        print("Network structure: {}".format(net_struc))
        net_struc = eval(net_struc)
    return net_struc


def mlp_train(model_file,x_train,y_train,x_test,y_test):
    # MLP molde given the best network structure of the trait

    ##X_train: training genotype data  (Numpy matrix - samples X viriants)
    # x_test: testing genotype data
    # y_train: training trait value data (Numpy vector - 1 X N)
    # x_test: testing trait value data
    # model_file: the file stores the top network structures selected in validation step
    # return the learned model, r and r2 performance of the model on testing data

    #Using 10% of the training samples as validation data to control convengence
    x_train, x_val, y_train, y_val = train_test_split(x_train, y_train, test_size=0.1, random_state=42)

    input_shape = x_train.shape[1]

    #read the optimal MLP network structure of a trait
    struc_param = read_net_structure(model_file)

    #complile the model with the given network strucure
    model = compile_model_mlp(struc_param, input_shape)

    model.fit(x_train, y_train,
              epochs=1200,
              verbose=2,
              validation_data=(x_val, y_val),
              callbacks=[early_stopper])

    r = pearsonr(y_test,model.predict(x_test).ravel())[0]
    r2 = r2_score(y_test, model.predict(x_test).ravel())

    return model,r2,r

