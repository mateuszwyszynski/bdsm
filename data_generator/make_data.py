import numpy as np
import csv

# This is a function to calculate new Xit based on prevoius values of X and Y, X_it = params_i[0]*Yi(t-1) + params_i[1]*X_1i(t-1) + params_i[2]*X_2i(t-1) + .... + c * ni + e_i
def multiply(params, Y, X, N, std_dev):
    X_new = params[0]*Y + np.random.normal(0, std_dev ,N)
    i = 0
    while i < len(params) - 1: 
        X_new = X_new + params[i+1]*X[:,i]
        i = i + 1
    return X_new

# This is a function to calculate new matrix X based on prevoius values of X and Y
def generate_X(N, Y_prev, X_prev, params_X, ni, c, std_devs):
    X_new = np.empty(X_prev.shape)
    X_new[:] = np.nan
    i = 0
    while i < X_prev.shape[1]:
        X_new[:,i] = multiply(params_X[i,:].T, Y_prev, X_prev, N, std_devs[i+1]) + c * ni # TODO need to think about indexes more
        i = i + 1
    return X_new

# N - number of countries
# T - numer of years 
# params_Y  - parameters which multiply variables in the equation for Y (a vector with first value being a multiplyer for Y_t-1, next one for X_t, etc)
# params_X  - parameters which multiply variables in the equation for X-es (a matrix with first row being the parameters which multiply variables in the equation for X_1 etc)
# ni_std_dev - standard deviation for ni
# eta_std_dev - standard deviation for eta
# c_mul - multiplyer in the eqation
# c_std_dev - standard deviation for c
# std_devs - a list of standard deviations for variables taken from the normal distribution that are being added to X-es

def generate_data(N, T, params_Y, params_X, ni_std_dev, eta_std_dev, c_mul, c_std_dev, std_devs):
    Y_prev = np.random.normal(0, 1 ,N)
    X_prev = np.random.normal(0, 1 ,[N, len(params_Y)-1])
    ni = np.random.normal(0, ni_std_dev, N)
    c_prev = np.random.normal(0, c_std_dev)
    for ti in range(T):
        c = c_mul * c_prev + np.random.normal(0, c_std_dev)
        eta = np.random.normal(0, eta_std_dev)
        X = generate_X(N, Y_prev, X_prev, params_X, ni, c, std_devs)
        Y = multiply(params_Y, Y_prev, X, N, std_devs[0]) + ni + np.array([eta] * N)
        yield X, Y, Y_prev
        Y_prev = Y
        X_prev = X
        c_prev = c

def write_data(file_name, data_generator):
    with open(file_name, "w") as f:
        writer = csv.writer(f)
        for i, (X,Y, Y_prev) in enumerate(data_generator):
            no_of_countries = X.shape[0]
            j = 0
            while j < no_of_countries:
                writer.writerow([j, i, Y[j], Y_prev[j]] + X[j,:].tolist())
                j = j+1