import numpy as np 
import make_data as md
 
# This is a list of parameters for function generate_data (in a correct order)
# N - number of countries
# T - numer of years 
# params_Y  - parameters which multiply variables in the equation for Y (a vector with first value being a multiplyer for Y_t-1, next one for X_t, etc)
# params_X  - parameters which multiply variables in the equation for X-es (a matrix with first row being the parameters which multiply variables in the equation for X_1 etc)
# ni_std_dev - standard deviation for ni
# eta_std_dev - standard deviation for eta
# c_mul - multiplyer in the eqation
# c_std_dev - standard deviation for c
# std_devs - a list of standard deviations for variables taken from the normal distribution that are being added to X-es

exogenous_data_generator = md.generate_data(4, 3, np.array([1, 1, 1.5, -1, 2, 0, 0]), np.matrix([[0, 1, 0, 0, 0, 0, 0], [0, 0, 1, 0, 0, 0, 0], [0, 0, 0, 1, 0, 0, 0], [0, 0, 0, 0, 1, 0, 0], [0, 0, 0, 0, 0, 1, 0], [0, 0, 0, 0, 0, 0, 1]]), 0.2, 0.3, 0.2, 0.1, np.array([0.1, 0.2, 0.3, 0.1, 0.2, 0.1, 0.1]))
relatively_exogenous_data_generator = md.generate_data(4, 3, np.array([1, 1, 1.5, -1, 2, 0, 0]), np.matrix([[0.5, 0, 0, 0, 0, 0, 0], [0, 0, 0, -0.5, 0, 0, 0], [0, 0, 0, 1, 0, 0, 0], [0, 0, 0, 0, 1, 0, 0], [0, 0, 0, 0, 0, 1, 0], [0, 0, 0, 0, 0, 0, 1]]), 0.2, 0.3, 0.2, 0.1, np.array([0.1, 0.2, 0.3, 0.1, 0.2, 0.1, 0.1]))
strongly_exogenous_data_generator = md.generate_data(4, 3, np.array([1, 1, 1.5, -1, 2, 0, 0]), np.matrix([[0.5, 1.0, 0, 0, 0, 2, 0], [0, 0, 0, -0.5, 1.5, 0, -1.0], [0.5, 2.0, 0, 0, 0, 1.0, 0], [0, 0, -0.5, 0, 0, 0, 2.0], [0, 0, 0, 0, 0, 1, 0], [0, 0, 0, 0, 0, 0, 1]]), 0.2, 0.3, 0.2, 0.1, np.array([0.1, 0.2, 0.3, 0.1, 0.2, 0.1, 0.1]))
test_data_generator = md.generate_data(4, 3, np.array([1, 1, 1.5, -1, 2, 0, 0]), np.matrix([[0, 1, 0, 0, 0, 0, 0], [0, 0, 2, 0, 0, 0, 0], [0, 0, 0, 1, 0, 0, 0], [0, 0, 0, 0, 1, 0, 0], [0, 0, 0, 0, 0, 1, 0], [0, 0, 0, 0, 0, 0, 1]]), 0, 0, 0, 0, np.array([0, 0, 0, 0, 0, 0, 0]))

# The below function takes 2 argiments 
# Name of the file you want to save your data in 
# The data generator you want to use (examples above)

md.write_data(
    "my_awesomest_data.csv",
    exogenous_data_generator
    )

# The result is in the form of Country, Timestamp, Y_t, Y_t-1, X_1_t, X_2_t, ... in each row