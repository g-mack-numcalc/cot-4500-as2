import numpy as np

np.set_printoptions(precision=7, suppress=True, linewidth=100)

def nevilles_method(x, y, point):
    """
    Calculates an approximation of the value of a function at a given point
    using Neville's method.

    Args:
    x (list): A list of x-coordinates of data points
    y (list): A list of y-coordinates of data points
    point (float): The point at which to approximate the value of the function

    Returns:
    float: The approximation of the function value at the given point
    """
    n = len(x)
    P = [[0] * n for i in range(n)]
    for i in range(n):
        P[i][0] = y[i]

    for j in range(1, n):
        for i in range(n - j):
            P[i][j] = ((point - x[i + j]) * P[i][j - 1] -
                       (point - x[i]) * P[i + 1][j - 1]) / (x[i] - x[i + j])

    return P[0][n - 1]

def divided_difference_table(x_points, y_points):
    # set up the matrix
    size: int = len(x_points)
    matrix: np.array = np.zeros((size, size))
    
    # fill the matrix
    for index, row in enumerate(matrix):
        row[0] = y_points[index]
        
    # populate the matrix
    for i in range(1, size):
        for j in range(1, i+1):
            # the numerator are the immediate left and diagonal left indices...
            numerator = matrix[i][j-1] - matrix[i-1][j-1]
            # the denominator is the X-SPAN...
            denominator = x_points[i] - x_points[i-j]
            operation = numerator / denominator
            matrix[i][j] = round(operation, 15)

    print(round(matrix[1][1], 15))
    print(round(matrix[2][2], 15))
    print(round(matrix[3][3], 15))
    return matrix

def get_approximate_result(matrix, x_points, value):
    # p0 is always y0 and we use a reoccuring x to avoid having to recalculate x 
    reoccuring_x_span = 1
    reoccuring_px_result = matrix[0][0]
    
    # we only need the diagonals...and that starts at the first row...
    for index in range(1, len(x_points)):
        polynomial_coefficient = matrix[index][index]
        # we use the previous index for x_points....
        reoccuring_x_span *= (value - x_points[index-1])
        
        # get a_of_x * the x_span
        mult_operation = polynomial_coefficient * reoccuring_x_span
        # add the reoccuring px result
        reoccuring_px_result += mult_operation
    
    # final result
    return reoccuring_px_result

def apply_div_dif(matrix: np.array):
    size = len(matrix)
    for i in range(2, size):
        for j in range(2, min(i+3, size)):
            # skip if value is prefilled (we dont want to accidentally recalculate...)
            if matrix[i][j] != 0:
                continue
            
            # get left cell entry
            left: float = matrix[i][j-1]
            # get diagonal left entry
            diagonal_left: float = matrix[i-1][j-1]
            # order of numerator is SPECIFIC.
            numerator: float = (left - diagonal_left)
            numerator2: float = (left - matrix[i-1][j-1])
            # denominator is current i's x_val minus the starting i's x_val....
            denominator = matrix[i][0] - matrix[i-j+1][0]
            if denominator == 0:
                matrix[i][j] = 0
            else:
                # something save into matrix
                operation = (numerator / denominator) if j == 2 else (numerator2 / denominator)
                matrix[i][j] = operation
    
    return matrix

def hermite_interpolation():
    x_points = [3.6, 3.8, 3.9]
    y_points = [1.675, 1.436, 1.318]
    slopes = [-1.195, -1.188, -1.182]
    # matrix size changes because of "doubling" up info for hermite 
    num_of_points = len(x_points)
    matrix = np.zeros((2*num_of_points, 2*num_of_points))
    # populate x values (make sure to fill every TWO rows)
    for i, x in enumerate(x_points):
        matrix[2*i][0] = x
        matrix[2*i+1][0] = x
    # prepopulate y values (make sure to fill every TWO rows)
    for i, y in enumerate(y_points):
        matrix[2*i][1] = y
        matrix[2*i+1][1] = y
    # prepopulate with derivatives (make sure to fill every TWO rows. starting row CHANGES.)
    for i, slope in enumerate(slopes):
        matrix[2*i+1][2] = slope
    filled_matrix = apply_div_dif(matrix)
    print(filled_matrix)

def cubic_spline_int():
    x = np.array([2, 5, 8, 10])
    y = np.array([3, 5, 7, 9])

    n = len(x)

    # Calculate the second derivatives of the function at the data points
    h = np.diff(x)

    A = np.zeros((n, n))
    A[0, 0] = 1
    A[n-1, n-1] = 1

    for i in range(1, n-1):
        A[i, i-1] = h[i-1]
        A[i, i] = 2 * (h[i-1] + h[i])
        A[i, i+1] = h[i]

    # Print the matrix A
    print(A)

    b = np.zeros(n)

    for i in range(1, n-1):
        b[i] = (3/h[i])*(y[i+1]-y[i]) - (3/h[i-1])*(y[i]-y[i-1])

    print(b)

    x = np.linalg.solve(A, b)

    print(x)

# #1
x = [3.6, 3.8, 3.9]
y = [1.675, 1.436, 1.318]
point = 3.7

result = nevilles_method(x, y, point)
print(result)

# #2
# point setup
x_points = [7.2, 7.4, 7.5, 7.6]
y_points = [23.5492, 25.3913, 26.8224, 27.4589]
divided_table = divided_difference_table(x_points, y_points)

# #3
approximating_x = 7.3
final_approximation = get_approximate_result(divided_table, x_points, approximating_x)
print(final_approximation)

# #4
hermite_interpolation()

# #5
cubic_spline_int()