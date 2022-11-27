import numpy as np
from LUFactorization import LU_Factorization
from LUFactorization import LU_solve
from QRFactorization import QR_Factorization
from QRFactorization import QR_solve
from HousholderReduction import Householder_Reduction
from HousholderReduction import Householder_solve
from GivensReduction import Givens_Reduction
from GivensReduction import Givens_solve
from URVFactorization import URV_Factorization
from URVFactorization import URV_solve
from MatrixExample import matrix_option, solve_option


def stop():
    print("Exit the program!")
    quit() 


def check_square(matrix):
    '''
    Check if the matrix is square

    Parameters:
    matrix: the input matrix

    Returns
    True, False
    '''
    matrix = np.array(matrix)
    row_len = matrix.shape[0]
    col_len = matrix.shape[1]

    if row_len == col_len:
        return True
    else:
        return False

 
def det(matrix):
    '''
    Calculate the determinant
    
    Parameters:
    matrix: the input matrix

    Returns:
    the determinant of the input matrix
    '''
    if len(matrix) == 1:
        return matrix[0][0]
        
    s = 0
    for i in range(len(matrix)):
        A = [matrix[j][1:] for j in range(len(matrix)) if j != i]
        
        if i % 2:
            s -= matrix[i][0] * det(A)
        else:
            s += matrix[i][0] * det(A)

    return s


def choose(method):
    factorization_method_option = {'1': LU_Factorization,\
                                   '2': QR_Factorization,\
                                   '3': Householder_Reduction,\
                                   '4': Givens_Reduction,\
                                   '5': URV_Factorization,\
                                   '6': stop,}
    
    solve_method_option = {'1': LU_solve,\
                           '2': QR_solve,\
                           '3': Householder_solve,\
                           '4': Givens_solve,\
                           '5': URV_solve}

    factorization_method = factorization_method_option.get(method)
    
    if factorization_method == stop:
        factorization_method()
    elif factorization_method:
        solve_method = solve_method_option.get(method)

        print("Please type in the index (1-7) of the matrix defined in the MatrixExample.py!\n")
        matrix = input()
        A, b = matrix_option.get(matrix)

        if A: 
            if check_square(A) == False:
                print("The input matrix is not square, and thus will not calsulate the determinant")
            else:
                print("The determinant of the input matrix is: ", det(A))

            if b == []:
                factorization_method(A)
            else:
                solve_method(A, b)
            
        else:
            print("Invalid input!\nPlease check and type in again!")


    else:
        print("Invalid input!\nPlease check and type in again!")


if __name__ == "__main__":
    while True:
        print("----------------------------------------------")
        print("Please choose one method:\n",\
          "1. LU Factorization\n",\
          "2. QR Factorization\n", \
          "3. Householder Reduction\n",\
          "4. Givens Reduction\n",\
          "5. URV Reduction\n",\
          "6. Exit the program!\n")
        method = input()
        choose(method)












