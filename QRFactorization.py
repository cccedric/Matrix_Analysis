import numpy as np
from numpy.linalg import matrix_rank
from numpy.linalg import norm


def check(matrix):
    '''
    Check if the number of the columns is equal to the rank

    Parameters:
    matrix: the input matrix

    Returns
    True, False
    '''
    rank = matrix_rank(matrix)
    n = matrix.shape[-1]
    if rank == n:
        return True
    else:
        return False

    
def check_unique_solution(A, b):
    '''
    Check if the Ax = b has unique solution

    Parameters:
    A, b

    Returns
    True, False
    '''
    Ab = np.c_[A, b.T]
    rank_A = matrix_rank(A)
    rank_Ab = matrix_rank(Ab)
    if rank_A == rank_Ab:
        return True
    else:
        return False

    
def QR_Factorization(matrix):
    '''
    QR Factorization

    Parameters:
    matrix: the input matrix

    Returns:
    Q, R
    '''
    np.set_printoptions(precision=3, suppress=True)

    matrix=np.array(matrix)
    print("The input matrix is: ")
    print(matrix)

    if check(matrix) == False:
        print("ERROR: The column vertor is not linearly independent.")
        return
    else:
        # Implement the Gram-Schmidt
        Q = np.zeros_like(matrix)

        col_index = 0
        for a in matrix.T:
            u = np.copy(a)
            for i in range(0, col_index):
                u = u - np.dot(np.dot(Q[:, i].T, a), Q[:, i])
            norm_factor = norm(u)
            Q[:, col_index] = u / norm_factor

            col_index += 1
        
        R = np.dot(Q.T, matrix)

        print("The result of the QR factorization is shown as follows: ")
        print("-----The Q matrix is: ")
        print(Q)
        print("-----The R matrix is: ")
        print(R)
        return (Q, R)


def QR_solve(A, b):
    '''
    Use LU Factorization to solve Ax = b
    
    Parameters:
    A, b

    Returns:
    x
    '''
    np.set_printoptions(precision=3, suppress=True)
    A = np.array(A, dtype=np.float32)
    b = np.array(b, dtype=np.float32)
    print("The input matrix A is: ")
    print(A)
    print("The input matrix b is: ")
    print(b)
    if check_unique_solution(A, b) == False:
        print("Error! The Ax = b does not have unique solution!")
    else:
        Q, R = QR_Factorization(A)
        n = len(A)
        # Solve Qy = b
        y = np.dot(Q.T, b)

        print("The solution of Qy = b is:")
        print(y)
        
        # Solve Rx = y
        x = np.zeros((n, 1))
        for i in range(len(A)-1, -1, -1):
            t = 0
            for j in range(i+1,len(A)):
                t += R[i][j] * x[j][0]
            t = y[i] - t
            if t != 0 and R[i][i] == 0:
                return 0
            x[i] = t / R[i][i]

        print("The solution of Rx = y is:")
        print(x)

        return x


if __name__ == "__main__":
    matrix = [[1.0, 2.0, 5.0], [4.0, 20.0, 6.0], [7.0, 1.0, 9.0]]
    QR_Factorization(matrix)
    print("----------------------------------------------")
    matrix = [[10.0, 2.0, 4.0], [6.0, 5.0, 2.0], [7.0, 15.0, 9.0]]
    QR_Factorization(matrix)
    print("----------------------------------------------")
    matrix = [[10.0, 2.0, 4.0, 8.0], [6.0, 5.0, 2.0, 1.0], [7.0, 15.0, 9.0, 6.0]]
    QR_Factorization(matrix)
    print("----------------------------------------------")
    matrix = [[10.0, 2.0, 4.0], [10.0, 2.0, 4.0], [7.0, 15.0, 9.0]]
    QR_Factorization(matrix)
    print("----------------------------------------------")
    A = [[1.0, 2.0, 3.0], [2.0, 2.0, 8.0], [-3.0, -10.0, -2.0]]
    b = [0.0, -4.0, -11.0]  
    QR_solve(A, b)
    print("----------------------------------------------")
    A = [[2.0, -3.0, 1.0], [1.0, -1.0, -2.0], [3.0, 1.0, -1.0]]
    b = [7.0, -2.0, 0.0]
    QR_solve(A, b)