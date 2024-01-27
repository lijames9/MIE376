import numpy as np
    
#===============================================================
# Q1 data
c1 = np.array([-3, 8, 0, 0])
b1 = np.array([12, 6])
A1 = np.array([[4, 2, 1, 0], [2, 3, 0, 1]])
slack_vars1 = 2
minimize1 = True
#===============================================================
#===============================================================
# Q2 data
c2 = np.array([3, 2, -1, -2, 1, 2, -1, 3, 4, -3, 0, 0, 0, 0, 0])
b2 = np.array([80, 50, 40, 90, 50])
A2 = np.array([[2, 1, 3, 1, 2, 1, 4, 1, -2, 3, 1, 0, 0, 0, 0],
              [1, -4, 1, 2, 3, 1, -1, 4, 1, 2, 0, 1, 0, 0, 0],
              [3, 2, -2, -1, 1, 3, 2, 1, 1, 1, 0, 0, 1, 0, 0],
              [2, 3, 1, 1, 1, 1, 1, 1, 3, 2, 0, 0, 0, 1, 0],
              [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 1]])
slack_vars2 = 5
minimize2 = False
#===============================================================

def revised_sim(A, b, c, minimize, slack_vars):
    # Iteration 1
    c_b = c[(-1*slack_vars):]
    c_n = c[:(-1*slack_vars)]
    B = A[:, (-1*slack_vars):]
    N = A[:, :(-1*slack_vars)]

    print("Basic coefficient vector:", c_b)
    print("Non-basic coefficient vector:", c_n)
    print("B Matrix: \n", B)
    print("N Matrix: \n", N)
    res_mat = c_b.dot(np.linalg.inv(B).dot(N)) - c_n
    iters = 1
    if minimize == True: 
        #Iteration 1
        if all(i <= 0 for i in res_mat) == True:
            x_b = np.linalg.inv(B).dot(b)
            print("The optimal solution is:", c_b.dot(x_b), ", obtained on iteration 1.")
            print(x_b)
        while all(i <= 0 for i in res_mat) == False: 
            # determining var entering basic set B
            res_mat = c_b.dot(np.linalg.inv(B).dot(N)) - c_n
            max_index = max(res_mat)
            res = [i for i, j in enumerate(res_mat) if j == max_index]

            print(res_mat)

            # determining var leaving non-basic set N
            a_entering = [row[res[0]] for row in A]

            # B(-1)*b/B(-1)*a 
            res_mat2 = np.divide(np.linalg.inv(B).dot(b), (np.linalg.inv(B).dot(a_entering)))
            print((np.linalg.inv(B).dot(a_entering)))

            max_index2 = max(res_mat2)
            res2 = [i for i, j in enumerate(res_mat2) if j == max_index2]

            enter_col = np.copy(N[:,res[0]])
            exit_col = np.copy(B[:,res2[0]])
            enter_cB = np.copy(c_n[res[0]])
            enter_cN = np.copy(c_b[res2[0]])    

            #if cost vectors are fine after iteration 1, compute z = cB*xB
            if all(i <= 0 for i in res_mat) == True:
               continue
            else:
                B[:,res2[0]] = enter_col
                N[:,res[0]] = exit_col
                c_b[res2[0]] = enter_cB
                c_n[res[0]] = enter_cN
                iters = iters+1
        # all other iterations
        if all(i <= 0 for i in res_mat) == True:
            x_b = np.linalg.inv(B).dot(b)
            print("The optimal solution is:", c_b.dot(x_b), ", obtained on iteration", iters)
            print(x_b)
    elif minimize == False:
        # Iteration 1
        if all(i >= 0 for i in res_mat) == True:
            x_b = np.linalg.inv(B).dot(b)
            print("The optimal solution is:", c_b.dot(x_b), ", obtained on iteration 1.")
        while all(i>= 0 for i in res_mat) == False:
            # determining var entering basic set B
            res_mat = c_b.dot(np.linalg.inv(B)).dot(N) - c_n

            min_index = min(res_mat)
            res = [i for i, j in enumerate(res_mat) if j == min_index]

            print("res_mat =", res_mat)
            print("The basic variable entering non-basic is index:", res[0] + 1)
            # determining var leaving non-basic set N

            a_entering = [row[res[0]] for row in A]

            # B(-1)*b/B(-1)*a 
            res_mat2 = np.divide(np.linalg.inv(B).dot(b), (np.linalg.inv(B).dot(a_entering)))

            min_index2 = min([n for n in res_mat2 if n>=0])
            
            res2 = [i for i, j in enumerate(res_mat2) if j == min_index2]
            print(res2)
            print("res_mat2 =", res_mat2)
            print("The non-basic variable entering is index:", res2[0] + 1)
            #print(c_n[res[0]])

            enter_col = np.copy(N[:,res[0]])
            exit_col = np.copy(B[:,res2[0]])
            enter_cB = np.copy(c_n[res[0]])
            enter_cN = np.copy(c_b[res2[0]])
            #if cost vectors are fine after iteration 1, compute z = cB*xB
            if all(i >= 0 for i in res_mat) == True:
                continue
            else:
                B[:,res2[0]] = enter_col
                N[:,res[0]] = exit_col
                c_b[res2[0]] = enter_cB
                c_n[res[0]] = enter_cN
                iters = iters + 1
        if all(i >= 0 for i in res_mat) == True:
            x_b = np.linalg.inv(B).dot(b)
            print("The optimal solution is:", c_b.dot(x_b), ", obtained on iteration", iters)

revised_sim(A1,b1,c1,minimize1,slack_vars1)
revised_sim(A2,b2,c2,minimize2,slack_vars2)