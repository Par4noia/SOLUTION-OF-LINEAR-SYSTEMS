import copy
import numpy as np

def fill_matrix1(a, b, file_name): #Заполнение матрицы из файла
    file = open(file_name, 'r')
    f = file.readlines()
    n = int(f[0])
    for i in range(1, n+1):
        a.append(f[i])
        a[i - 1] = a[i - 1].split(" ")
        for j in range(n - 1):
            a[i - 1][j] = int(a[i - 1][j])
        a[i - 1][n - 1] = int(a[i - 1][n - 1][:len(a[i - 1][n - 1]) - 1])
    f[n + 1] = f[n + 1].split(" ")
    for i in range(n):
        b.append(int(f[n + 1][i]))
    return n

def fill_matrix2(a, b):
    print("Enter matrix size")
    n = int(input())
    for i in range(n):
        print("Enter the", i + 1, "row of the matrix")
        a.append([])
        a[i] = list(map(int,input().split()))
    print("Enter the vector of solutions")
    s = input().split(" ")
    for i in range(n):
        b.append(int(s[i]))
    return n

def fill_formula(a, b):
    n = 30
    m = 9
    for i in range(n):
        a.append([])
        for j in range(n):
            if i != j:
                a[i].append((i + j + 2)/(m + n))
            else:
                a[i].append(n + m**2 + (j + 1)/m + (i + 1)/n)
    for i in range(1, n+1):
        b.append(i**2 - 100)

def swap_lines(a, ind1, ind2):#Функция меняющая строки в матрице местами
    a[ind1], a[ind2] = a[ind2], a[ind1]
    
def gauss_solution(a, b, c, n): #Решение СЛАУ методом Гаусса
    for i in range(0, n - 1): #Приводим к верхней треугольной форме
        for j in range(i + 1, n):
            if a[i][i] == 0:
                line = i
                for k in range(i + 1, n):
                    if (a[k][i] != 0):
                        line = k
                        break
                swap_lines(a, i, line)
                swap_lines(b, i, line)
            fact = a[j][i] / a[i][i]
            for k in range(i, n):
                a[j][k] -= fact * a[i][k]
            b[j] -= b[i] * fact
    #Ищем решения    
    if a[n - 1][n - 1] == 0:
        c[n - 1] = 0
    else:
        c[n - 1] = b[n - 1]/a[n-1][n-1]

    for i in range(n - 2, -1, -1):
        for j in range(n - 1, i, -1):
            b[i] -= a[i][j] * c[j]
            if b[i] == 0:
                c[i] = 0
            else:
                c[i] = b[i] / a[i][i]

def modified_gauss(a, b, c, n): #Решение СЛАУ методом Гаусса с выделеением главного элемента
    res = [0] * n
    res_num = list()
    for i in range(n):
        res_num.append(i)
    for i in range(n - 1):
        max_elem = abs(a[i][i])
        ind = i
        for k in range(i + 1, n):
            if (abs(a[i][k]) > max_elem):
                ind = k
                max_elem = abs(a[i][k])
        if i != ind:
            for k in range(0, n):
                swap_lines(a, i, ind)
                swap_lines(b, i, ind)
        res_num[i], res_num[ind] = res_num[ind], res_num[i]
        fact = 0
        for j in range(i + 1, n):
            fact = a[j][i]/a[i][i]
            for k in range(i, n):
                a[j][k] -= fact * a[i][k]
            b[j] -= fact * b[i]
    if a[n - 1][n -1] == 0:
        res[n - 1] = 0
    else:
        res[n - 1] = b[n - 1]/a[n - 1][n -1]
    for i in range(n - 2, -1, -1):
        fact = b[i]
        for j in range(n - 1, i, -1):
            fact -= a[i][j] * res[j]

        if fact == 0:
            res[i] = 0
        else:
            res[i] = fact / a[i][i]
    for i in range(n):
        if res_num[i] != i:
            res[i], res[res_num[i]] = res[res_num[i]], res[i]
            res_num[i], res_num[res_num[i]] = res_num[res_num[i]], res_num[i]
    for i in range(n):
        c[i] = res[i]
def det(a, n): # Функция считает определитель матрицы
    a2 = copy.deepcopy(a)
    sign = 1
    for i in range(0, n - 1): #Приводим к верхней треугольной форме
        for j in range(i + 1, n):
            if a2[i][i] == 0:
                line = i
                for k in range(i + 1, n):
                    if (a2[k][i] != 0):
                        line = k
                        break
                swap_lines(a2, i, line)
                sign *= -1
            fact = a2[j][i] / a2[i][i]
            for k in range(i, n):
                a2[j][k] -= fact * a2[i][k]
    ans = 1
    for i in range(n):
        ans *= a2[i][i]
    ans *= sign
    return ans

def inverse_matrix(a, n): # Поиск обратной матрицы
    a2 = copy.deepcopy(a)
    for i in range(n):
        for j in range(n):
            if j == i:
                a[i][j] = 1
            else:
                a[i][j] = 0
    for j in range(n):
        t = j
        while t < n and a2[t][j] == 0:
            t += 1
        if t != j:
            a2[j], a2[t] = a2[t], a2[j]
            a[j], a[t] = a[t], a[j]
        fact = a2[j][j]
        for k in range(n):
            a2[j][k] = a2[j][k] / fact
            a[j][k] = a[j][k] / fact
        for i in range(t + 1, n):
            if (a2[i][j] != 0):
                fact = a2[i][j]
                for k in range(n):
                    a2[i][k] -= a2[j][k] * fact
                    a[i][k] -= a[j][k] * fact
    for j in range(n - 1, 0, -1):
        for i in range(j - 1, -1, -1):
            for k in range(n):
                a[i][k] -= a2[i][j] * a[j][k]

def matrix_norm(a, n): # Первая матричная норма
    ans = 0
    for i in range(n):
        summ = 0
        for j in range(n):
            summ += abs(a[i][j])
        if summ > ans:
            ans = summ
    return ans

def conditionality_num(a, n):
    ss = matrix_norm(a, n)
    inverse_matrix(a, n)
    return ss * matrix_norm(a, n)

def relaxation(a, b, n, max_iters, omega, eps):
    c = np.zeros(n)
    counter = 0
    while counter < max_iters:
        c_new = np.copy(c)
        for i in range(n):
            s1 = 0
            s2 = 0
            for j in range(i):
                s1 += a[i][j] * c_new[j]
            for j in range(i, n):
                s2 += a[i][j] * c[j]
            c_new[i] = c[i] + ((b[i] - s1 - s2) * omega) / a[i][i]
        s = 0
        for i in range(n):
            s += (c_new[i] - c[i]) ** 2
        if np.sqrt(s) <= eps:
            break
        c = c_new
        counter += 1
    return counter
def gg(a):
    for i in range(len(a)):
        for j in range(len(a)):
            print("{0:.3f}".format(a[i][j]), end = " ")
        print()
def hh(a):
    for i in range(len(a)):
        print("{0:.3f}".format(a[i]), end = " ")
   
print("Enter 1 to fill the matrix from the keyboard")
print("or 2 to fill the matrix from the file")
print("or another number to fill by formula")
mode = int(input())
a = list()
b = list()
n = 30
c = [0]*n
if mode == 1:
    fill_matrix2(a, b)
elif mode == 2:
    print("Enter file name")
    #f = input()
    f = "input.txt"
    fill_matrix1(a, b, f)
else:
    fill_formula(a, b)

bb = 20000
ans = 0.1
for i in range(1, 20):
    q = relaxation(a, b, n, 20000, i * 0.1, 0.0001)
    print("w =", "{0:.1f}".format(0.1 * i), "iterations =", q, sep = " ")
    if (bb > q):
        bb = q
        ans = i
print("Best result with w =", "{0:.1f}".format(0.1 * ans))
print("Number of iterations =", bb)
