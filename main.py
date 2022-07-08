import sympy as sp
from sympy.utilities.lambdify import lambdify

class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[90m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'
    # Background colors:
    GREYBG = '\033[100m'
    REDBG = '\033[101m'
    GREENBG = '\033[102m'
    YELLOWBG = '\033[103m'
    BLUEBG = '\033[104m'
    PINKBG = '\033[105m'
    CYANBG = '\033[106m'

x = sp.symbols('x') # Symbolizes x as a parameter in which a value can be placed

def polynomial(p):
    """
    The function get a list of Coefficients and create Polynom
    :param p: list of Coefficients
    :return:polynom
    """
    n=len(p)-1
    poly=0
    for i in range(0,len(p)):
        if(n>=2):
            poly+=p[i]*x**n
        if(n==1):
            poly += p[i]*x
        if (n == 0):
            poly += p[i]
        n-=1
    return poly

def find_root(p,q, guess ):
    """
    The function find approximate value to root by using Newton Rapson method
    :param p: polynom
    :param q: Q(X) -polynom
    :param guess: the guess we want to approximate
    :return: approximate value to root
    """
    p_calc = lambdify(x, p)  # Conversion to a function that can be calculated when placing x
    q_calc = lambdify(x, q)
    x_n = guess - (p_calc(guess) / q_calc(guess)) # Newton Rapson approximate

    return x_n


def div_poly(p, guess):
    """
    synthetic division of p in x-guess binomial
    :param p: the origin polynom
    :param guess:the root of p that was found
    :return:new polynom (degree n-1)
    """
    x = sp.symbols('x')
    poly = sp.poly(p).all_coeffs() # get the Coefficients in a list
    result = []  # to contain all the new Coefficients
    result.append(poly[0])
    for i in range(1,len(poly)-1):
        x=poly[i]+result[i-1]*guess
        result.append(x)

    return polynomial(result)



def Horners_Method(p,n,guess):
    """
    find roots of polynom
    the function print to the screen and return the roots of the polynon (If there are any at all)
    :param p: the polynom
    :param n: the degree of the polynom
    :param guess: initial guess
    :return: roots
    """
    result = []  # contain all the roots that was found
    e = 0.000001  # The required level of accuracy
    p_calc = lambdify(x, p)
    print(bcolors.OKGREEN, 'p(x):', p, bcolors.ENDC)
    print('Guess: ', guess)
    for i in range(0, n):  # maxinum possible number of roots
        root = guess
        flag=True  # mark when there is no roots anymore
        j=1  # count the number of iterations
        print('Try Finding a root by approximation:')
        while abs(p_calc(root)) > e:
            Q_x = div_poly(p, root)
            root=find_root(p,Q_x,root)
            print(bcolors.FAIL, 'approx number',j,'= ', root, bcolors.ENDC)
            if j>100:  # if out of bound or above 100 iterations stop search roots
                flag=False
                print('No more roots ')
                break
            j += 1

        p=Q_x
        p_calc = lambdify(x, p)

        if flag:
            result.append(round(root, 5))  # round the result according to the selected epsilon
            print(bcolors.OKBLUE, "root=", root, bcolors.ENDC)
        else:
            break
        prt = 'p(x)='
        for i in range(0, len(result)):
            if -result[i] > 0:
                prt += '(x+'
            else:
                prt += '(x'
            prt += str(-result[i]) + ')'
        if i!=n-1:
            prt += '*(' + str(p) + ')'
        print(prt)
    if len(result)==0:
        print("No roots")
    else:
        print(bcolors.GREYBG, 'roots found:', result, bcolors.ENDC)
        return result


def main():
    x = sp.symbols('x')
    f = 5*x**5+4*x**3-2*x # change here if you want another function
    rank=5
    print("Please enter an initial guess: ")
    x0=float(input())
    f_calc = lambdify(x, f)
    if f_calc(x0)==0: # if x0 is a root
        x0 += 0.000000001
    Horners_Method(f, rank, x0)



main()

