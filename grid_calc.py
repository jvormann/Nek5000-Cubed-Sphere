def isprime(n):
#efficient way to test if n is a prime number
    n*=1.0
    if n%2==0 and n!=2 or n%3==0 and n!=3:
        return False
    for b in range(1,int((n**0.5+1)/6.0+1)):
        if n%(6*b-1)==0:
            return False
        if n %(6*b+1)==0:
           return False
    return True
#------------------------------------------------------
def aab_fact(n_el):
#Find a factorization n_el/6 = a*a*b with a,b integer and a approx b
  a_approx=int((n_el/6)**(1./3.))
  b_approx=a_approx
  a_list = range(1,int(2.*a_approx))
  b_list = a_list
  a_sol = 0
  b_sol = int(n_el/6)
  abdiff = abs(a_sol-b_sol)
  for b in b_list:
    for a in a_list:
      if ((a*a*b == n_el/6) and (abs(a-b) <= abdiff)):
	abdiff=abs(a-b)
	a_sol=a
	b_sol=b
  return [a_sol, b_sol] 
#------------------------------------------------------
def main():
  n_el=input("Minimal Number of Elements? ")
  n_p=input("Number of MPI ranks? ")
#  n_p=int(n_el/50.)
#  print "Max. # of ranks= ",n_p
  n_el_max=2*n_el
#Find a recommended number of ranks with at least 
#50 elements per rank
#  while (n_p % 96 != 0):
#    n_p -= 1
#  n_el=50*n_p
#  print "Recommended # of ranks= ", n_p

#Find a recommended number of elements where each rank has the
#same number of elements
  print '#El | #El/6 | #El_hor/4 | #El_rad'  
  while (n_el<n_el_max):
    decomp=[0, 0]
    if ((n_el % 6 == 0) and (isprime(n_el/6) == False) and (n_el % n_p == 0)):
      decomp=aab_fact(n_el)
      if (decomp[0] != 0):
	print n_el, n_el/6, decomp[0], decomp[1]
    n_el += 1

  

#------------------------------------------------------
main()
