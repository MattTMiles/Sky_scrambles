import numpy as np
import matplotlib.pyplot as plt
import csv 

with open('n_tries_100.csv', 'r') as file:
    csvreader = csv.reader(file)
    n_t = []    
    for row in csvreader:
        n_t.append(row)

n_tries = [int(item[0]) for item in n_t]
print(n_tries)


def count_zeros(arr):
    result = []
    count = 0
    for i in arr:
        if i == 0:
            count += 1
        else:
            result.append(count+1)
            count = 0
    return result

def cumulative(lists): 
    cu_list = [] 
    length = len(lists) 
    cu_list = [sum(lists[0:x:1]) for x in range(0, length+1)] 
    return cu_list[1:]


num_of_tries = count_zeros(n_tries)
#efficiency = np.reciprocal([float(i) for i in num_of_tries])
#cum_eff = cumulative(efficiency)
cum_num_of_tries = cumulative(num_of_tries)
cum_accepted = np.arange(len(num_of_tries))


plt.scatter(cum_accepted, cum_num_of_tries)
#plt.xscale('log') 
#plt.yscale('log')
plt.ylabel('number of tries before acceptance')
plt.xlabel('number of accepted scrambles');



plt.savefig('scrambles_val_100.png', bbox_inches='tight' )

