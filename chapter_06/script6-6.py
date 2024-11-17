import matplotlib.pyplot as plt

def gegenbauer(  n , x , alpha = -1/2 ):
    if n == 0 :
        return 1
    if n == 1 :
        return 2 *  alpha * x 
    if n >= 2 :
        return 1.0 / n * ( 2 * x * ( n + alpha - 1 ) * gegenbauer( n - 1 , x , alpha ) - ( n + 2 * alpha - 2 ) * gegenbauer( n - 2 , x , alpha ) )

y = []
for i in range(2,19):
    y.append( 4 * gegenbauer( i , 1/2 ) )

z = []
for i in range(0,len(y)):
    sum = 0
    for j in range(0,i+1):
        sum += y[j]
    z.append( sum )

x = range(2,len(z)+2)

plt.figure( figsize=(8,4) , dpi=960 )
plt.plot(x, z, marker='o')

plt.xlim(1, 19)
plt.locator_params(axis='x', nbins=9)
plt.ylim(1.45, 2.45)
plt.locator_params(axis='y', nbins=13)

ax = plt.subplot()
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['bottom'].set_position(('data',2))
ax.tick_params(direction='inout')
##plt.axhline(y=2.0, c="r", ls="--")


plt.xlabel(r'$M$')
plt.ylabel(r"$\mathrm{\sum_{n=2}^M E^{(n)}_0/\beta}$")



plt.savefig("./gegenbauer.png", format="png", transparent=True, bbox_inches='tight', pad_inches=0.0)