import math,json,time
p=13;Q=169;m=4*Q
wi=[math.comb(i+2*Q,2*Q)**3 for i in range(2*m)]
wseed=[math.comb(i+2,2)**3%p for i in range(8)]
def redfrac(n,d):
 while d%p==0:
  assert n%p==0;n//=p;d//=p
 return n*pow(d,-1,p)%p
records=[]
for a in [67,-70]:
 b=[0]
 for s in range(1,2*m):b.append((a*b[-1]+redfrac(Q**3,s**3))%p)
 bs=[0]
 for s in range(1,8):bs.append((a*bs[-1]+pow(s,-3,p))%p)
 count=0
 for i in range(m):
  for j in range(m,2*m):
   c=redfrac(Q*wi[i],j-i);d=-redfrac(Q*wi[j],j-i)%p
   ep=(-2*(c*b[i]+d*b[j]))%p
   kp=2*(c*pow(a,i,p)+d*pow(a,j,p))%p
   if i%Q==j%Q:
    u,v=i//Q,j//Q;r=i%Q
    cs=wseed[u]*pow(v-u,-1,p)%p;ds=-wseed[v]*pow(v-u,-1,p)%p
    exp_e=pow(a,r,p)*(-2*(cs*bs[u]+ds*bs[v]))%p
    exp_k=pow(a,r,p)*2*(cs*pow(a,u,p)+ds*pow(a,v,p))%p
   else:exp_e=exp_k=0
   assert (ep,kp)==(exp_e,exp_k),(a,i,j,ep,kp,exp_e,exp_k)
   count+=1
 records.append(dict(a=a,m=m,entries=count,checks=['Q4 rational-pole part of E','Q coefficient K3'],all_pass=True,polynomial_part='uniform proof; NOT an independent full-entry numerical test'))
print(json.dumps(records, indent=2))
