'''Independent exact oracle for the finite ellipse study.
Run from any directory with the default repository data path, or pass --csv.
Uses no nrcert imports and writes no files. JSON is emitted to stdout.
Checks a fixed 144-case grid; this is finite validation, not a general proof.
'''
import argparse,csv,hashlib,json
from fractions import Fraction as F
from pathlib import Path
def root(x,b):
    scale=1<<b; lo=0; hi=(x.numerator//x.denominator+1)*scale
    while lo<hi:
        mid=(lo+hi)//2
        if mid*mid*x.denominator>=x.numerator*scale*scale: hi=mid
        else: lo=mid+1
    return F(lo,scale)
def clip(poly,a,b,h):
    out=[]
    for i,current in enumerate(poly):
        previous=poly[i-1]; fp=a*previous[0]+b*previous[1]-h; fc=a*current[0]+b*current[1]-h
        if (fp<=0)!=(fc<=0):
            t=fp/(fp-fc); out.append(tuple(previous[j]+t*(current[j]-previous[j]) for j in (0,1)))
        if fc<=0: out.append(current)
    clean=[]
    for p in out:
        if not clean or p!=clean[-1]: clean.append(p)
    if len(clean)>1 and clean[0]==clean[-1]: clean.pop()
    changed=True
    while changed and len(clean)>2:
        changed=False
        for i,v in enumerate(clean):
            u,w=clean[i-1],clean[(i+1)%len(clean)]
            if (v[0]-u[0])*(w[1]-v[1])==(v[1]-u[1])*(w[0]-v[0]):
                clean.pop(i);changed=True;break
    return clean
parser=argparse.ArgumentParser(description='Independently recompute the exact ellipse-study CSV with standard-library arithmetic only.')
parser.add_argument('--csv',type=Path,default=(Path(__file__).resolve().parents[1]/'results'/'ellipse-error-budget.csv') if '__file__' in globals() else Path('results/ellipse-error-budget.csv'))
args=parser.parse_args()
raw=args.csv.read_bytes()
rows=list(csv.DictReader(raw.decode('utf-8').splitlines()))
expected_cases={(p,z,b,k,m) for p in ('analytic_ellipse','gershgorin') for z in ('z_squared','one_minus_z_squared') for b in (8,24,40) for k in (1,2,4,8) for m in (4,16,64)}
actual_cases=[(r['producer'],r['polynomial'],int(r['sqrt_bits']),int(r['resolution']),int(r['segments_per_edge'])) for r in rows]
if len(actual_cases)!=144 or len(set(actual_cases))!=144 or set(actual_cases)!=expected_cases:
    raise ArithmeticError('CSV does not contain exactly the prescribed 144 distinct parameter cases')
groups={}
for row in rows:
    bits,k,producer=int(row['sqrt_bits']),int(row['resolution']),row['producer']
    key=bits,k,producer
    if key not in groups:
        normals=sorted((a,b) for a in range(-k,k+1) for b in range(-k,k+1) if max(abs(a),abs(b))==k)
        supports=[];step=F(1,1<<bits)
        for a,b in normals:
            exact_root_upper=root(F(5*a*a+b*b,4),bits)
            g=exact_root_upper if producer=='analytic_ellipse' else abs(a)+root(F(a*a+b*b,4),bits)
            supports.append((a,b,g,exact_root_upper))
        hd={(a,b):g for a,b,g,_ in supports}
        xmin=-hd[-k,0]/k;xmax=hd[k,0]/k;ymin=-hd[0,-k]/k;ymax=hd[0,k]/k
        poly=[(xmin,ymin),(xmax,ymin),(xmax,ymax),(xmin,ymax)]
        for a,b,h,_ in supports: poly=clip(poly,a,b,h)
        eta=max((g-max(F(0),v-step))/k for a,b,g,v in supports)
        R=root(F(5,4),bits);rho=(R+eta)/(1-F(1,8*k*k));excess=eta+(R+rho)/(2*k)
        e2=max((v[0]-u[0])**2+(v[1]-u[1])**2 for u,v in zip(poly,poly[1:]+poly[:1]))
        if not all(x*x+y*y<=rho*rho for x,y in poly): raise ArithmeticError('independent vertex radius test failed')
        groups[key]=poly,eta,R,rho,excess,e2
    poly,eta,R,rho,excess,e2=groups[key]
    m=int(row['segments_per_edge']); max2=F(0);samples=0
    for u,v in zip(poly,poly[1:]+poly[:1]):
        for j in range(m+1):
            x=u[0]+(v[0]-u[0])*F(j,m);y=u[1]+(v[1]-u[1])*F(j,m)
            value2=(x*x+y*y)**2 if row['polynomial']=='z_squared' else (1-x*x+y*y)**2+4*x*x*y*y
            max2=max(max2,value2);samples+=1
    delta=root(e2,bits)/(2*m);D=2*rho;sample=root(max2,bits);U=sample+delta*D;budget=D*(excess+delta)+F(1,1<<bits)
    expected={'support_slack_upper':eta,'radius_upper':rho,'hausdorff_excess_upper':excess,'derivative_upper':D,'mesh_radius_upper':delta,'sample_upper':sample,'mesh_term':delta*D,'U_Q':U,'M_Q':F(5,4),'observed_excess':U-F(5,4),'excess_budget':budget}
    if any(F(row[name])!=value for name,value in expected.items()): raise ArithmeticError(('row mismatch',key,row['polynomial'],m))
    if int(row['vertices'])!=len(poly) or int(row['sample_count'])!=samples: raise ArithmeticError('count mismatch')
    if not F(0)<=U-F(5,4)<=budget: raise ArithmeticError('budget mismatch')
print(json.dumps({'rows_recomputed':len(rows),'geometric_parameter_groups':len(groups),'all_exact_fields_match':True,'complete_parameter_grid':True,'independent_components':['binary-search dyadic root','literal 2x2 Hermitian support formula','rational Sutherland-Hodgman half-plane clipping','closed-form polynomial modulus evaluation'],'shared_nrcert_helpers_used':False,'checks_survive_python_optimization':True,'scope':'independent arithmetic oracle for this ellipse study only; not general software proof','input_csv':str(args.csv),'input_csv_sha256':hashlib.sha256(raw).hexdigest()}))
