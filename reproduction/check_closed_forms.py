"""Independent 40-digit closed forms vs exported current Julia solutions.

Run scripts/reproduce_paper.jl first, then:
  python3 reproduction/check_closed_forms.py [--results=results/paper]
Only mpmath is required. No working manuscript folders are imported.
"""
from pathlib import Path
import argparse, csv, json
import mpmath as mp

mp.mp.dps = 40
LAM = mp.j

def power_pair(scale, shift, degree):
    C = mp.mpf(scale)**(mp.mpf(1)/degree)
    s = lambda z: C*(z+mp.mpf(shift))
    if degree == 2:
        q = mp.sqrt(1+4*LAM/C**2)
        return lambda z:s(z)**((-1+q)/2), lambda z:s(z)**((-1-q)/2)
    nu = mp.mpf(1-degree)/(2-degree)
    beta = lambda z:2*mp.sqrt(LAM)*s(z)**(mp.mpf(2-degree)/2)/(C*(2-degree))
    return (lambda z:s(z)**(mp.mpf(1-degree)/2)*mp.besseli(nu,beta(z)),
            lambda z:s(z)**(mp.mpf(1-degree)/2)*mp.besselk(nu,beta(z)))

def legendre(scale, left, right, typ):
    scale,left,right = map(mp.mpc,(scale,left,right))
    xi = lambda z:2*(z+left)/(left-right)-1
    nu = (-1+mp.sqrt(1+4*LAM/scale))/2
    return (lambda z:mp.legenp(nu,0,xi(z),type=typ),
            lambda z:mp.legenq(nu,0,xi(z),type=typ))

def pair(name):
    if name == 'affine':
        beta = lambda z:2*mp.sqrt(LAM*(mp.mpf('.3')*z+mp.mpf('.1')))/mp.mpf('.3')
        return lambda z:mp.besseli(0,beta(z)), lambda z:mp.besselk(0,beta(z))
    if name == 'double_root_quadratic':return power_pair('.25','.8',2)
    if name == 'cubic_power':return power_pair('.4','.6',3)
    if name == 'quartic_power':return power_pair('.3','.7',4)
    if name == 'general_quadratic':return legendre('.6','.5','1.3',3)
    if name == 'roots_straddle':return legendre('-.6','.2','-1.2',2)
    if name == 'roots_above':return legendre('.5','-1.4','-2.2',3)
    if name == 'roots_complex':
        return legendre('.05',mp.mpc('-.5','2'),mp.mpc('-.5','-2'),2)
    if name == 'restricted_cubic':
        k, r, s = map(mp.mpf,('.5','.6','1.4'))
        q = mp.sqrt(1+4*LAM/(k*(s-r)));a,b=(-1-q)/2,(-1+q)/2
        x = lambda z:-(z+r)/(s-r)
        return (lambda z:(z+r)**a*mp.hyp2f1(a,1-b,1+a-b,x(z)),
                lambda z:(z+r)**b*mp.hyp2f1(b,1-a,1+b-a,x(z)))
    if name == 'exponential':
        c=mp.mpf('1.8');K0=mp.mpf('.12')
        beta=lambda z:2*mp.sqrt(LAM*mp.exp(-c*z)/K0)/c
        return (lambda z:mp.exp(-c*z/2)*mp.besseli(1,beta(z)),
                lambda z:mp.exp(-c*z/2)*mp.besselk(1,beta(z)))
    raise ValueError(name)

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--results',type=Path,default=Path(__file__).resolve().parents[1]/'results/paper')
    args=parser.parse_args();groups={}
    with (args.results/'catalogue_velocities.csv').open() as stream:
        for row in csv.DictReader(stream):groups.setdefault((row['case'],int(row['layers'])),[]).append(row)
    checks=[]
    for (name,layers),rows in groups.items():
        f,g=pair(name)
        a,b=mp.lu_solve(mp.matrix([[f(0),g(0)],[f(1),g(1)]]),mp.matrix([-mp.j,0]))
        error=max(abs(mp.mpc(row['real'],row['imag'])-a*f(mp.mpf(row['z']))-b*g(mp.mpf(row['z']))) for row in rows)
        if error > mp.mpf('1e-11'):raise AssertionError((name,layers,str(error)))
        checks.append(dict(case=name,layers=layers,max_error=float(error)))
        print(f'{name:24s} N={layers}  {float(error):.3e}')
    expected={(n,N) for n in ('affine','double_root_quadratic','general_quadratic','cubic_power','restricted_cubic','quartic_power','exponential') for N in (1,4)}
    if not expected.issubset(groups):raise AssertionError('missing one of the fourteen catalogue runs')
    (args.results/'closed_form_checks.json').write_text(json.dumps(checks,indent=2)+'\n')

if __name__=='__main__':main()
