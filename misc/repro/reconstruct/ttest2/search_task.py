# Reconstructed driver, not the authors' script. It is meant to recover the command lines of the archived
# code/ttest2.py behind Table tab:tilt of Chapter ch:polya (ch/ch04.tex), which the book does not record.
# Each SLURM array task runs the archived ttest2.py unchanged (copied beside this file with code/tilted.py)
# for one candidate step h. The table fixes delta, the bits, the node count 2J+1 and the heights; J =
# ceil(xmax/h) in code/tilted.py, so for a given J the nodes x_j = j*h, |j| <= J, depend on h alone, and
# any xmax in ((J-1)h, Jh] gives the same run. The grid uses xmax = (J - 1/2)h, which is safe from the
# rounding of xmax/h, plus one check that another xmax with the same J gives the same output.
# Usage: python search_task.py TASK_ID   (writes out/search_TASK_ID.json)
#        python search_task.py --list    (prints the grid)
import sys, os, json, subprocess, time
from decimal import Decimal as D

CONF = {'A': dict(delta='0.10', prec='200', J=1000, heights=['200', '500', '1000']),
        'B': dict(delta='0.02', prec='220', J=4200, heights=['1000', '2000', '5000'])}

def grid():
    G = []
    hA = [D('0.0020') + D('0.0001') * k for k in range(41)]
    hA += [D(x) for x in ['0.00395', '0.00399', '0.003999', '0.004001', '0.00401', '0.00405']]
    for h in hA:
        G.append(('A', h, (CONF['A']['J'] - D('0.5')) * h))
    G.append(('A', D('0.004'), D('3.9965')))   # check: same J = 1000 as xmax = 3.998 and 4
    hB = [D('0.00070') + D('0.00001') * k for k in range(61)]
    hB += [D(x) for x in ['0.000995', '0.000999', '0.0009999', '0.0010001', '0.001001', '0.001005']]
    for h in hB:
        G.append(('B', h, (CONF['B']['J'] - D('0.5')) * h))
    G.append(('B', D('0.001'), D('4.1991')))   # check: same J = 4200 as xmax = 4.1995 and 4.2
    return G

if __name__ == '__main__':
    G = grid()
    if sys.argv[1] == '--list':
        for i, (c, h, x) in enumerate(G):
            print(i, c, h.normalize(), x.normalize())
        sys.exit(0)
    i = int(sys.argv[1])
    c, h, x = G[i]
    cf = CONF[c]
    args = [cf['delta'], cf['prec'], str(h.normalize()), str(x.normalize())] + cf['heights']
    t0 = time.time()
    p = subprocess.run([sys.executable, '-u', 'ttest2.py'] + args, capture_output=True, text=True)
    rec = dict(task=i, config=c, delta=cf['delta'], prec=cf['prec'], h=str(h.normalize()), xmax=str(x.normalize()),
               J=cf['J'], cmd='python ttest2.py ' + ' '.join(args), rc=p.returncode, stdout=p.stdout,
               stderr=p.stderr[-2000:], seconds=round(time.time() - t0, 2), host=os.uname().nodename)
    os.makedirs('out', exist_ok=True)
    json.dump(rec, open(f'out/search_{i}.json', 'w'), indent=1)
    print(json.dumps(rec))
