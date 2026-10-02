#!/usr/bin/env python3
# NOT an author script. Reconstruction helper written for the reproducibility audit.
# It is meant to reproduce verification_reciprocal_zeta/spot_checks_report.md.
#
# The archived spot_checks.py never opens a file: it prints its checks to standard
# output, and the archived spot_checks_report.md is that output saved verbatim (no
# header, no markdown). run_all.py therefore sent the output to logs/spot_checks.out,
# and verif/spot_checks_report.md was never written. This helper runs the archived
# script unchanged, with the same interpreter, and saves its standard output byte for
# byte as spot_checks_report.md, which is what
#     python3 -u spot_checks.py > spot_checks_report.md
# does. It also records when each line appeared, so that the time of every section is
# known, and the versions and the node, in spot_checks_run.json.
import hashlib
import json
import os
import platform
import socket
import subprocess
import sys
import time

import mpmath
import sympy

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPT = os.path.join(HERE, 'spot_checks.py')
REPORT = os.path.join(HERE, 'spot_checks_report.md')
ERR = os.path.join(HERE, 'spot_checks_stderr.txt')
INFO = os.path.join(HERE, 'spot_checks_run.json')


def cpu_model():
    try:
        with open('/proc/cpuinfo') as fh:
            for line in fh:
                if line.startswith('model name'):
                    return line.split(':', 1)[1].strip()
    except OSError:
        pass
    return platform.processor()


t0 = time.time()
lines, stamps = [], []
with open(ERR, 'wb') as err:
    p = subprocess.Popen([sys.executable, '-u', SCRIPT], cwd=HERE,
                         stdout=subprocess.PIPE, stderr=err)
    for raw in p.stdout:
        lines.append(raw)
        stamps.append(round(time.time() - t0, 2))
    rc = p.wait()
total = round(time.time() - t0, 1)
with open(REPORT, 'wb') as fh:
    fh.write(b''.join(lines))

# section times: from one section header to the next (the first starts at launch)
sections, start, name = [], 0.0, '(imports)'
for raw, t in zip(lines, stamps):
    s = raw.decode().rstrip('\n')
    if s.startswith('--- ') and s.endswith(' ---'):
        sections.append({'section': name, 'seconds': round(t - start, 1)})
        start, name = t, s.strip('- ')
sections.append({'section': name, 'seconds': round(total - start, 1)})

info = {
    'rc': rc,
    'total_seconds': total,
    'job': os.environ.get('SLURM_JOB_ID'),
    'host': socket.gethostname(),
    'cpu': cpu_model(),
    'python': sys.version.split()[0],
    'mpmath': mpmath.__version__,
    'mpmath_backend': mpmath.libmp.BACKEND,
    'sympy': sympy.__version__,
    'spot_checks_py_sha256': hashlib.sha256(open(SCRIPT, 'rb').read()).hexdigest(),
    'report_sha256': hashlib.sha256(b''.join(lines)).hexdigest(),
    'report_lines': len(lines),
    'report_bytes': sum(len(x) for x in lines),
    'stderr_bytes': os.path.getsize(ERR),
    'sections': sections,
}
with open(INFO, 'w') as fh:
    json.dump(info, fh, indent=1)
sys.stdout.write(b''.join(lines).decode())
print('\n=== run info ===')
print(json.dumps(info, indent=1))
sys.exit(rc)
