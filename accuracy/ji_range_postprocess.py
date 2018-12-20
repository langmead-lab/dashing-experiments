#!/usr/bin/env python

from __future__ import print_function
import re

with open('ji_range.latex_snippet', 'rt') as fh:
    for ln in fh:
        ln = re.sub(r'\\multicolumn\{1\}\{.\}\{([^\}]*)\}', r'\1', ln, flags=re.DOTALL)
        toks = ln.split('&')
        toks = list(map(lambda x: x.strip(), toks))
        if len(toks) > 8 and 'Mash' not in toks:
            num_toks = list(map(lambda x: float(x.strip()), toks[3:8]))
            num_sort = sorted(num_toks)
            for i, tok in enumerate(num_toks):
                if tok == num_sort[0] or tok == num_sort[1]:
                    num_toks[i] = r'\textbf{' + "{:,}".format(num_toks[i]) + r'}'
                else:
                    num_toks[i] = "{:,}".format(num_toks[i])
            toks = toks[:3] + list(map(str, num_toks)) + toks[8:]
            toks[8] = toks[8].replace('$', '')
            print(' & '.join(toks))
        else:
            print(ln, end='')
