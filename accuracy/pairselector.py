#!/usr/bin/env python
import sys
import numpy as np
import math


def main():
    import argparse
    p, a = argparse.ArgumentParser(), p.add_argument
    p.add_argument("--npairs", type=int, default=20,
                   help=("Number of pairs within bucket to store (to "
                         "account for potential inaccuracies"
                         " in our estimation."))
    p.add_argument("--linspace", "-l",
                   type=int, default=200,
                   help="Number of buckets between 0 and 1 to cover.")
    p.add_argument("-o", "--outpath", default="/dev/stdout")
    p.add_argument("--use-mash-distance", "-m", action='store_true')
    p.add_argument("-k", "--kmer-size", type=int, default=0)
    p.add_argument("distfile", help="Path to distance file")
    a = p.parse_args()
    p = a.npairs
    labels = None
    lsp = a.linspace
    idxpairs = [[] for i in range(lsp)]
    print("Len of idx pairs at start: %i" % len(idxpairs), file=sys.stderr)
    assert len(idxpairs) == lsp
    linenum = 0
    if a.use_mash_distance:
        assert a.kmer_size != 0
        ksinv = 1. / a.kmer_size

        def get_val(x):
            x = float(x)
            return (-math.log(2. * x / (1. + x)) * ksinv if x > 0. else 1.) if x != 1. else 0.
    else:
        def get_val(x):
            return float(x)

    def unfilled_bucket_count():
        return sum(map(lambda x: len(x) < p, idxpairs))
    names = next(open(a.distfile)).strip().split('\t')[1:]
    ln = 0
    for line in open(a.distfile):
        ln += 1
        if line[0] == '#':
            if not labels:
                labels = line[1:].strip().split()
            continue
        gen = (i for i in line.strip().split('\t')[1:] if i != '-')
        for ind, dist in enumerate(map(get_val, gen), 1):
            bucket = int(lsp * dist)
            if bucket < lsp:
                dest = idxpairs[bucket]
                if len(dest) < p:
                    dest.append((linenum, ind + linenum))
        linenum += 1
        if linenum % 100 == 0:
            ubc = unfilled_bucket_count()
            if ubc == 0:
                break
            print(f"Unfilled bucket count: {unfilled_bucket_count()}. "
                  f"Total buckets to fill: {len(idxpairs)}. Line num: {ln}",
                  file=sys.stderr)
    if unfilled_bucket_count():
        print(f"Unfilled bucket count: {unfilled_bucket_count()}. "
              "Emitting unfinished results.", file=sys.stderr)
    else:
        print("No unfilled bucket count any more, we got 'em all.",
              file=sys.stderr)
    with open(a.outpath, "w") as f:
        for pairlist in idxpairs:
            for x, y in pairlist:
                f.write("%s\t%s\n" % (names[x], names[y]))
        f.write("=" * 80)
        f.write("#Fraction\tIndices of pairs [0-based]...\n")
        for values, frac in zip(idxpairs, np.linspace(0., 1., lsp + 1)[:-1]):
            f.write("%f\t%s\n" % (frac, '\t'.join(map(str, values))))
    return 0


if __name__ == "__main__":
    sys.exit(main())
