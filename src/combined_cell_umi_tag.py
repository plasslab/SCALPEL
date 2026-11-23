#!/usr/bin/env python3
"""SAM stdin -> stdout with temporary XU=CB+UB; keep cells distinct at dedup."""
import json
import sys


def main():
    kept = skipped = 0
    for line in sys.stdin:
        if line.startswith('@'):
            sys.stdout.write(line)
            continue
        fields = line.rstrip('\n').split('\t')
        tags = {t[:2]: t[5:] for t in fields[11:] if t.startswith(('CB:Z:', 'UB:Z:'))}
        if not tags.get('CB') or not tags.get('UB'):
            skipped += 1
            continue
        fields = fields[:11] + [t for t in fields[11:] if not t.startswith('XU:')]
        fields.append('XU:Z:' + tags['CB'] + '_' + tags['UB'])
        sys.stdout.write('\t'.join(fields) + '\n')
        kept += 1
    print(json.dumps({'CB_UB_tagged_reads': kept, 'reads_without_CB_or_UB': skipped}), file=sys.stderr)


if __name__ == '__main__':
    main()
