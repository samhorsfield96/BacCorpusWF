#!/usr/bin/env python
import hashlib, os, stat, sys

## Return the hash of the contents of the specified file, as a hex string
def file_hash(name):
    h = hashlib.sha256()
    with open(name, "rb") as f:  # read binary
        while True:
            buf = f.read(16384)
            if not buf:
                break
            h.update(buf)
    return h.hexdigest()

## Traverse the specified path and update the hash with relative names + contents
def traverse(h, path, base):
    rs = os.lstat(path)
    relpath = os.path.relpath(path, base)  # relative path instead of absolute
    if stat.S_ISDIR(rs.st_mode):
        h.update(('dir ' + relpath + '\n').encode())
        for entry in sorted(os.listdir(path)):
            traverse(h, os.path.join(path, entry), base)
    elif stat.S_ISREG(rs.st_mode):
        h.update(('reg ' + relpath + ' ').encode())
        h.update((str(rs.st_size) + ' ').encode())
        h.update((file_hash(path) + '\n').encode())
    else:
        pass  # ignore symlinks and special files

def main():
    h = hashlib.sha256()
    for root in sys.argv[1:]:
        traverse(h, root, root)  # pass root as base
    h.update(('end\n').encode())
    print(h.hexdigest())

if __name__ == "__main__":
    main()
