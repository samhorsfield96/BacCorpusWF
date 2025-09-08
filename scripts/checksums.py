#! /usr/bin/env python
import hashlib, hmac, os, stat, sys

# solution from https://unix.stackexchange.com/questions/35832/how-do-i-get-the-md5-sum-of-a-directorys-contents-as-one-sum

## Return the hash of the contents of the specified file, as a hex string
def file_hash(name):
    h = hashlib.sha256()
    with open(name, "rb") as f:  # <-- "rb" means read binary
        while True:
            buf = f.read(16384)
            if not buf:
                break
            h.update(buf)
    return h.hexdigest()

## Traverse the specified path and update the hash with a description of its name and contents
def traverse(h, path):
    rs = os.lstat(path)
    quoted_name = repr(path)
    if stat.S_ISDIR(rs.st_mode):
        h.update(('dir ' + quoted_name + '\n').encode())
        for entry in sorted(os.listdir(path)):
            traverse(h, os.path.join(path, entry))
    elif stat.S_ISREG(rs.st_mode):
        h.update(('reg ' + quoted_name + ' ').encode())
        h.update((str(rs.st_size) + ' ').encode())
        h.update((file_hash(path) + '\n').encode())
    else: pass # silently symlinks and other special files

def main():
    h = hashlib.sha256()
    for root in sys.argv[1:]: traverse(h, root)
    h.update(('end\n').encode())
    print(h.hexdigest())

if __name__ == "__main__":
    main()