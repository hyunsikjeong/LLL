# Flatter Script

This is a script I use for [flatter](https://github.com/keeganryan/flatter).

I think the original source code itself is written by a CTF player. If you know who's the author of this short function, please let me know. I just want to put the name here...

## Script

```py
def flatter(M):
    from subprocess import check_output
    from re import findall
    z = "[[" + "]\n[".join(" ".join(map(str, row)) for row in M) + "]]"
    ret = check_output(["flatter"], input=z.encode())
    return matrix(M.nrows(), M.ncols(), map(int, findall(b"-?\\d+", ret)))
```