def printE(name, f=0, unit=""):
    if type(name) is not str:
        f = name
        name = None
    print(str(f"{name}: " if name else "") + format(float(f), "e") + (f" {unit}" if unit else ""))
    return f
