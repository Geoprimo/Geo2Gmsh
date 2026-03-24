#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import hashlib
from pathlib import Path

# -----------------------------
# SHA-256 function
# -----------------------------
def sha256sum(file_path):
    h = hashlib.sha256()
    with open(file_path, "rb") as f:
        for chunk in iter(lambda: f.read(8192), b""):
            h.update(chunk)
    return h.hexdigest()

# -----------------------------
# Reference checksums (embedded)
# -----------------------------
EXPECTED = {
    "examples/overview/outputs/overview.vtk":
        "fc16e53733b8145319893bafb4dcf2e7b7dc69f22c8d3ae91560312837d2ec3d",

    "examples/test_llanos/outputs/test_llanos.vtk":
        "9fb2da0775e8e4d47adfbe43442abe2fbd11bf3c8c3caf414d7ef4bee9a8666f",

    "examples/test_ringvent/outputs/test_ringvent.vtk":
        "231945cab11f10c2262ad6f29af4977820936c5ccebe53f3e4cc9fb4df106986",

    "examples/overview/outputs/overview.msh":
        "88026dde738fbcda328cb67e68e72cc55f56251adc609713918bdbe1af504381",

    "examples/test_llanos/outputs/test_llanos.msh":
        "99960cf09ad678e7acbaab19c5d136b1da912c8c9c31c1de296e7e7eb5c05f86",

    "examples/test_ringvent/outputs/test_ringvent.msh":
        "09b879a272ff6110b719fbdab647f64b61a9b676750ebfcbf6a970e9f124a31a",
}

# -----------------------------
# Validation
# -----------------------------
def main():
    print("\n🔍 Geo2Gmsh checksum validation\n")

    all_ok = True

    for path_str, expected_hash in EXPECTED.items():
        path = Path(path_str)

        if not path.exists():
            print(f"❌ MISSING FILE: {path}")
            all_ok = False
            continue

        actual_hash = sha256sum(path)

        if actual_hash == expected_hash:
            print(f"✅ OK      {path}")
        else:
            print(f"❌ FAIL    {path}")
            print(f"   expected: {expected_hash}")
            print(f"   got:      {actual_hash}")
            all_ok = False

    print("\n-----------------------------")
    if all_ok:
        print("🎉 ALL CHECKS PASSED")
        exit(0)
    else:
        print("⚠️ SOME CHECKS FAILED")
        exit(1)

if __name__ == "__main__":
    main()

