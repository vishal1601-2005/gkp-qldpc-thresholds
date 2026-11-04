"""Compile the C++ extension in place:  python cpp/build.py   (needs g++/clang++ and `pip install pybind11`).
Writes src/gkp_qldpc/_fastbp*.so (or .pyd). OpenMP is used if available."""
import subprocess, sys, sysconfig, pathlib, shutil
import pybind11
root = pathlib.Path(__file__).resolve().parent
dest = root.parent / "src" / "gkp_qldpc"
suffix = sysconfig.get_config_var("EXT_SUFFIX")
cxx = shutil.which("g++") or shutil.which("clang++") or "c++"
out = dest / f"_fastbp{suffix}"
base = [cxx, "-O3", "-std=c++17", "-shared", "-fPIC", "-DNDEBUG",
        f"-I{pybind11.get_include()}", f"-I{sysconfig.get_paths()['include']}", f"-I{root/'include'}", str(root / "src" / "pybind_module.cpp"), "-o", str(out)]
for flags in (["-fopenmp", "-march=native"], ["-fopenmp"], []):
    r = subprocess.run(base + flags, capture_output=True, text=True)
    if r.returncode == 0:
        print("built", out, "with", " ".join(flags) if flags else "no extra flags"); sys.exit(0)
print(r.stderr); sys.exit(1)
