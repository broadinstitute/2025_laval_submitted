import os
import re
import types
import glob

# Path to your conda environment
conda_env_path = '/home/shenrunx/software/anaconda3/envs/evo2'

# Directly find the path to the transformer_engine common init file
common_init_path = f"{conda_env_path}/lib/python3.12/site-packages/transformer_engine/common/__init__.py"
# Create a backup
backup_path = common_init_path + '.backup'

if os.path.exists(backup_path):
    print(f"Existing backup at {backup_path}, reading backup...")
    # Read the original file
    with open(backup_path, 'r') as f:
        content = f.read()
else:
    os.system(f'cp {common_init_path} {backup_path}')
    print(f"Created backup at {backup_path}")
    with open(common_init_path, 'r') as f:
        content = f.read()


# Create a new _load_nvrtc function
new_load_nvrtc = '''
def _load_nvrtc():
    """Load NVRTC shared library."""
    # Check conda environment first
    conda_env_path = '{}'
    extension = _get_sys_extension()
    
    # Look in conda environment
    conda_libs = glob.glob(f"{{conda_env_path}}/**/libnvrtc.{{extension}}*", recursive=True)
    conda_libs = list(filter(lambda x: not ("stub" in x or "libnvrtc-builtins" in x), conda_libs))
    if conda_libs:
        return ctypes.CDLL(conda_libs[0], mode=ctypes.RTLD_GLOBAL)
        
    # Continue with original search paths
    cuda_home = os.environ.get("CUDA_HOME") or os.environ.get("CUDA_PATH") or "/usr/local/cuda"
    libs = glob.glob(f"{{cuda_home}}/**/libnvrtc.{{extension}}*", recursive=True)
    libs = list(filter(lambda x: not ("stub" in x or "libnvrtc-builtins" in x), libs))
    libs.sort(reverse=True, key=os.path.basename)
    if libs:
        return ctypes.CDLL(libs[0], mode=ctypes.RTLD_GLOBAL)

    # Try ldconfig approach but don't fail if it doesn't work
    try:
        libs = subprocess.check_output("ldconfig -p | grep 'libnvrtc'", shell=True)
        libs = libs.decode("utf-8").split(\"\\n\")
        sos = []
        for lib in libs:
            if "stub" in lib or "libnvrtc-builtins" in lib:
                continue
            if "libnvrtc" in lib and "=>" in lib:
                sos.append(lib.split(">")[1].strip())
        if sos:
            return ctypes.CDLL(sos[0], mode=ctypes.RTLD_GLOBAL)
    except subprocess.CalledProcessError:
        pass  # Ignore errors from ldconfig

    # If all else fails, return a dummy object that won't crash the program
    print("Warning: Could not load libnvrtc. Some functionality may be limited.")
    class DummyLib:
        def __getattr__(self, name):
            return lambda *args, **kwargs: None
    return DummyLib()
'''.format(conda_env_path)


# Replace the problematic function
pattern = r'def _load_nvrtc\(\):.*?return ctypes\.CDLL\(.*?, mode=ctypes\.RTLD_GLOBAL\)'
content_new = re.sub(pattern, new_load_nvrtc, content, flags=re.DOTALL)

# Write the modified version
with open(common_init_path, 'w') as f:
    f.write(content_new)

print("Transformer Engine patched successfully!")