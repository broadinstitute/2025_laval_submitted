### Test: 
conda create -n evo2 python=3.12 -y
conda activate evo2
pip install ipykernel

### Install cudnn first
conda install nvidia::cudnn -y
<!-- conda install conda-forge::cudnn -->
<!-- conda install conda-forge::cudnn python=3.12 sqlite>=3.45 -->

### Install cuda-toolkit too
conda install nvidia::cuda-toolkit=12.6.0 -y
<!-- conda install nvidia::cuda-toolkit=12.6.0 -y -->
<!-- conda install conda-forge::cuda-toolkit=12.6.0 -->

### Pip
pip install .

## DONE!!!

## Try running test: python ./test/test_evo2.py --model_name evo2_7b
Traceback (most recent call last):
  File "/home/shenrunx/igvf/varchamp/2025_laval_submitted/4_compare_AI_scores/2_analysis/scripts/evo2/./test/test_evo2.py", line 9, in <module>
    from evo2 import Evo2
  File "/home/shenrunx/software/anaconda3/envs/evo2/lib/python3.12/site-packages/evo2/__init__.py", line 1, in <module>
    from .models import Evo2
  File "/home/shenrunx/software/anaconda3/envs/evo2/lib/python3.12/site-packages/evo2/models.py", line 12, in <module>
    from vortex.model.model import StripedHyena
  File "/home/shenrunx/igvf/varchamp/2025_laval_submitted/4_compare_AI_scores/2_analysis/scripts/evo2/vortex/vortex/model/model.py", line 14, in <module>
    from vortex.model.layers import (
  File "/home/shenrunx/igvf/varchamp/2025_laval_submitted/4_compare_AI_scores/2_analysis/scripts/evo2/vortex/vortex/model/layers.py", line 10, in <module>
    from transformer_engine.pytorch import Linear
  File "/home/shenrunx/software/anaconda3/envs/evo2/lib/python3.12/site-packages/transformer_engine/__init__.py", line 10, in <module>
    import transformer_engine.common
  File "/home/shenrunx/software/anaconda3/envs/evo2/lib/python3.12/site-packages/transformer_engine/common/__init__.py", line 128, in <module>
    _NVRTC_LIB_CTYPES = _load_nvrtc()
                        ^^^^^^^^^^^^^
  File "/home/shenrunx/software/anaconda3/envs/evo2/lib/python3.12/site-packages/transformer_engine/common/__init__.py", line 111, in _load_nvrtc
    libs = subprocess.check_output("ldconfig -p | grep 'libnvrtc'", shell=True)
           ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/home/shenrunx/software/anaconda3/envs/evo2/lib/python3.12/subprocess.py", line 468, in check_output
    return run(*popenargs, stdout=PIPE, timeout=timeout, check=True,
           ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/home/shenrunx/software/anaconda3/envs/evo2/lib/python3.12/subprocess.py", line 573, in run
    raise CalledProcessError(retcode, process.args,
subprocess.CalledProcessError: Command 'ldconfig -p | grep 'libnvrtc'' returned non-zero exit status 1.

## Library not found, try to fix:
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/shenrunx/software/anaconda3/envs/evo2/lib
## Didn't work. Try to modify the /home/shenrunx/software/anaconda3/envs/evo2/lib/python3.12/site-packages/transformer_engine/common/__init__.py to bypass the import issue
python patch_transformer_engine.py
## Worked!!!

## New error:
pip3 install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cu126

### Original Installation guide from official git repo:
### =============================================================================================================
## Installation guide: https://github.com/ArcInstitute/evo2/tree/main
conda create -n evo2 python=3.12
conda activate evo2

pip install .
## Failed due to missing cudatoolkit
Error message:
Collecting transformer_engine_torch==2.1.0 (from transformer_engine[pytorch])
        Using cached transformer_engine_torch-2.1.0.tar.gz (147 kB)
        Preparing metadata (setup.py): started
        Preparing metadata (setup.py): finished with status 'error'
        error: subprocess-exited-with-error

        × python setup.py egg_info did not run successfully.
        │ exit code: 1
        ╰─> [17 lines of output]
            Traceback (most recent call last):
              File "<string>", line 2, in <module>
              File "<pip-setuptools-caller>", line 34, in <module>
              File "/tmp/pip-install-1wgujc59/transformer-engine-torch_40b66c5a95f14f2db911862943bf7a21/setup.py", line 46, in <module>
                setup_pytorch_extension(
              File "/tmp/pip-install-1wgujc59/transformer-engine-torch_40b66c5a95f14f2db911862943bf7a21/build_tools/pytorch.py", line 58, in setup_pytorch_extension
                cuda_architectures = cuda_archs()
                                     ^^^^^^^^^^^^
              File "/tmp/pip-install-1wgujc59/transformer-engine-torch_40b66c5a95f14f2db911862943bf7a21/build_tools/utils.py", line 193, in cuda_archs
                version = cuda_version()
                          ^^^^^^^^^^^^^^
              File "/tmp/pip-install-1wgujc59/transformer-engine-torch_40b66c5a95f14f2db911862943bf7a21/build_tools/utils.py", line 204, in cuda_version
                _, nvcc_bin = cuda_path()
                              ^^^^^^^^^^^
              File "/tmp/pip-install-1wgujc59/transformer-engine-torch_40b66c5a95f14f2db911862943bf7a21/build_tools/utils.py", line 186, in cuda_path
                raise FileNotFoundError(f"Could not find NVCC at {nvcc_bin}")
            FileNotFoundError: Could not find NVCC at /usr/local/cuda/bin/nvcc
            [end of output]

        note: This error originates from a subprocess, and is likely not a problem with pip.
      error: metadata-generation-failed

      × Encountered error while generating package metadata.
      ╰─> See above for output.

conda install -c nvidia cuda-toolkit
## Check nvcc is correctly installed
nvcc --version
nvcc: NVIDIA (R) Cuda compiler driver
Copyright (c) 2005-2025 NVIDIA Corporation
Built on Fri_Feb_21_20:23:50_PST_2025
Cuda compilation tools, release 12.8, V12.8.93
Build cuda_12.8.r12.8/compiler.35583870_0

## Output paths
export CUDA_HOME=/home/shenrunx/software/anaconda3/envs/evo2
export PATH=$PATH:$CUDA_HOME/bin
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$CUDA_HOME/lib

pip install .
## Still failed with the same error

## Try to install failed package directly
pip install transformer_engine[pytorch]==2.1.0
## More errors show up about new dependency issues
Missing nvtx3/nvToolsExt.h: This is the NVIDIA Tools Extension header file
Missing cudnn.h: This is the CUDA Deep Neural Network library header file

## Install missing dependencies
conda install -c nvidia cuda-nvtx cudnn
Error:
Channels:
 - defaults
 - conda-forge
 - nvidia
Platform: linux-64
Collecting package metadata (repodata.json): done
Solving environment: failed

InvalidSpec: The package "nvidia/linux-64::cuda-toolkit==12.8.1=0" is not available for the specified platform
### =============================================================================================================
### So many issues that couldn't be resolved, retry with a different order of installing packages to resolve them.
### Check out the final approach at the beginning of this installation document.