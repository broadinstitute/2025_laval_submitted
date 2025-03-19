conda create -y -n cytoself python=3.9

conda activate cytoself

# (Optional: Install pytorch GPU following the official instruction)
# For our machine, we install Pytorch with GPU support and cuda with version cu126 
pip3 install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cu126

## Install cytoself dependencies
pip install -e .

