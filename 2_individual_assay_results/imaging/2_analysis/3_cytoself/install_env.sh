conda create -y -n cytoself python=3.9

conda activate cytoself

# (Optional: Install pytorch GPU following the official instruction)
# For our machine, we install Pytorch with GPU support and cuda with version cu126 
pip3 install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cu126

## Install cytoself dependencies
pip install -e .

## Install tools for easy downloading and process data from CellPainting Gallery
## https://github.com/broadinstitute/cpg/tree/main/cpgdata
pip install cpgdata
cpg index sync -o "cp_gallery_index_files" ## cp_gallery_index_files here is a symbolic link to the directory in which the CP gallery index files are stored
