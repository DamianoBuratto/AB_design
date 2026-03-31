#######################################################################
#
# Script with step-by-step try and error installaiton for RFantibody
#
#########################################################################

source /public/software/profile.d/cuda12.2-env.sh
# module load cuda/12.2 # this has the same result as source /public/software/profile.d/cuda12.2-env.sh
source /public/software/profile.d/compiler_intel-2021.3.0.sh

conda create -n RFantibody_HPC python=3.10 -y
source /public/home/damiano/Programs/miniconda3/bin/activate RFantibody_HPC

conda install -c conda-forge numpy scipy pandas -y
conda install -c conda-forge make cmake ninja gcc_linux-64 gxx_linux-64 -y
conda install -y pytorch=2.3.* torchvision torchaudio torchdata pytorch-cuda=12.1 -c pytorch -c nvidia

# pip install dgl-2.4.0+cu121-cp310-cp310-manylinux1_x86_64.whl # searching for connection and failing
# conda install -c conda-forge setuptools dgl -y # downgraiding pytorch to CPU-ONLY

# pip install -i https://mirrors.zju.edu.cn/pypi/web/simple dgl # DGL import failed: No module named 'setuptools.extern'
# conda install setuptools==65.5.0  # solution to setuptools.extern (downgrade from v80) -> didn't work

# download dgl-2.5.0+cu121-cp310-cp310-manylinux1_x86_64.whl from https://data.dgl.ai/wheels/torch-2.3/cu121/repo.html
pip install dgl-2.5.0+cu121-cp310-cp310-manylinux1_x86_64.whl

# Now you can test that torch-cuda and dgl are correctly installed
# Define the Python launcher with the GLIBC shim so RFantibody finds the right runtime
PYTHON="/public/software/mathlib/glibc/2.34/lib/ld-linux-x86-64.so.2 --library-path \"/public/software/mathlib/glibc/2.34/lib:/public/home/jiangyw/miniconda3/envs/AF3v3/lib:$LD_LIBRARY_PATH\" /public/home/xuziyi/miniconda/envs/RFantibody_HPC/bin/python"
# You use the test.slurm script to run the test
# check the install folder at the WS (/data/RFantibody/install)
# Example: eval "$PYTHON analyze_rf2_pdb.py"

# compile the USalign like in the setup.sh script
cd USalign &&
make &&

# Install all Python dependencies from poetry.lock
conda install -c conda-forge annotated-types antlr4-python3-runtime asttokens certifi charset-normalizer colorama cuda-python cython e3nn exceptiongroup executing filelock fsspec hydra-core icecream idna iniconfig jinja2 markupsafe mpmath networkx omegaconf opt-einsum packaging pillow pluggy psutil pydantic pydantic-core pygments pyrsistent pytest python-dateutil pytz pyyaml requests six sympy tomli tqdm typing-extensions tzdata  urllib3 -y
# conda install -c conda-forge opt-einsum-fx # This cannot be found in the ZJU mirrors, we can try without it

# Install MKL dependencies 
conda install -c conda-forge mkl intel-openmp tbb -y

# Change the way to start your RFantibody script in a similar way to the test by calling $PYTHON



# ziyi test:
# on HPC, use slurm debug to test
# rfdiffusion_inference.py: No module named 'hydra.core'
pip install --upgrade --force-reinstall hydra-core==1.3.2 omegaconf==2.3.0 
# rfdiffusion_inference.py: /lib64/libm.so.6: version `GLIBC_2.27' not found (required by /public/home/xuziyi/miniconda/envs/RFantibody_HPC/lib/python3.10/site-packages/dgl/libdgl.so)
# File "/public/home/xuziyi/miniconda/envs/RFantibody_HPC/lib/python3.10/site-packages/torch/cuda/nvtx.py", line 12, in _fail     raise RuntimeError(untimeError: NVTX functions not installed. Are you sure you have a CUDA build?
pip install nvidia-nvtx-cu12==12.1.105 # after this still reported NVTX error
# Node diagnostics show PyTorch (2.3.0) is present but CUDA is not available (torch.cuda.is_available() = False).
# torch.cuda.nvtx can be imported, but calling nvtx APIs (e.g., nvtx.range or range_push) raises RuntimeError: "NVTX functions not installed. Are you sure you have a CUDA build?" — this causes SE3Transformer's with nvtx_range(...) to crash RFdiffusion.

# Pre-install a no-op torch.cuda.nvtx module into sys.modules before importing torch. This ensures imports like from torch.cuda.nvtx import range return a safe no-op implementation rather than a backend that will raise at call time.
# After importing torch, perform a runtime probe; if the real torch.cuda.nvtx is present but raises on call, replace it with the same no-op module so downstream calls are safe.
# Added an opt-out via environment variable: set RFANTIBODY_ALLOW_NVTX=1 to skip the shim and use the real NVTX (useful for profiling when the node truly supports it).


# Professor fixed:
# torch.cuda was failing before, pytorch was running the cpu version from the default channel, I had to force the installation from pytorch channel to have the cuda version

# ziyi test:
'''
[2025-11-09 14:33:48,258][rfantibody.rfdiffusion.inference.utils][INFO] - Sampled motif RMSD: 3.83
Error executing job with overrides: ['antibody.target_pdb=/public/home/xuziyi/RFantibody/data/p53_R175H_cut.pdb', 'antibody.framework_pdb=/public/home/xuziyi/RFantibody/scripts/examples/example_inputs/hu-4D5-8_Fv.pdb', 'inference.ckpt_override_path=/public/home/xuziyi/RFantibody/weights/RFdiffusion_Ab.pt', 'antibody.design_loops=[L1:9-12,L2:6-8,L3:9-12,H1:6-8,H2:5-7,H3:8-15]', 'inference.num_designs=100', 'inference.output_prefix=/public/home/xuziyi/RFantibody/outputs/rfdiffusion/design', 'ppi.hotspot_res=[C7,C8,C9,A73,A74,A76,A77,A80,A84,A143,A146,A147,A150]', 'diffuser.T=200', 'inference.cautious=False']
Traceback (most recent call last):
  File "<string>", line 147, in main
  File "/public/home/xuziyi/RFantibody/src/rfantibody/rfdiffusion/inference/model_runners.py", line 642, in sample_step
    x_t_1, seq_t_1, tors_t_1, px0 = self.denoiser.get_next_pose(
  File "/public/home/xuziyi/RFantibody/src/rfantibody/rfdiffusion/inference/utils.py", line 861, in get_next_pose
    frames_next = get_next_frames(xt, px0, t, diffuser=self.diffuser,
  File "/public/home/xuziyi/RFantibody/src/rfantibody/rfdiffusion/inference/utils.py", line 119, in get_next_frames
    R_0 = scipy_R.from_matrix(R_0.squeeze().numpy()).as_matrix()
  File "_rotation.pyx", line 1136, in scipy.spatial.transform._rotation.Rotation.from_matrix
ValueError: Non-positive determinant (left-handed or null coordinate frame) in rotation matrix 237: [[0. 0. 0.]
 [0. 0. 0.]
 [0. 0. 0.]].

Set the environment variable HYDRA_FULL_ERROR=1 for a complete stack trace.
[2025-11-09 22:33:49] Error: RFdiffusion execution failed
'''






 
