import torch
import os

import hydra
from hydra.core.hydra_config import HydraConfig

import rfantibody.rf2.modules.util as util
import rfantibody.rf2.modules.pose_util as pu
from rfantibody.rf2.modules.model_runner import AbPredictor
from rfantibody.rf2.modules.preprocess import pose_to_inference_RFinput, Preprocess

# Dynamically determine config path based on installation location
def get_config_path():
    """Find RF2 config directory relative to this script or from environment"""
    # Try environment variable first (for flexible deployment)
    if 'RF2_CONFIG_PATH' in os.environ:
        return os.environ['RF2_CONFIG_PATH']
    
    # Try relative to script location
    script_dir = os.path.dirname(os.path.abspath(__file__))
    project_root = os.path.dirname(script_dir)  # Go up from scripts/ to project root
    relative_config = os.path.join(project_root, 'src/rfantibody/rf2/config')
    
    if os.path.exists(relative_config):
        return relative_config
    
    # Fallback to Docker path
    return '/home/src/rfantibody/rf2/config'

@hydra.main(version_base=None, config_path=get_config_path(), config_name='base')
def main(conf: HydraConfig) -> None:
    """
    Main function
    """
    print(f'Running RF2 with the following configs: {conf}')
    done_list=util.get_done_list(conf)
    device = torch.device('cuda:0' if torch.cuda.is_available() else 'cpu')
    preprocessor=Preprocess(pose_to_inference_RFinput, conf)
    predictor=AbPredictor(conf, preprocess_fn=preprocessor, device=device)
    for pose, tag in pu.pose_generator(conf):
        if tag in done_list and conf.inference.cautious:
            print(f'Skipping {tag} as output already exists')
            continue
        predictor(pose, tag)

if __name__ == '__main__':
    main()
