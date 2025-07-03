process TestGPU {
    label 'gpu_process'
    container 'docker://pytorch/pytorch:2.7.1-cuda11.8-cudnn9-runtime'
    script:
    """
    echo CUDA_VISIBLE_DEVICES=\${CUDA_VISIBLE_DEVICES:-}
    python -c "import torch; print('torch.cuda.is_available():', torch.cuda.is_available())"
    """
}

workflow {
    TestGPU()
}
