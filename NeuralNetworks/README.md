# Neural Network codes

:heavy_exclamation_mark: *Warning: don't call `cmsenv`, as this will conflict with user-installed packages...*

## Setup

*(Instructions adapted from David Walter's DeepPotato code.)*

- First, install conda (anaconda or the lightweight miniconda):

Download the installer (for anaconda, check the newest version on [www.anaconda.com](https://www.anaconda.com))
```
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh DIRPATH
```
or
```
wget https://repo.anaconda.com/archive/Anaconda3-2018.12-Linux-x86_64.sh DIRPATH
```

Install it.
```
bash DIRPATH/Miniconda3-latest-Linux-x86_64.sh
```
or
```
bash DIRPATH/Anaconda3-2018.12-Linux-x86_64.sh
```

- Update your .bashrc (or equivalent) file:

Set paths to your conda.
```
source ~/.conda_init
```
with `conda_init` sourcing/adding to path the relevant files:
```
# >>> conda initialize >>>
# !! Contents within this block are managed by 'conda init' !!
__conda_setup="$('/nfs/dust/cms/user/ntonon/miniconda3/bin/conda' 'shell.zsh' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f "/nfs/dust/cms/user/ntonon/miniconda3/etc/profile.d/conda.sh" ]; then
        . "/nfs/dust/cms/user/ntonon/miniconda3/etc/profile.d/conda.sh"
    else
        export PATH="/nfs/dust/cms/user/ntonon/miniconda3/bin:$PATH"
    fi
fi
unset __conda_setup
# <<< conda initialize <<<
```

*NB: there might be issues if the PYTHONPATH is set.
If needed, execute these lines whenever you use conda:*
```
export PYTHONPATH=
export PATH=/your/path/to/conda/bin:$PATH
```

- Setup conda, install tensorflow and all the necessary packages. In order to use GPU support, the cuda environment has to be set beforehand.
This is the case e.g. on 'naf-cms-gpu01.desy.de'.
```
conda install matplotlib pandas pytables
conda install -c conda-forge root_numpy

conda install tensorflow-gpu
conda install keras
conda install scikit-learn
conda install shap
conda install seaborn
```
*NB1: Python tells you if any other package is missing.*
*NB2: order may matter.*
*NB3: basic ROOT should be installed from root_numpy package.*

## Train, test, evaluate

- The main code `Train_Neural_Network.py` calls all the main functions. The user is able to tune a large choice of options there, and then shall simply run the code:
```
python Train_Neural_Network.py
```

- For SM vs EFT strategies, you may run the code `Train_Neural_Network.py` on a previously-trained NN in order to only re-evaluate its performance over user-defined data (options must be set directly in the code):
```
python StandaloneValidation.py
```

## HTCondor

*Not yet possible...*
