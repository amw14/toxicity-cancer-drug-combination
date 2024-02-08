#!/bin/bash

# install, create, and activate virtualenv called "env"
python3 -m pip install virtualenv
python3 -m venv env
module load python/3.7.4
module load cuda/11.3.1
module load cudnn/8.2.0
source env/bin/activate

# update pip and install required packages
pip install -U pip
pip install -r requirements.txt

echo created venv
