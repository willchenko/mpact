Quickly find, curate, store, and retrieve relevant information on metabolic pathways to be compatible as an attachment to a a Genome-Scale Model (GEM). 

# mpact

## Setting up virtual environment
$MPACT_PATH should be set to the location of mpact folder
~~~
virtualenv  ~/.envs/mpact-env
source ~/.envs/mpact-env/bin/activate
cd $MPACT_PATH
pip install -r requirements.txt
~~~

## starting mpact
~~~
python3 
import mpact
from mpact import *
~~~

## Tutorial
Basic functions and workflow can be seen in "master_script.py".
Edit this file to suit your pathway and copy/paste functions into python terminal.
 