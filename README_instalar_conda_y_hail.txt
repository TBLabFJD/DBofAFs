# todo desde nodo login uam
cd /home/cserrano

#descargar instalador
curl -O https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh

#instalar en /home/cserrano/miniconda3
bash ~/Miniconda3-latest-Linux-x86_64.sh

# source bash rc
source ~/.bashrc

# si hago cat .bashrc puedo ver que se ha metido el bloque de conda de inicializacion

# check conda list (paquetes de conda) -> poner todo el path
/home/cserrano/miniconda3/bin/conda list

## aceptar unos permisos 
/home/cserrano/miniconda3/bin/conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main
/home/cserrano/miniconda3/bin/conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r

# create the hail environment:
# he visto que para el hail que yo instale en su dia me pedian python=3.8 y java 8 
#en: cat /home/graciela/anaconda3/envs/hail/conda-meta/history

/home/cserrano/miniconda3/bin/conda create -y -n hail python=3.8 openjdk=8 -c defaults

# activar el environment
source /home/cserrano/miniconda3/bin/activate hail

# instalar hail dentro del environment
pip install hail==0.2.120
pip show hail

#cerrar entorno
conda deactivate
