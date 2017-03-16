apt-get update -y
apt-get install -y build-essential python-pip  python-dev  python-numpy python-scipy wget python-networkx

# install python2.7
pip install pythonbrew
pythonbrew_install && source "/root/.pythonbrew/etc/bashrc" && pythonbrew install 2.7.3
