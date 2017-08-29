#!/usr/bin/env bash


# Define help function
function help(){
    echo "Anaconda3 Install - User install anaconda";
    echo "Usage example:";
    echo "anaconda3-install [(-h|--help)] [(-a|--anaconda3_base) string]";
    echo "Options:";
    echo "-h or --help: Display help information.";
    echo "-a or --anaconda3_base string: Base for anaconda3 installation.";
    exit 1;
}

# Declare vars. Flags initalizing to 0.
anaconda3_base="${HOME}/.local/";

# Execute getopt
ARGS=$(getopt -o "ha:" -l "help,anaconda3_base:" -n "Anaconda3 Install" -- "$@");

#Bad arguments
if [ $? -ne 0 ];
then
    help;
fi

eval set -- "$ARGS";

while true; do
    case "$1" in
        -h|--help)
            shift;
            help;
            ;;
        -a|--anaconda3_base)
            shift;
                    if [ -n "$1" ];
                    then
                        anaconda3_base="$1";
                        shift;
                    fi
            ;;

        --)
            shift;
            break;
            ;;
    esac
done

# Check required arguments
if [ -z "$anaconda3_base" ]
then
    echo "anaconda3_base is required";
    help;
fi


#########################################################
# 2. Install the base anaconda
#########################################################

echo "===================================================================="
echo ""
echo "Installing Anaconda3. This could take some time..."
export ANACONDA3_BASE="${anaconda3_base}/anaconda3"

curl -s -L https://repo.continuum.io/miniconda/Miniconda3-4.2.12-Linux-x86_64.sh  > miniconda.sh && \
  openssl md5 miniconda.sh | grep d0c7c71cc5659e54ab51f2005a8d96f3

chmod 777 miniconda.sh ; bash miniconda.sh -b -p $ANACONDA3_BASE
rm miniconda.sh && \
export PATH=${ANACONDA3_BASE}/bin:$PATH && \
conda config --set show_channel_urls True && \
#Add default packages
conda config --add create_default_packages setuptools && \
conda config --add create_default_packages ipython && \
#Add default channels
conda config --add channel nyuad-cgsb && \
conda config --add channels conda-forge && \
conda config --add channels defaults && \
conda config --add channels r && \
conda config --add channels bioconda && \
conda update --all --yes && \
conda clean -tipy

conda install -y anaconda-client

#########################################################
# 2. Configure Anaconda Environment
#########################################################

echo "===================================================================="
echo "Finished installing Anaconda3."
echo "It is located in ${ANACONDA3_BASE}."
echo ""

echo "The following are helpful environmental variables that you probably want to include in your ~/.bashrc or other shell file."
echo ""
echo "export PATH=${ANACONDA3_BASE}/bin:${PATH}"
echo ""

echo "===================================================================="
echo "To see the full list of nyuad-cgsb analysis modules - please head over to: "
echo "https://anaconda.org/nyuad-cgsb/environments"
echo "Or"
echo "https://jerowe.gitbooks.io/nyuad-gencore-hpc-modules/content/"

echo "Additionally, use the anaconda-client in order to search for packages."
echo "anaconda-client show nyuad-cgsb"

echo "To install a gencore env the command is:"
echo "export PATH=${ANACONDA3_BASE}/bin:${PATH} #Only if this is not in your ~/.bashrc or equivalent"
echo "conda env create -p nyuad-cgsb/gencore_variant_detection"
