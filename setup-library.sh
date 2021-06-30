if [ ! -d "./git-actions-practice/" ]; then
    git clone https://github.com/KazegamiKuon/git-actions-practice.git    
fi
cd ./git-actions-practice/
git pull
cd ../
read -n 1 -p "will create environment? [Y/n]: " create_env
echo ""
if [ "$create_env" != "${create_env#[Yy]}" ] ;then
source ./git-actions-practice/conda-environment/automatically_initialize_environment.sh
fi