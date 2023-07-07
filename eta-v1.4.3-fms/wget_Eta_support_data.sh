wget -c -t 10  http://ftp1.cptec.inpe.br/pesquisa/grpeta/VII-WorkEta/model/Eta_support_data_worketa.tgz
tar -zxvf Eta_support_data_worketa.tgz
mv Eta_support_data_worketa Eta_support_data
cd Eta_support_data/static/topo
gzip -d U*
cd ../../..
rm -f Eta_support_data_worketa.tgz
