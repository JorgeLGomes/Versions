############################################################################################
Para o uso do ETA no cluster Una, foi criada uma estrutura capaz de gerar, em
poucos passos, todo o modelo para execucao.

Existem 2 scripts principais (gera_set_parmeta.sh e buildall) que sao usados,
respectivamente, na criacao de um set_parmeta configurado para a resolucao
199X249X38 e na compilacao do modelo.

Para criar um set_parmeta execute:

./set_parmeta 2 10

Note que 2 e 10 sao, respectivamente, os valores de INPES e JNPES

O arquivo gerado sera: set_parmeta_199X249X38_4proc

Para compilar o modelo definido no set_parmeta acima use:

./buildall 199X249X38_4proc

Um diretorio ../199X249X38_4proc sera criado com os executaveis e demais
arquivos necessarios para a rodada.


Duvidas ou sugestoes:

Daniel M. Lamosa PAD/CPTEP/INPE 
email:lamosa@cptec.inpe.br
07/2007
############################################################################################

Diretorio eta/install

Alteracoes:

foi inserido nos script "configure" e "buildall" um argumento que define o diretorio onde os
executaveis serao gerados.
Com esta alteracao o script "configure" chama o arquivo de configuracao set_parmeta_{ARGUMENTO}

Exemplo:
Foi configurado dois modelos, o primeiro com 7 processadores e o segundo com 5. Para tanto foi gerado
dois arquivos de configuracao: set_parmeta_7proc e set_parmeta_5proc.
onde a diferenca entre os dois arquivos estah na definicao do numero de processadores utilizados para
a rodada.
set_parmeta_7proc: 
INPES=1
JNPES=7
#
set_parmeta_5proc: 
INPES=1
JNPES=5
#


O passo segunte serah compilar o modelo executando o comando abaixo:

buildall 7proc
o script irah criar um diretorio ./worketa_all_tmp/eta/7proc e ./worketa_all_tmp/eta/7proc/exe
e os executaveis gerados na compilacao serao armazenados no diretorio ./worketa_all_tmp/eta/7proc/exe

