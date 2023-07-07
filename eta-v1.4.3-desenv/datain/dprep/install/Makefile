
.SILENT:

ROOT   = /home/jorge/EtaLab/eta-v1.4.3/datain/dprep
SRC    = $(ROOT)/src
EXEDIR = $(ROOT)/exe


EXEDIRS   =	$(SRC)/dgetacpt_eta15   \
		$(SRC)/dgetacpt_eta40   \
		$(SRC)/dglobal_cptec    \
		$(SRC)/dggposeta        \
		$(SRC)/dglobal_bamclim  \
		$(SRC)/dggfs2gr0.5      \
		$(SRC)/dggfs2gr0.25     \
		$(SRC)/dgcfs            \
		$(SRC)/dgERA5 
exe:
	@for dir in $(EXEDIRS) ; \
	  do \
	  echo " "; \
	  echo Making executable from directory $$dir ;\
	  (cd $$dir; \
	  $(MAKE) exe "EXEDIR=$(EXEDIR)"; if [ $$? != 0 ] ; then \
	      echo "Exit status from make was $$?" ; exit 1 ; fi ;) ;\
	done
