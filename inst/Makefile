
fig:
	(cd .. && Rscript --vanilla inst/figures/components-overview.R > components-overview.out)

dependencies:
	(cd .. && Rscript --vanilla inst/code/dependencies.R > dependencies.out)

mimic:
	(cd .. && Rscript --vanilla inst/code/mimic.R > mimic.out)

multiomics:
	(cd .. && Rscript --vanilla inst/code/multiomics.R > multiomics.out)

ccle:
	(cd .. && Rscript --vanilla inst/code/ccle.R > ccle.out)

dep:
	bash dispatch.sh run-install.sh 1200

mim:
	bash dispatch.sh run-mimic.sh 21600

mul:
	bash dispatch.sh run-multiomics.sh 3600

sig:
	bash dispatch.sh run-ccle.sh 3600

clean:
	rm ../Rplots.pdf ../*.out

all: dependencies fig ccle multiomics mimic
