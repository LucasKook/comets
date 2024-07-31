
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

timing:
	(cd .. && Rscript --vanilla inst/code/timing.R > ccle.out)

clean:
	rm ../Rplots.pdf ../*.out *.out

all: dependencies fig ccle multiomics mimic
